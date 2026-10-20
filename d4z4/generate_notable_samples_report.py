#!/usr/bin/env python3
"""
Generate Notable Samples Report

Creates paginated SVG/PDF reports for notable D4Z4 samples by category:
- FSHD1-like samples (short array + β-satellite terminal)
- Samples without any arrays
- Multi-array samples
- Truncated arrays

Each page contains 4-5 samples with legend repeated.

Usage:
    python generate_notable_samples_report.py
    python generate_notable_samples_report.py --output-dir results/notable_samples
    python generate_notable_samples_report.py --samples-per-page 5 --format pdf

For PDF output on macOS with homebrew cairo:
    DYLD_LIBRARY_PATH=/opt/homebrew/lib python generate_notable_samples_report.py --format pdf
"""

import argparse
import json
import os
import subprocess
import sys
from collections import defaultdict
from io import BytesIO

# Fix library path for cairo on macOS (homebrew)
if sys.platform == 'darwin':
    homebrew_lib = '/opt/homebrew/lib'
    if os.path.exists(homebrew_lib):
        current_path = os.environ.get('DYLD_LIBRARY_PATH', '')
        if homebrew_lib not in current_path:
            os.environ['DYLD_LIBRARY_PATH'] = f"{homebrew_lib}:{current_path}"

import pandas as pd

# Import the plotting module
import plot_d4z4_arrays as plotter


def parse_args():
    parser = argparse.ArgumentParser(description="Generate notable samples report")
    parser.add_argument("--analysis-file",
                        default="results/d4z4_array_analysis.tsv",
                        help="D4Z4 analysis results TSV")
    parser.add_argument("--output-dir", default="results/notable_samples",
                        help="Output directory for reports")
    parser.add_argument("--samples-per-page", type=int, default=4,
                        help="Samples per page (default: 4)")
    parser.add_argument("--format", choices=["svg", "pdf", "both"], default="svg",
                        help="Output format (default: svg)")
    parser.add_argument("--background", choices=["white", "black"], default="white",
                        help="Background color (default: white)")
    parser.add_argument("--width", type=int, default=1000,
                        help="Page width in pixels (default: 1000)")
    return parser.parse_args()


def get_notable_samples(analysis_file):
    """Extract notable samples by category and chromosome from analysis results."""

    df = pd.read_csv(analysis_file, sep='\t')

    categories = {}

    # Helper functions for complex categories
    def has_mixed_arrays(row):
        if pd.isna(row.get('arrays_info')) or row['num_arrays'] <= 1:
            return False
        try:
            arrays = json.loads(row['arrays_info'])
            types = set(a.get('type') for a in arrays)
            return 'canonical' in types and 'degraded' in types
        except:
            return False

    def has_only_degraded(row):
        if pd.isna(row.get('arrays_info')) or row['num_arrays'] == 0:
            return False
        try:
            arrays = json.loads(row['arrays_info'])
            types = set(a.get('type') for a in arrays)
            return types == {'degraded'}
        except:
            return False

    def has_solitary_plam(row):
        if pd.isna(row.get('clusters_info')):
            return False
        try:
            clusters = json.loads(row['clusters_info'])
            return any(c.get('type') == 'solitary_plam' for c in clusters)
        except:
            return False

    # Define category masks
    category_masks = {
        'FSHD1-like': (
            (~df['is_truncated']) &
            (df['repeat_count_main'] <= 10) &
            (df['terminal_feature'] == 'β-satellite terminal')
        ),
        'No arrays': df['num_arrays'] == 0,
        'Multi-array': df['num_arrays'] > 1,
        'Mixed': df.apply(has_mixed_arrays, axis=1),
        'Only degraded': df.apply(has_only_degraded, axis=1),
        'Solitary pLAM': df.apply(has_solitary_plam, axis=1),
        'Previously truncated': df['array_length'] > 250000,
    }

    # Separate each category by chromosome (chr4 vs chr10)
    for category_name, mask in category_masks.items():
        category_df = df[mask]

        # chr4 samples
        chr4_haps = category_df[category_df['chromosome'] == 'chr4']['haplotype'].tolist()
        if chr4_haps:
            categories[f'{category_name} (chr4)'] = chr4_haps

        # chr10 samples
        chr10_haps = category_df[category_df['chromosome'] == 'chr10']['haplotype'].tolist()
        if chr10_haps:
            categories[f'{category_name} (chr10)'] = chr10_haps

    # C- on chr10 (10qB - rare, consistent with 4q->10q subtelomeric transfer)
    cminus_chr10 = df[
        (df['chromosome'] == 'chr10') &
        (df['terminal_feature'] == 'No β-satellite terminal')
    ]['haplotype'].tolist()
    if cminus_chr10:
        categories['10qB C- terminal (chr10)'] = cminus_chr10

    # Functional pLAM on chr10 (rare - potential 4q/10q translocation)
    plam_file = os.path.join(os.path.dirname(analysis_file), 'd4z4_array_coordinates.plam.tsv')
    if os.path.exists(plam_file):
        plam_df = pd.read_csv(plam_file, sep='\t')
        chr10_func = plam_df[(plam_df['chromosome'] == 'chr10') & (plam_df['plam_functional'] == 'func_Plam')]
        if len(chr10_func) > 0:
            func_haps = []
            for _, row in chr10_func.iterrows():
                match = df[(df['sample'] == row['sample']) &
                           (df['chromosome'] == 'chr10') &
                           (df['hap_num'] == row['hap'])]
                if len(match) > 0:
                    func_haps.append(match.iloc[0]['haplotype'])
            if func_haps:
                categories['Functional pLAM (chr10)'] = func_haps

    return categories


def paginate(items, per_page):
    """Split list into pages."""
    for i in range(0, len(items), per_page):
        yield items[i:i + per_page]


def generate_category_report(category_name, haplotypes, output_dir, args):
    """Generate paginated report for a single category."""

    if not haplotypes:
        print(f"  {category_name}: No samples")
        return []

    pages = list(paginate(haplotypes, args.samples_per_page))
    output_files = []

    print(f"  {category_name}: {len(haplotypes)} samples -> {len(pages)} page(s)")

    # Sanitize category name for filename
    safe_name = category_name.lower().replace(' ', '_').replace('-', '_').replace('(', '').replace(')', '')

    for page_num, page_haplotypes in enumerate(pages, 1):
        # Build output filename
        if len(pages) == 1:
            svg_file = os.path.join(output_dir, f"{safe_name}.svg")
        else:
            svg_file = os.path.join(output_dir, f"{safe_name}_page{page_num}.svg")

        # Build command to run plot_d4z4_arrays.py
        cmd = [
            sys.executable, 'plot_d4z4_arrays.py',
            '--haplotypes'] + page_haplotypes + [
            '--output', svg_file,
            '--background', args.background,
            '--width', str(args.width),
            '--adaptive-window'
        ]

        # Run the plotting command
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    Error on page {page_num}: {result.stderr}")
            continue

        output_files.append(svg_file)
        print(f"    Page {page_num}: {svg_file}")

    return output_files


def convert_svgs_to_pdf(svg_files, pdf_output):
    """Convert multiple SVGs to a single multi-page PDF using cairosvg."""
    # Try to import cairosvg
    try:
        import cairosvg
    except (ImportError, OSError) as e:
        print(f"\nCannot load cairosvg: {e}")
        print("To generate PDF, ensure cairo is installed: brew install cairo")
        return False

    # Try to import a PDF merger
    PdfWriter = None
    PdfMerger = None
    try:
        from pypdf import PdfWriter, PdfReader
    except ImportError:
        try:
            from PyPDF2 import PdfMerger, PdfReader
        except ImportError:
            print("\nTo merge PDFs, install: pip install pypdf")
            return False

    pdf_pages = []

    for svg_file in svg_files:
        # Convert SVG to PDF
        try:
            pdf_bytes = cairosvg.svg2pdf(url=svg_file)
            pdf_pages.append(BytesIO(pdf_bytes))
        except Exception as e:
            print(f"  Error converting {svg_file}: {e}")
            continue

    if not pdf_pages:
        print("  No pages converted successfully")
        return False

    # Merge PDFs
    if PdfWriter is not None:
        writer = PdfWriter()
        for pdf_io in pdf_pages:
            reader = PdfReader(pdf_io)
            for page in reader.pages:
                writer.add_page(page)
        with open(pdf_output, 'wb') as f:
            writer.write(f)
    else:
        merger = PdfMerger()
        for pdf_io in pdf_pages:
            merger.append(pdf_io)
        merger.write(pdf_output)
        merger.close()

    print(f"\nPDF saved: {pdf_output}")
    return True


def main():
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("Notable Samples Report Generator")
    print("=" * 60)
    print(f"Analysis file: {args.analysis_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Samples per page: {args.samples_per_page}")
    print(f"Background: {args.background}")
    print()

    # Get notable samples by category
    print("Extracting notable samples...")
    categories = get_notable_samples(args.analysis_file)

    # Generate reports for each category
    print("\nGenerating reports:")
    all_svg_files = {}

    for category_name, haplotypes in categories.items():
        svg_files = generate_category_report(
            category_name, haplotypes, args.output_dir, args
        )
        if svg_files:
            all_svg_files[category_name] = svg_files

    # Convert to PDF if requested
    if args.format in ("pdf", "both"):
        print("\nConverting to PDF...")
        for category_name, svg_files in all_svg_files.items():
            safe_name = category_name.lower().replace(' ', '_').replace('-', '_').replace('(', '').replace(')', '')
            pdf_output = os.path.join(args.output_dir, f"{safe_name}.pdf")
            convert_svgs_to_pdf(svg_files, pdf_output)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for category_name, haplotypes in categories.items():
        print(f"  {category_name}: {len(haplotypes)} samples")

    total_pages = sum(
        (len(haps) + args.samples_per_page - 1) // args.samples_per_page
        for haps in categories.values() if haps
    )
    print(f"\nTotal pages: {total_pages}")
    print(f"Output directory: {args.output_dir}")


if __name__ == '__main__':
    main()
