#!/usr/bin/env python3
"""
Figure 5: D4Z4 macrosatellite structural variation resolved by KaryoScope.

Panel (a): Locus schematic with real KaryoScope track example (HG01358 hap1)
Panel (b): Repeat count distribution (violin) + structural config frequencies (waffle/bar)
Panel (c): Multi-array examples (HG01109 chr4, HG02717 chr10)
Panel (d): Degraded/missing examples (HG01928 chr10, HG00290 chr4)
Panel (e): Haplotype asymmetry scatter + exchange signatures bar

Output: results/karyoscope_companion/
"""

import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results' / 'karyoscope_companion'

# Register Barthel fonts
FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
if FONT_DIR.exists():
    for font_file in FONT_DIR.glob("BasicSans-*.otf"):
        fm.fontManager.addfont(str(font_file))
    plt.rcParams['font.family'] = 'Arial'
plt.rcParams['svg.fonttype'] = 'none'

# Barthel brand colors
C_GREEN = '#40D392'
C_BLUE = '#60A5FA'
C_CORAL = '#F07167'
C_YELLOW = '#FBBF24'
C_EMERALD = '#10B981'
C_ROYAL = '#3B82F6'
C_GRAY = '#545454'
C_LAVENDER = '#C4A9E8'

# Array type colors (consistent with plot_d4z4_arrays.py)
C_CPLUS = '#90EE90'   # C+ canonical with β-satellite
C_CMINUS = '#87CEEB'  # C- canonical without
C_DEGRADED = '#FF6B6B'  # Degraded

# Config category colors — based on original Barthel palette, more diverging
# Normal alleles (shown at 50% alpha in waffle, except 10qB on chr10)
CONFIG_COLORS = {
    'A-type (C+)': '#40D392',      # Green (Accent 1) — normal A-type
    'B-type (C-)': '#60A5FA',      # Blue (Accent 2) — normal B-type
    # Structural abnormalities (full opacity, maximally diverging)
    '10qB (C-)': '#3B82F6',        # Royal Blue (Accent 6) — abnormal on chr10
    'Functional pLAM': '#10B981',   # Emerald (Accent 5)
    'Multi-array': '#C4A9E8',      # Lavender (Light 2)
    'Mixed': '#FBBF24',            # Yellow (Accent 4)
    'Degraded-only': '#F07167',    # Coral (Accent 3)
    'No arrays': '#545454',        # Gray (Dark 2)
}

# Alpha values per category — normal alleles are subdued
CONFIG_ALPHA = {
    'A-type (C+)': 0.25,
    'B-type (C-)': 0.25,
    '10qB (C-)': 1.0,        # abnormal on chr10 = full visibility
    'Functional pLAM': 1.0,
    'Multi-array': 1.0,
    'Mixed': 1.0,
    'Degraded-only': 1.0,
    'No arrays': 1.0,
}

# Categories in display order
CATEGORIES_CHR4 = ['A-type (C+)', 'B-type (C-)', 'Multi-array', 'Mixed', 'Degraded-only', 'No arrays']
CATEGORIES_CHR10 = ['A-type (C+)', 'B-type (C-)', '10qB (C-)', 'Functional pLAM', 'Multi-array', 'Mixed', 'Degraded-only', 'No arrays']
# Kept for backward compat
CATEGORIES_BASE = CATEGORIES_CHR4


def categorize_haplotype(row, plam_lookup=None):
    """Categorize a haplotype into structural config type."""
    if row['num_arrays'] == 0:
        return 'No arrays'
    if row['terminal_feature'] == 'degraded (no terminal)':
        return 'Degraded-only'

    # Check for functional pLAM on chr10 (takes priority)
    if row['chromosome'] == 'chr10' and plam_lookup is not None:
        hap = row['haplotype']
        if hap in plam_lookup:
            return 'Functional pLAM'

    # C- on chr10 (10qB) — abnormal, gets its own category
    if row['chromosome'] == 'chr10' and row['terminal_feature'] == 'No β-satellite terminal':
        return '10qB (C-)'

    if row['num_arrays'] > 1:
        if isinstance(row['arrays_info'], str):
            try:
                arrays = json.loads(row['arrays_info'])
                types = set(a.get('type') for a in arrays)
                if 'canonical' in types and 'degraded' in types:
                    return 'Mixed'
            except:
                pass
        return 'Multi-array'

    # Single-array canonical: split into A-type and B-type
    if row['terminal_feature'] == 'β-satellite terminal':
        return 'A-type (C+)'
    else:
        return 'B-type (C-)'


def build_plam_lookup():
    """Build lookup of haplotypes with functional pLAM on chr10."""
    plam_file = BASE_DIR / 'results' / 'd4z4_array_coordinates.plam.tsv'
    if not plam_file.exists():
        return set()
    plam_df = pd.read_csv(plam_file, sep='\t')
    chr10_func = plam_df[
        (plam_df['chromosome'] == 'chr10') &
        (plam_df['plam_functional'] == 'func_Plam')
    ]
    # Map back to haplotype names
    analysis_df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    func_haps = set()
    for _, row in chr10_func.iterrows():
        match = analysis_df[
            (analysis_df['sample'] == row['sample']) &
            (analysis_df['chromosome'] == 'chr10') &
            (analysis_df['hap_num'] == row['hap'])
        ]
        if len(match) > 0:
            func_haps.add(match.iloc[0]['haplotype'])
    return func_haps


def generate_panel_a():
    """Locus schematic with real track example."""
    output_path = OUTPUT_DIR / 'fig5_panel_a_schematic.svg'
    # Generate a single-haplotype visualization as basis
    cmd = [
        sys.executable, str(BASE_DIR / 'plot_d4z4_arrays.py'),
        '--haplotypes', 'chr4_hap1_HG01358#1#CM089053.1',
        '--chromosome', 'chr4',
        '--output', str(output_path),
        '--background', 'white',
        '--width', '1200',
        '--adaptive-window',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
    else:
        print(f"  Saved: {output_path}")
    print("  Note: Add locus schematic above the track in Illustrator/Inkscape")


def generate_panel_b():
    """Repeat count distribution + structural config frequencies."""
    df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')

    from scipy.stats import gaussian_kde

    fig, ax1 = plt.subplots(1, 1, figsize=(8, 5))

    # --- Ridge plot of repeat counts (4 ridges) ---
    ridges = [
        ('chr4', 'C+', 'β-satellite terminal', C_CPLUS, '#2E8B57'),
        ('chr4', 'C-', 'No β-satellite terminal', C_CMINUS, '#4682B4'),
        ('chr10', 'C+', 'β-satellite terminal', C_EMERALD, '#0E7C61'),
        ('chr10', 'C-', 'No β-satellite terminal', C_ROYAL, '#2B5EA7'),
    ]

    x_range = np.linspace(0, 145, 300)
    ridge_height = 0.8
    overlap = 0.3

    for i, (chrom, label, terminal, fill_color, line_color) in enumerate(ridges):
        chrom_df = df[(df['chromosome'] == chrom) & (df['terminal_feature'] == terminal)]
        counts = chrom_df['repeat_count_main'].dropna()
        counts = counts[counts > 0]

        y_offset = (len(ridges) - 1 - i) * (ridge_height - overlap)

        if len(counts) >= 3:
            kde = gaussian_kde(counts, bw_method=0.25)
            density = kde(x_range)
            density = density / density.max() * ridge_height  # normalize height

            ax1.fill_between(x_range, y_offset, y_offset + density,
                            alpha=0.6, color=fill_color, zorder=len(ridges) - i)
            ax1.plot(x_range, y_offset + density,
                    color=line_color, linewidth=1.2, zorder=len(ridges) - i + 1)

        # Baseline
        ax1.axhline(y=y_offset, color=C_GRAY, linewidth=0.3, alpha=0.3,
                   xmin=0, xmax=1, zorder=0)

        # Label
        n = len(counts)
        med = int(counts.median()) if len(counts) > 0 else 0
        ax1.text(-2, y_offset + ridge_height * 0.3,
                f'{chrom} {label}\n(n={n}, med={med})',
                ha='right', va='center', fontsize=9, fontweight='bold',
                color=line_color)

    ax1.set_xlim(0, 145)
    top = len(ridges) * (ridge_height - overlap) + ridge_height + 0.1
    ax1.set_ylim(-0.1, top)
    ax1.set_xlabel('Repeat units', fontsize=12)
    ax1.set_title('D4Z4 repeat count distribution', fontsize=13, fontweight='bold')
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'fig5_panel_b_ridge'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")

    # --- Structural config frequencies (separate plot) ---
    generate_config_bar(df)

    # --- Patient-level stacked bar ---
    generate_patient_stacked_bar(df)


def generate_config_bar(df):
    """Stacked bar of structural configurations per chromosome (separate plot)."""
    plam_lookup = build_plam_lookup()
    df = df.copy()
    df['config'] = df.apply(lambda row: categorize_haplotype(row, plam_lookup), axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    for i, chrom in enumerate(['chr4', 'chr10']):
        chrom_df = df[df['chromosome'] == chrom]
        total = len(chrom_df)
        categories = CATEGORIES_CHR10 if chrom == 'chr10' else CATEGORIES_BASE
        bottom = 0
        for cat in categories:
            count = len(chrom_df[chrom_df['config'] == cat])
            if count == 0:
                continue
            pct = count / total * 100
            ax.bar(i, pct, bottom=bottom, color=CONFIG_COLORS[cat],
                   edgecolor='white', linewidth=0.5)
            if pct > 3:
                ax.text(i, bottom + pct/2, f'{pct:.1f}%',
                       ha='center', va='center', fontsize=8,
                       color='black' if cat in ('A-type (C+)', 'B-type (C-)') else 'white',
                       fontweight='bold')
            bottom += pct

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['chr4\n(n=249)', 'chr10\n(n=333)'], fontsize=11)
    ax.set_ylabel('Percentage of haplotypes', fontsize=12)
    ax.set_title('Structural configurations', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    all_cats = list(dict.fromkeys(CATEGORIES_BASE + CATEGORIES_CHR10))
    config_patches = [mpatches.Patch(color=CONFIG_COLORS[c], label=c)
                      for c in all_cats if c not in ('A-type (C+)', 'B-type (C-)')]
    ax.legend(handles=config_patches, loc='center right', fontsize=7, framealpha=0.9)

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'fig5_panel_b_config_bar'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def generate_patient_stacked_bar(df):
    """Stacked bar at patient level: % of patients with each config type, separated by chromosome."""
    plam_lookup = build_plam_lookup()
    df = df.copy()
    df['config'] = df.apply(lambda row: categorize_haplotype(row, plam_lookup), axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(8, 5), sharey=True)

    for ax_idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[ax_idx]
        chrom_df = df[df['chromosome'] == chrom]
        # Get unique patients for this chromosome
        patients = chrom_df['sample'].unique()
        n_patients = len(patients)

        cats = CATEGORIES_CHR10 if chrom == 'chr10' else CATEGORIES_BASE
        non_standard = [c for c in cats if c not in ('A-type (C+)', 'B-type (C-)')]

        # For each category, count how many patients have at least one haplotype of that type
        cat_pcts = []
        cat_labels = []
        cat_colors = []
        for cat in non_standard:
            patients_with = chrom_df[chrom_df['config'] == cat]['sample'].nunique()
            if patients_with > 0:
                cat_pcts.append(patients_with / n_patients * 100)
                cat_labels.append(cat)
                cat_colors.append(CONFIG_COLORS[cat])

        if not cat_pcts:
            ax.set_title(f'{chrom} (n={n_patients})', fontsize=13, fontweight='bold')
            ax.axis('off')
            continue

        # Stacked horizontal bar (single bar, stacked by category)
        bottom = 0
        for pct, label, color in zip(cat_pcts, cat_labels, cat_colors):
            bar = ax.barh(0, pct, left=bottom, color=color,
                         edgecolor='white', linewidth=0.5, height=0.5)
            if pct > 2:
                ax.text(bottom + pct / 2, 0, f'{pct:.1f}%',
                       ha='center', va='center', fontsize=8,
                       color='white', fontweight='bold')
            bottom += pct

        ax.set_xlim(0, max(bottom * 1.15, 15))
        ax.set_xlabel('% of patients', fontsize=11)
        ax.set_yticks([])
        ax.set_title(f'{chrom} (n={n_patients} patients)', fontsize=13, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # Shared legend
    all_cats = list(dict.fromkeys(CATEGORIES_BASE + CATEGORIES_CHR10))
    non_standard_all = [c for c in all_cats if c not in ('A-type (C+)', 'B-type (C-)')]
    # Only include categories that actually appear
    all_configs = df['config'].unique()
    legend_cats = [c for c in non_standard_all if c in all_configs]
    legend_patches = [mpatches.Patch(color=CONFIG_COLORS[c], label=c) for c in legend_cats]
    fig.legend(handles=legend_patches, loc='lower center', ncol=min(len(legend_cats), 4),
              fontsize=8, framealpha=0.9, bbox_to_anchor=(0.5, -0.08))

    plt.suptitle('Patients with non-standard D4Z4 configurations',
                fontsize=13, fontweight='bold')
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'fig5_panel_b_patient_bar'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def generate_panel_b_waffle():
    """Alternative: waffle chart version of structural config frequencies, separated by chromosome."""
    df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')

    plam_lookup = build_plam_lookup()
    df['config'] = df.apply(lambda row: categorize_haplotype(row, plam_lookup), axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(10, 6))

    for ax_idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[ax_idx]
        chrom_df = df[df['chromosome'] == chrom]
        total = len(chrom_df)

        # Waffle: each square = 1 haplotype
        # Determine grid dimensions
        ncols = int(np.ceil(np.sqrt(total * 1.3)))  # slightly wider than tall
        nrows = int(np.ceil(total / ncols))

        # Build flat list of (color, alpha) in category order
        cats = CATEGORIES_CHR10 if chrom == 'chr10' else CATEGORIES_CHR4
        cell_data = []  # list of (color, alpha)
        for cat in cats:
            count = len(chrom_df[chrom_df['config'] == cat])
            alpha = CONFIG_ALPHA.get(cat, 1.0)
            cell_data.extend([(CONFIG_COLORS[cat], alpha)] * count)

        # Pad to fill grid
        while len(cell_data) < nrows * ncols:
            cell_data.append(('white', 0.0))

        # Draw waffle
        for i, (color, alpha) in enumerate(cell_data):
            row = i // ncols
            col = i % ncols
            rect = plt.Rectangle((col, nrows - 1 - row), 0.9, 0.9,
                                  facecolor=color, alpha=alpha,
                                  edgecolor='#E5E7EB', linewidth=0.3)
            ax.add_patch(rect)

        ax.set_xlim(-0.5, ncols + 0.5)
        ax.set_ylim(-0.5, nrows + 0.5)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'{chrom} (n={total})', fontsize=13, fontweight='bold', pad=10)

        # Annotation: counts for non-normal categories
        y_text = -1.5
        normal_cats = {'A-type (C+)', 'B-type (C-)'}
        for cat in cats:
            if cat in normal_cats:
                continue
            count = len(chrom_df[chrom_df['config'] == cat])
            if count > 0:
                ax.text(ncols / 2, y_text,
                       f'{cat}: {count} ({count/total*100:.1f}%)',
                       ha='center', fontsize=9, color=C_GRAY,
                       fontweight='bold')
                y_text -= 1.2

    # Shared legend — all categories, with alpha applied
    all_cats = list(dict.fromkeys(CATEGORIES_CHR4 + CATEGORIES_CHR10))
    legend_patches = []
    for c in all_cats:
        import matplotlib.colors as mcolors
        rgba = list(mcolors.to_rgba(CONFIG_COLORS[c]))
        rgba[3] = CONFIG_ALPHA.get(c, 1.0)
        legend_patches.append(mpatches.Patch(color=rgba, label=c))
    fig.legend(handles=legend_patches, loc='lower center', ncol=min(len(all_cats), 4),
              fontsize=7, framealpha=0.9, bbox_to_anchor=(0.5, -0.05))

    plt.suptitle('D4Z4 structural configurations', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = OUTPUT_DIR / 'fig5_panel_b_waffle'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def generate_panel_c():
    """Multi-array examples: HG01109 chr4 + HG02717 chr10."""
    output_path = OUTPUT_DIR / 'fig5_panel_c_multi_array.svg'
    cmd = [
        sys.executable, str(BASE_DIR / 'plot_d4z4_arrays.py'),
        '--haplotypes',
        'chr4_hap1_HG01109#1#CM092086.1',
        'chr10_hap2_HG02717#2#CM088975.1',
        '--output', str(output_path),
        '--background', 'white',
        '--width', '1200',
        '--adaptive-window',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
    else:
        print(f"  Saved: {output_path}")


def generate_panel_d():
    """Degraded + missing: HG01928 chr10 + HG00290 chr4."""
    output_path = OUTPUT_DIR / 'fig5_panel_d_degraded_missing.svg'
    cmd = [
        sys.executable, str(BASE_DIR / 'plot_d4z4_arrays.py'),
        '--haplotypes',
        'chr10_hap2_HG01928#2#CM089008.1',
        'chr4_hap1_HG00290#1#CM089980.1',
        '--output', str(output_path),
        '--background', 'white',
        '--width', '1200',
        '--adaptive-window',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
    else:
        print(f"  Saved: {output_path}")


def generate_panel_e():
    """Haplotype asymmetry scatter (1x1 only) + rare chr10 configurations."""
    asym_df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_variability' / 'haplotype_asymmetry.tsv', sep='\t')

    # Filter to 1x1 samples only (both haplotypes have exactly 1 array)
    analysis_df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    single_array_samples = set()
    for sample in analysis_df['sample'].unique():
        for chrom in ['chr4', 'chr10']:
            sdf = analysis_df[(analysis_df['sample'] == sample) & (analysis_df['chromosome'] == chrom)]
            hap1 = sdf[sdf['hap_num'] == 'hap1']
            hap2 = sdf[sdf['hap_num'] == 'hap2']
            if len(hap1) == 1 and len(hap2) == 1:
                if hap1.iloc[0]['num_arrays'] == 1 and hap2.iloc[0]['num_arrays'] == 1:
                    single_array_samples.add((sample, chrom))

    asym_df = asym_df[asym_df.apply(lambda r: (r['sample'], r['chromosome']) in single_array_samples, axis=1)]

    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6))

    # Asymmetry scatter (1x1 only)
    for chrom, color, marker in [('chr4', C_BLUE, 'o'), ('chr10', C_EMERALD, 's')]:
        subset = asym_df[asym_df['chromosome'] == chrom]
        ax1.scatter(subset['hap1_total_ru'], subset['hap2_total_ru'],
                   s=25, alpha=0.5, color=color, marker=marker,
                   edgecolor='white', linewidth=0.3, label=chrom, zorder=3)

    # Identity line
    max_val = max(asym_df['hap1_total_ru'].max(), asym_df['hap2_total_ru'].max()) + 5
    ax1.plot([0, max_val], [0, max_val], 'k-', alpha=0.2, linewidth=1)

    # Fold-difference lines
    x = np.linspace(1, max_val, 100)
    for fold in [3, 5]:
        ax1.plot(x, x * fold, '--', color=C_GRAY, alpha=0.3, linewidth=0.8)
        ax1.plot(x * fold, x, '--', color=C_GRAY, alpha=0.3, linewidth=0.8)
        ax1.text(max_val * 0.7, max_val * 0.7 / fold, f'{fold}×',
                fontsize=8, color=C_GRAY, alpha=0.6, rotation=0)

    # Label extreme asymmetry cases
    extremes = asym_df.nlargest(3, 'fold_difference')
    for _, row in extremes.iterrows():
        ax1.annotate(f"{row['sample']}\n({row['fold_difference']:.1f}×)",
                    xy=(row['hap1_total_ru'], row['hap2_total_ru']),
                    fontsize=7, ha='left', va='bottom',
                    xytext=(5, 5), textcoords='offset points',
                    color=C_GRAY)

    ax1.set_xlabel('Haplotype 1 repeat units', fontsize=12)
    ax1.set_ylabel('Haplotype 2 repeat units', fontsize=12)
    n_samples = len(asym_df)
    ax1.set_title(f'Haplotype asymmetry (single-array, n={n_samples})', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_xlim(0, max_val)
    ax1.set_ylim(0, max_val)
    ax1.set_aspect('equal')

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'fig5_panel_e_asymmetry'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def generate_panel_g():
    """Haplotype config pairing heatmap (triangle) — separate for 4q and 10q.

    Since hap1/hap2 are randomly assigned, we collapse unordered pairs
    into a lower triangle: count(A,B) = count where one hap is A and other is B.
    Overwrites the original haplotype_config_crosstab output.
    """
    from matplotlib.colors import LinearSegmentedColormap

    df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')
    plam_lookup = build_plam_lookup()
    df['config'] = df.apply(lambda row: categorize_haplotype(row, plam_lookup), axis=1)

    # Blue colormap matching the original crosstab
    cmap_blues = LinearSegmentedColormap.from_list('custom_blues',
        ['#F0F4FA', '#A8C4E0', '#4A7FB5', '#1A3A5C'])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    for ax_idx, (chrom, cats) in enumerate([
        ('chr4', CATEGORIES_CHR4),
        ('chr10', CATEGORIES_CHR10),
    ]):
        ax = axes[ax_idx]
        chrom_df = df[df['chromosome'] == chrom]

        # Filter to categories that actually appear; exclude 'No arrays' on chr4
        exclude = {'No arrays'} if chrom == 'chr4' else set()
        present_cats = [c for c in cats if c in chrom_df['config'].values and c not in exclude]

        # Build per-sample pairs
        pair_counts = {}
        total = 0
        for sample in chrom_df['sample'].unique():
            sdf = chrom_df[chrom_df['sample'] == sample]
            hap1 = sdf[sdf['hap_num'] == 'hap1']
            hap2 = sdf[sdf['hap_num'] == 'hap2']
            if len(hap1) != 1 or len(hap2) != 1:
                continue
            c1 = hap1.iloc[0]['config']
            c2 = hap2.iloc[0]['config']
            a, b = sorted([c1, c2], key=lambda x: present_cats.index(x) if x in present_cats else 999)
            pair_counts[(a, b)] = pair_counts.get((a, b), 0) + 1
            total += 1

        # Build lower-triangle matrix
        n = len(present_cats)
        matrix = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(i + 1):
                matrix[i, j] = pair_counts.get((present_cats[j], present_cats[i]), 0)

        # Plot
        ax.imshow(matrix, cmap=cmap_blues, aspect='equal', origin='upper')
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(present_cats, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(present_cats, fontsize=8)

        # Annotate cells
        max_val = np.nanmax(matrix) if np.any(~np.isnan(matrix)) else 1
        for i in range(n):
            for j in range(i + 1):
                val = int(matrix[i, j])
                if val > 0:
                    color = 'white' if val > max_val * 0.5 else 'black'
                    ax.text(j, i, str(val), ha='center', va='center',
                            fontsize=9, fontweight='bold', color=color)

        # Mask upper triangle
        for i in range(n):
            for j in range(i + 1, n):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor='white', edgecolor='white', zorder=2))

        ax.set_title(f'{chrom} (n={total} patients)',
                     fontsize=12, fontweight='bold')

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'haplotype_config_crosstab'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def generate_panel_h():
    """Repeat count comparison: shorter vs longer haplotype (sorted pair).

    Since hap1/hap2 are random, we sort each pair so the shorter array
    is always on x and longer on y. 4q and 10q shown with different colors.
    Overwrites the original fig5_panel_e_asymmetry output.
    Labels the top 3 most asymmetric samples as in the original panel E.
    """
    df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')

    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    all_pairs = []  # collect for annotation

    for chrom, color, marker, label in [
        ('chr4', C_BLUE, 'o', 'chr4'),
        ('chr10', C_EMERALD, 's', 'chr10'),
    ]:
        chrom_df = df[df['chromosome'] == chrom]
        pairs = []
        for sample in chrom_df['sample'].unique():
            sdf = chrom_df[chrom_df['sample'] == sample]
            hap1 = sdf[sdf['hap_num'] == 'hap1']
            hap2 = sdf[sdf['hap_num'] == 'hap2']
            if len(hap1) != 1 or len(hap2) != 1:
                continue
            if hap1.iloc[0]['num_arrays'] != 1 or hap2.iloc[0]['num_arrays'] != 1:
                continue
            r1 = hap1.iloc[0]['repeat_count_main']
            r2 = hap2.iloc[0]['repeat_count_main']
            if pd.isna(r1) or pd.isna(r2) or r1 == 0 or r2 == 0:
                continue
            shorter, longer = sorted([r1, r2])
            fold = longer / shorter if shorter > 0 else float('inf')
            pairs.append((shorter, longer))
            all_pairs.append({
                'sample': sample, 'chromosome': chrom,
                'shorter': shorter, 'longer': longer,
                'fold': fold,
            })

        if pairs:
            xs, ys = zip(*pairs)
            ax.scatter(xs, ys, s=25, alpha=0.5, color=color, marker=marker,
                       edgecolor='white', linewidth=0.3, label=f'{label} (n={len(pairs)})',
                       zorder=3)

    # Identity line
    all_vals = df['repeat_count_main'].dropna()
    max_val = all_vals.max() + 5
    ax.plot([0, max_val], [0, max_val], 'k-', alpha=0.2, linewidth=1)

    # Fold-difference guide lines
    x = np.linspace(1, max_val, 100)
    for fold in [3, 5]:
        ax.plot(x, x * fold, '--', color=C_GRAY, alpha=0.3, linewidth=0.8)
        ax.text(max_val * 0.7, max_val * 0.7 / fold, f'{fold}\u00d7',
                fontsize=8, color=C_GRAY, alpha=0.6)

    # Annotate top 3 most asymmetric samples
    pairs_df = pd.DataFrame(all_pairs)
    if len(pairs_df) > 0:
        extremes = pairs_df.nlargest(3, 'fold')
        for _, row in extremes.iterrows():
            ax.annotate(f"{row['sample']}\n({row['fold']:.1f}\u00d7)",
                        xy=(row['shorter'], row['longer']),
                        fontsize=7, ha='left', va='bottom',
                        xytext=(5, 5), textcoords='offset points',
                        color=C_GRAY)

    n_total = len(all_pairs)
    ax.set_xlabel('Shorter haplotype (repeat units)', fontsize=12)
    ax.set_ylabel('Longer haplotype (repeat units)', fontsize=12)
    ax.set_title(f'Haplotype asymmetry (single-array, n={n_total})',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')

    plt.tight_layout()
    output_path = OUTPUT_DIR / 'fig5_panel_e_asymmetry'
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Figure 5: D4Z4 Structural Variation Resolved by KaryoScope")
    print("=" * 60)

    print("\nPanel (a): Locus schematic with real track")
    generate_panel_a()

    print("\nPanel (b): Repeat count distribution + config frequencies")
    generate_panel_b()

    print("\nPanel (b) alt: Waffle chart version")
    generate_panel_b_waffle()

    print("\nPanel (c): Multi-array examples")
    generate_panel_c()

    print("\nPanel (d): Degraded + missing examples")
    generate_panel_d()

    print("\nPanel (e): Haplotype asymmetry + exchange signatures")
    generate_panel_e()

    print("\nPanel (g): Diplotype config pairing heatmap (triangle)")
    generate_panel_g()

    print("\nPanel (h): Repeat count pairing (shorter vs longer)")
    generate_panel_h()

    print(f"\nAll panels saved to {OUTPUT_DIR}/fig5_panel_*.svg")
    print("Assemble panels a-h into final figure in Illustrator/Inkscape.")


if __name__ == '__main__':
    main()
