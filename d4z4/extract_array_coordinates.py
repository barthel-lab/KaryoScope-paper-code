#!/usr/bin/env python3
"""
Extract D4Z4 array coordinates from analysis TSV.

Outputs a TSV with array coordinates, classification, and β-satellite terminal positions.
Array index is numbered from telomere (rightmost = #1) to match plot convention.
"""

import json
import pandas as pd


def extract_coordinates(input_file, output_file):
    """Extract array coordinates from analysis TSV."""
    df = pd.read_csv(input_file, sep='\t')

    rows = []
    for _, row in df.iterrows():
        chromosome = row['chromosome']
        sample = row['sample']
        hap = row['hap_num']

        # Extract contig from haplotype name (e.g., "chr4_hap1_HG02129#1#CM085848.1" -> "CM085848.1")
        haplotype = row['haplotype']
        contig = haplotype.split('#')[-1] if '#' in haplotype else haplotype

        # Parse arrays_info JSON
        arrays_info_str = row.get('arrays_info')
        if pd.isna(arrays_info_str) or not arrays_info_str:
            continue

        try:
            arrays_info = json.loads(arrays_info_str)
        except (json.JSONDecodeError, TypeError):
            continue

        if not arrays_info:
            continue

        # Sort arrays by position and re-number from telomere (rightmost = 1)
        sorted_arrays = sorted(arrays_info, key=lambda a: a.get('start', 0))
        total_arrays = len(sorted_arrays)

        for i, arr in enumerate(sorted_arrays):
            # Number from telomere: rightmost = #1
            array_index = total_arrays - i

            # Classify array
            arr_type = arr.get('type', 'canonical')
            has_plam = arr.get('has_plam', False)
            has_ct = arr.get('has_ct_terminal', False)

            if arr_type == 'degraded':
                classification = 'D'
            elif has_plam:
                classification = 'C+'
            elif has_ct:
                classification = 'C-'
            else:
                classification = 'C-'  # Canonical without terminal

            # Array coordinates
            array_start = int(arr.get('start', 0))
            array_end = int(arr.get('end', 0))
            element_count = arr.get('element_count', 0)

            # β-satellite terminal coordinates (when present)
            bsat_start = ''
            bsat_end = ''
            if has_plam:
                plam_end = arr.get('plam_end')
                if plam_end:
                    bsat_start = array_end  # β-sat starts after array
                    bsat_end = int(plam_end)

            rows.append({
                'chromosome': chromosome,
                'sample': sample,
                'hap': hap,
                'contig': contig,
                'array_index': array_index,
                'array_type': classification,
                'element_count': element_count,
                'array_start': array_start,
                'array_end': array_end,
                'bsat_start': bsat_start,
                'bsat_end': bsat_end
            })

    # Create output dataframe and save
    output_df = pd.DataFrame(rows)
    output_df = output_df.sort_values(['chromosome', 'sample', 'hap', 'array_index'],
                                       ascending=[True, True, True, False])
    output_df.to_csv(output_file, sep='\t', index=False)

    print(f"Extracted {len(output_df)} arrays from {len(df)} haplotypes")
    print(f"Output: {output_file}")

    return output_df


def merge_plam(coords_file, plam_file, output_file, window_size=500000):
    """Merge pLAM annotations into array coordinates.

    The pLAM file uses telomere-0-based coordinates. Convert to the same
    coordinate system as the array coordinates by subtracting from window_size.
    """
    coords = pd.read_csv(coords_file, sep='\t')
    plam = pd.read_csv(plam_file, sep='\t')

    # Convert telomere-based pLAM coordinates to match array coordinate system
    plam_cols = plam[['chromosome', 'sample', 'hap', 'array_index',
                       'plam_start', 'plam_end', 'plam_status']].copy()

    # Convert: 500k - telo_coord, swapping start/end since orientation flips
    has_plam = plam_cols['plam_start'].notna() & (plam_cols['plam_start'] > 0)
    plam_cols.loc[has_plam, ['plam_start', 'plam_end']] = (
        window_size - plam_cols.loc[has_plam, ['plam_end', 'plam_start']].values
    )
    # Zero out rows without pLAM data
    plam_cols.loc[~has_plam, ['plam_start', 'plam_end']] = pd.NA

    plam_cols = plam_cols.rename(columns={'plam_status': 'plam_functional'})

    # Left join
    merged = coords.merge(plam_cols, on=['chromosome', 'sample', 'hap', 'array_index'], how='left')

    # Convert plam_start/plam_end to Int64 (nullable integer) for clean output
    merged['plam_start'] = pd.array(merged['plam_start'], dtype='Int64')
    merged['plam_end'] = pd.array(merged['plam_end'], dtype='Int64')

    merged.to_csv(output_file, sep='\t', index=False)

    n_plam = merged['plam_start'].notna().sum()
    print(f"Merged pLAM annotations: {n_plam}/{len(merged)} rows with pLAM data")
    print(f"Functional: {(merged['plam_functional'] == 'func_Plam').sum()}, "
          f"Non-functional: {(merged['plam_functional'] == 'nonfunc_Plam').sum()}")
    print(f"Output: {output_file}")

    return merged


if __name__ == '__main__':
    input_file = 'results/d4z4_array_analysis.tsv'
    output_file = 'results/d4z4_array_coordinates.tsv'
    extract_coordinates(input_file, output_file)

    # Merge pLAM annotations
    plam_file = 'data/plam/d4z4_array_coordinates.plam.array_coords.tsv'
    plam_output = 'results/d4z4_array_coordinates.plam.tsv'
    merge_plam(output_file, plam_file, plam_output)
