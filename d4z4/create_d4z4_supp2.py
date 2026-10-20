#!/usr/bin/env python3
"""
Supplementary Figure 2: Additional examples of non-standard D4Z4 configurations.

Generates individual panel SVGs for a 5-row × 2-column grid:
  Row 1: Multi-array chr4 (HG03942, HG06807)
  Row 2: Multi-array chr10 (HG06807, HG02818)
  Row 3: Mixed canonical + degraded chr4 (NA18570, HG00133)
  Row 4: Degraded-only chr10 (HG00438, HG01361)
  Row 5: Functional pLAM on chr10 (HG02572, HG02976)

Output: results/karyoscope_companion/supp2_row{N}_{category}.svg
"""

import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).parent
SCRIPT = BASE_DIR / 'plot_d4z4_arrays.py'
OUTPUT_DIR = BASE_DIR / 'results' / 'karyoscope_companion'

# Each row: (filename_suffix, chromosome, [haplotype_ids])
ROWS = [
    ("row1_multi_array_chr4", "chr4", [
        "chr4_hap1_HG03942#1#CM088391.1",
        "chr4_hap2_HG06807#2#CM101432.1",
    ]),
    ("row2_multi_array_chr10", "chr10", [
        "chr10_hap2_HG06807#2#CM101438.1",
        "chr10_hap2_HG02818#2#CM091991.1",
    ]),
    ("row3_mixed_chr4", "chr4", [
        "chr4_hap1_NA18570#1#CM087657.1",
        "chr4_hap2_HG00133#2#CM090063.1",
    ]),
    ("row4_degraded_only_chr10", "chr10", [
        "chr10_hap1_HG00438#1#CM089176.1",
        "chr10_hap1_HG01361#1#CM089027.1",
    ]),
    ("row5_functional_plam_chr10", "chr10", [
        "chr10_hap1_HG02572#1#CM088159.1",
        "chr10_hap2_HG02976#2#CM087126.1",
    ]),
]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for suffix, chrom, haplotypes in ROWS:
        output_path = OUTPUT_DIR / f'supp2_{suffix}.svg'
        print(f"\nGenerating {suffix}...")
        print(f"  Haplotypes: {', '.join(h.split('_', 2)[-1] for h in haplotypes)}")

        cmd = [
            sys.executable, str(SCRIPT),
            '--haplotypes', *haplotypes,
            '--chromosome', chrom,
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

    print(f"\nAll panels saved to {OUTPUT_DIR}/supp2_row*.svg")
    print("Assemble rows into final supplementary figure in Illustrator/Inkscape.")


if __name__ == '__main__':
    main()
