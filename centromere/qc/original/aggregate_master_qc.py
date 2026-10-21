# PROVENANCE: these 4 comment lines were added in this snapshot and are not part of the original.
# KaryoScope-archive workflow/scripts/aggregate_master_qc.py: UNTRACKED (in no commit), from the working tree, dated 2025-11-19.
# Original sha256 6de17992c64103e2a4139b5d3d78b0d716a7d04c8afd212843ed99e67c5330f8; drop the first 4 lines to get it.
# Snakemake script of rule aggregate_master_qc_table; build_master_qc_table.py is the runnable version.
import sys
import os
from collections import defaultdict

# --- Globals from Snakemake ---
SAMPLES = snakemake.params.samples
HAPS = snakemake.params.haps
CHROMS = snakemake.params.chroms
TOOLS = snakemake.params.tools
# This list matches the {zoom} wildcard from pangenome.smk
REGIONS = snakemake.params.regions 

# --- Step 1: Build the T2T Map ---
# This map stores (sample, hap, chrom) tuples that have a T2T contig.
print("Building T2T map from scaffold stats files...")
t2t_map = set()
sample_hap_chrom_to_contig = {}
for stats_file in snakemake.input.scaffold_stats:
    try:
        with open(stats_file, "r") as f:
            for line in f:
                try:
                    # Logic from get_regions.py and your clarification
                    contig_id = line.strip().split()[0]
                    contig = contig_id.split("_")[2]
                    chromosome, haplotype, sample = contig_id.split("#")[0].split("_")
                    simple_name = contig_id.split("#")[2]
                    
                    # This is your new combined T2T logic
                    if simple_name.startswith("CM"):
                        t2t_map.add((sample, haplotype, chromosome))
                        sample_hap_chrom_to_contig[(sample, haplotype, chromosome)] = contig
                except Exception:
                    # Skip malformed lines in stats file
                    continue
    except FileNotFoundError:
        # This sample might not have a stats file
        print(f"Warning: Could not find stats file {stats_file}")
        continue

print(f"Found {len(t2t_map)} T2T (sample, hap, chrom) combinations.")

# --- Step 2: Build Intersection File Map ---
# This maps (sample, hap, chrom, region, tool) -> "path/to/intersection.bed"
print("Mapping intersection files...")
intersection_file_map = {}
for f in snakemake.input.intersections:
    # Re-build the key from the file path
    # Path: results/pangenome/errors/{sample}/{hap}/{chrom}/{sample}.{hap}.{chrom}.{zoom}.{qc}.errors.bed
    try:
        parts = f.split("/")
        basename = parts[-1]
        chrom = parts[-2]
        hap = parts[-3]
        sample = parts[-4]
        
        basename_parts = basename.split(".")
        qc = basename_parts[-3]
        region = basename_parts[-4] # This is the {zoom} wildcard

        key = (sample, hap, chrom, region, qc)
        intersection_file_map[key] = f
    except Exception as e:
        print(f"Warning: Could not parse file path {f}: {e}")
        
# --- Step 3: Iterate and Write Master Table ---
print("Generating master table...")
with open(snakemake.output[0], "w") as out_f:
    out_f.write("sample\thaplotype\tchromosome\tregion\ttool\tcontig\tstart\tstop\tflag\n")
    
    for sample in SAMPLES:
        for hap in HAPS:
            for chrom in CHROMS:
                # Check if this (s, h, c) combo has a T2T contig
                is_t2t = (sample, hap, chrom) in t2t_map
                
                for region in REGIONS:
                    for tool in TOOLS:
                        key = (sample, hap, chrom, region, tool)
                        
                        if not is_t2t:
                            # Not a T2T contig, write N/A
                            out_f.write(f"{sample}\t{hap}\t{chrom}\t{region}\t{tool}\tN/A\tN/A\tN/A\tN/A\n")
                        else:
                            # Is a T2T contig, check for errors
                            intersection_file = intersection_file_map.get(key)
                            
                            if not intersection_file:
                                print(f"Error: Missing intersection file for {key}. Writing '.'")
                                out_f.write(f"{sample}\t{hap}\t{chrom}\t{region}\t{tool}\t.\t.\t.\t.\n")
                                continue
                            
                            found_errors = []
                            try:
                                with open(intersection_file, "r") as in_f:
                                    for line in in_f:
                                        if line.strip():
                                            found_errors.append(line.strip())
                                            
                            except FileNotFoundError:
                                # This can happen if a job failed.
                                # We treat it as "empty" and write "."
                                print(f"Warning: File not found {intersection_file}. Writing '.'")
                                # No need to do anything, found_errors is still []

                            if not found_errors:
                                # File is empty, write "."
                                contig = sample_hap_chrom_to_contig[(sample, hap, chrom)]
                                out_f.write(f"{sample}\t{hap}\t{chrom}\t{region}\t{tool}\t{contig}\t.\t.\t.\n")
                            else:
                                # File has errors, write all of them
                                for error_line in found_errors:
                                    # Assuming intersection file is: contig \t start \t stop \t flag ...
                                    parts = error_line.split("\t")
                                    # Take first 4 columns to be safe
                                    contig, start, stop, flag = parts[0:4]
                                    out_f.write(f"{sample}\t{hap}\t{chrom}\t{region}\t{tool}\t{contig}\t{start}\t{stop}\t{flag}\n")
                                    
print("Master QC table generation complete.")
