# PROVENANCE: these 4 comment lines were added in this snapshot and are not part of the original.
# KaryoScope-archive workflow/scripts/get_regions.py at commit 0bf3a7d ("master pangenome QC table").
# Original sha256 ddaba0d356f0449ec4bf1d90ec9ac699b35e9c542c560be0231add6cd47a37c7; drop the first 4 lines to get it.
# Snakemake script of rule get_regions; its centromere input (1 Mb annotation) no longer exists.
import sys
import os

wanted_chromosome = snakemake.wildcards.chrom
wanted_haplotype = snakemake.wildcards.hap
wanted_sample = snakemake.wildcards.sample

is_t2t = False
contig_for_fai = None
full_contig_id_for_cen = None # The centromere file parsing is different

# --- Step 1: Check for T2T status ---
try:
    with open(snakemake.input[0], "r") as input_file:
      for line in input_file:
        full_contig_id = line.strip().split()[0]
        try:
            # Parse components based on original script and your logic
            chromosome, haplotype, sample = full_contig_id.split("#")[0].split("_")
            simple_name = full_contig_id.split("#")[2]
        except IndexError:
            continue # Skip malformed line in stats file

        if (chromosome, haplotype, sample) == (wanted_chromosome, wanted_haplotype, wanted_sample):
          # If we find the combo, check if it's a "CM" contig
          if simple_name.startswith("CM"):
            is_t2t = True
            # The original script uses a strange parsing for the contig name
            #
            # We will replicate it.
            contig_for_fai = full_contig_id.split("_")[2]
            full_contig_id_for_cen = full_contig_id
          break # Found our S/H/C combo, no need to search more
except FileNotFoundError:
    print(f"Warning: Stats file not found: {snakemake.input[0]}", file=sys.stderr)
    # is_t2t remains False


# --- Step 2: Write files based on T2T status ---

# If not a T2T combo, write empty files and exit successfully
if not is_t2t:
  try:
    with open(snakemake.output.subtelo_p, "w") as f_p, \
         open(snakemake.output.subtelo_q, "w") as f_q, \
         open(snakemake.output.centromere, "w") as f_c:
      pass # Write empty files
  except AttributeError:
      # Fallback if output wildcards aren't named in the rule
      with open(snakemake.output[0], "w") as f_p, \
           open(snakemake.output[1], "w") as f_q, \
           open(snakemake.output[2], "w") as f_c:
        pass # Write empty files
  sys.exit(0) # Exit successfully

# --- If we get here, it IS a T2T combo ---

# --- Step 3: Process FAI file (p- and q-subtelomeres) ---
subtelo_boundary = 500000 #
found_in_fai = False
try:
    with open(snakemake.input[1], "r") as input_file, \
         open(snakemake.output[0], "w") as output_file_p, \
         open(snakemake.output[1], "w") as output_file_q:
      for line in input_file:
        contig_from_fai, length = line.strip().split()[:2]
        length = int(length)
        
        # Using the original script's logic
        if contig_from_fai == contig_for_fai:
          output_file_p.write(f"{contig_from_fai}\t0\t{subtelo_boundary}\n")
          output_file_q.write(f"{contig_from_fai}\t{length-subtelo_boundary}\t{length}\n")
          found_in_fai = True
          break
except FileNotFoundError:
    print(f"Warning: FAI file not found: {snakemake.input[1]}", file=sys.stderr)
    # Let empty files be created

if not found_in_fai:
    print(f"Warning: T2T contig {contig_for_fai} not found in FAI {snakemake.input[1]}", file=sys.stderr)
    # Create empty files if they weren't already
    if not os.path.exists(snakemake.output[0]):
      with open(snakemake.output[0], "w"): pass
    if not os.path.exists(snakemake.output[1]):
      with open(snakemake.output[1], "w"): pass

# --- Step 4: Process centromeres bed file ---
is_first_line = True
found_in_cen = False
start, stop = "0", "0" # initialize
try:
    with open(snakemake.input[2], "r") as input_file, \
         open(snakemake.output[2], "w") as output_file_cen:
      
      # NOTE: This logic from the original script is complex and assumes
      # the centromere file is sorted and contiguous.
      for line in input_file:
        contig_from_cen = line.strip().split()[0] # Full ID from cen file
        
        # Check if we've passed the block for our contig
        if not is_first_line and contig_from_cen != full_contig_id_for_cen:
          output_file_cen.write(f"{contig_for_fai}\t{start}\t{stop}\n")
          found_in_cen = True # We found it and wrote it
          break
          
        if contig_from_cen == full_contig_id_for_cen:
          if is_first_line:
            start, stop = line.strip().split()[1:3]
            is_first_line = False
          else:
            stop = line.strip().split()[2] # Keep extending the end
      
      # Handle case where it's the last contig in the file
      if not is_first_line and not found_in_cen:
          output_file_cen.write(f"{contig_for_fai}\t{start}\t{stop}\n")
          found_in_cen = True
      
      if not found_in_cen:
          # Contig wasn't in the centromere file at all, write empty file
          pass # File is already open for writing, just close it empty

except FileNotFoundError:
    print(f"Warning: Centromere file not found: {snakemake.input[2]}", file=sys.stderr)
    with open(snakemake.output[2], "w"): pass # Create empty file
