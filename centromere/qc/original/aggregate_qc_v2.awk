# PROVENANCE: these 4 comment lines were added in this snapshot and are not part of the original.
# aggregate_qc_v2.awk: UNTRACKED file from a KaryoScope-paper working copy (in no commit). Run as
# gawk -f aggregate_qc_v2.awk master_qc_table.full.tsv > aggregate_qc_v4.tsv. Original sha256
# 266dc64890868605a3abdc5c7fd52438188bce7f9b371ba7e94e1b999b9a4914; drop the first 4 lines to get it.
#!/usr/bin/awk -f

# BEGIN block: Runs once before processing any lines
BEGIN {
    # Set the Input (FS) and Output (OFS) Field Separators to a tab character
    FS = OFS = "\t"

    # Define the full list of error flags you want to count, in the order
    # you want them to appear in the output.
    all_flags = "Col Dup Err NNN COLLAPSE COLLAPSE_OTHER COLLAPSE_VAR MISJOIN HET"
    
    # Split this string into an array called 'flag_list'
    n_flags = split(all_flags, flag_list, " ")

    # --- Print the Header Row ---
    # NOTE: The 'tool' column ($5) has been removed from the key columns.
    printf "sample\thaplotype\tchromosome\tregion\tcontig"
    for (i = 1; i <= n_flags; i++) {
        printf "%s_n_%s", OFS, flag_list[i]
    }
    printf "\n"
    
    # Initialize previous_key to an empty string to track the "group"
    previous_key = ""
}

# Main processing block: Runs for every line *except* the first (NR > 1)
# This neatly skips the header row of the input file.
NR > 1 {
    # Get the flag value from the 9th column
    flag = $9

    # Check for "N/A" rows, which should be skipped entirely
    if (flag == "N/A") {
        next # Skip to the next line of input
    }

    # Create the unique aggregation key from the first 4 columns ($5, tool, is ignored) and $6 (contig)
    key = $1 FS $2 FS $3 FS $4 FS $6

    # --- Control Break Logic ---
    
    # If this is the very first data line, just set the previous_key
    # and wait for the next line.
    if (previous_key == "") {
        previous_key = key
    }

    # If the key has changed, it's time to print the results
    # for the *previous* key group we were tracking.
    if (key != previous_key) {
        # Print the key itself
        printf "%s", previous_key

        # Now, loop through our ordered list of flags
        for (i = 1; i <= n_flags; i++) {
            current_flag = flag_list[i]

            # Get the count for this key and this flag.
            # We add 0 to force a numeric value (0 instead of empty string)
            count = counts[previous_key][current_flag] + 0
            
            # Print the output separator (a tab) and the count
            printf "%s%d", OFS, count
        }
        printf "\n" # Finish the row
        
        # Free the memory for the key we just printed
        delete counts[previous_key]
        
        # Update the previous_key to the new key
        previous_key = key
    }

    # --- Aggregation Logic ---
    # Always increment the count for the *current* key (unless flag is ".")
    # This will now correctly count "HET" and other flags.
    if (flag != ".") {
        counts[key][flag]++
    }
}

# END block: Runs once after all lines have been processed
END {
    # We must print the very last group, as the "key != previous_key"
    # check won't have run for the final set of records.
    # We must check if previous_key is not empty (in case of empty/filtered input)
    if (previous_key != "") {
        printf "%s", previous_key
        for (i = 1; i <= n_flags; i++) {
            current_flag = flag_list[i]
            count = counts[previous_key][current_flag] + 0
            printf "%s%d", OFS, count
        }
        printf "\n"
        
        # Clean up the last entry
        delete counts[previous_key]
    }
}