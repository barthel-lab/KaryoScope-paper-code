# VENDORED, NOT OURS. Pinned copy of KaryoScope-archive:workflow/scripts/bin_features.py
# as of commit 18dccc6 (2026-02-11) -- the version that produced the PUBLISHED
# pangenome inputs. Kept verbatim here because the binning must match exactly and
# the archive working tree has since moved on: commit ea5c2fe (2026-03-23) added the
# ">50% novel" rule and f15c29b (2026-04-09) the v2 vocabulary, neither of which was
# in effect for the published D4Z4 / centromere BEDs.
# Retrieve the original with:
#   git -C <KaryoScope-archive> show 18dccc6:workflow/scripts/bin_features.py
import sys
import argparse
import gzip  # Added to support .gz input
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser(
        description="Bin a sorted BED file and assign the feature with max overlap to each bin."
    )
    parser.add_argument(
        "-b", "--binsize", 
        type=int, 
        required=True, 
        help="Size of the bins in base pairs."
    )
    parser.add_argument(
        "-i", "--input", 
        default="-", 
        help="Input BED file. Use '-' for stdin (default). Supports .gz extensions."
    )
    parser.add_argument(
        "-o", "--output", 
        default=None, 
        help="Output BED file. Defaults to stdout."
    )
    parser.add_argument(
        "--specific", 
        action="store_true", 
        help="Prioritize features ending in '_specific'. Falls back to max overlap of any feature if no specific feature exists."
    )
    return parser.parse_args()

class BedMerger:
    def __init__(self, out_stream):
        self.out_stream = out_stream
        self.last_chrom = None
        self.last_start = -1
        self.last_end = -1
        self.last_feature = None

    def add(self, chrom, start, end, feature):
        # Check if we can merge with the previous bin
        if (chrom == self.last_chrom and 
            start == self.last_end and 
            feature == self.last_feature):
            # Extend the end coordinate
            self.last_end = end
        else:
            # Flush the previous buffer if it exists
            self.flush_buffer()
            # Start a new buffer
            self.last_chrom = chrom
            self.last_start = start
            self.last_end = end
            self.last_feature = feature

    def flush_buffer(self):
        if self.last_chrom is not None:
            self.out_stream.write(f"{self.last_chrom}\t{self.last_start}\t{self.last_end}\t{self.last_feature}\n")
            self.last_chrom = None

def get_best_feature(counts, prioritize_specific):
    if not prioritize_specific:
        return max(counts, key=counts.get)
    
    # Try to find the best "specific" feature (ends with '_specific' or equals 'novel')
    # We use a generator to filter keys, then find the max based on values (overlap length)
    specific_winner = max((k for k in counts if k.endswith("_specific")), key=counts.get, default=None)
    
    if specific_winner is not None:
        return specific_winner
    
    # Fallback: No specific feature found, return standard max overlap
    return max(counts, key=counts.get)

def flush_bins(active_bins, current_chrom, limit_idx, bin_size, merger, max_end_coord, prioritize_specific):
    """
    Outputs and deletes bins from the dictionary that are strictly 
    before the limit_idx.
    """
    # We sort keys to ensure output is sorted by position
    # We only look at keys that are ready to be flushed
    sorted_keys = sorted([k for k in active_bins.keys() if k < limit_idx])
    
    for bin_idx in sorted_keys:
        counts = active_bins[bin_idx]
        
        best_feature = get_best_feature(counts, prioritize_specific)
        
        start = bin_idx * bin_size
        end = start + bin_size
        
        # Truncate the bin if it extends beyond the maximum coordinate seen so far.
        final_end = min(end, max_end_coord)
        
        if final_end > start:
            merger.add(current_chrom, start, final_end, best_feature)
        
        # Delete to free memory
        del active_bins[bin_idx]

def flush_all(active_bins, current_chrom, bin_size, merger, max_end_coord, prioritize_specific):
    """
    Flushes everything remaining in the buffer (used at end of chrom or file).
    """
    if not active_bins:
        return

    for bin_idx in sorted(active_bins.keys()):
        counts = active_bins[bin_idx]
        best_feature = get_best_feature(counts, prioritize_specific)
        
        start = bin_idx * bin_size
        end = start + bin_size
        
        # Truncate the bin if it extends beyond the maximum coordinate seen so far.
        final_end = min(end, max_end_coord)
        
        if final_end > start:
            merger.add(current_chrom, start, final_end, best_feature)
    
    active_bins.clear()

def main():
    args = get_args()
    bin_size = args.binsize
    prioritize_specific = args.specific
    
    # Handle Input Stream (Updated for gzip support)
    if args.input == '-':
        in_stream = sys.stdin
    elif args.input.endswith('.gz'):
        in_stream = gzip.open(args.input, 'rt') # Read as text
    else:
        in_stream = open(args.input, 'r')

    # Handle Output Stream
    if args.output is None or args.output == '-':
        out_stream = sys.stdout
    else:
        out_stream = open(args.output, 'w')
    
    # Helper to merge adjacent bins
    merger = BedMerger(out_stream)
    
    # Data structure: { bin_index: { feature_name: total_overlap_len } }
    active_bins = defaultdict(lambda: defaultdict(int))
    
    current_chrom = None
    current_max_end = 0
    
    for line in in_stream:
        # Fast parsing
        parts = line.rstrip().split('\t')
        
        # Safety check for empty lines
        if not parts: 
            continue
            
        chrom = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
            feature = parts[3]
        except (IndexError, ValueError):
            continue # Skip malformed lines
            
        # Handle Chromosome Switch
        if chrom != current_chrom:
            if current_chrom is not None:
                # Flush using the max_end of the OLD chromosome
                flush_all(active_bins, current_chrom, bin_size, merger, current_max_end, prioritize_specific)
                merger.flush_buffer()
            current_chrom = chrom
            current_max_end = 0
            
        # Update the max coordinate seen for this chromosome.
        current_max_end = max(current_max_end, end)
        
        start_bin_idx = start // bin_size
        # (end - 1) handles cases where end aligns exactly with bin boundary
        end_bin_idx = (end - 1) // bin_size
        
        # Flush older bins to keep memory footprint low
        if active_bins:
            # We can safely print any bin index less than the current feature's start bin
            flush_bins(active_bins, current_chrom, start_bin_idx, bin_size, merger, current_max_end, prioritize_specific)
        
        # Add overlaps to relevant bins
        for b_idx in range(start_bin_idx, end_bin_idx + 1):
            bin_start_coord = b_idx * bin_size
            bin_end_coord = bin_start_coord + bin_size
            
            # Calculate overlap
            overlap_start = max(start, bin_start_coord)
            overlap_end = min(end, bin_end_coord)
            overlap_len = overlap_end - overlap_start
            
            if overlap_len > 0:
                active_bins[b_idx][feature] += overlap_len

    # End of file: flush whatever is left
    if current_chrom is not None and active_bins:
        flush_all(active_bins, current_chrom, bin_size, merger, current_max_end, prioritize_specific)
    merger.flush_buffer()

    # Cleanup file handles (streams)
    if in_stream is not sys.stdin:
        in_stream.close()
    if out_stream is not sys.stdout:
        out_stream.close()

if __name__ == "__main__":
    main()
