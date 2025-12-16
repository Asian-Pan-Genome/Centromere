import os
import sys
import subprocess
import pandas as pd

'''
This script firstly adds rDNA annotation if rDNA file exists,
and changed color for different satellite,
and merges intervals related to specific satellite types.
'''

# Satellite color mapping
satellite_hexcolor = {
    'ASat': '#891640',
    'HSat1A': '#18988B',
    'HSat1B': '#18988B',
    'HSat2': '#323366',
    'HSat3': '#77A0BA',
    'BSat': '#DCAED0',
    'GSat': '#758660',
    'rDNA': '#000000',
    'Censat': '#959595',
    'ct': '#E0E0E0',
}


def add_rDNA(f_rdna):
    """Add rDNA annotation from file if it exists"""
    if f_rdna and os.path.exists(f_rdna) and os.path.getsize(f_rdna) > 0:
        rdna = pd.read_csv(f_rdna, sep="\t", header=None, names=['chrom', 'start', 'end'])
        rdna['type'] = 'rDNA'
        rdna['score'] = '.'
        rdna['strand'] = '+'
        rdna['s'] = rdna['start']
        rdna['e'] = rdna['end']
        rdna['color'] = "#000000"
        print(f"Loaded rDNA annotation from {f_rdna}")
        print(rdna.head())
        return rdna
    else:
        print("No valid rDNA file provided or file is empty. Proceeding without rDNA annotation.")
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'type', 'score', 'strand', 's', 'e', 'color'])


def rename_type(row):
    """Rename satellite types"""
    if not any(x in row for x in ['HSat', 'BSat', 'GSat', 'rDNA']):
        return 'ASat'
    return row


def merge_intervals_with_distance(bed_file, rdna_df, output_file):
    """
    Merge intervals with specified distances, optionally including rDNA annotation
    """
    # Read the main annotation bed file
    annobed = pd.read_csv(bed_file, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'type', 'score', 'strand', 's', 'e', 'color'])
    
    # Rename types
    annobed['type'] = annobed['type'].astype(str)
    annobed['type'] = annobed['type'].apply(rename_type)
    
    # Combine with rDNA if available
    if not rdna_df.empty:
        bed = pd.concat([annobed, rdna_df], axis=0, ignore_index=True).copy()
        print(f"Added {len(rdna_df)} rDNA intervals")
    else:
        bed = annobed.copy()
        print("Processing without rDNA intervals")
    
    # Apply colors from the satellite_hexcolor dictionary
    bed['color'] = bed['type'].map(satellite_hexcolor)
    
    # Sort for merging
    bed = bed.sort_values(by=['chrom', 'start', 'type']).reset_index(drop=True)
    
    # Define merge distances for different satellite types
    merge_distance_dict = {
        'ASat': 10000,
        'BSat': 5000,
        'GSat': 5000,
        'rDNA': 0 
    }
    
    merged_intervals = []
    
    if len(bed) == 0:
        print("Warning: No intervals to process")
        # Create empty output file
        pd.DataFrame().to_csv(output_file, sep='\t', header=False, index=False)
        return pd.DataFrame()
    
    current_interval = bed.iloc[0].copy()
    
    for i in range(1, len(bed)):
        next_interval = bed.iloc[i]
        merge_distance = merge_distance_dict.get(current_interval['type'], 0)
        
        # Check if intervals are from the same chromosome and type
        # and if they are within the specified merge distance
        if (current_interval['chrom'] == next_interval['chrom'] and
            current_interval['type'] == next_interval['type'] and
            current_interval['strand'] == next_interval['strand'] and
            current_interval['end'] + merge_distance >= next_interval['start']):
            # Merge intervals by updating the end position
            current_interval['end'] = max(current_interval['end'], next_interval['end'])
            current_interval['e'] = max(current_interval['end'], next_interval['end'])
        else:
            merged_intervals.append(current_interval)
            current_interval = next_interval.copy()
    
    merged_intervals.append(current_interval)
    
    merged_df = pd.DataFrame(merged_intervals)
    merged_df['score'] = '0'
    
    print(f"Merged from {len(bed)} to {len(merged_df)} intervals")
    print(f"Output shape: {merged_df.shape}")
    
    # Save merged results
    merged_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    return merged_df


def compress_bed(bed_file):
    """Compress BED file with bgzip and index with tabix"""
    if os.path.exists(bed_file) and os.path.getsize(bed_file) > 0:
        cmd0 = ["bgzip", "-k", "-f", bed_file]
        cmd1 = ["tabix", "-p", "bed", bed_file + ".gz"]
        subprocess.run(cmd0, check=True)
        subprocess.run(cmd1, check=True)
    else:
        print(f"Warning: {bed_file} is empty or doesn't exist. Skipping compression.")


def main():
    """Main function to handle command line arguments"""
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_bed> <output_bed> [rDNA_file]")
        print("  input_bed: Input BED file")
        print("  output_bed: Output BED file")
        print("  rDNA_file: Optional rDNA annotation file (can be empty or non-existent)")
        sys.exit(1)
    
    inbed = sys.argv[1]
    outbed = sys.argv[2]
    
    # Check if rDNA file argument is provided
    if len(sys.argv) > 3:
        f_rdna = sys.argv[3]
        # Check if the file exists and is not empty
        if f_rdna and f_rdna.lower() != 'none' and os.path.exists(f_rdna) and os.path.getsize(f_rdna) > 0:
            print(f"Using rDNA file: {f_rdna}")
            rdna = add_rDNA(f_rdna)
        else:
            print(f"rDNA file '{f_rdna}' is empty or doesn't exist. Proceeding without rDNA.")
            rdna = pd.DataFrame(columns=['chrom', 'start', 'end', 'type', 'score', 'strand', 's', 'e', 'color'])
    else:
        print("No rDNA file provided. Proceeding without rDNA annotation.")
        rdna = pd.DataFrame(columns=['chrom', 'start', 'end', 'type', 'score', 'strand', 's', 'e', 'color'])
    
    # Merge intervals
    merged_df = merge_intervals_with_distance(inbed, rdna, outbed)
    
    # Compress the output BED file
    compress_bed(outbed)
    
    print(f"Processing completed. Output saved to {outbed}")


if __name__ == "__main__":
    main()
