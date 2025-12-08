import os
import sys
import subprocess
import pandas as pd
from color_anno import  *
from pybedtools import BedTool

'''
This script firstly adds rDNA annotation,
and changed color for different satellite,
and merges invervals related to specific satellite types,
and finally convert into bigbed format.
'''


def add_rDNA(f_rdna):
    rdna = pd.read_csv(f_rdna, sep = "\t", header = None, names = ['chrom', 'start', 'end'])
    rdna['type'] = 'rDNA'
    rdna['score']  = '.'
    rdna['strand'] = '+'
    rdna['s'] = rdna['start']
    rdna['e'] = rdna['end']
    rdna['color'] = "#000000"

    print(rdna.head())
    return rdna



def merge_intervals_with_distance(bed_file, rdna, output_file):

    def rename_type(row):
        if not any(x in row for x in ['HSat', 'BSat', 'GSat', 'rDNA']):
            return 'ASat'
        return row

    annobed = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'type', 'score', 'strand', 's',  'e', 'color'])
    annobed['type'] = annobed['type'].apply(rename_type)
    bed = pd.concat([annobed, rdna], axis=0, ignore_index=True).copy()
    bed['color'] = bed['type'].map(satellite_hexcolor)
    
    bed = bed.sort_values(by=['chrom', 'start', 'type']).reset_index(drop=True)

    print(bed.head(30))

    merge_distance_dict = {
    'ASat': 10000,
    'BSat': 5000,
    'GSat': 5000,
    }
    merged_intervals = []

    current_interval = bed.iloc[0].copy()

    for i in range(1, len(bed)):
        next_interval = bed.iloc[i]
        merge_distance = merge_distance_dict.get(current_interval['type'], 0)

        # Check if intervals are from the same chromosome and type
        # and if they are within the specified merge distance
        if (current_interval['chrom'] == next_interval['chrom'] and
            current_interval['type'] == next_interval['type'] and
            current_interval['end'] + merge_distance >= next_interval['start']):
            # Merge intervals by updating the end position
            current_interval['end'] = max(current_interval['end'], next_interval['end'])
        else:
            merged_intervals.append(current_interval)
            current_interval = next_interval.copy()

    merged_intervals.append(current_interval)

    merged_df = pd.DataFrame(merged_intervals)
    merged_df['score'] = '0'
    print(merged_df.shape)
    
    merged_df.to_csv(output_file, sep='\t', header=False, index=False)

    return merged_df

def convert_bed_to_bigbed(bed_file, bigbed_file, ffai):
    base_name = os.path.basename(ffai)
    prefix_name, _ = os.path.splitext(base_name)
    chrom_size_file = prefix_name + ".chrom.sizes" 

    df = pd.read_csv(ffai, sep = "\t", names = ['chrom', 'length', 'col3', 'col4', 'col5'])
    subdf = df[['chrom', 'length']]
    subdf.to_csv(chrom_size_file, sep = "\t", header = False, index = False)
    
    cmd = ['bedToBigBed', bed_file, chrom_size_file , bigbed_file]
    subprocess.run(cmd, check=True)

def compress_bed(bed_file):
    cmd0 = ["bgzip", "-k", bed_file]
    cmd1 = ["tabix", "-p", "bed", bed_file+".gz" ]
    subprocess.run(cmd0, check=True)
    subprocess.run(cmd1, check=True)


if __name__ == "__main__":
    inbed = sys.argv[1]
    outbed = sys.argv[2]
    f_rdna = sys.argv[3]
    #ffai = sys.argv[4]
    #bigbed_file = outbed + ".bb"
    rdna =  add_rDNA(f_rdna)
    merged_df = merge_intervals_with_distance(inbed, rdna, outbed)
    #convert_bed_to_bigbed(outbed, bigbed_file, ffai)
    compress_bed(outbed)
