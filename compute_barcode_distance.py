#!/usr/bin/env python
### adapted from https://github.com/lh3/biofast/blob/master/fqcnt/fqcnt_py7x_pysam.py
import sys
import pysam
import os
import sys
import re 
import argparse 
from datetime import datetime as dt
##########
def logging_statement(string_to_print):
    date_time_obj = dt.now()
    timestamp_str = date_time_obj.strftime("%Y/%b/%d %H:%M:%S:%f")
    #############
    final_str = f"[ {timestamp_str} ] {string_to_print}"
    return print(f"{final_str}")
#############
def hamming_distance(string1, string2):
    if len(string1) != len(string2):
        raise ValueError("Strings must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(string1, string2))

def compute_hamming(barcode_of_interest,barcode_list):
    barcode_dist = dict()
    ### initialize
    for i in barcode_list:
        #print(f"Comparing {i} to {barcode_of_interest}")
        barcode_dist[i] = hamming_distance(i,barcode_of_interest)
    return barcode_dist
# read in barcodes summary to get whitelist
def get_whitelist(barcodes_file):
    barcode_whitelist = dict()
    key = None
    with open(barcodes_file,"r") as f:
        for line in f.readlines():
            line_cleaned = line.rstrip()
            if re.search("^#",line_cleaned):
                print("Found key")
                line_split = line_cleaned.split("-")
                key = line_split[len(line_split)-1]
            else:
                if key in list(barcode_whitelist.keys()):
                    ### append to array
                    barcode_whitelist[key] = barcode_whitelistp[key].append(line_cleaned)
                else:
                    ## initialize array
                    barcode_whitelist[key] = [line_cleaned]
    return barcode_whitelist

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_file',default=None, required=True, type=str, help="plain text whitelist barcode file")
    parser.add_argument('--fastq', default=None, type=str, help="FASTQ file path")
    parser.add_argument('--scrna_barcode',default="0_7+11_16+20_25+31_38",type=str,help="string containing barcode index")
    args, extras = parser.parse_known_args()
    ################################
    barcode_file = args.barcode_file
    fastq = args.fastq
    whitelist_lookup = get_whitelist(barcode_file)
    ##############################################
    scrna_barcode_str = args.scrna_barcode
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    print(f"barcode_sequence,block1_dist,block2_dist,block3_dist,block4_dist,total_dist")
    with pysam.FastxFile(fastq) as fastx:
        for record in fastx:
            sequence_split = [x for x in str(record.sequence)]
            barcode_str = ""
            for idx,val in enumerate(scrna_barcode_starts):
                barcode_subset = sequence_split[int(scrna_barcode_starts[idx]):(int(scrna_barcode_ends[idx])+1)]
                barcode_str += "".join(barcode_subset)
                block1_dists = compute_hamming(barcode_subset[0],whitelist_lookup["Block1"])
                block1_dist = min(block1_dist)
                block2_dists = compute_hamming(barcode_subset[1],whitelist_lookup["Block2"])
                block2_dist = min(block2_dist)
                block3_dists = compute_hamming(barcode_subset[2],whitelist_lookup["Block3"])
                block3_dist = min(block3_dist)
                block4_dists = compute_hamming(barcode_subset[3],whitelist_lookup["Block4"])
                block4_dist = min(block4_dist)
                total_dist = block1_dist + block2_dist + block3_dist + block4_dist
            #qlen += len(record.quality)
            #print(f"{barcode_str},{min(block1_dist)},{min(block2_dist)},{min(block3_dist)},{min(block4_dist)},{total_dist}"})

if __name__ == "__main__":
    main()    