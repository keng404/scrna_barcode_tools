#!/usr/bin/env python
### adapted from https://github.com/lh3/biofast/blob/master/fqcnt/fqcnt_py7x_pysam.py
import sys
import pysam
import os
import sys
import re 
import argparse 
from datetime import datetime as dt
import collections
##########
def logging_statement(string_to_print):
    date_time_obj = dt.now()
    timestamp_str = date_time_obj.strftime("%Y/%b/%d %H:%M:%S:%f")
    #############
    final_str = f"[ {timestamp_str} ] {string_to_print}"
    return sys.stderr.write(f"{final_str}\n")
#############
def whitelist_component_lookup(whitelist):
    whitelist_dict = dict()
    for barcode in whitelist:
        d = collections.defaultdict(list)
        for k, v in enumerate(barcode):
            d[v].append(k)
        whitelist_dict[barcode] = d
    return whitelist_dict

def barcode_to_dict(barcode):
    d = collections.defaultdict(list)
    for k, v in enumerate(barcode):
        d[v].append(k)  
    return d  

def hamming_distance_by_dict(dict1,dict2):
    ### taking in 2 default dicts
    keys_to_check = list(set(list(dict1.keys()) + list(dict2.keys())))
    hamming_collection = []
    for ktc in keys_to_check:
        l1 = dict1[ktc]
        l2 = dict2[ktc]
        if len(l2) > len(l1):
            for i in l2:
                if i not in l1:
                    hamming_collection.append(i)
        elif len(l1) >= len(l2):
            for i in l1:
                if i not in l2:
                    hamming_collection.append(i)
    return len(list(set(hamming_collection)))

def compute_hamming_v2(barcode_of_interest,barcode_defaultdict):
    barcode_dist = []
    barcode_dict = barcode_to_dict(barcode_of_interest)
    barcode_list = list(barcode_defaultdict.keys())
    ### initialize
    for i in barcode_list:
        #sys.stderr.write(f"Comparing {i} to {barcode_of_interest}")
        barcode_dist.append(hamming_distance_by_dict(barcode_dict,barcode_defaultdict[i]))
    return barcode_dist

def hamming_distance(string1, string2):
    if len(string1) != len(string2):
        raise ValueError("Strings must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(string1, string2))

def compute_hamming(barcode_of_interest,barcode_list):
    barcode_dist = []
    ### initialize
    for i in barcode_list:
        #sys.stderr.write(f"Comparing {i} to {barcode_of_interest}")
        barcode_dist.append(hamming_distance(i,barcode_of_interest))
    return barcode_dist
# read in barcodes summary to get whitelist
def get_whitelist(barcodes_file):
    barcode_whitelist = dict()
    key = None
    line_num = 0
    with open(barcodes_file,"r", encoding='utf-8') as f:
        content = f.readlines()
        for line in content:
            line_num += 1
            line_cleaned = line.strip()
            #print(f"{line}")
            if re.search("#",line_cleaned) is not None:
                line_split = line_cleaned.split("-")
                key = line_split[len(line_split)-1]
                sys.stderr.write(f"Found key {key}\n")
            else:
                sys.stderr.write(f"LINE_NUM: {line_num},LINE: {line_cleaned},key: {key}\n")
                if len(list(barcode_whitelist.keys())) > 0:
                    if key in list(barcode_whitelist.keys()):
                        ### append to array
                        barcode_whitelist[key].append(line_cleaned)
                    else:
                        barcode_whitelist[key] = [line_cleaned]
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
    output_file = None
    if output_file is None:
        output_dir = os.path.dirname(barcode_file)
        initial_file_name_split = os.path.basename(barcode_file).split(".")
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["R1_barcode_hamming_distance",initial_file_name_split[-1]]
        output_filename = ".".join(output_filename_components)
        output_file = output_dir + "/"+ output_filename
        logging_statement(f"No output file provided. Will write out hamming distance table to {output_file}")
    ##############################################
    scrna_barcode_str = args.scrna_barcode
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    header_line_str = f"barcode_sequence,block1_dist,block2_dist,block3_dist,block4_dist,total_dist"
    print(f"{header_line_str}")
    #with open(output_file,"w") as outfile:
    #    outfile.write(header_line_str + "\n")
    num_records = 0
    with pysam.FastxFile(fastq) as fastx:
        for record in fastx:
            num_records += 1
            if num_records % 10000 == 0:
                logging_statement(f"Processed {num_records} reads in {fastq}")
            sequence_split = [x for x in str(record.sequence)]
            barcode_str = ""
            barcode_blocks = []
            for idx,val in enumerate(scrna_barcode_starts):
                barcode_subset = sequence_split[int(scrna_barcode_starts[idx]):(int(scrna_barcode_ends[idx])+1)]
                barcode_subset_str = "".join(barcode_subset)
                barcode_blocks.append(barcode_subset_str)
                barcode_str += "".join(barcode_subset)
            ############
            block1_defaultdict = whitelist_component_lookup(whitelist_lookup['Block1'])
            block1_dists = compute_hamming_v2(barcode_blocks[0],block1_defaultdict)
            block1_dist = min(block1_dists)
            
            block2_defaultdict = whitelist_component_lookup(whitelist_lookup['Block2'])
            block2_dists = compute_hamming_v2(barcode_blocks[1],block2_defaultdict)
            block2_dist = min(block2_dists)
            
            block3_defaultdict = whitelist_component_lookup(whitelist_lookup['Block3'])
            block3_dists = compute_hamming_v2(barcode_blocks[2],block3_defaultdict)
            block3_dist = min(block3_dists)
            
            block4_defaultdict = whitelist_component_lookup(whitelist_lookup['Block4'])
            block4_dists = compute_hamming_v2(barcode_blocks[3],block4_defaultdict)
            block4_dist = min(block4_dists)
            total_dist = block1_dist + block2_dist + block3_dist + block4_dist
            
            output_line_str = f"{barcode_str},{block1_dist},{block2_dist},{block3_dist},{block4_dist},{total_dist}" 
            print(f"{output_line_str}")
            #with open(output_file,"a+") as outfile:
            #    outfile.write(output_line_str + "\n")

if __name__ == "__main__":
    main()    