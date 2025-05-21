#!/usr/bin/env python
### adapted from https://github.com/lh3/biofast/blob/master/fqcnt/fqcnt_py7x_pysam.py
import sys
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
    return sys.stderr.write(f"{final_str}\n")
#############
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
##########################
def break_barcode_blocks(barcode):
    s1 = "".join(barcode[0:8])
    s2 = "".join(barcode[8:14])
    s3 = "".join(barcode[14:20])
    s4 = "".join(barcode[20:28])
    barcode_components = [s1,s2,s3,s4]
    return barcode_components

def get_tags(line_parsed):
    parsed_dict = dict()
    tags_of_interest = ["XB","CR","RX"]
    for i in line_parsed:
        if re.search(":",i) is not None:
            i_split = i.split(":")
            if i_split[0] in tags_of_interest:
                parsed_dict[i_split[0]] = i_split[len(i_split)-1]
    return parsed_dict



def main():
    tags_of_interest = ["XB","CR","RX"]
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_file',default=None, required=True, type=str, help="plain text whitelist barcode file")
    parser.add_argument('--bam_tag_records',default=None,type=str,help="SAM File containing bam tag records (XB,CR tags)")
    args, extras = parser.parse_known_args()
    ################################
    barcode_file = args.barcode_file
    bam_tag_records = args.bam_tag_records
    whitelist_lookup = get_whitelist(barcode_file)
    ### output_file1 is for valid barcodes
    output_file1 = None
    if output_file1 is None:
        output_dir = os.path.dirname(bam_tag_records)
        initial_file_name_split = os.path.basename(bam_tag_records).split(".")
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["barcode_hamming_distance","csv"]
        output_filename = ".".join(output_filename_components)
        output_file1 = output_dir + "/"+ output_filename
        logging_statement(f"No output file provided. Will write out hamming distance of valid barcodes to {output_file1}")
    #### output_file2 is for invalid barcodes
    #output_file2 = None
    #if output_file2 is None:
    #    output_dir = os.path.dirname(bam_tag_records)
    #    initial_file_name_split = os.path.basename(bam_tag_records).split(".")
    #    output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["filtered_barcode_hamming_distance","csv"]
    #    output_filename = ".".join(output_filename_components)
    #    output_file2 = output_dir + "/"+ output_filename
    #    logging_statement(f"No output file provided. Will write out hamming distance of filtered barcodes to {output_file2}")
    ##############################################
    header_line_str1 = f"read_name,barcode_sequence,corrected_sequence,binning_index,total_dist"
    #header_line_str2 = f"read_name,barcode_sequence,block1_dist,block2_dist,block3_dist,block4_dist,total_dist"
    ### write output file headers
    with open(output_file1,"w") as outfile:
        outfile.write(header_line_str1 + "\n")
    #with open(output_file2,"w") as outfile:
    #    outfile.write(header_line_str2 + "\n")
    ###############################
    num_records = 0
    reads_skipped = 0
    corrected_sequence_compare = 0
    with open(bam_tag_records,"r", encoding='utf-8') as f:
        for line in f:
            num_records += 1
            if num_records % 1000000 == 0:
                logging_statement(f"Processed {num_records} lines from {bam_tag_records}")
            line_cleaned = line.strip()
            line_split = line_cleaned.split("\t")
            #print(f"{line_split}")
            read_name = line_split[0]
            mapping_quality = int(line_split[4])
            if mapping_quality >= 1:
                line_dict = get_tags(line_split)
                binning_index = None
                if 'RX' in list(line_dict.keys()):
                    binning_index = line_dict['RX']
                #if 'XB' not in list(line_dict.keys()):
                #    print(f"{line_split}")
                if 'CR' in list(line_dict.keys()) and 'XB' not in list(line_dict.keys()):
                    line_dict['XB'] = line_dict['CR']
                if 'XB' not in list(line_dict.keys()):
                    reads_skipped += 1
                    #sys.exit(1)
                else:    
                    barcode_sequence = line_dict['XB']

                    if "CR" in list(line_dict.keys()):
                        corrected_sequence_compare += 1
                        #print(f"{line_split}")
                        #logging_statement(f"Found valid barcode : {line_split[0]}\t{line_dict}")
                        corrected_sequence = line_dict['CR']
                        #print(f"Comparing {line_dict['XB']} and {line_dict['CR']}")
                        total_dist = hamming_distance(line_dict['XB'],line_dict['CR'])
                        #print(f"Total_dist: {total_dist}")
                        #if corrected_sequence_compare == 10:
                        #    sys.exit(0)
                        output_str = f"{read_name},{barcode_sequence},{corrected_sequence},{binning_index},{total_dist}"
                        with open(output_file1,"a+") as outfile:
                            outfile.write(output_str + "\n")
                    else:
                        #logging_statement(f"Found filtered barcode : {line_split[0]}\t{line_dict}")
                        barcode_blocks = break_barcode_blocks(line_dict["XB"])
                        block1_dists = compute_hamming(barcode_blocks[0],whitelist_lookup["Block1"])
                        block1_dist = min(block1_dists)
                        block1_str = None
                        for i,v in enumerate(whitelist_lookup["Block1"]):
                            if block1_dists[i] == block1_dist:
                                block1_str = v
                        block2_dists = compute_hamming(barcode_blocks[1],whitelist_lookup["Block2"])
                        block2_dist = min(block2_dists)
                        block2_str = None
                        for i,v in enumerate(whitelist_lookup["Block2"]):
                            if block2_dists[i] == block2_dist:
                                block2_str = v
                        block3_dists = compute_hamming(barcode_blocks[2],whitelist_lookup["Block3"])
                        block3_dist = min(block3_dists)
                        block3_str = None
                        for i,v in enumerate(whitelist_lookup["Block3"]):
                            if block3_dists[i] == block3_dist:
                                block3_str = v
                        block4_dists = compute_hamming(barcode_blocks[3],whitelist_lookup["Block4"])
                        block4_dist = min(block4_dists)
                        block4_str = None
                        for i,v in enumerate(whitelist_lookup["Block4"]):
                            if block4_dists[i] == block4_dist:
                                block4_str = v
                        total_dist = block1_dist + block2_dist + block3_dist + block4_dist
                        ########################################
                        corrected_sequence = block1_str + block2_str + block3_str + block4_str
                        output_line_str = f"{read_name},{barcode_sequence},{corrected_sequence},{binning_index},{total_dist}" 
                        with open(output_file1,"a+") as outfile:
                            outfile.write(output_line_str + "\n")     
            else:
                reads_skipped += 1       
    logging_statement(f"Skipped {reads_skipped} out of {num_records} reads")
if __name__ == "__main__":
    main()    