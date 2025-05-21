#!/usr/bin/env python
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

def compute_hamming(barcode_list,barcode_of_interest):
    barcode_dist = dict()
    ### initialize
    for i in barcode_list:
        #print(f"Comparing {i} to {barcode_of_interest}")
        barcode_dist[i] = hamming_distance(i,barcode_of_interest)
    return barcode_dist

# read in barcodes summary to get whitelist
def get_whitelist(barcodes_summary_file):
    barcode_whitelist = []
    with open(barcodes_summary_file,"r") as f:
        for line in f.readlines():
            line_cleaned = line.rstrip()
            original_line = line_cleaned.split('\t')
            original_line = [str(x) for x in original_line]
            if original_line[0] == "ID":
                print("Found Header")
            else:
                barcode_whitelist.append(original_line[0])
    return barcode_whitelist

def get_barcode_size(barcodes_file):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            lines_processed += 1
    return lines_processed

def get_barcode(barcodes_file,idx):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            lines_processed += 1
            if lines_processed == idx:
                return str(line.strip("\n"))
#####################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_summary_file',default=None, required=True, type=str, help="barcodeSummary.tsv file from DRAGEN")
    parser.add_argument('--barcode_file',default=None, required=True, type=str, help="plain text barcode list extracted from FASTQ")
    parser.add_argument('--output_file', default=None, type=str, help="output file path")
    args, extras = parser.parse_known_args()
    ################################
    barcode_summary_file = args.barcode_summary_file
    barcode_file = args.barcode_file
    output_file = args.output_file

    if output_file is None:
        output_dir = os.path.dirname(barcode_summary_file)
        initial_file_name_split = os.path.basename(barcode_summary_file).split(".")
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["R1_barcode_hamming_distance",initial_file_name_split[-1]]
        output_filename = ".".join(output_filename_components)
        output_file = output_dir + "/"+ output_filename
        logging_statement(f"No output file provided. Will write out hamming distance table to {output_file}")
    logging_statement(f"Reading in Barcode Summary Table {barcode_summary_file}")
    whitelisted_barcodes = get_whitelist(barcode_summary_file)
    logging_statement(f"Reading in FASTQ barcodes file {barcode_file}")
    fastq_barcodes_size = get_barcode_size(barcode_file)
    logging_statement("Finished reading in FASTQ barcodes line")
    output_lines = []
    barcodes_processed = 0
    ## create header line for barcode distance table
    logging_statement("Creating header line")
    header_line = ["read1_barcode"]
    for w in whitelisted_barcodes:
        header_line.append(w)
    header_line_str = "\t".join(header_line)
    output_lines.append(header_line_str)
    with open(output_file,"w") as outfile:
        outfile.write(header_line_str + "\n")
 
    #### create lines -- each row contains the hamming distances of each read1 barcode against the white list barcodes
    barcodes_processed = 0
    while barcodes_processed < fastq_barcodes_size:
        output_line = []
        barcodes_processed += 1
        if barcodes_processed % 10000 == 0:
            logging_statement(f"Computed Hamming Distance for {barcodes_processed} lines")
        ##############    
        fastq_barcode = get_barcode(barcode_file,barcodes_processed)
        output_line.append(fastq_barcode)
        #print(f"Compute hamming distances of whitelist against {fastq_barcode}")
        hamming_distances_dict = compute_hamming(whitelisted_barcodes,fastq_barcode)
        for whitelisted_barcode in whitelisted_barcodes:
            val_to_add = hamming_distances_dict[whitelisted_barcode]
            output_line.append(str(val_to_add))
        output_line_str = "\t".join(output_line)
        with open(output_file,"a+") as outfile:
            outfile.write(output_line_str + "\n")

    logging_statement(f"Finished writing {output_file}")


if __name__ == "__main__":
    main()    