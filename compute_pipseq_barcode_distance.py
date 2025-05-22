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
#### 
def get_readname_list_len(read_names_file):
    line_num = 0
    with open(read_names_file,"r", encoding='utf-8') as f:
        content = f.readlines()
        for line in content:
            line_num += 1
    return line_num
# read in read names  to get whitelist of readnames to extract barcodes for
def read_name_whitelist(read_names_file,line_number_to_start = 0):
    read_name_whitelist = []
    key = 1
    num_lines_to_block = 1000000
    line_num_end = line_number_to_start + num_lines_to_block
    line_num = 0
    with open(read_names_file,"r", encoding='utf-8') as f:
        content = f.readlines()
        for line in content:
            line_num += 1
            line_cleaned = line.strip()
            if line_num >= line_number_to_start and line_num < line_num_end:
                read_name_whitelist.append(line_cleaned)
            if line_num >= line_num_end:
                break
    whitelist_obj = dict()
    whitelist_obj['read_name_whitelist'] = read_name_whitelist
    whitelist_obj['line_num'] = line_num
    return whitelist_obj

############################################################
def get_str(dists,barcodes,min_dist):
    my_str = None
    for idx,val in enumerate(dists):
        if val == min_dist:
            my_str = barcodes[idx]
            break
    return my_str

def compute_hamming_by_tier(barcode_blocks,whitelist_lookup):
    block1_dists = compute_hamming(barcode_blocks[0],whitelist_lookup["Block1"])
    block1_dist = min(block1_dists)
    block1_str = get_str(block1_dists,whitelist_lookup["Block1"],block1_dist)
    block2_dists = compute_hamming(barcode_blocks[1],whitelist_lookup["Block2"])
    block2_dist = min(block2_dists)
    block2_str = get_str(block2_dists,whitelist_lookup["Block2"],block2_dist)
    block3_dists = compute_hamming(barcode_blocks[2],whitelist_lookup["Block3"])
    block3_dist = min(block3_dists)
    block3_str = get_str(block3_dists,whitelist_lookup["Block3"],block3_dist)
    block4_dists = compute_hamming(barcode_blocks[3],whitelist_lookup["Block4"])
    block4_dist = min(block4_dists)
    block4_str = get_str(block4_dists,whitelist_lookup["Block4"],block4_dist)
    hamming_dist_dict = dict()
    hamming_dist_dict['block1_dists'] = block1_dists
    hamming_dist_dict['block1_dist'] = block1_dist
    hamming_dist_dict['block1_str'] = block1_str
    hamming_dist_dict['block2_dists'] = block2_dists
    hamming_dist_dict['block2_dist'] = block2_dist
    hamming_dist_dict['block2_str'] = block1_str
    hamming_dist_dict['block3_dists'] = block3_dists
    hamming_dist_dict['block3_dist'] = block3_dist
    hamming_dist_dict['block3_str'] = block1_str
    hamming_dist_dict['block4_dists'] = block4_dists
    hamming_dist_dict['block4_dist'] = block4_dist
    hamming_dist_dict['block4_str'] = block1_str
    return hamming_dist_dict
####
def get_barcode_sequences(sequence_of_interest, scrna_barcode_str,pad=0):
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    reverse_complement = False
    if len(scrna_barcode_intervals) < 2:
        scrna_barcode_intervals = scrna_barcode_str.split("-")
        reverse_complement = True # Not sure if this is needed for this specifc use-case
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    barcode_str = ""
    barcode_blocks = []
    for idx,val in enumerate(scrna_barcode_starts):
        start_idx = int(scrna_barcode_starts[idx]) + pad
        end_idx = (int(scrna_barcode_ends[idx])+1) + pad
        barcode_subset = sequence_of_interest[start_idx:end_idx]
        barcode_subset_str = "".join(barcode_subset)
        barcode_blocks.append(barcode_subset_str)
        barcode_str += "".join(barcode_subset)
    ##############################    
    results_dict = dict()
    results_dict['barcode_str'] = barcode_str
    results_dict['barcode_blocks'] = barcode_blocks
    return results_dict

def get_linker_sequences(sequence_of_interest, scrna_barcode_str,pad=0):
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    reverse_complement = False
    if len(scrna_barcode_intervals) < 2:
        scrna_barcode_intervals = scrna_barcode_str.split("-")
        reverse_complement = True # Not sure if this is needed for this specifc use-case
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    scrna_linker_starts = []
    scrna_linker_ends = []
    for idx,val in enumerate(scrna_barcode_starts):
        if idx < len(scrna_barcode_starts)-1:
            end_idx = int(scrna_barcode_starts[idx+1]) -1
            start_idx = int(scrna_barcode_ends[idx]) + 1
            scrna_linker_ends.append(end_idx)
            scrna_linker_starts.append(start_idx)
    linker_blocks = []
    for idx,val in enumerate(scrna_linker_starts):
        start_idx = int(scrna_linker_starts[idx]) + pad
        end_idx = (int(scrna_linker_ends[idx])+1) + pad
        linker_subset = sequence_of_interest[start_idx:end_idx]
        linker_subset_str = "".join(linker_subset)
        linker_blocks.append(linker_subset_str)
    ##############################    
    results_dict = dict()
    results_dict['linker_blocks'] = linker_blocks
    return results_dict

def score_candidate(candidate_scores):
    candidate_scores1 = [int(x) for x in candidate_scores]
    candidate_scores = candidate_scores1
    num_hits = sum([x == 0 for x in candidate_scores])
    total_diff = sum(candidate_scores)
    if num_hits == len(candidate_scores):
        return num_hits
    elif num_hits == 0:
        return float(0.00001/total_diff)
    else:
        return float(num_hits/total_diff)

#### can variable phase 0-3 bases from beginning of read when trying to extract barcode
def phased_check(original_sequence,barcode_string,lookup_list):
    linker_sequences = ['ATG','CAG','TCGAG']
    ### initialize scores
    my_best_hit = dict()
    dummy_barcode = "N" * len(original_sequence)
    my_best_hit['barcode_seq'] = dummy_barcode
    my_best_hit['best_barcode'] = dummy_barcode
    my_best_hit['best_dists'] = [100000,100000,100000,100000]
    my_best_hit['best_tiered_seqs'] = ["N" * 8, "N" * 6,"N" * 6,"N" * 8]
    # CURRRENTLY will not check for linker sequences, but will return barcode with the lowest hamming distance
    #############################
    phase_lag = [0 , 1, 2, 3]
    keys_to_check = ['block1_dist','block2_dist','block3_dist','block4_dist']
    for i in phase_lag:
        candidate_barcodes = get_barcode_sequences(original_sequence, barcode_string,pad = i)
        candidate_linkers = get_linker_sequences(original_sequence, barcode_string,pad = i)
        linkers_match = candidate_linkers['linker_blocks'] == linker_sequences
        candidate_dists = compute_hamming_by_tier(candidate_barcodes['barcode_blocks'],lookup_list)
        #####################################################
        candidate_scores = []
        for k in keys_to_check:
            candidate_scores.append(candidate_dists[k])
        #################################################  
        logging_statement(f"BEST_HIT:\t{my_best_hit['best_dists']}\tCANDIDATE_HIT:\t{candidate_scores}\t{i}")  
        if score_candidate(candidate_scores) > score_candidate(my_best_hit['best_dists']) or linkers_match is True:
            #for idx,j in enumerate(keys_to_check):
            my_best_hit['best_dists'] = candidate_scores
            my_best_hit['best_tiered_seqs'] = [candidate_dists['block1_str'],candidate_dists['block2_str'],candidate_dists['block3_str'],candidate_dists['block4_str']]
            my_best_hit['barcode_seq'] = candidate_barcodes['barcode_str']
            my_best_hit['best_barcode'] = "".join( my_best_hit['best_tiered_seqs'])
    return my_best_hit
    
def found_val(my_list,my_val):
    found_match = False
    if my_val in my_list:
        return True
    return found_match

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_names_file',default=None, type=str, help="plain text whitelist read_names file")
    parser.add_argument('--barcode_file',default=None, required=True, type=str, help="plain text whitelist barcode file")
    parser.add_argument('--fastq', default=None, required=True,type=str, help="FASTQ file path")
    parser.add_argument('--scrna_barcode',default="0_7+11_16+20_25+31_38",type=str,help="string containing barcode index")
    args, extras = parser.parse_known_args()
    ################################
    barcode_file = args.barcode_file
    read_names_file = args.read_names_file
    fastq = args.fastq
    whitelist_lookup = get_whitelist(barcode_file)
    get_read_names = False
    read_name_len = 0
    read_name_start = 0
    #### read in and return dict of read names to extract and compute barcode distance for
    if read_names_file is not None:
        get_read_names = True
        read_name_len = get_readname_list_len(read_names_file)
        logging_statement(f"Found {read_name_len} reads to check")
    else:
        read_names = []
    output_file = None
    if output_file is None:
        if read_names_file is not None:
            output_dir = os.path.dirname(read_names_file)
            initial_file_name_split = os.path.basename(read_names_file).split(".")
        else:
            output_dir = os.path.dirname(fastq)
            initial_file_name_split = os.path.basename(fastq).split(".")            
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["R1_barcode_hamming_distance","csv"]
        output_filename = ".".join(output_filename_components)
        output_file = output_dir + "/"+ output_filename
        logging_statement(f"No output file provided. Will write out hamming distance table to {output_file}")
    ##############################################
    scrna_barcode_str = args.scrna_barcode
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    header_line_str = f"read_name,barcode_sequence_extracted,closest_whitelist_barcode,block1_dist,block1_str,block2_dist,block2_str,block3_dist,block3_str,block4_dist,block4_dist,total_dist"
    #print(f"{header_line_str}")
    with open(output_file,"w") as outfile:
        outfile.write(header_line_str + "\n")
    num_records = 0
    line_num = 0
    if read_names_file is not None:
        with open(read_names_file,"r", encoding='utf-8') as f:
            content = f.readlines()
            for line in content:
                line_num += 1
                line_cleaned = line.strip()
                read_name = line_cleaned
                with pysam.FastxFile(fastq) as fastx:
                    for record in fastx:
                        if num_records % 10000 == 0 and num_records > 0:
                            logging_statement(f"Processed {num_records} reads in {read_names_file}")
                        if record.name == read_name:
                            num_records += 1
                            sequence_split = [x for x in str(record.sequence)]
                            extracted_barcode = phased_check(sequence_split,scrna_barcode_str,whitelist_lookup)
                            ### barcode_seq is extracted sequence
                            barcode_str = extracted_barcode['barcode_seq']
                            ### best_barcode is best matching barcode
                            best_barcode = extracted_barcode['best_barcode']
                            ### best_dists is array of the hamming distances for each tier
                            dists = extracted_barcode['best_dists']
                            block1_dist = dists[0]
                            block2_dist = dists[1]
                            block3_dist = dists[2]
                            block4_dist = dists[3]
                            total_dist = block1_dist + block2_dist + block3_dist + block4_dist
                            ### best_tiered_seqs is array of sequences for each barcode block
                            tiered_seqs  = extracted_barcode['best_tiered_seqs']
                            
                            output_line_str = f"{read_name},{barcode_str},{best_barcode},{block1_dist},{tiered_seqs[0]},{block2_dist},{tiered_seqs[1]},{block3_dist},{tiered_seqs[2]},{block4_dist},{tiered_seqs[3]},{total_dist}" 
                            with open(output_file,"a+") as outfile:
                                outfile.write(output_line_str + "\n")
                            #print(f"{output_line_str}")
                            break
    else:
        num_records = 0
        with pysam.FastxFile(fastq) as fastx:
            for record in fastx:
                num_records += 1
                if num_records % 10000 == 0 and num_records > 0:
                    logging_statement(f"Processed {num_records} reads in {fastq}")

                sequence_split = [x for x in str(record.sequence)]
                extracted_barcode = phased_check(sequence_split,scrna_barcode_str,whitelist_lookup)
                ### barcode_seq is extracted sequence
                barcode_str = extracted_barcode['barcode_seq']
                ### best_barcode is best matching barcode
                best_barcode = extracted_barcode['best_barcode']
                ### best_dists is array of the hamming distances for each tier
                dists = extracted_barcode['best_dists']
                block1_dist = dists[0]
                block2_dist = dists[1]
                block3_dist = dists[2]
                block4_dist = dists[3]
                total_dist = block1_dist + block2_dist + block3_dist + block4_dist
                ### best_tiered_seqs is array of sequences for each barcode block
                tiered_seqs  = extracted_barcode['best_tiered_seqs']
                            
                output_line_str = f"{record.name},{barcode_str},{best_barcode},{block1_dist},{tiered_seqs[0]},{block2_dist},{tiered_seqs[1]},{block3_dist},{tiered_seqs[2]},{block4_dist},{tiered_seqs[3]},{total_dist}" 
                with open(output_file,"a+") as outfile:
                    outfile.write(output_line_str + "\n")
                #print(f"{output_line_str}")

if __name__ == "__main__":
    main()    