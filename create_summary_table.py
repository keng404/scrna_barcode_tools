import pandas as pd
from scipy import stats
import argparse
import os
from datetime import datetime as dt
##########
def logging_statement(string_to_print):
    date_time_obj = dt.now()
    timestamp_str = date_time_obj.strftime("%Y/%b/%d %H:%M:%S:%f")
    #############
    final_str = f"[ {timestamp_str} ] {string_to_print}"
    return print(f"{final_str}")
## uses pandas to compute descriptive statistics and scipy to compute normal test p-value
def compute_stats(data_arr):
    df = pd.DataFrame(data_arr)
    stats_dict = df.describe().to_dict()
    stats_dict = stats_dict[0]
    ## coefficient of variation
    stats_dict['cv'] = stats_dict['std']/stats_dict['mean']
    ## kurtosis
    kurtosis = df.kurtosis().to_dict()
    kurtosis = kurtosis[0]
    stats_dict["kurtosis"] = kurtosis
    ## skewness
    skew = 0
    try:
        skew = df.skew().to_dict()
        skew = skew[0]
    except:
        skew = 0
    stats_dict["skewness"] = skew
    ## normal p-value
    try:
        normal_pvalue = stats.normaltest(df).pvalue[0]
    except:
        normal_pvalue = -1
    stats_dict["normal_p"] = normal_pvalue
    return stats_dict

def get_summary_stats(line_of_interest):
    output_line = []
    full_res_stats = compute_stats(line_of_interest)
    if len(full_res_stats) > 0:
        fields_of_interest = ["min","max","mean","25%","50%","75%","std","cv"] 
        for f in fields_of_interest:
            if f in full_res_stats.keys():
                output_line.append(str(full_res_stats[f]))
            else:
                output_line.append("UNKNOWN")
    return output_line


def get_table_size(barcodes_file):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            lines_processed += 1
    return lines_processed

def get_barcode_line(barcodes_file,idx):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            if lines_processed == idx:
                line_cleaned = line.strip("\n")
                line_cleaned = line_cleaned.split("\t")
                all_metrics = []
                for idx,val in enumerate(line_cleaned):
                    if idx > 0:
                        all_metrics.append(int(val))
                return  all_metrics
            lines_processed += 1
        

def get_barcode(barcodes_file,idx):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            if lines_processed == idx:
                line_cleaned = line.strip("\n")
                line_cleaned = line_cleaned.split("\t")
                return  str(line_cleaned[0]) 
            lines_processed += 1
                 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_distance_file',default=None, required=True, type=str, help="R1_barcode_hamming_distance.tsv file from DRAGEN")
    parser.add_argument('--output_file', default=None, type=str, help="output file path")
    args, extras = parser.parse_known_args()
####################
    output_file = args.output_file
    barcode_distance_file = args.barcode_distance_file
    if output_file is None:
        output_dir = os.path.dirname(barcode_distance_file)
        initial_file_name_split = os.path.basename(barcode_distance_file).split(".")
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["R1_barcode_hamming_distance_summary",initial_file_name_split[-1]]
        output_filename = ".".join(output_filename_components)
        output_file = output_dir + "/"+ output_filename
        logging_statement(f"No output file provided. Will write out hamming distance summary table to {output_file}")
    ###############
    num_total_lines = get_table_size(barcode_distance_file)
    lines_processed = 0
    fields_of_interest = ["min","max","mean","25%","50%","75%","std","cv"] 
    ### write header
    header_line = ["read1_barcode"]
    for w in fields_of_interest:
        header_line.append(w)
    header_line_str = "\t".join(header_line)
    with open(output_file,"w") as outfile:
        outfile.write(header_line_str + "\n")

    while lines_processed < num_total_lines:
        lines_processed += 1
        if lines_processed % 10000 == 0:
            logging_statement(f"Computed Hamming Distance for {lines_processed} lines")
        my_line = get_barcode_line(barcode_distance_file,lines_processed)
        my_line_stats = get_summary_stats(my_line)
        current_barcode = get_barcode(barcode_distance_file,lines_processed)
        full_line_components = [current_barcode]
        for statistic in my_line_stats:
            full_line_components.append(str(statistic))
        full_line_str = "\t".join(full_line_components)
        # update output file line by line
        with open(output_file,"a+") as outfile:
            outfile.write(full_line_str + "\n")
    logging_statement(f"Finished writing {output_file}")
if __name__ == "__main__":
    main()    