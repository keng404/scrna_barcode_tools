import argparse
import os
from datetime import datetime as dt
import math
##########
def logging_statement(string_to_print):
    date_time_obj = dt.now()
    timestamp_str = date_time_obj.strftime("%Y/%b/%d %H:%M:%S:%f")
    #############
    final_str = f"[ {timestamp_str} ] {string_to_print}"
    return print(f"{final_str}")

def get_table_size(barcodes_file):
    lines_processed = 0
    with open(barcodes_file,"r") as f:
        for line in f:
            lines_processed += 1
    return lines_processed

def get_barcode_lines(barcodes_file,start_idx,end_idx):
    lines_processed = 0
    lines_to_keep = []
    with open(barcodes_file,"r") as f:
        for line in f:
            if lines_processed >= start_idx and lines_processed <= end_idx :
                line_cleaned = line.strip("\n")
                lines_to_keep.append(line_cleaned)
            lines_processed += 1
    return lines_to_keep

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--barcode_file',default=None, required=True, type=str, help="barcodes.txt file from FASTQ")
    parser.add_argument('--chunk_size', default=1000000, type=int, help="chunk_size")
    args, extras = parser.parse_known_args()
####################
    chunk_size = args.chunk_size
    barcode_file = args.barcode_file
    logging_statement(f"Getting number of lines in {barcode_file}")
    num_lines = get_table_size(barcode_file)
    files_to_generate = math.ceil(num_lines/chunk_size)
    chunks = list(range(1,files_to_generate + 1 ))
    logging_statement(f"Creating {len(chunks)} chunks")
    for i in chunks:
        output_dir = os.path.dirname(barcode_file)
        initial_file_name_split = os.path.basename(barcode_file).split(".")
        output_filename_components = initial_file_name_split[0:len(initial_file_name_split)-1] + ["chunk",str(i),initial_file_name_split[-1]]
        output_filename = ".".join(output_filename_components)
        output_file = output_dir + "/"+ output_filename
        start_idx = (i-1) * chunk_size
        end_idx = (i * chunk_size) -1
        if end_idx > num_lines:
            end_idx = num_lines
        lines_to_write = get_barcode_lines(barcode_file,start_idx,end_idx)
        logging_statement(f"Writing out to {output_file}")
        with open(output_file,"w") as outfile:
            outfile.write("\n".join(lines_to_write))
if __name__ == "__main__":
    main()    