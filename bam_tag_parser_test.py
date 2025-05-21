import sys
import re

def break_barcode_blocks(barcode):
    s1 = "".join(barcode[0:8])
    s2 = "".join(barcode[8:14])
    s3 = "".join(barcode[14:20])
    s4 = "".join(barcode[20:28])
    barcode_components = [s1,s2,s3,s4]
    return barcode_components

def get_tags(line_parsed):
    parsed_dict = dict()
    tags_of_interest = ["XB","CR"]
    for i in line_parsed:
        i_split = i.split(":")
        if i_split[0] in tags_of_interest:
            parsed_dict[i_split[0]] = i_split[2]
    return parsed_dict

tags_of_interest = ["XB","CR"]

bam_tag_file = "XRQ-DM-032725.bam_tags.csv"
with open(bam_tag_file,"r", encoding='utf-8') as f:
    for line in f:
        line_cleaned = line.strip()
        line_split = line_cleaned.split(",")
        #print(f"{line_split}")
        line_dict = get_tags(line_split)
        
        if "CR" in list(line_dict.keys()):
            sys.stderr.write(f"Found valid barcode : {line_split[0]}\t{line_dict}\n")
            sys.exit(0)
        else:
            sys.stderr.write(f"Found filtered barcode : {line_split[0]}\t{line_dict}\n")
            print(f"{break_barcode_blocks(line_dict["XB"])}")
            