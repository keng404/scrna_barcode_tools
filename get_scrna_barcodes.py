#!/usr/bin/env python
import sys
import pysam
### adapted from https://github.com/lh3/biofast/blob/master/fqcnt/fqcnt_py7x_pysam.py
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: {} <in.fq.gz>".format(sys.argv[0]))
        sys.exit(0)
    scrna_barcode_str = "0_7+11_16+20_25+31_38"
    scrna_barcode_intervals = scrna_barcode_str.split("+")
    scrna_barcode_starts = [x.split("_")[0] for x in scrna_barcode_intervals ]
    scrna_barcode_ends = [x.split("_")[1] for x in scrna_barcode_intervals ]
    n, slen, qlen = 0, 0, 0
    with pysam.FastxFile(sys.argv[1]) as fastx:
        for record in fastx:
            sequence_split = [x for x in str(record.sequence)]
            barcode_str = ""
            for idx,val in enumerate(scrna_barcode_starts):
                barcode_subset = sequence_split[int(scrna_barcode_starts[idx]):(int(scrna_barcode_ends[idx])+1)]
                barcode_str += "".join(barcode_subset)
            n += 1
            slen += len(record.sequence)
            qlen += len(record.quality)
            print(barcode_str)
    sys.stderr.write("Processed {} Reads\t Number of Bases {}\t{}".format(n, slen, qlen))