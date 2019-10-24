#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    if not args.outname:
        outname = args.input_gbk.split(".")[0]+".embl"
    else:
        outname = args.outname

    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    print seq_records
    SeqIO.write(seq_records, outname, "embl")
