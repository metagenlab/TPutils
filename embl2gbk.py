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
    parser.add_argument("-i", '--input_embl', type=str, help="input embl file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    if not args.outname:
        outname = args.input_embl.split(".")[0]+".gbk"
    else:
        outname = args.outname

    input_handle = open(args.input_embl, "rU")
    seq_records = list(SeqIO.parse(input_handle, "embl"))
    SeqIO.write(seq_records, outname, "genbank")
