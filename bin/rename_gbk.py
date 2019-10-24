#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert rename gbk with accession
# Date: 2015
# ---------------------------------------------------------------------------



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))

    if not args.outname:
        outname = seq_records[0].name+".gbk"

    else:
        outname = args.outname

    out_handle = open(outname, "w")
    SeqIO.write(seq_records, out_handle, "genbank")
