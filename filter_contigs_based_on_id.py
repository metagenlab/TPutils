#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def filter_contigs(seq_records, filter_list, outname, format = False):
    import re
    output_handle = open(outname, "w")

    keep = []
    for record in seq_records:
        if record.name not in filter_list:
            keep.append(record)
    SeqIO.write(keep, output_handle, 'fasta')



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-f", '--filter', type=str, help="filter list (tablulated file with contig ids)", default=False)
    parser.add_argument("-o", '--outname', type=str, help="outname", default=False)
    args = parser.parse_args()



    input_handle = open(args.input_fasta, "rU")
    seq_records = list(SeqIO.parse(input_handle, "fasta"))
    if not args.outname:
        # use input file to rename the file
        #outname = args.input_gbk.split(".")[0]+".faa"

        # use record id to rename the file, remove version number using split
        outname = seq_records[0].id.split('.')[0] + "_filter.faa"
    else:
        outname = args.outname


    with open(args.filter) as f:
        rm_list = [line.rstrip() for line in f]

    filter_contigs(seq_records, rm_list, outname)
