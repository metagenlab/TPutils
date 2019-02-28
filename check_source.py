#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# check if all records in a file have identical sources
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2017
# ---------------------------------------------------------------------------

def check_sources(records):
    from Bio.SeqUtils import GC

    source_list = []
    for n, record in enumerate(records):

        source = record.annotations['source']
        if n == 0:
            source_list.append(source)
        else:
            if source in source_list:
                continue
            else:
                print 'problem:', source, record.id

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", nargs='+')
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    all_records = []
    for genbank in args.input_gbk:
        input_handle = open(genbank, "rU")
        seq_records = list(SeqIO.parse(input_handle, "genbank"))
        all_records+=seq_records



    check_sources(all_records)
