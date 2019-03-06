#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# check if all records in a file have identical sources
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2017
# ---------------------------------------------------------------------------

def check_plasmids(records, outname):
    from Bio.SeqUtils import GC

    source_list = []
    plasmid_count = 0
    chromosome_count = 0
    check=False
    for n, record in enumerate(records):

        source = record.annotations['source']
        if 'plasmid' in record.features[0].qualifiers:
            plasmid_count+=1
            plasmid_name = record.features[0].qualifiers['plasmid'][0]
            if 'plasmid' not in record.description:
                check=True
                print 'adding plasmid...'
                record.description+=' plasmid'
            if plasmid_name not in record.description:
                check=True
                print 'adding plasmid name... %s -- %s' % (record.name, plasmid_name)
                record.description+=' %s' % plasmid_name
                print record.description
        else:
            if chromosome_count == 0:
                chromosome_count+=1
            else:
                chromosome_count+=1
                print 'multiple chromosomes!'
                print record.description
    if check:
        with open(outname, 'w') as f:
            SeqIO.write(records, f, 'genbank')

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
        out_name = genbank.split('.')[0] + '_plasmids.gbk'
        seq_records = list(SeqIO.parse(input_handle, "genbank"))
        all_records+=seq_records



    check_plasmids(all_records, out_name)
