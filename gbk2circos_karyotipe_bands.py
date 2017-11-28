#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

def gbk2bands(seq_records, chromosome_id, add=200):
    
    start = 0
    end = 0
    colc = 0
    for n, record in enumerate(seq_records):
        if n == 0:
            pass
        else:
            start += add
        end = start + len(record.seq)
        if colc == 4:
            colc=0
        print 'band %s band%s band%s %s %s greys-3-seq-%s' % (chromosome_id, n+1, n+1, start, end, colc)
        start = end
        colc+=1
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()
    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "fasta"))

    gbk2bands(seq_records, "KpGe")
