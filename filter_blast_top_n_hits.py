#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def filter_fasta(blast_files, max_n_hits=100):
    for one_blast_file in blast_files:
        keep = []

        with open(one_blast_file, 'r') as f:
            previous_locus = ''
            locus_count = 0
            for n, row in enumerate(f):
                data = row.rstrip().split('\t')
                if n == 0:
                    keep.append(data)
                    locus_count+=1
                    previous_locus = data[0]
                else:
                    if data[0] == previous_locus:
                        if locus_count < max_n_hits:
                            keep.append(data)
                            locus_count+=1
                        else:
                            locus_count+=1
                            continue
                    else:
                        # print first hits new locus
                        previous_locus = data[0]
                        keep.append(data)
                        locus_count=1
        out_name = one_blast_file.split('.')[0] + '_filtered_top_%s.tab' % max_n_hits
        with open(out_name, 'w') as f:
            for row in keep:
                f.write('\t'.join(row)+'\n')





if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast_tab', type=str, help="input blast tab file", nargs='+')

    args = parser.parse_args()

    filter_fasta(args.input_blast_tab, 100)

