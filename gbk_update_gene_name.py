#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def gbk_update(gbk_file, annot_table):
    from Bio import SeqIO

    out = gbk_file.split('.')[0] + '_annot.gbk'

    locus2gene = {}
    locus2product = {}
    with open(annot_table, 'r') as f:
        for line in f:
            data = line.rstrip().split('\t')
            if len(data)>2:
                locus2product[data[0]] = data[2]
            locus2gene[data[0]] = data[1]
    with open(gbk_file, 'r') as f:
        records = [i for i in SeqIO.parse(f, 'genbank')]

        for record in records:
            for feature in record.features:
                if 'locus_tag' in feature.qualifiers:
                    if feature.qualifiers['locus_tag'][0] in locus2gene.keys():
                        feature.qualifiers['gene'] = locus2gene[feature.qualifiers['locus_tag'][0]]
                    if feature.qualifiers['locus_tag'][0] in locus2product.keys():
                        feature.qualifiers['product'] = locus2product[feature.qualifiers['locus_tag'][0]]
    with open(out, 'w') as f:
        SeqIO.write(records, f, 'genbank')

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-a", '--annot_table', help="annot table", default=False)

    args = parser.parse_args()

    gbk_update(args.input_gbk,
               args.annot_table)