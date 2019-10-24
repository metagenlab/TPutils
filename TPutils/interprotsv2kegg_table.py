#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def interpro2kegg(interpro_file):
    import re
    pattern = re.compile(".*KEGG:*")
    locus2kegg = {}
    with open(interpro_file) as f:
        data = [line.rstrip().split('\t') for line in f]
        for line in data:
            if pattern.match(line[-1]):
                locus = line[0]
                kegg = line[-1].split('|')
                for one_entry in kegg:
                    details =  one_entry.split(':')
                    if details[0] == 'KEGG':
                        ko = 'K'+details[1].split('+')[0][1:]
                        ec = details[1].split('+')[1]
                        #print locus, ko
                        if locus not in locus2kegg:
                            locus2kegg[locus] = [ko]
                        else:
                            if ko not in locus2kegg[locus]:
                                locus2kegg[locus].append(ko)
    return locus2kegg




if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_interpro', type=str, help="input interpro tsv file")

    args = parser.parse_args()
    locus2kegg = interpro2kegg(args.input_interpro)
    for locus in locus2kegg:
        for ko in locus2kegg[locus]:
            print locus + '\t' + ko