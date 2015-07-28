#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to ffn
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------



def locus2EC(priam_file, count_locus=False):

    locus_tag2EC_dico = {}
    with open(priam_file, "r") as f:
        lines = [i.rstrip() for i in f]
        for line in lines:
            if len(line) == 0:
                continue
            elif line[0] == '#':
                continue
            elif line[0] == '>':
                locus_tag = line.split(' ')[0][1:]
                if locus_tag not in locus_tag2EC_dico:
                    locus_tag2EC_dico[locus_tag] = []

            else:
                # valid EC
                locus_tag2EC_dico[locus_tag].append(line.split('\t'))

    if count_locus:
        i = 0
        for locus in locus_tag2EC_dico:
            if len(locus_tag2EC_dico[locus]) > 0:
                i+=1
        print i
    else:
        for locus in locus_tag2EC_dico:
            if len(locus_tag2EC_dico[locus]) > 0:
                for i in locus_tag2EC_dico[locus]:
                    print '%s\t%s' % (locus, i[0])



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_priam', type=str, help="input priam EC file")
    parser.add_argument("-c", '--count_locus', action='store_true', help="count locus with minimum 1 EC")

    args = parser.parse_args()


    locus2EC(args.input_priam, args.count_locus)
