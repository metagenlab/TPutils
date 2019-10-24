#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def get_code2locus_dico(seq_file):
    code2locus = {}
    with open(seq_file, 'r') as f:
        for line in f:
            data = line.rstrip().split(': ')
            code2locus[data[0]] = data[1].split(' ')[0]
    return code2locus

def rename_blast(blast_file, seq_file):

    code2locus = get_code2locus_dico(seq_file)

    o = open("renamed_blast.tab", 'w') 
    with open(blast_file, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            data[0] = code2locus[data[0]]
            data[1] = code2locus[data[1]]
            o.write('\t'.join(data) + '\n')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", '--merged_orthofinder_blast_files', type=str, help="merged orthofinder blast file")
    parser.add_argument("-i", '--SequenceIDs', type=str, help="SequenceIDs.txt")


    args = parser.parse_args()
    rename_blast(args.merged_orthofinder_blast_files,args.SequenceIDs)
