#! /usr/bin/env python



def fasta2tab(fasta_file):
    from Bio import SeqIO
    import re

    try:
        file_description = re.findall('G[A-Za-z0-9_]+\.[0-9]+',fasta_file)[0]
    except:
        file_description = fasta_file.split('.')[0]

    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            print "%s\t%s" % (file_description, record.name)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input fasta file")

    args = parser.parse_args()
    fasta2tab(args.input)