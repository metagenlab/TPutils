#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Contactenate contigs with 200 (default) NNN between each contig
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: septembre 2015
# ---------------------------------------------------------------------------

def merge_fasta(fasta_records):
    seq = ''
    for one_seq in fasta_records:
        seq+=one_seq.seq + 200*'N'
    return seq[0:-200]



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--input_fna', type=str, help="input fna file")
    parser.add_argument("-n", '--header', type=str, help="Header")



    args = parser.parse_args()
    fasta_records = list(SeqIO.parse(args.input_fna, "fasta"))
    merged_seq = merge_fasta(fasta_records)

    output_name = args.input_fna.split('.')[0] + '_concat.fna'
    if not args.header:
        args.header = args.input_fna.split('.')[0]
    
    with open(output_name, 'w') as f:
        f.write('>%s\n%s' % (args.header, merged_seq))
