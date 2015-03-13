#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# input: fasta
# input: sequence id(s)
# return seq record matching sequence id(s) in fasta
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# ---------------------------------------------------------------------------


def locate_sequence(fasta_handle, seq_ids):
    import sys
    for record in SeqIO.parse(fasta_handle, "fasta"):
        if record.id in seq_ids:
            SeqIO.write(record, sys.stdout, "fasta")
            #sys.stdout.write(record.seq)


if __name__ == '__main__':
    
    from Bio import SeqIO
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--input_fasta',type=str,help="input file ")
    parser.add_argument("-i",'--sequence_ids',type=str,help="sequence id(s)", nargs='+')
    args = parser.parse_args()

    handle = open(args.input_fasta, "rU")

    locate_sequence(handle, args.sequence_ids)