#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# input: fasta
# input: sequence id(s)
# return seq record matching sequence id(s) in fasta
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# ---------------------------------------------------------------------------


def locate_sequence(fasta_handle, id):
    import sys
    for record in SeqIO.parse(handle, "fasta"):
        #print record.id
        if record.id in args.sequence_ids:
            sys.stdout.write(">" + record.description + "\n")
            sys.stdout.write(record.seq + "\n")


if __name__ == '__main__':
    
    from Bio import SeqIO
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--input_fasta',type=str,help="input file ")
    parser.add_argument("-i",'--sequence_ids',type=str,help="sequence id(s)", nargs='+')
    args = parser.parse_args()


    #seq = SeqIO.read(file(ars.input_fasta), "fasta").seq
    handle = open(args.input_fasta, "rU")

