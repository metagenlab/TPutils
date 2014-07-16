#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-






if __name__ == '__main__':
    
    from Bio import SeqIO
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--input_fasta',type=str,help="input file ")
    parser.add_argument("-i",'--sequence_ids',type=str,help="sequence id(s)", nargs='+')
    args = parser.parse_args()


    #seq = SeqIO.read(file(ars.input_fasta), "fasta").seq
    handle = open(args.input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        #print record.id
        if record.id in args.sequence_ids:
            print ">", record.description
            print record.seq

