#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk',type=str,help="gbk file")
    parser.add_argument("-n", '--fna',type=str,help="fna file")
    
    args = parser.parse_args()
    genbank = [i for i in SeqIO.parse(args.gbk, "genbank")]

    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna, generic_protein
    fna = [i for i in SeqIO.parse(args.fna, "fasta")]
    print fna[0]
    for record_gbk, record_fna in zip(genbank, fna):

        record_gbk.seq = Seq(str(record_fna.seq), generic_dna)

    output_handle = open("%s_seq.gbk" % args.gbk.split('.')[0], 'w')
    SeqIO.write(genbank,output_handle, "genbank")
