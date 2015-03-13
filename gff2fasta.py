#! /usr/bin/env python


# TODO utile?
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO

def examine(gff_file, fasta_file):
    gff_handle = open(gff_file)
    
    fasta_handle = open(fasta_file)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))
    
    for rec in GFF.parse(gff_handle):
        #print rec.id
        for one_feature in rec.features:
            #print  one_feature.qualifiers.keys()
            locus = one_feature.qualifiers["Name"][0]
            product = one_feature.qualifiers["product"][0]
            seq =  one_feature.extract(fasta_dict[rec.id].seq)
            out = ">%s|%s|%s\n%s" % (rec.id, locus, product, seq)
            print out
            
    gff_handle.close()
    fasta_handle.close()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fasta',type=str,help="fasta file")
    parser.add_argument("-g",'--gff',type=str,help="gff file")

    
    args = parser.parse_args()
    examine(args.gff, args.fasta)
