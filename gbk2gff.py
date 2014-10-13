#! /usr/bin/env python

"""
Genbank to gff conversion using the GFF python module
"""


from Bio import SeqIO
from BCBio import GFF
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="genbank file")
parser.add_argument("-o", "--outfile", help="output name", default = False)
 
args = parser.parse_args()

if not args.outfile:
    outname = args.input.split(".")[0]+".gff"
else:
    outname = args.outfile
    
in_handle = open(args.input)
out_handle = open(outname, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()
