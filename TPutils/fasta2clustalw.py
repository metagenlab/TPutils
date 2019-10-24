#! /usr/bin/env python


# Genbank to gff conversion using the GFF python module
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


from Bio import AlignIO
from BCBio import GFF
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input fasta file")
parser.add_argument("-o", "--outfile", help="output name", default = False)
 
args = parser.parse_args()

if not args.outfile:
    outname = args.input.split(".")[0]+".clw"
else:
    outname = args.outfile
    
in_handle = open(args.input)
out_handle = open(outname, "w")
 
AlignIO.write(AlignIO.parse(in_handle, "fasta"), out_handle, 'clustal')
 
in_handle.close()
out_handle.close()
