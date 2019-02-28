#! /usr/bin/env python


# Genbank to gff conversion using the GFF python module
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


from Bio import SeqIO, SeqRecord, Seq
from BCBio import GFF
import argparse
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="genbank file")
parser.add_argument("-o", "--outfile", help="output name", default = False)
 
args = parser.parse_args()

if not args.outfile:
    outname = args.input.split(".")[0]+".gbk"
else:
    outname = args.outfile
    
in_handle = open(args.input)
out_handle = open(outname, "w")



#SeqIO.write(GFF.parse(in_handle, "gff"), out_handle)
#for rec in GFF.parse(in_handle):
gff_records = [i for i in GFF.parse(in_handle)]
for i,record in enumerate(gff_records):
    print i, record
    my_seq = Seq.Seq(str(record.seq), generic_dna)
    record.seq = my_seq
    SeqIO.write(record, out_handle, 'genbank')
    print len(record.seq)
#rec = SeqRecord.SeqRecord(seq, id=id, name=name, description=description)
#SeqIO.write(, out_handle, 'genbank')

#in_handle.close()
#out_handle.close()
