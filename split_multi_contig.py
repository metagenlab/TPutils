#! /usr/bin/python

# produce one fasta file / record

from Bio import SeqIO
from optparse import OptionParser
import re


parser = OptionParser()

parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="fasta file", metavar="FILE")

(options, args) = parser.parse_args()










for seq_record in SeqIO.parse(options.input_file, "fasta"):
    locus= seq_record.id.split("|")[3]
    out_name= locus.split(".")[0]+".fna"
    print locus
    #print repr(seq_record.seq)
    print len(seq_record)
    SeqIO.write(seq_record, out_name, "fasta")
