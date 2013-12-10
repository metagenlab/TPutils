#!/usr/bin/env python

from Bio import SeqIO
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--file",dest="filename",action="store", type="string", help="input fasta file", metavar="FILE")
parser.add_option("-o", "--outfile",dest="outfile",action="store", type="string", help="output fasta file", metavar="FILE")

(options, args) = parser.parse_args()

handle = open(options.filename)

allrecords=list(SeqIO.parse(handle, "fasta"))
for i in range(0,len(allrecords)):
        #print seq_record

        if str(allrecords[i].seq[0:3])=="atg":
            seq=allrecords[i].seq
        elif str(allrecords[i].seq[1:4])=="atg":
            print "OK!"
            seq=allrecords[i].seq[1:-2]
        elif str(allrecords[i].seq[2:5])=="atg":
            seq=allrecords[i].seq[2:-1]
        else:
            seq="atg"+allrecords[i].seq[3:]
        allrecords[i].seq=seq

        

output_handle = open(options.outfile, "w")
SeqIO.write(allrecords, output_handle, "fasta")
