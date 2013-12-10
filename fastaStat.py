#!/usr/bin/env python


from Bio import SeqIO
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--file",dest="filename",action="store", type="string", help="input fasta file", metavar="FILE")
parser.add_option("-c", "--contigs",dest="contigs",action="store_true",default=False, help="print length of all individual contigs")
parser.add_option("-l", "--length",dest="length",action="store_true",default=False, help="return total length of the concatenated sequences")
parser.add_option("-o", "--concat",dest="concat",action="store_true",default=False, help="return concatenated sequence")
parser.add_option("-p", "--percent",dest="percent",action="store_true",default=False, help="return % of each base")

(options, args) = parser.parse_args()




def baseCount(seq):
    dico={}
    for base in seq:
        if base in dico.keys():
            dico[base]+=1
        else:
            dico[base]=1
    return dico


handle = open(options.filename)

concatSeq=""
concatLen=0
for seq_record in SeqIO.parse(handle, "fasta"):
        #print seq_record
        if options.contigs:
                print len(seq_record.seq)
        concatLen+=len(seq_record.seq)
        if options.percent or  options.concat:
                concatSeq+=seq_record.seq
handle.close()

if options.length:
        print concatLen

if options.concat:
        print concatSeq

if options.percent:
        basecount=baseCount(concatSeq)
        for base in basecount.keys():
                print base, basecount[base]
