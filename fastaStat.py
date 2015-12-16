#!/usr/bin/env python

# commpute basic statistics of fasta files
# print length of all individual records
# return total length of the concatenated sequences
# return base frequency
# return concatenated sequence

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2013
# ---------------------------------------------------------------------------

from Bio import SeqIO
from optparse import OptionParser
import sys


def baseCount(seq):
    dico={}
    for base in seq:
        if base in dico.keys():
            dico[base]+=1
        else:
            dico[base]=1
    return dico


if __name__ == '__main__':
    parser = OptionParser()

    parser.add_option("-f", "--file",dest="filename",action="store", type="string", help="input fasta file", metavar="FILE")
    parser.add_option("-c", "--contigs",dest="contigs",action="store_true",default=False, help="print length of all individual records")
    parser.add_option("-l", "--length",dest="length",action="store_true",default=False, help="return total length of the concatenated sequences")
    parser.add_option("-o", "--concat",dest="concat",action="store_true",default=False, help="return concatenated sequence")
    parser.add_option("-p", "--percent",dest="percent",action="store_true",default=False, help="return base frequency")

    (options, args) = parser.parse_args()

    handle = open(options.filename)

    concatSeq=""
    concatLen=0
    for seq_record in SeqIO.parse(handle, "fasta"):
            #print seq_record
            if options.contigs:
                sys.stdout.write(seq_record.id +"\t" + str( len(seq_record.seq)) + "\n")
            concatLen+=len(seq_record.seq)
            if options.percent or options.concat:
                    concatSeq += seq_record.seq
    handle.close()

    if options.length:
            sys.stdout.write(str(concatLen) + "\n")

    if options.concat:
            sys.stdout.write(concatSeq + "\n")

    if options.percent:
            basecount=baseCount(concatSeq)
            for base in basecount.keys():
                    sys.stdout.write(base + "\t" + str(basecount[base]) + "\n")
