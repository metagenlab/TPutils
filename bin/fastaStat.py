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
    parser.add_option("-c", "--concat",dest="concat",action="store_true",default=False, help="concat records (need either -l or -p to be set)")
    parser.add_option("-l", "--length",dest="length",action="store_true",default=False, help="return length of the sequence(s)")
    parser.add_option("-o", "--concatseq",dest="concatseq",action="store_true",default=False, help="return concatenated sequence")
    parser.add_option("-p", "--percent",dest="percent",action="store_true",default=False, help="return base frequency")

    (options, args) = parser.parse_args()

    handle = open(options.filename)

    if options.concat:
        concatSeq=""
        concatLen=0
        for seq_record in SeqIO.parse(handle, "fasta"):
            #print seq_record
            #if options.contigs:
            #    sys.stdout.write(seq_record.id +"\t" + str( len(seq_record.seq)) + "\n")
            concatLen+=len(seq_record.seq)
            if (options.percent and options.concat) or options.concatseq:
                concatSeq += seq_record.seq
    

    if options.length and options.concat:
            sys.stdout.write(str(concatLen) + "\n")
    elif options.length:
        for seq_record in SeqIO.parse(handle, "fasta"):
            sys.stdout.write("%s\t%s\n" % (seq_record.name, len(seq_record.seq)))
    elif options.concatseq:
            sys.stdout.write(concatSeq + "\n")
    elif options.percent:
        if options.concat:
            basecount=baseCount(concatSeq)
            for base in basecount.keys():
                    sys.stdout.write(base + "\t" + str(basecount[base]) + "\n")
        else:
            rec = list(SeqIO.parse(handle, "fasta"))
            all_names = [seq_record.name for seq_record in rec]
            all_basecount = [baseCount(seq_record.seq) for seq_record in rec]
            sys.stdout.write('\t' + '\t'.join(all_names) + '\n')
            for base in ['A','T','G','C','N', '-']:
                sys.stdout.write("%s\t" % base)
                for i, seq_record in enumerate(rec):
                    #basecount=baseCount(seq_record.seq)
                    try:
                        sys.stdout.write("%s (%s %%)\t" % (all_basecount[i][base], round(float(all_basecount[i][base])/len(seq_record.seq)*100, 2)))
                    except:
                        sys.stdout.write("-")
                sys.stdout.write("\n")
    handle.close()
