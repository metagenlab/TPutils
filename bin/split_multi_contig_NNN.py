#! /usr/bin/env python

# produce one fasta file / record

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = OptionParser()

parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="fasta file", metavar="FILE")

(options, args) = parser.parse_args()

print options.input_file


out_name = options.input_file.split('.')[0] + '_NNN.fna'

new = []
for seq_record in SeqIO.parse(options.input_file, "fasta"):
    

    contigs = str(seq_record.seq).split('N'*200)

    print 'n contigs', len(contigs)
    
    for i, contig in enumerate(contigs):
        new.append(SeqRecord(Seq(contig,
                           seq_record.seq.alphabet),
                           id='CONTIG_%s' % i, name='CONTIG_%s',
                           description=seq_record.description))
    
SeqIO.write(new, out_name, "fasta")
