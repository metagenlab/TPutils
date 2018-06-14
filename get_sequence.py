#! /usr/bin/env python

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqRecord import SeqRecord
import os
from time import time
from StringIO import StringIO

def get_sequence(file_name, file_format, start=False, stop=False, contig=False, protein_id=False, flanking_region_size_bp=0, locus_tag = False, revcomp=False,translate=False):

    import sys
    #print "File:", file_name
    #print "Format:", file_format
    #print "Protein ID:", protein_id

    contig_name2record =  SeqIO.to_dict(SeqIO.parse(file_name, file_format))

    # account for python count system
    start = start-1
    #print 'start, end', start, stop

    tag = False
    if not contig:
        #print 'AAA'
        for record_name in contig_name2record:
            #print record_name
            record = contig_name2record[record_name]
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                else:
                    
                    if protein_id in feature.qualifiers["protein_id"]:
                        tag = True
                    elif protein_id in feature.qualifiers["locus_tag"]:
                        tag = True
                if tag:
                    out_handle = StringIO()
                    start = int(feature.location.start-1) - flanking_region_size_bp
                    end = int(feature.location.end-1) + flanking_region_size_bp
                    sub_record = record[start:end]
                    if revcomp:
                        sub_record.seq = sub_record.seq.reverse_complement()
                    print len(sub_record)
                    SeqIO.write(sub_record, out_handle, "fasta")
                    new_gb_output = out_handle.getvalue() 
                    print new_gb_output
                    quit()
    else:

        if revcomp:
            if start is not False and stop is not False:
                if translate:
                    print contig_name2record[contig][start-flanking_region_size_bp:stop+flanking_region_size_bp].reverse_complement().seq.translate()
                    #sys.stdout.write(contig_name2record[contig][start-flanking_region_size_bp:stop+flanking_region_size_bp].reverse_complement().seq.translate())
                else:
                    print "start-flanking_region_size_bp:stop+flanking_region_size_bp: %s:%s" % (start-flanking_region_size_bp,stop+flanking_region_size_bp)
                    SeqIO.write(contig_name2record[contig][start-flanking_region_size_bp:stop+flanking_region_size_bp].reverse_complement(), sys.stdout ,"fasta")
            else:
                if translate:
                    SeqIO.write(contig_name2record[contig].reverse_complement().translate(), sys.stdout ,"fasta")
                else:
                    SeqIO.write(contig_name2record[contig].reverse_complement(), sys.stdout ,"fasta")
        else:
            if start is not False and stop is not False:
		print "start-flanking_region_size_bp:stop+flanking_region_size_bp: %s:%s" % (start-flanking_region_size_bp,stop+flanking_region_size_bp)
                SeqIO.write(contig_name2record[contig][start-flanking_region_size_bp:stop+flanking_region_size_bp], sys.stdout ,"fasta")
            else:
                SeqIO.write(contig_name2record[contig], sys.stdout ,"fasta")
    if protein_id and file_format=="fasta":
        records = []
        for one_prot in protein_id:
            record = contig_name2record[one_prot]
            records.append(record)
        out_handle = StringIO()
        SeqIO.write(records, out_handle, "fasta")
        new_gb_output = out_handle.getvalue()
        sys.stdout.write(new_gb_output)




if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",'--genbank_file',type=str,help="input genbank file", default=False)
    parser.add_argument("-f",'--fasta_file',type=str,help="input FASTA file", default=False)
    parser.add_argument("-s",'--start',type=int,help="start", default=False)
    parser.add_argument("-e",'--end',type=int,help="end", default=False)
    parser.add_argument("-c",'--contig',type=str,help="contig_name", default=False)
    parser.add_argument("-p",'--protein_id',type=str,help="protein id", default=False, nargs='+')
    parser.add_argument("-r",'--region',type=int,help="extract flanking region (bp)", default=0)
    parser.add_argument("-v",'--revcomp',action="store_true",help="reverse complement", default=False)
    parser.add_argument("-t",'--translate',action="store_true",help="translate", default=False)

    args = parser.parse_args()

    if args.genbank_file:
            get_sequence(file_name = args.genbank_file, file_format="genbank", protein_id = args.protein_id, flanking_region_size_bp = args.region, revcomp=args.revcomp)
    if args.fasta_file:
            get_sequence(file_name = args.fasta_file, file_format="fasta", protein_id = args.protein_id, contig=args.contig, start=args.start, stop=args.end, revcomp=args.revcomp, flanking_region_size_bp=args.region, translate=args.translate)
            

