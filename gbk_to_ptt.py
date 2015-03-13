#! /usr/bin/env python

# produce one ptt / record present in the genbank file

# produce basic ppt file from genbank files
# TODO not all ptt fields supported, should be rewritten
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# ---------------------------------------------------------------------------

from Bio import SeqIO
from optparse import OptionParser
import re


parser = OptionParser()

parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="genbank file", metavar="FILE")
parser.add_option("-o", "--output",dest="output_file",action="store",default="record", type="string", help="genbank file", metavar="FILE")
(options, args) = parser.parse_args()



gb_file = options.input_file
record_list=[]
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
    # now do something with the record
    record_list.append(gb_record)
    #print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
    #print repr(gb_record.seq)

print "this gbk file contain %s records" % len(record_list)

def find_index(pattern, seq):
  """Return first item in sequence where f(item) == True."""
  for item in seq:
    if re.match(pattern,item): 
      return seq.index(item)




def format_ptt(one_record):
    out_name= one_record.id.split(".")[0]+".ptt"
    f = open(out_name,'w')
    # write header
    f.write(one_record.description+"\n")
    f.write(str(len(one_record.features))+" proteins\n")
    f.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
    #print "Location\t Strand\t gi\t geneID\t proteinID\t locus\t proteinID\t product"
    for i in range(0,len(one_record.features)):
        if one_record.features[i].type == "CDS":
            #print one_record.features[i]

            location= str(one_record.features[i].location.start)+".."+str(one_record.features[i].location.end)
            location=re.sub(">|<","",location)
            if one_record.features[i].strand ==1:
                strand="+"
            else:
                strand="-"
            #print one_record.features[i].qualifiers
            length= "-"
            COG="-"
            synonym="-"
            code="-"
            try:
                gene=one_record.features[i].qualifiers['gene'][0]
            except:
                gene="-"
            try:
                gi_position=find_index("GO*",one_record.features[i].qualifiers['db_xref'])
                gi= one_record.features[i].qualifiers['db_xref'][gi_position][3:]
            except:
                gi="-"
            
            #geneID= one_record.features[i].qualifiers['db_xref'][1][7:]
            locus= one_record.features[i].qualifiers['locus_tag'][0]
            try:
                proteinID= one_record.features[i].qualifiers['protein_id'][0]
            except:
                proteinID="-"
            product= one_record.features[i].qualifiers['product'][0]
            line=location+"\t"+strand+"\t"+length+"\t"+gi+"\t"+gene+"\t"+synonym+"\t"+code+"\t"+COG+"\t"+product+"\n"
            #print line
            f.write(line)


for i in range(0,len(record_list)):
    print "record",i
    one_record=record_list[i]
    format_ptt(one_record)
    
