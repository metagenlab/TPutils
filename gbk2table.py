#! /usr/bin/env python

# produce one ptt / record present in the genbank file

from Bio import SeqIO
from optparse import OptionParser
import re
from Bio import SeqUtils

def find_index(pattern, seq):
  """Return first item in sequence where f(item) == True."""
  for item in seq:
    if re.match(pattern,item): 
      return seq.index(item)


class Feature:
    def __init__(self):
        
        self.contig = "-" 
        self.type = "-" 
        self.start = "-" 
        self.stop = "-" 
        self.length = "-" 
        self.GC = "-" 
        self.sequence = "-" 
        self.strand = "-" 
        self.gene = "-" 
        self.function = "-" 
        self.gi = "-"
        self.locus = "-"
        self.protein_id = "-"
        self.inference = "-"
        self.translation = "-"
        self.seq = "-"
        self.id = "-"
        self.product = "-"

class Record:
    def __init__(self, record):
        self.seq = record.seq
        #print "length", len(self.seq)
        self.contig =  record.name
        #print "name", self.contig
        self.features = self.get_one_record_features(record)

    def get_one_record_features(self, one_record):
            
        feature_list = [None] * len(one_record.features)
        for i in range(0,len(one_record.features)):
            #print one_record.features[i]
            if one_record.features[i].type == "misc_feature":
                continue
            feature_list[i] = Feature() 
            #print  one_record.features[i]
            feature_list[i].type = one_record.features[i].type
            feature_list[i].contig = one_record.name
            feature_list[i].start = one_record.features[i].location.start
            feature_list[i].stop = one_record.features[i].location.end
            feature_list[i].length = len(one_record.features[i].location)
            feature_list[i].strand = one_record.features[i].strand
            try:
                feature_list[i].gene = one_record.features[i].qualifiers['gene'][0]
            except:
                pass
            try:
                gi_position=find_index("GO*",one_record.features[i].qualifiers['db_xref'])
                feature_list[i].gi = one_record.features[i].qualifiers['db_xref'][gi_position][3:]
            except:
                pass
         
            #geneID= one_record.features[i].qualifiers['db_xref'][1][7:]
            try:
                feature_list[i].locus = one_record.features[i].qualifiers['locus_tag'][0]
            except:
                pass
            try:
                feature_list[i].protein_id = one_record.features[i].qualifiers['protein_id'][0]
            except:
                pass
            try:
                feature_list[i].product = one_record.features[i].qualifiers['product'][0]
            except:
                pass
            try:
              
              feature_list[i].inference = one_record.features[i].qualifiers['inference']
            except:
              pass
            try:
                feature_list[i].translation = one_record.features[i].qualifiers['translation'][0]
            except:
                pass
            feature_list[i].seq = one_record.features[i].extract(self.seq)
          
            feature_list[i].GC = SeqUtils.GC(feature_list[i].seq)
            

        return feature_list




if __name__ == '__main__':
    
    
    parser = OptionParser()

    parser.add_option("-i", "--input",dest="input_file",action="store", type="string", help="genbank file", metavar="FILE")
    parser.add_option("-o", "--output",dest="output_file",action="store",default="record", type="string", help="genbank file", metavar="FILE")
    (options, args) = parser.parse_args()



    gb_file = options.input_file
    # get list of all records present in the gbk file
    record_list=[]
   
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        #print gb_record
        record_list.append(Record(gb_record))
    

    # format output tab delimited table
    print "contig\ttype\tstart\tstop\tlength\tGC\tstrand\tgene\tfunction\tinference\tgi\tlocus\ttranslation\tsequence"
    for record in record_list:
      for feature in record.features:
        if feature is None:
          continue
        if feature.type == "source":
          pass
          #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t \t " % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus)
        else:
          print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (feature.contig, feature.type, feature.start, feature.stop, feature.length, feature.GC, feature.strand, feature.gene, feature.product, feature.inference, feature.gi, feature.locus, feature.translation, feature.seq)
        
