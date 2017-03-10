#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def split_gbk(seq_records, outname, format = False):
    import re
    from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
    output_handle = open(outname, "w")

    merged_record = ''
    fasta_record = False
    for i, record in enumerate(seq_records):
        print i

        for feature in record.features:
            if feature.type == "fasta_record":
                fasta_record = True
                merged_record+=record[feature.location.start:feature.location.end]
                merged_record += "N" * 200
                my_start_pos = ExactPosition(len(merged_record)-200)
                my_end_pos = ExactPosition(len(merged_record))
                my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
                my_feature = SeqFeature(my_feature_location, type="assembly_gap")
                merged_record.features.append(my_feature)
               
            elif feature.type == 'source' and fasta_record == False:
                merged_record+=record[feature.location.start:feature.location.end]
                merged_record += "N" * 200
                my_start_pos = ExactPosition(len(merged_record)-200)
                my_end_pos = ExactPosition(len(merged_record))
                my_feature_location = FeatureLocation(my_start_pos,my_end_pos)
                my_feature = SeqFeature(my_feature_location, type="assembly_gap")
                merged_record.features.append(my_feature)
            

    to_remove = []
    for n, feature in enumerate(merged_record.features):
        if feature.type == 'source' or feature.type == "fasta_record":
           to_remove.append(n)
           

    for index in sorted(to_remove, reverse=True):
        if index != 0:
            #print index
            del merged_record.features[index]

    merged_record.id = seq_records[0].annotations["accessions"][-1]
    try:
        merged_record.description = "%s (merged contigs)" % seq_records[0].annotations["organism"]
    except:
        merged_record.description = 'Unkown bacteria'
    merged_record.annotations = seq_records[0].annotations
    merged_record.name = seq_records[0].annotations["accessions"][-1]
    return merged_record[0:-200]
    

        
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-f", '--format', action='store_true', help="format header: >locus description", default=False)

    args = parser.parse_args()
    


    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    if not args.outname:
        # use input file to rename the file
        #outname = args.input_gbk.split(".")[0]+".faa"

        # use record id to rename the file, remove version number using split
        outname = seq_records[0].id.split('.')[0] + ".faa"
    else:
        outname = args.outname


    merged_record = split_gbk(seq_records, outname, args.format)
    SeqIO.write(merged_record, sys.stdout, "genbank")
    
