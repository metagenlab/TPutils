#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys

def merge_gbk(gbk_records, filter_size=0, gi=False):

    from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
    from Bio.SeqRecord import SeqRecord
    n=0
    if len(gbk_records) == 1:
        merged_rec = gbk_records[0]
    else:
        for i, rec in enumerate(gbk_records):
            # remove source feature
            rec.features.pop(0)
            # filter small contigs
            if len(rec) > filter_size:
                if n == 0:
                    n+=1
                    merged_rec = rec
                else:
                    merged_rec+=rec
                # you could insert a spacer if needed
                # do not add spacer after the last contig
                if i != len(gbk_records)-1:
                    merged_rec += "N" * 200

                    my_start_pos = ExactPosition(len(merged_rec)-200)
                    my_end_pos = ExactPosition(len(merged_rec))

                    my_feature_location = FeatureLocation(my_start_pos,my_end_pos)

                    my_feature = SeqFeature(my_feature_location, type="assembly_gap")

                    merged_rec.features.append(my_feature)
    print "gbk_records[0].annotations", gbk_records[0].annotations
    try:
        merged_rec.id = gbk_records[0].annotations["accessions"][-1]
    except KeyError:
        merged_rec.id = gbk_records[0].id
        
    if gi:
        merged_rec.annotations["gi"] = gi

    merged_rec.description = "%s (merged contigs)" % gbk_records[0].annotations["organism"]
    merged_rec.annotations = gbk_records[0].annotations
    try:
        merged_rec.name = gbk_records[0].annotations["accessions"][-1]
    except KeyError:
        merged_rec.name = gbk_records[0].id
    return merged_rec


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file(s)", nargs='+')
    parser.add_argument("-g", '--gi', default=False, type=str, help="gi (optional)", nargs='+')
    args = parser.parse_args()

    infile = "multi-rec.gbk"

    records = []
    for gbk in args.input_gbk:
        records+=SeqIO.parse(open(gbk,"r"), "genbank")
    merged_record = merge_gbk(records, args.gi)
    SeqIO.write(merged_record, sys.stdout, "genbank")
