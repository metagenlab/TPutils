#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys

def merge_gbk(gbk_records):

    from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

    merged_rec = ''


    for i, rec in enumerate(gbk_records):
        print rec
        merged_rec += rec
        # you could insert a spacer if needed
        if i != len(gbk_records):
            merged_rec += rec + ("N" * 200)

            my_start_pos = ExactPosition(len(merged_rec)-200)
            my_end_pos = ExactPosition(len(merged_rec))

            my_feature_location = FeatureLocation(my_start_pos,my_end_pos)

            my_feature = SeqFeature(my_feature_location, type="assembly_gap")

            merged_rec.features.append(my_feature)

    merged_rec.id = gbk_records[0].annotations["accessions"][-1]
    merged_rec.description = "Merged record from %s" % gbk_records[0].annotations["organism"]
    print merged_rec
    SeqIO.write(merged_rec, sys.stdout, "genbank")


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file(s)", nargs='+')

    args = parser.parse_args()

    infile = "multi-rec.gbk"

    records = []
    for gbk in args.input_gbk:
        records+=SeqIO.parse(open(gbk,"r"), "genbank")
    merge_gbk(records)