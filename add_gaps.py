#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys

def add_gaps(gbk_record, start_end_list):

    from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

    merged_rec = ''


    for start, end in start_end_list:
            #print start, end
            my_start_pos = ExactPosition(start)
            my_end_pos = ExactPosition(end)
            my_feature_location = FeatureLocation(my_start_pos, my_end_pos)
            my_feature = SeqFeature(my_feature_location, type="assembly_gap")
            gbk_record.features.append(my_feature)

    #print gbk_record[40000:50000].features
    SeqIO.write(gbk_record, sys.stdout, "genbank")


def parse_search_input(match_file):
    #MATCH: 44324 44336 CTAGCTAGCTAG
    all_matches_positions = []
    for line in open(match_file, 'r'):
        cols = line.split(' ')
        start = cols[1]
        stop = cols[2]
        all_matches_positions.append([int(start), int(stop)])
    return all_matches_positions

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-m", '--match_file', type=str, help="input match file")

    args = parser.parse_args()

    infile = "multi-rec.gbk"

    record = SeqIO.read(open(args.input_gbk,"r"), "genbank")
    gaps = parse_search_input(args.match_file)
    #print gaps
    add_gaps(record, gaps)