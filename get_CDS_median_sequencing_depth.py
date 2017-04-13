#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def gbk2median_CDS_cov(gbk_file, depth_file):

    from Bio import SeqIO
    record2locus2location = {}
    with open(gbk_file, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            record2locus2location[record.name] = {}
        for feature in record.features:
            if feature.type == 'CDS':
                record2locus2location[record.name][feature.qualifiers['locus_tag']] = [feature.location.start, feature.location.end]

def parse_smatools_depth(samtools_depth):
    import pandas

    with open(samtools_depth, 'r') as f:
        table = pandas.read_csv(f, sep='\t', header=None, index_col=0)
    return table



def get_coverage_target_coordinates(samtools_depth, coordinate_file):

    '''

    coordinate file format:
    gene, record.name, start, end

    :param samtools_depth:
    :param coordinate_file:
    :return:
    '''

    import numpy

    record2gene2coord = {}
    with open(coordinate_file) as f:
        for row in f:
            data = row.rstrip().split('\t')
            if data[1] not in record2gene2coord:
                record2gene2coord[data[1]] = {}
                record2gene2coord[data[1]][data[0]] = [data[2], data[3]]
            else:
                record2gene2coord[data[1]][data[0]] = [data[2], data[3]]

    samtools_dataframe = parse_smatools_depth(samtools_depth)
    #print samtools_dataframe
    for record in record2gene2coord:
        for gene in record2gene2coord[record]:
            # attention range
            # index commence à 0
            # range n'inclu pas le dernier chiffre
            start_pos = int(record2gene2coord[record][gene][0])-1
            end_pos = int(record2gene2coord[record][gene][1])
            if start_pos > end_pos:
                start_pos = int(record2gene2coord[record][gene][1])-1
                end_pos = int(record2gene2coord[record][gene][0])
            #median_coverage = numpy.median(samtools_dataframe.loc[0:10,2])
            subset_table = samtools_dataframe.loc[record]
            print "%s\t%s\t%s\t%s\t%s" % (record, gene, start_pos, end_pos, numpy.median(subset_table.iloc[start_pos:end_pos,1]))
    print 'all_assembly\t-\t1\t%s\t%s' % (len(samtools_dataframe.iloc[:,1]),
                                                     numpy.median(samtools_dataframe.iloc[:,1]))










if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--input_depth', type=str, help="input samtools depth file")
    parser.add_argument("-g", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-c", '--coord_file', type=str, help="coord file")

    args = parser.parse_args()

    #parse_smatools_depth(args.input_depth)
    get_coverage_target_coordinates(args.input_depth, args.coord_file)