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



def get_coverage_target_coordinates(samtools_depth, coordinate_file, ref_fasta):

    '''

    coordinate file format:
    gene, record.name, start, end

    :param samtools_depth:
    :param coordinate_file:
    :return:
    '''

    import numpy
    import contig_coverage_from_depth  
    
    id2l = contig_coverage_from_depth.contig_id2contig_length(ref_fasta)
    id2cov, id2depth = contig_coverage_from_depth.contig_name2contig_coverage(samtools_depth)
    contig2med = contig_coverage_from_depth.contig_id2median_id(id2cov, id2depth)    

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
    median_depth = numpy.median(samtools_dataframe.iloc[:,1])
    #print samtools_dataframe
    print ("Contig\tgene\tstart\tend\tdepth\tratio_assembly\tcontig_depth\tcontig_ratio_depth\tcontig_length")
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
            gene_median = numpy.median(subset_table.iloc[start_pos:end_pos,1])
            print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (record, gene, start_pos, end_pos, gene_median, round(gene_median/median_depth, 2), contig2med[record], round(contig2med[record]/median_depth,2),  id2l[record]))
    print ('all_assembly\t-\t1\t%s\t%s\t-\t-\t-' % (len(samtools_dataframe.iloc[:,1]),
                                                     numpy.median(samtools_dataframe.iloc[:,1])))










if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--input_depth', type=str, help="input samtools depth file")
    parser.add_argument("-g", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-c", '--coord_file', type=str, help="coord file")
    parser.add_argument("-f", '--ref_fasta', type=str, help="fasta file")

    args = parser.parse_args()

    #parse_smatools_depth(args.input_depth)
    get_coverage_target_coordinates(args.input_depth, args.coord_file, args.ref_fasta)
