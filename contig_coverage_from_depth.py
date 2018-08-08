#!/usr/bin/env python

import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
rpy2.robjects.numpy2ri.activate()



def contig_name2contig_coverage(samtool_depth_file):

    contig2coverage = {}
    contig2depth = {}
    for line in open(samtool_depth_file, 'r'):
        data = line.rstrip().split('\t')
        if data[0] not in contig2coverage:
            contig2coverage[data[0]] = 1
            contig2depth[data[0]] = [int(data[2])]

        else:
            contig2coverage[data[0]] += 1
            contig2depth[data[0]].append(int(data[2]))
    return contig2coverage, contig2depth

def contig_id2contig_length(contigs_file):

    from Bio import SeqIO

    fasta_handle = open(contigs_file, 'r')
    id2length = {}
    for record in SeqIO.parse(fasta_handle, "fasta"):
        id2length[record.name] = len(record)
    return id2length


def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
        return None
    if len(lst) % 2 == 1:
        return lst[((len(lst)+1)/2)-1]
    else:
        return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

def contig_id2median_id(id2cov, id2depth):
    contig2med = {}
    for id in id2cov:
        m = median(id2depth[id])
	contig2med[id] = m
    return contig2med


if __name__ == '__main__':
    import argparse
    import numpy
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_samtools_depth', type=str, help="input samtools depth files")
    parser.add_argument("-m", '--input_contigs', type=str, help="input contigs")


    args = parser.parse_args()

    id2l = contig_id2contig_length(args.input_contigs)
    id2cov, id2depth = contig_name2contig_coverage(args.input_samtools_depth)
    contig2med = contig_id2median_id(id2cov, id2depth)
    for id in contig2med:
        print '%s\t%s\t%s\t%s\t%s' % (id, id2cov[id], id2l[id], round(float(id2cov[id])/id2l[id]*100), contig2med[id])
