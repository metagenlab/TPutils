#! /usr/bin/env python

# calculate identity matrix based on multiple sequence alignment
# input: fasta MSA
# identity calculation: identical_sites/aligned_sites (gaps are not taken into account)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 02.2015
# ---------------------------------------------------------------------------

from Bio import AlignIO
import numpy as np


def pairewise_identity(seq1, seq2):
    '''
    :param seq1:
    :param seq2:
    :return: sequences identity (%)
    '''
    A=list(seq1)
    B=list(seq2)
    identical_sites = 0
    aligned_sites = 0
    gaps=0
    for n in range(0, len(A)):
        if A[n] != "-" and B[n] != "-":
            aligned_sites+=1
        else:
            continue
        if A[n]==B[n]:
            identical_sites+=1
    if aligned_sites == 0:
        # return false if align length = 0
        return False
    else:
        return 100*(identical_sites/float(aligned_sites))


def get_identity_matrix_from_multiple_alignment(alignment): 
    identity_matrix = np.empty([len(alignment), len(alignment)])
    for x in range(0,len(alignment)):
        for y in range(0, len(alignment)):
            print type(alignment[x])
            identity_matrix[x,y] = pairewise_identity(alignment[x], alignment[y])
    return identity_matrix

def write_id_table(align, array):
    ids = [i.name for i in align]
    print "\t" + "\t".join(ids)
    for i in range(0,len(array)):
        name = ids[i]
        row = array[i,:]
        row = [str(round(i,2)) for i in row]
        print name + "\t" + "\t".join(row) #+ "\n"

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--align',type=str,help="input alignment (fasta FORMAT)")

    args = parser.parse_args()
    align = AlignIO.read(args.align, "fasta")
    id_table = get_identity_matrix_from_multiple_alignment(align)
    write_id_table(align, id_table)
