#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def compare_sequences(seq1, seq2, print_var):
    if len(seq1) != len(seq2):
        raise IOError('Input sequences should have the same length!')
    ident = 0
    diff = 0
    for i in range(0, len(seq1)):
        
        if seq1[i].lower() == 'n' or  seq2[i].lower() == 'n':
            continue
        elif seq1[i].lower() == '-' or  seq2[i].lower() == '-':
            continue
        else:
            if seq1[i] == seq2[i]:
                ident+=1
            else:
                if print_var:
                    print "%s\t%s\t%s" % (i+1, seq1[i], seq2[i])
                diff+=1
    return ident, diff
def count_diff(records, print_var):
    diff_count = {}
    for i in range(0, len(records)):
        for y in range(0, len(records)):
            if records[i].name not in diff_count:
                diff_count[records[i].name] = {}
            diff_count[records[i].name][records[y].name] = list(compare_sequences(records[i].seq, records[y].seq, print_var))
    return diff_count
                




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--multifasta', type=str, help="input multifasta")
    parser.add_argument("-p",'--print_var',action='store_true',help="print variants")
    args = parser.parse_args()
    from Bio import SeqIO
    diff_count = count_diff([i for i in SeqIO.parse(args.multifasta, "fasta")], args.print_var)

    print '\t' + '\t'.join(diff_count.keys())
    for record1 in diff_count:
        row = "%s\t" % record1
        for record2 in diff_count:
            row+='%s\t' % diff_count[record1][record2][1]
        print row
