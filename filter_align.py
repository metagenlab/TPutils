#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# remove invariant sites from fasta alignment
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------


def remove_invar(align_record):
    from Bio.Align import MultipleSeqAlignment
    #new_align = MultipleSeqAlignment()
    from Bio.Alphabet import generic_dna
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    keep = []

    for column in range(0, len(align_record[0, :])):
        col_char = align[:, column]
        lst = list(set([str(i[0]) for i in col_char]))
        try:
            lst.pop(lst.index('-'))
        except:
            pass
        try:
            lst.pop(lst.index('N'))
        except:
            pass
        if len(lst) > 1:
            keep.append(column)
    print 'len keep', len(keep)
    new_sequences = []
    for seq in align_record:
        new_seq_str = "".join([seq[i] for i in keep])
        #print "new_seq_str", len(new_seq_str), new_seq_str, type(new_seq_str)
        new_seq = SeqRecord(Seq(new_seq_str, generic_dna), id=str(seq.id))
        new_sequences.append(new_seq)

    return MultipleSeqAlignment(new_sequences)

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    from Bio import AlignIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    args = parser.parse_args()



    #input_handle = open(args.input_gbk, "rU")
    align = AlignIO.read(args.input_fasta, "fasta")
    if not args.outname:
        # use input file to rename the file
        #outname = args.input_gbk.split(".")[0]+".faa"

        # use record id to rename the file, remove version number using split
        outname = args.input_fasta.split('.')[0] + "_no_invar.fa"
    else:
        outname = args.outname


    new_align = remove_invar(align)

    #print type(new_align), new_align

    handle = open(outname, 'w')

    AlignIO.write(new_align, outname, "fasta")