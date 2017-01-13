#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# get 6 frames tranlastion of nucl sequence
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------


def find_orfs_with_trans(seq, trans_table=11, min_protein_length=30):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

def format_seqquences(orf_list, outname, header='orf'):
    with open(outname, 'w') as f:
        # start, end, strand, pro
        for i, cds in enumerate(orf_list):
            if str(cds[2]) == '-1':
                f.write(">%s:%s: complement(%s..%s)\n%s\n" % (header,
                                                    i,
                                                    cds[0],
                                                    cds[1],
                                                    cds[3]))
            else:
                f.write(">%s:%s: %s..%s\n%s\n" % (header,
                                                    i,
                                                    cds[0],
                                                    cds[1],
                                                    cds[3]))




if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fna', type=str, help="input gbk file")
    parser.add_argument("-l", '--min_len', type=int, help="minimum translated sequence length (default 30)", default=30)
    parser.add_argument("-p", '--trans_table', type=int, help="translation table (default: 11)", default=11)

    args = parser.parse_args()
    handle = open(args.input_fna, 'r')
    record = SeqIO.read(handle, 'fasta')

    all_cds = find_orfs_with_trans(record.seq, trans_table=args.trans_table, min_protein_length=args.min_len)
    format_seqquences(all_cds)
