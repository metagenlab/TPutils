#!/usr/bin/env python

def get_size(records, size_cutoff=0):
    cumul = 0
    for record in records:
        if len(record.seq) > size_cutoff:
            cumul += len(record.seq)
    return cumul



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--input_fna', type=str, help="input fna file")
    parser.add_argument("-n", '--contig_size_limit', type=int, help="Min contig size to consider", default=0)



    args = parser.parse_args()
    fasta_records = list(SeqIO.parse(args.input_fna, "fasta"))
    cumul_size = get_size(fasta_records, args.contig_size_limit)

    print(cumul_size)