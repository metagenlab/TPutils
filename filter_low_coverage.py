#!/usr/bin/python


def filter_coverage(fasta_file, outname):
    from Bio import SeqIO
    o = open(outname, 'w')
    keep=[]
    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, "fasta")
        for record in records:
            coverage= float(record.id.split('_')[-3])
            if coverage>1:
                keep.append(record)
    SeqIO.write(keep, o, 'fasta')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-o", '--output_fasta', type=str, help="output fasta file", default="out.fa")
    args = parser.parse_args()
    filter_coverage(args.input_fasta, args.output_fasta)
