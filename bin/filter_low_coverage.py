#!/usr/bin/python


def filter_coverage(fasta_file, outname, cov_limit):
    print cov_limit
    from Bio import SeqIO
    o = open(outname, 'w')
    keep=[]
    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, "fasta")
        for record in records:
            print record.id.split('_')[-1]
            coverage= float(record.id.split('_')[-1])
            if coverage>(int(cov_limit)):
                keep.append(record)
    SeqIO.write(keep, o, 'fasta')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-o", '--output_fasta', type=str, help="output fasta file", default="out.fa")
    parser.add_argument("-c", '--cov_limit', type=int,help='minimal spades coverage, default =1', default=1)
    args = parser.parse_args()
    filter_coverage(args.input_fasta, args.output_fasta, cov_limit=args.cov_limit)
