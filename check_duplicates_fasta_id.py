#!/usr/bin/python

def verify_alphabet(records, rm_dupli=False, outname=False):
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.Data import IUPACData
    from Bio.Alphabet import _verify_alphabet
    from Bio.Alphabet import IUPAC

    letters = list(IUPACData.extended_protein_letters)
    #print letters
    #print IUPAC.extended_protein.contains
    locus2count = {}
    nr_records = []
    for record in records:
        locus = record.name
        if locus not in locus2count: 
            locus2count[locus] = 1
            nr_records.append(record)
        else:
            locus2count[locus] += 1
    for locus in locus2count:
        if locus2count[locus] > 1:
            print '%s\t%s' % (locus, locus2count[locus])
            
    if rm_dupli:
        with open(outname, 'w') as f:
            SeqIO.write(nr_records, f, 'fasta')

if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input', type=str,help="input fasta")
    parser.add_argument("-r", '--rm_dupli', action='store_true',help="remove duplicate")
    args = parser.parse_args()
    records = [i for i in SeqIO.parse(args.input, 'fasta')]
    out_name = args.input.split('.')[0]+'_NR.faa'
    verify_alphabet(records, rm_dupli=args.rm_dupli, outname=out_name)
