#!/usr/bin/python

def verify_alphabet(records):
    from Bio.Seq import Seq
    from Bio.Data import IUPACData
    from Bio.Alphabet import _verify_alphabet
    from Bio.Alphabet import IUPAC

    letters = list(IUPACData.extended_protein_letters)
    #print letters
    #print IUPAC.extended_protein.contains
    
    for record in records:
        my_seq = Seq(str(record.seq),
                     IUPAC.extended_protein)
        if  _verify_alphabet(my_seq) is True:
            continue
        else:
            illegal_char = []
            for letter in list(str(record.seq)):
                if letter not in letters:
                    illegal_char.append(letter)
            print '%s\t%s' % (record.name, '-'.join(illegal_char))


if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input', type=str,help="input fasta")
    args = parser.parse_args()
    records = [i for i in SeqIO.parse(args.input, 'fasta')]
    verify_alphabet(records)
