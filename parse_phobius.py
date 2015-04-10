#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# parse phobius short output files
# Date: 2015
# ---------------------------------------------------------------------------


def parse_short_phobius(*input_files):
    import re
    phobius_result = {}
    for input_file in input_files:
        with open(input_file) as f:
            for i, row in enumerate(f):
                if i == 0:
                    assert (row.strip() == 'SEQENCE ID                     TM SP PREDICTION'), "Unexpected Phobis header"

                data = re.split("\s+", row.rstrip())

                if "|" in data[0]:
                    data[0] = data[0].split("|")[3]
                one_protein = {}
                one_protein['TM'] = data[1]
                one_protein['SP'] = data[2]
                one_protein['prediction'] = data[3]
                phobius_result[data[0]] = one_protein
    return phobius_result


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_phobius', type=str, help="input phobius file", nargs ='+')

    #f = SeqIO.parse('FR872580.faa', "fasta")
    #for i in f:
    #    print i

    args = parser.parse_args()

    print len(parse_short_phobius(*args.input_phobius).values())


