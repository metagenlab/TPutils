#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk',type=str,help="gbk file")
    parser.add_argument("-f", '--faa',type=str,help="faa file")

    args = parser.parse_args()
    genbank = [i for i in SeqIO.parse(args.gbk, "genbank")]
    faa = SeqIO.parse(args.faa, "fasta")
    print faa
    locus2seq = {}
    for record in faa:
        locus2seq[record.name] = record.seq
    #print locus2seq

    for record in genbank:
        for feature in record.features:
            if feature.type == 'CDS':
                               
                feature.qualifiers['translation'] = [str(locus2seq[feature.qualifiers['locus_tag'][0]])]
    output_handle = open("test.gbk", 'w')            
    SeqIO.write(genbank,output_handle, "genbank")
