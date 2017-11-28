#! /usr/bin/env python



def fasta2tab(fasta_file):
    from Bio import SeqIO
    import re

    try:
        file_description = re.findall('G[A-Za-z0-9_]+\.[0-9]+',fasta_file)[0]
    except:
        file_description = fasta_file.split('.')[0]

    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            print "%s\t%s" % (file_description, record.name)

fasta2tab('GCF_000721275.1_ASM72127v1_protein.faa')