#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# given input fasta protein and gen omes sequences
# get heatmap oh hit/non hit in each genome
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------


def tblan2heatmap(fasta_file, genome_files):
    from Bio.Blast.Applications import NcbitblastnCommandline
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
    import shell_command
    import sys

    out_list = []

    for genome in genome_files:
        sys.stdout.write('formatting %s...\n' % genome)
        cmd = 'formatdb -i %s -p F' % genome
        #a, b, c = shell_command.shell_command(cmd)
        sys.stdout.write('blasting %s vs %s\n' % (fasta_file, genome))
        out = genome.split('.')[0] + '_blast.xml'

        tblastn_cline = NcbitblastnCommandline(query=fasta_file,
                                             db=genome,
                                             evalue=0.001,
                                             outfmt=5,
                                             out=out)
        out_list.append(out)
        stdout, stderr = tblastn_cline()

    blasted_protein2seq_length = {}
    with open(fasta_file, 'r') as f:
        reference_fasta = SeqIO.parse(f, 'fasta')
        for seq in reference_fasta:
            blasted_protein2seq_length[seq.name] = len(seq.seq)
    sys.stdout.write('Protein ID list:\n')

    genome2protein2hit = {}
    for one_genome in out_list:
        with open(one_genome, 'r') as result_handle:
            genome2protein2hit[one_genome] = {}
            blast_records = [i for i in NCBIXML.parse(result_handle)]
            for one_protein in blast_records:
                try:
                    best_hit = one_protein.alignments[0].title.split(' ')[1] # .hsps[0]
                    #print dir(one_protein.alignments[0].hsps[0])
                    query_align_length = one_protein.alignments[0].hsps[0].query_end-one_protein.alignments[0].hsps[0].query_start

                    e_value = one_protein.alignments[0].hsps[0].expect
                    n_idenitical = one_protein.alignments[0].hsps[0].identities
                    align_length = one_protein.alignments[0].hsps[0].align_length
                    #print dir(one_protein)
                    query_name =  one_protein.query.split(' ')[0]
                    query_coverage = float(query_align_length)/blasted_protein2seq_length[query_name]
                    percent_identity = float(n_idenitical)/align_length

                    #print "% id", percent_identity, "% cov", query_coverage

                    if e_value < 0.05 and percent_identity > 0.8 and query_coverage > 0.5:
                        genome2protein2hit[one_genome][query_name] = percent_identity
                    else:
                        genome2protein2hit[one_genome][query_name] = 0
                except IndexError:
                    query_name =  one_protein.query.split(' ')[0]
                    genome2protein2hit[one_genome][query_name] = 0


    with open('results.tab', 'w') as w:
        proteins = genome2protein2hit[genome2protein2hit.keys()[0]].keys()
        header = '\t' + '\t'.join(proteins) + '\n'
        w.write(header)
        for genome in genome2protein2hit:
            values = [str(genome2protein2hit[genome][i]) for i in proteins]
            row = genome + '\t' + '\t'.join(values) +  '\n'
            w.write(row)


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-g", '--genomes', type=str, help="genomes files", nargs="+")


    args = parser.parse_args()


    tblan2heatmap(args.input_fasta, args.genomes)
