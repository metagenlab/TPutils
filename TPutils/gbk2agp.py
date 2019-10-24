#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------


def get_assembly_concatenated_length(gbk_file, unaligned_contig_list):

    concatenated_length = 0

    input_handle = open(gbk_file, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    for i, record in enumerate(seq_records):
        if record.id not in unaligned_contig_list:
            concatenated_length+=len(record.seq)
    return concatenated_length

def gbk2agp(gbk_file,
            out_agp="kpge.agp",
            unaligned_contig_list = [],
            prefix="scaffold",
            chromosome=False,
            fasta=False):

    from Bio import SeqIO
    from optparse import OptionParser
    import re
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    '''

    by default:
        - all sequences are part of the same chromosome
        - and ordered based on genome alignment
        ==> will be scaffolded with gaps of unknown size

        ATTENTION scaffold with NNN are not yet considered

    :param gbk_file:
    :return:
    '''

    header = '''##agp-version	2.0
            # ORGANISM: Klebsiella pneumoniae KpGe
            # TAX_ID: 573
            # ASSEMBLY DATE: March-2017
            # DESCRIPTION: assembly of scaffolds from WGS contigs assembled with SPAdes v3.10.1
            # Scaffolds reordered based on the reference genome: Klebsiella pneumoniae subsp. pneumoniae HS11286 (NC_016845)'''
    print header
    if chromosome:
        concat_assembly_length_excluding_unaligned = get_assembly_concatenated_length(gbk_file, unaligned_contig_list)

    input_handle = open(gbk_file, "rU")
    if not fasta:
        seq_records = list(SeqIO.parse(input_handle, "genbank"))
    else:
        print 'fasta'
        seq_records = list(SeqIO.parse(input_handle, "fasta"))
    scaffold_count = 1
    contig_count = 1

    print seq_records
    contig_list = []
    for i, record in enumerate(seq_records):
        if not chromosome:
            print 'splitting sequence...'
            # split scaffolds into contigs
            seq = str(record.seq)
            contig_sequences = re.split("N{200,}", seq)
            gap_list = re.findall("N{200,}", seq)
            #if record.accession not in unaligned_contig_list:
            start = 1
            part_count = 1
            for n, one_contig in enumerate(contig_sequences):
                #print record.id, len(contig_sequences)

                contig_list.append(SeqRecord(Seq(one_contig,
                                           record.seq.alphabet),
                                           id='contig_%s' % contig_count, name='contig_%s' % contig_count,
                                           description=record.description))

                end=start+len(one_contig)-1
                contig_row='%s\t%s\t%s\t%s\tD\t%s_%s\t1\t%s\t+' % ("Chr1",
                                                             start,
                                                             end,
                                                             part_count,
                                                                   prefix,
                                                             contig_count,
                                                             len(one_contig)
                                                             )
                print contig_row
                part_count+=1
                if len(contig_sequences) > 1 and n < len(contig_sequences)-1:
                    #print 'start, end',start, end
                    gap_start = end+1
                    gap_end = gap_start + 99#len(gap_list[n-1])-1
                    gap_row='%s\t%s\t%s\t%s\tU\t%s\t%s\t%s\tna' % ("Chr1",
                                                                 gap_start,
                                                                 gap_end,
                                                                 part_count,
                                                                 100,
                                                                 "scaffold", #gap_type
                                                                 "no" #linkage==> yes
                                                                 )
                    part_count+=1
                    start = end+101# unkown gap size == 100   len(gap_list[n-1])+1
                    print gap_row
                #print contig_row
                contig_count+=1


        if chromosome:
            entry_length = concat_assembly_length_excluding_unaligned


        #else:
        #    scaffold_count+=1
        #    new_row='chr_scaffold%s\t1\t%s\t1\tW\t%s\t1\t%s\t?' % (scaffold_count,
        #                                             len(record.seq),
        #                                             record.accession,
        #                                                        len(record.seq))

        #if i != len(seq_records):
        #    gap_row = 'chr1\t'
        out_contig = gbk_file.split('.')[0] + '_contigs.fa'
        SeqIO.write(contig_list, out_contig, 'fasta')


def gbk2unlocated_file_list(gbk_file,
                            chromosome_file="chromosomes.txt",
                            unlocated="unlocated.txt"):

    '''
    create one chromosome file
    create one
    :param gbk_file:
    :return:
    '''






if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-f", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-p", '--prefix_accession', type=str, help="prefix for the accession of sub-elements of the agp")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)

    args = parser.parse_args()

    if args.input_fasta:
        gbk2agp(args.input_fasta, prefix=args.prefix_accession,fasta=True)
    elif args.input_gbk:
        gbk2agp(args.input_gbk, prefix=args.prefix_accession)
