#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import locate_origin_of_replication



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk', type=str, help="input gbk")
    parser.add_argument("-c", '--input_contigs', type=str, help="input contigs")
    parser.add_argument("-r", '--reference_genbank', type=str, help="input reference gbk")
    args = parser.parse_args()
    import glob
    import os
    from Bio import SeqIO

    locate_origin_of_replication.reorder_contigs_with_mauve(args.reference_genbank,
                                                            args.input_contigs,
                                                            output_folder="mauve_reorder")
    print 'ok'


    working_dir = os.getcwd()
    mauve_dir = os.path.join(working_dir,"mauve_reorder/")

    mauve_dir_list = glob.glob(mauve_dir+'*')
    print mauve_dir_list
    ndir= len(mauve_dir_list)


    accession2record = SeqIO.to_dict(SeqIO.parse(open(args.gbk), 'genbank'))

    fasta = os.path.join(mauve_dir, 'alignment%s/%s.fas' % (ndir, args.input_contigs))
    with open(fasta, 'r') as f:
        for n, record in enumerate(SeqIO.parse(f, 'fasta')):
            if n == 0:
                ordered_rtecord = [accession2record[record.id]]
            else:
                ordered_rtecord.append(accession2record[record.id])
    with open(args.input_contigs.split('.')[0]+'_reordered.gbk', 'w') as f:
        SeqIO.write(ordered_rtecord, f,'genbank')