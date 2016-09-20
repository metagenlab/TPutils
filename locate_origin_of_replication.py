#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def locate_origin(contig_file, reference_dnaa=False, base_add=517):

    '''

    perform a tblastn to identify location od dnaa
    return chromosome name and position of the split

    :param contig_file:
    :param reference_dnaa: seqrecord object with reference dnaa sequence
    :return:
    '''

    '''
    :param contig_file:
    :param reference_dnaa:
    :return:
    '''

    import sys
    from Bio import SeqIO, Seq, SeqRecord
    from Bio.Alphabet import generic_dna, generic_protein
    from Bio.Blast.Applications import NcbitblastnCommandline
    import shell_command
    from Bio.Blast import NCBIXML

    if not reference_dnaa:

        seq = Seq.Seq('MSEKEIWEKVLEIAQEKLSAVSYSTFLKDTELYTIKDGEAIVLSS'
                      ' IPFNANWLNQQYAEIIQAILFDVVGYEVKPHFITTEELANYSNNETATPKETTKPSTET'
                      'TEDNHVLGREQFNAHNTFDTFVIGPGNRFPHAASLAVAEAPAKAYNPLFIYGGVGLGKT'
                      'HLMHAIGHHVLDNNPDAKVIYTSSEKFTNEFIKSIRDNEGEAFRERYRNIDVLLIDDIQ'
                      'FIQNKVQTQEEFFYTFNELHQNNKQIVISSDRPPKEIAQLEDRLRSRFEWGLIVDITPP'
                      'DYETRMAILQKKIEEEKLDIPPEALNYIANQIQSNIRELEGALTRLLAYSQLLGKPITT'
                      'ELTAEALKDIIQAPKSKKITIQDIQKIVGQYYNVRIEDFSAKKRTKSIAYPRQIAMYLS'
                      'RELTDFSLPKIGEEFGGRDHTTVIHAHEKISKDLKEDPIFKQEVENLEKEIRNV', generic_protein)

        reference_dnaa = SeqRecord.SeqRecord(seq,
                                 id="ADC36215.1",
                                 name="DnaA",
                                 description="Chromosomal replication initiator protein DnaA")

    cmd = 'formatdb -i %s -p F' % contig_file
    out, err, code = shell_command.shell_command(cmd)
    print out, err, code
    if code != 0:
        sys.stdout.write('problem with command: \n %s' % cmd)
    #path = os.path.abspath(contig_file)
    #blast_dir = os.path.dirname(path)
    #out, err, code = shell_command.shell_command('export BLASTDB=$BLASTDB:%s' % (blast_dir))
    #print out, err, code


    handle = open('dnaa.temp', 'w')
    SeqIO.write(reference_dnaa, handle, 'fasta')
    handle.close()

    tblastn_cline = NcbitblastnCommandline(query='dnaa.temp',
                                         db=contig_file,
                                         evalue=0.001,
                                         outfmt=5,
                                         out="dnaa_blast2.xml")

    stdout, stderr = tblastn_cline()

    result_handle = open("dnaa_blast2.xml", 'r')
    blast_records = [i for i in NCBIXML.parse(result_handle)]
    best_hit = blast_records[0].alignments[0].hsps[0]

    #print blast_records[0].alignments[0].hit_id
    contig = blast_records[0].alignments[0].hit_def
    identity_percentage = (best_hit.identities / best_hit.align_length)*100

    sbjct_end = best_hit.sbjct_end
    sbjct_start = best_hit.sbjct_start

    if sbjct_start > sbjct_end:
        #print 'strand 1', sbjct_start, sbjct_end
        split_location = sbjct_start - base_add
    else:
        #print 'strand -1', sbjct_start, sbjct_end
        split_location = sbjct_end + base_add

    sys.stdout.write('contig: %s --> hit with evalue of %s and identity of %s\n' % (contig, best_hit.expect, identity_percentage))

    return (contig, split_location)

def split_origin(fasta_contigs, contig_name, split_position, base_add=517):

    from Bio import SeqIO

    handle = open(fasta_contigs, 'r')
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    part1 = record_dict[contig_name][0:split_position]
    part1.name = contig_name + '_part_a'
    part1.id = contig_name + '_part_a'
    part1.description = contig_name + '_part_a'

    part2 = record_dict[contig_name][split_position:]
    part2.name = contig_name + '_part_b'
    part2.id = contig_name + '_part_b'
    part2.description = contig_name + '_part_b'


    record_dict[contig_name + '_part_a'] = part1
    record_dict[contig_name + '_part_b'] = part2
    del record_dict[contig_name]

    return [record_dict[i] for i in record_dict]


def reorder_contigs_with_mauve(reference, contigs, output_folder="mauve_ordering"):
    import shell_command

    cmd = 'link=$(readlink -f `which Mauve.jar`); java -Xmx500m -cp $link org.gel.mauve.contigs.ContigOrderer -output %s -ref %s -draft %s' % (output_folder,
                                                                                                                                               reference,
                                                                                                                                               contigs)
    out, err, code = shell_command.shell_command(cmd)

    if code != 0:
        sys.stdout.write('problem with command: \n %s\n' % cmd)



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO, SeqRecord, Seq
    import os
    import sys
    from Bio.Alphabet import generic_dna, generic_protein

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_contigs', type=str, help="input contig file")
    parser.add_argument("-g", '--reference_genbank', type=str, help="reference genbank file")
    parser.add_argument("-r", '--only_reorder', action='store_true', help="only reorder (do not cut any contig)")

    args = parser.parse_args()

    genbank_handle = open(args.reference_genbank)

    genbank_record = SeqIO.read(genbank_handle, 'genbank')

    for feature in genbank_record.features:
        if feature.type == 'CDS':

            try:
                name=feature.qualifiers['gene'][0]
            except:
                name= ''
            try:
                description=feature.qualifiers['product'][0]
            except:
                description='unknown product'
            try:
                id=feature.qualifiers['protein_id'][0]

            except:
                id='nooid'
            try:
                seq = Seq.Seq(feature.qualifiers['translation'][0], generic_protein)
            except:
                import sys
                print 'no translation for firt orf of the reference, please provide reference with translated sequences'
                sys.exit(255)


            reference = SeqRecord.SeqRecord(seq, id=id, name=name, description=description)

            reference_start = int(feature.location.start)
            break
    if not args.only_reorder:
        contig, position = locate_origin(args.input_contigs, reference, reference_start)

        sys.stdout.write('Spliting fasta using reference %s ...\n' % args.reference_genbank)

        new_record = split_origin(args.input_contigs, contig, position)

        out_name_prefix = os.path.basename(args.input_contigs).split('.')[0]


        output_handle = open(out_name_prefix + '_split.fa', 'w')
        sys.stdout.write('Writing edited contig file to %s ...\n' % (out_name_prefix + '_split.fa'))
        SeqIO.write(new_record, output_handle, 'fasta')

        sys.stdout.write('Reordering contigs with mauve\n')
        reorder_contigs_with_mauve(args.reference_genbank, out_name_prefix + '_split.fa')
        sys.stdout.write('Reordered contigs can be found in the folder "mauve_reorder"\n')
    else:
        sys.stdout.write('Reordering contigs with mauve\n')
        out_folder = 'mauve_%s' % (args.input_contigs.split('.')[0])
        reorder_contigs_with_mauve(args.reference_genbank, args.input_contigs,output_folder=out_folder)
        sys.stdout.write('Reordered contigs can be found in the folder "mauve_reorder"\n')
