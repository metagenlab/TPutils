#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# reannotate genomes unising prokka
# reannotation with prokka and transfer of the locus tag and protein ID when identical ORF are found
# input: either 1. genbank file only (homogenisation of the annotation with a single annotation pipeline)
#               2. 1x genbank file 1x fna file with new version of the genome (transfer of the annotation to a new version
#                  of the genome)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: août 2015
# ---------------------------------------------------------------------------



def prokka_reannotation(*genbank_files):
    import re
    import datetime
    all_locus = []
    l = open('reannotation_prokka_log.txt', 'w')
    l.write('reference_name\tnew_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical\n')
    print 'reference_name\tnew_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical'

    for genbank_file in genbank_files:
        from Bio import SeqIO
        import shell_command
        seq_records = list(SeqIO.parse(genbank_file, "genbank"))
        if len(seq_records) > 1:
            raise IOError('Wrong input, only scaffolded genomes should be reannotated with this script')
        else:
            record = seq_records[0]
        #print dir(record)
        #print record.annotations
        record_annotations = record.annotations
        record_annotations['comment'] = 'Genome reannotated using PROKKA version 1.1'
        record_name = record.name
        record_id = record.id
        record_description = record.description
        record_dbxrefs = record.dbxrefs

        # create locus tag based on genus and species name
        organism = re.sub('\'','', record_annotations['source']).split(' ')

        if len(organism) > 2:
            locus_tag = "P%s%s%s" % (organism[0][0], organism[1][0:2], organism[2][0])
        elif len(organism) == 2:
            locus_tag = "P%s%s" % (organism[0][0:2], organism[1][0])
        else:
            print record

        # check if the new locus_tag is unique, otherwise add a count
        i = 2
        no_match=False
        if locus_tag in all_locus:
            #print record
            while no_match == False:
                locus_tag = locus_tag+str(i)
                #print locus_tag
                if locus_tag in all_locus:
                    i+=1
                else:
                    no_match = True
                    all_locus.append(locus_tag)
        else:
            all_locus.append(locus_tag)

        # wite fasta and annotate it using prokka
        with open('temp_genome.fna', 'w') as f:
            f.write('>temp_seq\n%s' % record.seq)

        cmd = 'prokka --kingdom Bacteria --compliant --centre CHUV --locustag %s --outdir %s -genus temp_genus -strain temp_strain temp_genome.fna' % (locus_tag, locus_tag)


        today = datetime.date.today()
        date = today.strftime('%m%d%Y')
        prokka_genbank = '%s/%s_%s.gbk' %(locus_tag, locus_tag, date)
        #print prokka_genbank


        cmd2 = 'prokka --kingdom Bacteria --compliant --proteins proteins.faa --locustag Citr -genus Citronella -strain virus citro_spades_123_1000.fa'


        out, err, n = shell_command.shell_command(cmd)
        reanotated_gbk = list(SeqIO.parse(prokka_genbank, "genbank"))

        reanotated_gbk[0].id = record_id
        reanotated_gbk[0].name = record_name
        reanotated_gbk[0].annotations = record_annotations
        reanotated_gbk[0].description = record_description
        reanotated_gbk[0].dbxrefs = record_dbxrefs
        reanotated_gbk[0].features[0] = record.features[0]

        #count number of identical ORF
        ref_CDS = 0
        for feature in record.features:
            if feature.type=='CDS':
                if ref_CDS == 0:
                    ref_locus_tag = feature.qualifiers['locus_tag'][0].split('_')[0]
                ref_CDS+=1

        new_CDS = 0
        identical_CDS = 0
        # count number of identical features (exact same location)
        for new_feature in reanotated_gbk[0].features:
            if new_feature.type == 'CDS':
                new_CDS+=1
                for ref_feature in record.features:
                    if ref_feature.type == 'CDS':
                        if ref_feature.location.start == new_feature.location.start and ref_feature.location.end == new_feature.location.end:
                            identical_CDS +=1
                            break

        accession = reanotated_gbk[0].annotations["accessions"][0]

        out_name = '%s_prokka_reannot.gbk' % accession
        with open(out_name, 'w') as f:
            SeqIO.write(reanotated_gbk[0], f, 'genbank')

        #print 'Ref CDS', ref_CDS, 'New CDS', new_CDS, 'Identical CDS', identical_CDS
        #'reference_file_name\tnew_file_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical'
        l.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (genbank_file, out_name, accession, ref_locus_tag, locus_tag, ref_CDS, new_CDS, identical_CDS))
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (genbank_file, out_name, accession, ref_locus_tag, locus_tag, ref_CDS, new_CDS, identical_CDS)


def prokka_annotation_transfer(genbank_file, new_fna_file):
    import re
    all_locus = []
    l = open('reannotation_prokka_log.txt', 'w')
    l.write('reference_name\tnew_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical\n')
    print 'reference_name\tnew_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical'


    from Bio import SeqIO
    import shell_command
    seq_records = list(SeqIO.parse(genbank_file, "genbank"))
    if len(seq_records) > 1:
        raise IOError('Wrong input, only scaffolded genomes should be reannotated with this script')
    else:
        record = seq_records[0]

    record_annotations = record.annotations
    record_annotations['comment'] = 'Genome reannotated using PROKKA version 1.1'
    record_name = record.name
    record_id = record.id
    record_description = record.description
    record_dbxrefs = record.dbxrefs

    organism = re.sub('\'','', record_annotations['source']).split(' ')

    if len(organism) > 2:
        locus_tag = "P%s%s%s" % (organism[0][0], organism[1][0:2], organism[2][0])
    elif len(organism) == 2:
        locus_tag = "P%s%s" % (organism[0][0:2], organism[1][0])
    else:
        print record

    '''
    cmd = 'prokka --kingdom Bacteria --compliant --locustag %s --centre CHU --outdir %s -genus temp_genus -strain temp_strain %s' % (locus_tag, locus_tag, new_fna_file)
    '''
    import datetime
    today = datetime.date.today()
    date = today.strftime('%m%d%Y')
    prokka_genbank = '%s/%s_%s.gbk' %(locus_tag, locus_tag, date)
    '''
    out, err, n = shell_command.shell_command(cmd)
    print cmd
    print out
    print err
    '''
    reanotated_gbk = list(SeqIO.parse(prokka_genbank, "genbank"))

    reanotated_gbk[0].id = record_id
    reanotated_gbk[0].name = record_name
    reanotated_gbk[0].annotations = record_annotations
    reanotated_gbk[0].description = record_description
    reanotated_gbk[0].dbxrefs = record_dbxrefs
    #reanotated_gbk[0].features[0] = record.features[0]

    # count reference ORF
    ref_CDS = 0
    for feature in record.features:
        if feature.type=='CDS':
            if ref_CDS == 0:
                ref_locus_tag = feature.qualifiers['locus_tag'][0].split('_')[0]
            ref_CDS+=1

    new_CDS = 0
    identical_CDS = 0

    # count number of identical features (exact same translation)
    for new_feature in reanotated_gbk[0].features:
        if new_feature.type == 'CDS':
            new_CDS+=1
            match = 0
            locus = []
            for ref_feature in record.features:
                if ref_feature.type == 'CDS':
                    if ref_feature.qualifiers['translation'][0] == new_feature.qualifiers['translation'][0]:
                        identical_CDS +=1
                        match +=1
                        locus.append(ref_feature.qualifiers['locus_tag'][0])
            if match == 0:

                print '%s\t%s' % (new_feature.qualifiers['locus_tag'][0], '-')
            else:
                all_match = ', '.join(locus)
                print '%s\t%s' % (new_feature.qualifiers['locus_tag'][0], all_match)

    print '\nmulticopy locus'
    multicopy_list = []


    # find identical translations
    for new_feature1 in reanotated_gbk[0].features:
        if new_feature1.type == 'CDS':
            multicopy_list_temp = []
            multicopy_annotation_temp = ''
            for new_feature2 in reanotated_gbk[0].features:
                if new_feature2.type == 'CDS':
                    if new_feature1.qualifiers['translation'][0] == new_feature2.qualifiers['translation'][0]:
                        if new_feature2.qualifiers['locus_tag'][0] in multicopy_list:
                            continue
                        elif new_feature1.qualifiers['locus_tag'][0] == new_feature2.qualifiers['locus_tag'][0]:
                            continue
                        if new_feature1.qualifiers['locus_tag'][0] not in multicopy_list_temp:
                            multicopy_list_temp.append(new_feature1.qualifiers['locus_tag'][0])
                        multicopy_annotation_temp = new_feature2.qualifiers['product'][0]
                        multicopy_list_temp.append(new_feature2.qualifiers['locus_tag'][0])
            if len(multicopy_list_temp) > 0:
                multicopy_list+=(multicopy_list_temp)
                print '\t'.join(multicopy_list_temp) + '\t%s' % multicopy_annotation_temp

    # find missing translation in new assembly version
    print '\nNo match:'
    for ref_feature in record.features:
        if ref_feature.type == 'CDS':
            match = 0
            for new_feature in reanotated_gbk[0].features:
                if new_feature.type == 'CDS':
                    if new_feature.qualifiers['translation'][0] == ref_feature.qualifiers['translation'][0]:
                        match+=1
            if match == 0:
                print ref_feature.qualifiers['locus_tag'][0] + '\t%s' % ref_feature.qualifiers['product'][0]

    accession = reanotated_gbk[0].annotations["accessions"][0]

    out_name = '%s_prokka_reannot.gbk' % accession
    with open(out_name, 'w') as f:
        SeqIO.write(reanotated_gbk[0], f, 'genbank')

    #print 'Ref CDS', ref_CDS, 'New CDS', new_CDS, 'Identical CDS', identical_CDS
    #'reference_file_name\tnew_file_name\taccession\tref_locus_tag\tnew_locus_tag\tref_n_CDS\tnew_n_CDS\tn_CDS_identical'
    #l.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (genbank_file, out_name, accession, ref_locus_tag, locus_tag, ref_CDS, new_CDS, identical_CDS))
    #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (genbank_file, out_name, accession, ref_locus_tag, locus_tag, ref_CDS, new_CDS, identical_CDS)


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", nargs='+')
    parser.add_argument("-f", '--input_fna', type=str, help="input fna file (optionnal)")
    parser.add_argument("-t",'--anotation_transfer',action='store_true', help="transfer annotation to a new fna")


    args = parser.parse_args()

    if not args.anotation_transfer:
        prokka_reannotation(*args.input_gbk)
    else:
        prokka_annotation_transfer(args.input_gbk[0], args.input_fna)