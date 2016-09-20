#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys


def filter_plasmid(record_list):
    '''
    Return 2 lists: one with plasmid record(s), one without plasmids
    :param record_list:
    :return:
    '''

    plasmid_record_list = []
    chromosome_record_list = []

    for record in record_list:
        if 'plasmid' in record.description:
            plasmid_record_list.append(record)
        else:
            chromosome_record_list.append(record)
    return (chromosome_record_list, plasmid_record_list)


def is_annotated(gbk_record):
    if len(gbk_record.features) == 1 and gbk_record.features[0].type == 'source':
        return False
    else:
        return True


def count_missing_locus_tags(gbk_record):
    count_CDS = 0
    count_no_locus = 0
    for feature in gbk_record.features:
        if feature.type=='CDS':
            count_CDS += 1
            try:
                test = feature.qualifiers['locus_tag']
                #print 'locus ok'
            except:
                # missing locus case
                count_no_locus += 1
        pass
    return count_no_locus, count_CDS






def remove_record_taxon_id(record):
    if record.features[0].type == 'source':
        # delete evendual taxon_id (taxon id will be reattributed in the db based on the organism name)
        #try:
        if 'db_xref' in record.features[0].qualifiers:
            for item in record.features[0].qualifiers['db_xref']:
                if 'taxon' in item:
                    index = record.features[0].qualifiers['db_xref'].index(item)
                    record.features[0].qualifiers['db_xref'].pop(index)
                    if 'strain' in record.features[0].qualifiers:
                        if ';' in record.features[0].qualifiers['strain'][0]:
                            print 'ACHRTUNG: record %s has 2 strain names! check and edit source manually' % record.name
                            # put everythink lower size
                            strain = record.features[0].qualifiers['strain'][0].split(';')[1]
                        else:
                            strain = record.features[0].qualifiers['strain'][0]
                        if strain.lower() not in record.annotations['source'].lower():
                            msg = 'changing source %s' % record.annotations['source']
                            record.annotations['source'] += ' %s' % strain
                            print "ACHRTUNG\t" + msg + "\t--> %s " % record.annotations['source']
                            record.annotations['organism'] = record.annotations['source']

                    else:
                        print 'ACHTUNG\t no strain  record \t%s, source uniqueness should be checked' % record.name
    else:
        print 'ACHRTUNG\t no source for record \t%s' % record.name
    return record
    
def clean_description(description):
    import re
    description = re.sub(", complete genome\.", "", description)
    description = re.sub(", complete genome", "", description)
    description = re.sub(", complete sequence\.", "", description)
    description = re.sub("strain ", "", description)
    description = re.sub("str\. ", "", description)
    description = re.sub(" complete genome sequence\.", "", description)
    description = re.sub(" complete genome\.", "", description)
    description = re.sub(" chromosome", "", description)
    description = re.sub(" DNA", "", description)
    description = re.sub("Merged record from ", "", description)
    description = re.sub(", wgs", "", description)
    description = re.sub("Candidatus ", "", description)
    description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
    description = re.sub("complete genome, isolate", "", description)
    description = re.sub(" complete", "", description)
    #description = re.sub(", : 1.", "", description)
    
    return description

def check_gbk(gbff_files):
    from Bio import SeqIO
    import reannotate_genomes
    reannotation_list = []

    for gbff_file in gbff_files:
        print gbff_file
        records = list(SeqIO.parse(open(gbff_file, "r"), "genbank"))

        for record in records:
            n_missing, total = count_missing_locus_tags(record)
            if n_missing > 0:
                print '%s/%s missing locus tag for record %s' % (n_missing, total, record.name)

        chromosome, plasmids = filter_plasmid(records)

        reannotated_records = []

        cleaned_records = []
        plasmid_reannot = False
        chromosome_reannot = False

        if len(plasmids) > 0:
            #print '######### plasmid(s) ############'
            for plasmid in plasmids:
                if is_annotated(plasmid):
                    out_name = plasmid.name + '.gbk'
                    plasmid.description = clean_description(plasmid.description)
                    plasmid = remove_record_taxon_id(plasmid)
                    with open(out_name, 'w') as f:
                        SeqIO.write(plasmid, f, 'genbank')
                    cleaned_records.append(plasmid)
                else:
                    plasmid_reannot = True
                    reannotation_list.append([plasmid])
                    #reannotated_record = reannotate_genomes.prokka_reannotation([[plasmid]])
                    #accession = reannotated_record[0].annotations["accessions"][0]
                    #out_name = '%s_prokka_reannot.gbk' % accession
                    #with open(out_name, 'w') as f:
                    #    SeqIO.write(reannotated_record[0], f, 'genbank')



        if len(chromosome) > 0:

            '''

            Assume single chromosome bacteria.
            If multiple record founds, consider those as contigs.
            Contigs contatenation with 200 N between each (labelled as assembly_gap feature)

            '''

            #print '########## chromosome ###########'
            import concat_gbk

            if  chromosome[0].seq == 'N'*len(chromosome[0].seq):
                #print 'No sequences for %s, skipping! #################' % gbff_file
                continue
            if is_annotated(chromosome[0]):
                #print '## %s annotated (file: %s), %s contigs' % (chromosome[0].name, gbff_file, len(chromosome))
                print 'number of chromosomes:', len(chromosome)
                if len(chromosome) > 1:
                    merged_record = concat_gbk.merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]



                merged_record.description = clean_description(merged_record.description)
                merged_record = remove_record_taxon_id(merged_record)

                out_name = chromosome[0].name + '.gbk'
                with open(out_name, 'w') as f:
                    SeqIO.write(merged_record, f, 'genbank')
                cleaned_records.append(merged_record)
            else:
                chromosome_reannot = True
                #print '## %s NOT annotated (file: %s), %s contigs' % (chromosome[0].name, gbff_file, len(chromosome))
                #print 'Concatenating and reannotating gbk files...'
                merged_record = concat_gbk.merge_gbk(chromosome, filter_size=1000)
                #print '########### merged record ##############'
                #print merged_record
                reannotation_list.append(merged_record)
                #reannotated_record = reannotate_genomes.prokka_reannotation([[merged_record]])[0]
                #print '########### reannotated record #########'
                #print reannotated_record

                #out_name = reannotated_record.id
                #print '############ final out name ###############'
                #print out_name
                #with open(out_name, 'w') as f:
                #    SeqIO.write(reannotated_record, f, 'genbank')
            #cleaned_records.append(merged_record)
        if plasmid_reannot and not chromosome_reannot and len(chromosome) > 0:
            print "plasmid", plasmid_reannot
            print "chr", chromosome_reannot
            raise TypeError('combination of unannotated plasmid(s) and annotated chromosome')
        elif not plasmid_reannot and chromosome_reannot and len(plasmids) > 0:
            print "plasmid", plasmid_reannot
            print "chr", chromosome_reannot
            raise TypeError('combination of annotated plasmid(s) and unannotated chromosome')
        elif plasmid_reannot or chromosome_reannot:
            print 'record will be reannotated'

        else:
            out_name = gbff_file.split('.')[0] + '_merged.gbk'
            with open(out_name, 'w') as f:
                SeqIO.write(cleaned_records, f, 'genbank')

    #print '############ reannotation list #####################', len(reannotation_list)
    print "reannotating %s genomes" % len(reannotation_list)
    reannotated_genomes = reannotate_genomes.prokka_reannotation(reannotation_list)
    for reannotated_record in reannotated_genomes:
        #print reannotated_record
        out_name = reannotated_record.id + '_merged.gbk'
        out_name2 = reannotated_record.id + '.gbk'
        with open(out_name, 'w') as f:
            SeqIO.write(reannotated_record, f, 'genbank')
        with open(out_name2, 'w') as f:
            SeqIO.write(reannotated_record, f, 'genbank')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--gbff_files', type=str, help="input gbff file(s)", nargs="+")
    args = parser.parse_args()
    check_gbk(args.gbff_files)
