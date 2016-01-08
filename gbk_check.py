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
                             
                            record.annotations['source'] += ' %s' % strain
                            record.annotations['organism'] = record.annotations['source']
                            print record.annotations['source']
                    else:
                        print 'ACHTUNG: no strain for %s, source uniqueness should be checked' % record.name
    else:
        print 'ACHRTUNG: no source for record %s' % record.name
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
    description = re.sub(" DNA", "S.", description)
    description = re.sub("Merged record from ", "", description)
    description = re.sub(", wgs", "", description)
    description = re.sub("Candidatus ", "", description)
    description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
    return description

def check_gbk(gbff_files):
    from Bio import SeqIO
    import reannotate_genomes
    reannotation_list = []

    for gbff_file in gbff_files:
        records = SeqIO.parse(open(gbff_file, "r"), "genbank")

        chromosome, plasmids = filter_plasmid(records)

        reannotated_records = []

        if len(plasmids) > 0:
            print '######### plasmid(s) ############'
            for plasmid in plasmids:
                if is_annotated(plasmid):
                    out_name = plasmid.name + '.gbk'
                    plasmid.description = clean_description(plasmid.description)
                    plasmid = remove_record_taxon_id(plasmid)
                    with open(out_name, 'w') as f:
                        SeqIO.write(plasmid, f, 'genbank')
                else:
                    reannotation_list.append([plasmid])
                    #reannotated_record = reannotate_genomes.prokka_reannotation([[plasmid]])
                    #accession = reannotated_record[0].annotations["accessions"][0]
                    #out_name = '%s_prokka_reannot.gbk' % accession
                    #with open(out_name, 'w') as f:
                    #    SeqIO.write(reannotated_record[0], f, 'genbank')

        if len(chromosome) > 0:
            print '########## chromosome ###########'
            import concat_gbk

            if  chromosome[0].seq == 'N'*len(chromosome[0].seq):
                print 'No sequences for %s, skipping! #################' % gbff_file
                continue
            if is_annotated(chromosome[0]):
                print '## %s annotated (file: %s), %s contigs' % (chromosome[0].name, gbff_file, len(chromosome))

                if len(chromosome) > 1:
                    merged_record = concat_gbk.merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]



                merged_record.description = clean_description(merged_record.description)
                merged_record = remove_record_taxon_id(merged_record)

                out_name = chromosome[0].name + '.gbk'
                with open(out_name, 'w') as f:
                    SeqIO.write(merged_record, f, 'genbank')
            else:
                print '## %s NOT annotated (file: %s), %s contigs' % (chromosome[0].name, gbff_file, len(chromosome))
                print 'Concatenating and reannotating gbk files...'
                merged_record = concat_gbk.merge_gbk(chromosome, filter_size=1000)
                print '########### merged record ##############'
                print merged_record
                reannotation_list.append(merged_record)
                #reannotated_record = reannotate_genomes.prokka_reannotation([[merged_record]])[0]
                #print '########### reannotated record #########'
                #print reannotated_record

                #out_name = reannotated_record.id
                #print '############ final out name ###############'
                #print out_name
                #with open(out_name, 'w') as f:
                #    SeqIO.write(reannotated_record, f, 'genbank')
    print '############ reannotation list #####################', len(reannotation_list)
    print reannotation_list
    reannotated_genomes = reannotate_genomes.prokka_reannotation(reannotation_list)
    for reannotated_record in reannotated_genomes:
        print reannotated_record
        out_name = reannotated_record.id + '.gbk'
        with open(out_name, 'w') as f:
            SeqIO.write(reannotated_record, f, 'genbank')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--gbff_files', type=str, help="input gbff file(s)", nargs="+")
    args = parser.parse_args()
    check_gbk(args.gbff_files)
