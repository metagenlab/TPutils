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




def rename_source(record):

    if 'strain' in record.features[0].qualifiers:
        if ';' in record.features[0].qualifiers['strain'][0]:
            print 'ACHRTUNG: record has 2 strain names! \t%s\t --> check and edit source manually' % record.name
            # put everythink lower size
            strain = record.features[0].qualifiers['strain'][0].split(';')[1]
        else:
            strain = record.features[0].qualifiers['strain'][0]
        if strain == 'strain':
            return False
        if strain.lower() not in record.annotations['source'].lower():
            msg = '%s' % record.annotations['source']
            print "ACHTUNG changing source\t%s\t--> %s " % (msg, record.annotations['source'] + strain)

        return strain, "%s %s" % (record.annotations['source'], strain)
    else:
        return False, False

def remove_record_taxon_id(record):

    if record.features[0].type == 'source':
        # delete evendual taxon_id (taxon id will be reattributed in the db based on the organism name)
        # we need unique SOURCE to have unique taxon_id in the biosql database
        # check
        #try:

        if 'db_xref' in record.features[0].qualifiers:
            for item in record.features[0].qualifiers['db_xref']:
                if 'taxon' in item:
                    index = record.features[0].qualifiers['db_xref'].index(item)
                    record.features[0].qualifiers['db_xref'].pop(index)
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

    # count the number opf identical source names
    source2count = {}
    accession2count = {}
    for gbff_file in gbff_files:
        #print gbff_file
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
                    strain, new_source = rename_source(plasmid)
                    print 'new source:', new_source
                    if new_source:
                        if strain.lower() not in record.annotations['source'].lower():
                            record.description = new_source
                        if strain.lower() not in record.annotations['organism'].lower():
                            record.annotations['organism'] = new_source
                        if strain.lower() not in record.annotations['source'].lower():
                            record.annotations['source'] = new_source
                    else:
                        print 'ACHTUNG\t no strain name for \t%s\t, SOUCE uniqueness should be checked manually' % gbff_file
                    # check if accession is meaningful
                    if 'NODE_' in record.id or 'NODE_' in record.name:
                        print 'ACHTUNG\t accession probably not unique (%s) for \t%s\t --> should be checked manually' % (merged_record.id,
                                                                                                                        gbff_file)

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

            if chromosome[0].seq == 'N'*len(chromosome[0].seq):
                #print 'No sequences for %s, skipping! #################' % gbff_file
                continue
            if is_annotated(chromosome[0]):
                #print '## %s annotated (file: %s), %s contigs' % (chromosome[0].name, gbff_file, len(chromosome))
                #print 'number of chromosomes:', len(chromosome)
                if len(chromosome) > 1:
                    merged_record = concat_gbk.merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]



                merged_record.description = clean_description(merged_record.description)

                # check if source is unique (include strain name)
                merged_record = remove_record_taxon_id(merged_record)
                strain, new_source = rename_source(merged_record)
                print 'new source:', new_source
                if new_source:
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.description = new_source
                    if strain.lower() not in merged_record.annotations['organism'].lower():
                        merged_record.annotations['organism'] = new_source
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.annotations['source'] = new_source
                else:
                    print 'ACHTUNG\t no strain name for \t%s\t --> SOUCE uniqueness should be checked manually' % gbff_file
                # check if accession is meaningful
                if 'NODE_' in merged_record.id or 'NODE_' in merged_record.name:
                    print 'ACHTUNG\t accession probably not unique (%s) for \t%s\t --> should be checked manually' % (merged_record.id,
                                                                                                                    gbff_file)
                out_name = chromosome[0].name + '.gbk'
                with open(out_name, 'w') as f:
                    SeqIO.write(merged_record, f, 'genbank')
                cleaned_records.append(merged_record)
            else:
                # unannotated record: merge genbank and keep it in memory
                chromosome_reannot = True
                merged_record = concat_gbk.merge_gbk(chromosome, filter_size=1000)
                reannotation_list.append(merged_record)

            # count the source to identifiy redundant sources
            if merged_record.annotations['source'] not in source2count:
                source2count[merged_record.annotations['source']] = 1
            else:
                source2count[merged_record.annotations['source']] += 1


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

    # check if source list is non redundant
    for one_source in source2count:
        if source2count[one_source] > 1:
            print '%s\tsource is redundant (found %sx), add a unique source' % (one_source, source2count[one_source])



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--gbff_files', type=str, help="input gbff file(s)", nargs="+")
    args = parser.parse_args()
    check_gbk(args.gbff_files)
