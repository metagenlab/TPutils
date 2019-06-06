#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------

from Bio import Entrez

Entrez.email = "trestan.pillonel@unil.ch"



def taxon_id2taxonomy(taxon_id):
    # term=txid813[Organism:exp]
    handle = Entrez.esearch(db="nucleotide", term="txid%s[Organism:exp]" % taxon_id, retmax=1)
    search_record = Entrez.read(handle)
    match_id = search_record['IdList'][0]
    print ('match!', match_id)
    handle_record = Entrez.efetch(db="nucleotide", id=match_id, rettype="gb", retmode="text")
    print (handle_record)
    records = [i for i in SeqIO.parse(handle_record, "genbank")]
    print (dir(records[0]))
    print (records[0].annotations)
    return [records[0].annotations['source'], records[0].annotations['taxonomy'], records[0].annotations['organism']]

def reformat_gbk(gbk_file, study,
                 publication_title,
                 publication_authors,
                 publication_journal,
                 locus_tag_prefix,
                 taxon_id,
                 scaffold_prefix,
                 strain,
                 plasmid=False,
                 locus_count_start=1):

    '''

    - remove protein_id
    - split scaffolds into contigs ==> name contigs contig_XXX
    - generate agp file

    :param gbk_file:
    :param study:
    :param publication:
    :param locus_tag_prefix:
    :param plasmid:
    :return:
    '''

    source, taxonomy, organism = taxon_id2taxonomy(taxon_id)

    print(source)
    print()
    print(taxonomy)
    print()

    new_records = []
    from Bio import SeqIO
    import copy
    import copy
    from Bio.SeqFeature import Reference
    from Bio.SeqFeature import FeatureLocation
    with open(gbk_file, 'r') as f:



        records = [i for i in SeqIO.parse(f, 'genbank')]
        #locus_count=1

        contig_records = []
        contig_count = 1

        for new_record in records:
            start = 0
            end = len(new_record.seq)
            print (dir(new_record))
            for feature in new_record.features:
                '''
                if feature.type == 'assembly_gap':
                    print 'GAP-------'
                    print feature
                    contig = new_record[start:int(feature.location.start)]
                    # update start location
                    start = int(feature.location.end)

                    # rename contig record LOCUS

                    contig.id = "contig_%s" % contig_count
                    contig.name = "contig_%s" % contig_count

                    contig_records.append(contig)
                    contig_count += 1
                '''
            contig = new_record[start:end]

            contig.id = "%s_%02d" % (scaffold_prefix, contig_count)
            contig.name = "%s_%02d" % (scaffold_prefix, contig_count)
            contig_records.append(contig)
            contig_count += 1

        for n, record in enumerate(contig_records):

            ref = Reference()
            ref.authors = publication_authors
            ref.journal = publication_journal
            ref.title = publication_title

            '''
            ref_seq = Refserence()
            ref.authors = "Trestan Pillonel"
            ref.journal = "RL   Submitted (09-APRIL-2019) to the INSDC."
            '''
            
            #print record
            #print dir(record)
            #print "id", record.id
            #print "name", record.name
            #print record.annotations
            #print record.description
            #print record.dbxrefs
            #record.id = ''
            record.annotations['source'] = source
            record.annotations['taxonomy'] = taxonomy
            record.annotations['organism'] = organism
            record.description = '%s %s scaffold_%s' % (organism,
                                                        strain,
                                                        n+1)

            if record.features[0].type != 'source':

                print ('NOT SOURCE-------------------')
                record.features = [copy.copy(record.features[0])] + record.features
                record.features[0].qualifiers = {}
                record.features[0].type = 'source'
                record.features[0].location = FeatureLocation(0, len(record.seq))
            else:
                print ('SOURCE!!!!!!!!!!!!!!!!')
            record.features[0].qualifiers['db_xref'] = ["taxon:%s" % taxon_id]
            record.features[0].qualifiers['mol_type'] = ["genomic DNA"]
            record.features[0].qualifiers['organism'] = ["%s" % organism]
            record.features[0].qualifiers['strain'] = ["%s" % strain]

            if plasmid:
                #     /mol_type="genomic DNA"
                #     /organism="Klebsiella pneumoniae"
                #     /strain="KpGe"
                #record.features[0].type = "source"
                #record.features[0].qualifiers['organism'] = ["Klebsiella pneumoniae"]
                #record.features[0].qualifiers['strain'] = ["KpGe"]
                record.features[0].qualifiers['plasmid'] = ["p%s" % strain]


            record.annotations['mol_type'] = ["genomic DNA"]
            ref.location = [record.features[0].location]
            #print 'location!', ref.location
            record.annotations['references'] = [ref]
            record.dbxrefs = ['BioProject:%s' % study]
            for i, feature in enumerate(record.features):
                if "protein_id" in feature.qualifiers:
                    del feature.qualifiers['protein_id']
                if feature.type == 'gene':

                    '''
                    if not plasmid:
                        locus = "%s_%05d" % (locus_tag_prefix, locus_count)
                    else:
                        print 'rename locus!', locus_tag_prefix
                        locus = "%s_p%04d" % (locus_tag_prefix, locus_count)
                    '''
                    locus = "%s_%05d" % (locus_tag_prefix, locus_count_start)
                    locus_count_start+=1
                    feature.qualifiers['locus_tag'] = locus
                    record.features[i+1].qualifiers['locus_tag'] = locus
            new_records.append(record)





    return new_records



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    import shell_command

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", required=True)
    parser.add_argument("-b", '--bioproject', type=str, help="bioproject", required=True)
    parser.add_argument("-t", '--title', type=str, help="title", required=True)
    parser.add_argument("-l", '--locus_tag', type=str, help="putput_name", required=True)
    parser.add_argument("-a", '--authors', type=str, help="authors", required=True)
    parser.add_argument("-j", '--journal', type=str, help="authors", required=True)
    parser.add_argument("-tx", '--taxon_id', type=str, help="taxon_id", required=True)
    parser.add_argument("-s", '--scaffold_prefix', type=str, help="prefix for scaffolds", required=True)
    parser.add_argument("-st", '--strain_name', type=str, help="strain name", required=True)
    parser.add_argument("-p", '--plasmid', action="store_true", help="plasmid?", required=False, default=False)
    parser.add_argument("-c", '--count_start', type=int,  help="count start for locus_tag ", required=False, default=1)

    args = parser.parse_args()


    new_record =reformat_gbk(gbk_file=args.input_gbk,
                             study=args.bioproject,
                             publication_title=args.title,
                             locus_tag_prefix=args.locus_tag,
                             publication_authors=args.authors,
                             publication_journal=args.journal,
                             taxon_id=args.taxon_id,
                             scaffold_prefix=args.scaffold_prefix,
                             strain=args.strain_name,
                             plasmid=args.plasmid,
                             locus_count_start=args.count_start)


    print()
    print()
    for record in new_record:
        print (record.annotations['references'][0].location)
    out = args.input_gbk.split('.')[0] + '_edit.embl'
    h = open(out, 'w')
    out2 = args.input_gbk.split('.')[0] + '_edit.gbk'
    h2 = open(out2, 'w')
    SeqIO.write(new_record, h, 'embl')
    SeqIO.write(new_record, h2, 'genbank')
    h.close()
    h2.close()
    # update header line
    a,b,c = shell_command.shell_command("sed -i 's/^ID.*/ID   XXX; XXX; linear; XXX; XXX; XXX; XXX./' %s" % out)
    # add star
    #a,b,c = shell_command.shell_command("sed -i 's/AC   /AC * /' %s" % out)
    #print a,b,c
