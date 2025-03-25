#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------



def parse_blast_outfmt6(*blast_result_files):
    blast_result = {}
    for one_blast_file in blast_result_files:
        #print one_blast_file
        input_handle = open(one_blast_file, "rU")
        #qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
        for line in input_handle:
            line = line.rstrip().split('\t')
            if len(line) == 1:
                return None
            data = {}

            query_accession = line[1].split("|")[1]
            subject_accession = line[3]

            data["query_gi"] = line[0]
            data["subject_gi"] = line[2]

            data["subject_scientific_name"] = line[4]
            data["subject_kingdom"] = line[5]
            data["subject_taxid"] = line[6]
            data["evalue"] = line[7]
            data["n_identical"] = line[8]
            data["percent_identity"] = line[9]
            data["positive"] = line[10]
            data["gaps"] = line[11]
            data["length"] = line[12]
            data["query_start"] = line[13]
            data["query_end"] = line[14]
            data["query_cov"] = line[15]
            data["subject_start"] = line[16]
            data["subject_end"] = line[17]
            data["subject_strand"] = line[18]
            data["subject_title"] = line[19]

            if query_accession not in blast_result:
                blast_result[query_accession] = {}

            if subject_accession not in blast_result[query_accession]:

                blast_result[query_accession][subject_accession] = [data]
            else:
                blast_result[query_accession][subject_accession].append(data)

    return blast_result

def blast_dico2n_blast_hits(blast_dico):
    accession2n_blast_hits = {}
    for n, accession in enumerate(blast_dico.keys()):
        print (accession, n)
        accession2n_blast_hits[accession] = len(blast_dico[accession])
    return accession2n_blast_hits

def blast_dico2n_eukaryote_hits(blast_dico):
    accession2n_eukaryote_hits = {}
    for n, accession in enumerate(blast_dico.keys()):

        # for each hit, check if part of the eukaryote kingdom, count
        euk_count = 0
        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            if blast_dico[accession][one_hit_accession][0]['subject_kingdom'] == 'Eukaryota':
                euk_count+=1
        accession2n_eukaryote_hits[accession] = euk_count
    return accession2n_eukaryote_hits

def blast_dico2n_bacterial_hits(blast_dico):
    accession2n_bacterial_hits = {}
    for n, accession in enumerate(blast_dico.keys()):

        # for each hit, check if part of the eukaryote kingdom, count
        bact_count = 0
        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            if blast_dico[accession][one_hit_accession][0]['subject_kingdom'] == 'Bacteria':
                bact_count+=1
        accession2n_bacterial_hits[accession] = bact_count
    return accession2n_bacterial_hits


def blast_dico2classification(blast_dico):
    accession2n_bacterial_hits = {}
    for n, accession in enumerate(blast_dico.keys()):

        # for each hit, check if part of the eukaryote kingdom, count
        bact_count = 0
        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            if blast_dico[accession][one_hit_accession][0]['subject_kingdom'] == 'Bacteria':
                bact_count+=1
        accession2n_bacterial_hits[accession] = bact_count
    return accession2n_bacterial_hits


def blast_dico2taxonomy_info(blast_dico):
    import sequence_id2scientific_classification

    accession2hits_classification = {}
    for n, accession in enumerate(blast_dico.keys()):
        print (n, accession)
        accession2hits_classification[accession] = {}
        all_taxon_ids = []
        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            taxon_ids = blast_dico[accession][one_hit_accession][0]['subject_taxid'].split(";")
            all_taxon_ids+=taxon_ids


        print ("ids:", all_taxon_ids)

        taxon_id2classification = sequence_id2scientific_classification.taxon_id2scientific_classification(all_taxon_ids)

        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            accession2hits_classification[accession][one_hit_accession] = []
            taxon_ids = blast_dico[accession][one_hit_accession][0]['subject_taxid'].split(";")
            for taxon in taxon_ids:
                try:
                    accession2hits_classification[accession][one_hit_accession].append(taxon_id2classification[taxon])
                except KeyError:
                    accession2hits_classification[accession][one_hit_accession].append(None)

    return accession2hits_classification

def get_complete_classification(taxon, id2parent, taxid2rank, id2scientific, target_rank="order"):
    """
    superkingdom
    tribe
    subgenus
    family
    species subgroup
    species group
    phylum
    superclass
    subspecies
    species
    no rank
    superorder
    infraorder
    subclass
    superphylum
    kingdom
    subtribe
    subphylum
    subkingdom
    forma
    infraclass
    varietas
    subfamily
    class
    superfamily
    parvorder
    suborder
    genus
    order
    """


    import parseTaxonomy

    taxons = str(taxon).split(";")


    n = 1
    for one_taxon in taxons:

        rank = taxid2rank[one_taxon]
        print (n, "Taxon rank:", rank)
        print (n, "Name", id2scientific[one_taxon])


        target_ids = parseTaxonomy.find_rank_recursive(one_taxon, target_rank, id2parent, taxid2rank)
        print (n, "target IDs", target_ids)
        n+=1



def investigate_classification(blast_dico, classification_dico):
    order_count = {}
    sum_order_count = {}
    chlamydiales_count = {}
    non_chlamydiales_count = {}
    for n, accession in enumerate(blast_dico.keys()):
        order_count[accession] = {}
        for n, one_hit_accession in enumerate(blast_dico[accession].keys()):
            taxon_ids = blast_dico[accession][one_hit_accession][0]['subject_taxid'].split(";")
            #for taxon in taxon_ids:
            #print classification_dico[accession][one_hit_accession][0]
            try:
                if classification_dico[accession][one_hit_accession][0]['superkingdom'] == 'Bacteria':
                    try:
                        order = classification_dico[accession][one_hit_accession][0]['order']

                        # check if chlamydiales
                        #print order
                        if order == 'Chlamydiales':
                            if accession not in chlamydiales_count:
                                chlamydiales_count[accession] = 1
                            else:
                                chlamydiales_count[accession] +=1
                        else:
                            if accession not in non_chlamydiales_count:
                                non_chlamydiales_count[accession] = 1
                            else:
                                non_chlamydiales_count[accession] +=1


                        # contabilisation totale des ordres
                        if order not in sum_order_count:
                            sum_order_count[order] =1
                        else:
                            sum_order_count[order] +=1

                        # comptabilisation pour chaque accession
                        if order not in order_count[accession]:
                            order_count[accession][order] = 1
                        else:
                            order_count[accession][order] += 1
                    # si rank order pas disponible, ajouter à no rank
                    except KeyError:
                        # comptabilisation totale
                        if 'no rank' not in sum_order_count:
                            sum_order_count['no rank'] = 1
                        else:
                            sum_order_count['no rank'] += 1

                        # comptabilisation par accession
                        if 'no rank' not in order_count[accession]:
                            order_count[accession]['no rank'] = 1
                        else:
                            order_count[accession]['no rank'] += 1
                    except TypeError:
                        print ("NA")
            # si superkingdom pas disponible, ajouter no rank
            except:

                if 'no rank' not in sum_order_count:
                    sum_order_count['no rank'] =1
                else:
                    sum_order_count['no rank'] +=1

                if 'no rank' not in order_count[accession]:
                    order_count[accession]['no rank'] = 1
                else:
                    order_count[accession]['no rank'] += 1



    return order_count, sum_order_count, chlamydiales_count, non_chlamydiales_count



if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input fasta tab file", nargs='+')
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-t", '--taxonomy_dico', type=str, help="taxonomy dico", default=False)
    parser.add_argument("-g", '--get_all_hit_order_classif', type=str, help="get_all_hit_order_classif", default=False)
    args = parser.parse_args()


    all_blast_results = parse_blast_outfmt6(*args.input_tab)

    '''
    accession2hits_classification = blast_dico2taxonomy_info(all_blast_results)


    with open('accession2hits_classification.json', 'w') as f:
        json.dump(accession2hits_classification, f)
    '''

    with open(args.taxonomy_dico) as f:
        taxid2classification=json.load(f)

    classif, sum_order_count, chlamydiales_count, non_chlamydiales_count = investigate_classification(all_blast_results, taxid2classification)

    print (chlamydiales_count)
    print (non_chlamydiales_count)
    '''
    for accession in classif:
        for hit in classif[accession]:
            print accession, hit, classif[accession][hit]

    '''

    if args.get_all_hit_order_classif:
        for n, order in enumerate(sum_order_count):
            print ('%s\t%s\t%s' %(n, order, sum_order_count[order]))






    #print all_blast_results.keys()
    #print all_blast_results[all_blast_results.keys()[0]]

    #print blast_dico2n_blast_hits(all_blast_results)
    #blast_dico2n_eukaryote_hits(all_blast_results)

