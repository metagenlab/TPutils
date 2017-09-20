#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def convert_leaf_labels(input_tree, biodb_name, accession2taxon=False, taxon2accession=False):
    import shell_command
    import manipulate_trees
    import os

    print accession2taxon, taxon2accession

    cmd = 'newick2phyloxml.pl -i %s' % input_tree

    file_name = input_tree.split('.')[0]
    
    a,b,c = shell_command.shell_command(cmd)
    print a, b, c
    
    dirpath = os.getcwd()

    input_file = os.path.join(dirpath, file_name + '.phyloxml')
    output_file = os.path.join(dirpath, file_name + '_renamed.tree')

    print input_file
    print output_file

    if accession2taxon:
        manipulate_trees.convert_tree_accession2taxon_id(biodb_name,
                                               input_file,
                                               output_file)

    elif taxon2accession:
        manipulate_trees.convert_tree_taxon_id2accession(biodb_name,
                                               input_file,
                                               output_file)



    else:
        manipulate_trees.convert_tree_taxon2genome(biodb_name,
                                               input_file,
                                               output_file)

def convert_leaf_labels_from_genbank(input_tree,
                                     input_gbk_list,
                                     show_rank=False,
                                     use_gbk_file_names=False):
    import gbk2accessiontodefinition
    import parse_newick_tree


    id2description = gbk2accessiontodefinition.get_coressp(input_gbk_list, use_gbk_file_names=use_gbk_file_names)
    if show_rank:
        for id in id2description:
            print 'searching rank for %s...' % id
            try:

                id2description[id] = id2description[id] + ' (%s)' % accession2taxon_rank(id, 'phylum')
            except:
                print 'no phylum for %s' % id
                try:

                    id2description[id] = id2description[id] + ' (order: %s)' % accession2taxon_rank(id, 'order')
                except:
                    print 'no order for %s' % id
                    id2description[id] = id2description[id] + ' (?)'

    new_tree = parse_newick_tree.convert_terminal_node_names(input_tree, id2description, 1)

    return new_tree


def accession2taxon_rank(accession, rank='phylum'):
    from Bio import Entrez
    import sequence_id2scientific_classification

    Entrez.email = "trestan.pillonel@unil.ch"

    handle1 = Entrez.esearch(db="nuccore", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][0]

    handle2 = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=ncbi_id)
    record2 = Entrez.read(handle2)

    id = record2[0]['LinkSetDb'][0]['Link'][0]['Id']

    taxo_data = sequence_id2scientific_classification.taxon_id2scientific_classification([id])
    try:
        data = taxo_data[id][rank]
    except:
        print accession, ncbi_id,taxo_data[id]
        return False
    return data


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    from Bio import Phylo
    import os
    import ete2
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tree', type=str, help="input tree (newick)")
    parser.add_argument("-d", '--database_name', type=str, help="corresponding database name")
    parser.add_argument("-g", '--genbank_files', type=str, help="genbank files of all leaf nodes", nargs='+')
    parser.add_argument("-a", '--accession2taxon', action='store_true', help="convert accession to taxon_id used in biosqldb")
    parser.add_argument("-t", '--taxon2accession', action='store_true', help="convert taxon id to accessions used in biosqldb")
    parser.add_argument("-r", '--show_rank', action='store_true', help="show rank (phylum)")
    parser.add_argument("-f", '--tab_file', type=str, help="tabulated table with accession\tdescription", default=False)
    parser.add_argument("-gp", '--gbk_prefix', action='store_true', help="gbk file name were used as label for the newick tree")

    args = parser.parse_args()

    if args.database_name:
        print args.accession2taxon, args.taxon2accession
        convert_leaf_labels(args.input_tree,
                            args.database_name,
                            accession2taxon=args.accession2taxon,
                            taxon2accession=args.taxon2accession)
    if args.tab_file:
        import parse_newick_tree
        id2description = {}
        with open(args.tab_file, 'r') as f:
            for row in f:
                data = row.rstrip().split('\t')
                id2description[data[0]] = data[1]
        new_tree = parse_newick_tree.convert_terminal_node_names(args.input_tree, id2description, 1)
        file_name = args.input_tree.split('.')[0]
        output_file = file_name + '_renamed_tab.tree'
        new_tree.write(format=1, outfile=output_file)

    if args.genbank_files:



        new_tree = convert_leaf_labels_from_genbank(args.input_tree,
                                                    args.genbank_files,
                                                    show_rank = args.show_rank,
                                                    use_gbk_file_names=args.gbk_prefix)
        file_name = args.input_tree.split('.')[0]
        output_file = file_name + '_renamed.tree'

        new_tree.write(format=1, outfile=output_file)

        #with open(output_file, 'w') as f:
            #Phylo.write(new_tree, f, 'newick')
        #import ete2
        #t = ete2.Tree(output_file)
        #print t.write()

