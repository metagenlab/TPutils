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
    
    cmd = 'newick2phyloxml.pl -i %s' % input_tree

    file_name = input_tree.split('.')[0]
    
    a,b,c = shell_command.shell_command(cmd)
    print a,b,c
    
    dirpath = os.getcwd()

    input_file = os.path.join(dirpath, file_name + '.phyloxml')
    output_file = os.path.join(dirpath, file_name + '_renamed.tree')

    print input_file
    print output_file

    if accession2taxon:
        manipulate_trees.convert_tree_accession2taxon_id(biodb_name,
                                               input_file,
                                               output_file)

    if taxon2accession:
        manipulate_trees.convert_tree_taxon_id2accession(biodb_name,
                                               input_file,
                                               output_file)



    else:
        manipulate_trees.convert_tree_taxon2genome(biodb_name,
                                               input_file,
                                               output_file)

def convert_leaf_labels_from_genbank(input_tree, input_gbk_list):
    import gbk2accessiontodefinition
    import parse_newick_tree


    id2description = gbk2accessiontodefinition.get_coressp(input_gbk_list)
    new_tree = parse_newick_tree.convert_terminal_node_names(input_tree, id2description, 'newick')

    return new_tree
    
    

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    from Bio import Phylo
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tree', type=str, help="input tree (newick)")
    parser.add_argument("-d", '--database_name', type=str, help="corresponding database name")
    parser.add_argument("-g", '--genbank_files', type=str, help="genbank files of all leaf nodes", nargs='+')
    parser.add_argument("-a", '--accession2taxon', action='store_true', help="convert accession to taxon_id used in biosqldb")
    parser.add_argument("-t", '--taxon2accession', action='store_true', help="convert taxon id to accessions used in biosqldb")

    args = parser.parse_args()
    if args.database_name:
        convert_leaf_labels(args.input_tree,
                            args.database_name,
                            accession2taxon=args.accession2taxon,
                            taxon2accession=args.taxon2accession)
    if args.genbank_files:
        new_tree = convert_leaf_labels_from_genbank(args.input_tree, args.genbank_files)
        file_name = args.input_tree.split('.')[0]
        output_file = file_name + '_renamed.tree'
        Phylo.write(new_tree, output_file, 'newick')
