#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def convert_table_labels(input_tree, biodb_name, accession2taxon=False, taxon2accession=False):
    import shell_command
    import manipulate_trees
    import os

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




def convert_table_labels_from_genbank(input_table, input_gbk_list):
    import gbk2accessiontodefinition
    import parse_newick_tree


    id2description = gbk2accessiontodefinition.get_coressp(input_gbk_list)

    new_table = []
    with open(input_table, 'r') as f:
	for i, row in enumerate(f):
	    if i == 0:
                header = row.rstrip().split()
		new_header = []
                for column in header:
		    try:
                    	new_header.append(id2description[column])
                    except KeyError:
			new_header.append(column)
                new_table.append('\t'.join(new_header)+'\n')
	    else:
		data = row.rstrip().split()
                try:
		    data[0] = id2description[data[0]]
                except KeyError:
                    pass
		new_table.append('\t'.join(data)+'\n')
    return new_table
    
    

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    from Bio import Phylo
    import os
    import ete2
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_table', type=str, help="input table")
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
        new_table = convert_table_labels_from_genbank(args.input_table, args.genbank_files)
        file_name = args.input_table.split('.')[0]
        output_file = file_name + '_renamed.tab'

        with open(output_file, 'w') as f:
	   for row in new_table:
              f.write(row ) 

        #with open(output_file, 'w') as f:
            #Phylo.write(new_tree, f, 'newick')
        #import ete2
        #t = ete2.Tree(output_file)
        #print t.write()

