#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create html table from blast tabulated file
# headers: accession	size	gi	n proteins	n contigs	gc 	description
# add 4 columns with links of the form /assets/chlamdb/ffn/ for gbk/faa/ffm/fna
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2015
# ---------------------------------------------------------------------------


def blast_tab2htm(input_file):

    complete_file = ''

    with open(input_file, "r") as f:
        i = 0
        complete_file += '<html>\n'
        complete_file +='<body>\n'
        complete_file +='{%% include "chlamdb/user.html" %%}\n'

        complete_file +='<h2>Blast nr hit of <a href="/chlamdb/locusx/chlamydia_03_15/%s">%s<a></h2>\n'

        complete_file += '<table cellspacing="0" border="0">\n'
        complete_file += '<tr>\n'
        complete_file += '    <th>Subject</th>\n'

        complete_file += '    <th>Kingdom</th>\n'
        complete_file += '    <th>e-value</th>\n'
        complete_file += '    <th>n identical</th>\n'
        complete_file += '    <th>percent identity</th>\n'
        complete_file += '    <th>n positive</th>\n'
        complete_file += '    <th>n gaps</th>\n'
        complete_file += '    <th>length</th>\n'
        complete_file += '    <th>query start</th>\n'
        complete_file += '    <th>query end</th>\n'
        complete_file += '    <th>query cov.</th>\n'
        complete_file += '    <th>subject start</th>\n'
        complete_file += '    <th>subject end</th>\n'
        complete_file += '    <th>subject title</th>\n'
        complete_file += '    <th>Taxonomy</th>\n'
        complete_file += '    </tr>\n'

        for blast_hit in f:
            line = blast_hit.rstrip().split('\t')
            # #qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
            query_accession = line[1].split("|")[1]
            subject_accession = line[3]
            query_gi = line[0]
            subject_gi = line[2]
            subject_scientific_names = line[4].split(';')
            subject_kingdom = line[5]
            subject_taxids = line[6].split(';')
            evalue = line[7]
            n_identical = line[8]
            percent_identity = line[9]
            positive = line[10]
            gaps = line[11]
            length = line[12]
            query_start = line[13]
            query_end = line[14]
            query_cov = line[15]
            subject_start = line[16]
            subject_end = line[17]
            subject_strand = line[18]
            subject_title = line[19]

            complete_file += '<tr>\n'
            complete_file += '    <td><a href="http://www.ncbi.nlm.nih.gov/protein/%s">%s<a></td>\n' % (subject_accession, subject_accession)

            all_taxonomy = ''
            for taxon, name in zip(subject_taxids, subject_scientific_names):
                all_taxonomy += '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s">%s<a>,' % (taxon, name)

            complete_file += '    <td>%s</td>\n' % subject_kingdom
            complete_file += '    <td>%s</td>\n' % evalue
            complete_file += '    <td>%s</td>\n' % n_identical
            complete_file += '    <td>%s</td>\n' % percent_identity
            complete_file += '    <td>%s</td>\n' % positive
            complete_file += '    <td>%s</td>\n' % gaps
            complete_file += '    <td>%s</td>\n' % length
            complete_file += '    <td>%s</td>\n' % query_start
            complete_file += '    <td>%s</td>\n' % query_end
            complete_file += '    <td>%s</td>\n' % query_cov
            complete_file += '    <td>%s</td>\n' % subject_start
            complete_file += '    <td>%s</td>\n' % subject_end
            complete_file += '    <td>%s</td>\n' % subject_title.split('[')[0]
            complete_file += '    <td>%s</td>\n' % all_taxonomy
            complete_file += '</tr>\n'
        complete_file += '</table>\n'


        with open(query_accession  + '.html', 'w') as f:
            f.write(complete_file % (query_accession, query_accession))







if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input tabulated file")


    args = parser.parse_args()

    blast_tab2htm(args.input_tab)


