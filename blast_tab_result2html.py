#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create html table from blast tabulated file
# headers: accession	size	gi	n proteins	n contigs	gc 	description
# add 4 columns with links of the form /assets/chlamdb/ffn/ for gbk/faa/ffm/fna
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2015
# ---------------------------------------------------------------------------


def blast_tab2htm(input_file):



    with open(input_file, "r") as f:

        header_file = ''
        header_file += '<html>\n'
        header_file +='<body>\n'
        header_file +='<h2>Blast nr hit of <a href="{%% url \'locusx\' chlamydia_03_15 %s True %%}">%s<a></h2>\n'

        header_file += '<table id=blast_nr_table class="sortable">\n'
        header_file += '<tr>\n'
        header_file += '    <th id=blast_rank></th>\n'
        header_file += '    <th id=blast_subject>Subject</th>\n'

        header_file += '    <th id=blast_kingdom>Kingdom</th>\n'
        header_file += '    <th id=blast_evalue>eval.</th>\n'
        header_file += '    <th id=blast_percent_id>ID(%%)</th>\n'
        header_file += '    <th id=blast_n_id>N id</th>\n'
        header_file += '    <th id=blast_n_positive>N pos.</th>\n'
        header_file += '    <th id=blast_gaps>N gaps</th>\n'
        header_file += '    <th id=blast_length>Len.</th>\n'
        header_file += '    <th id=blast_query_start>Q. start</th>\n'
        header_file += '    <th id=blast_query_end>Q. end</th>\n'
        header_file += '    <th id=blast_query_cov>Q. cov.</th>\n'
        header_file += '    <th id=blast_subject_start>S. start</th>\n'
        header_file += '    <th id=blast_subject_end>S. end</th>\n'
        header_file += '    <th id=blast_subject_title>S. title</th>\n'
        header_file += '    <th id=blast_taxonomy>Taxonomy</th>\n'
        header_file += '    </tr>\n'


        input_file = [i.rstrip().split('\t') for i in f]

        complete_file = ''
        i = 0
        for n, line in enumerate(input_file):
            # #qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle

            query_accession = line[1].split("|")[3]

            if n == 0:
                i+=1
                pass
            elif line[1] != input_file[n-1][1]:
                out_accession = input_file[n-1][1].split("|")[3]
                f = open(out_accession + '.html', 'w')
                out_file = header_file + complete_file + '</table>\n'

                f.write(out_file % (out_accession, out_accession))
                f.close()
                complete_file = ''
                i = 1
            else:
                if line[3] != input_file[n-1][3]:
                    i+=1
                pass

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
            complete_file += '    <td>%s</td>\n' % i
            complete_file += '    <td><a href="http://www.ncbi.nlm.nih.gov/protein/%s">%s<a></td>\n' % (subject_accession, subject_accession)

            all_taxonomy = ''
            for taxon, name in zip(subject_taxids, subject_scientific_names):
                all_taxonomy += '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s">%s<a>,' % (taxon, name)

            complete_file += '    <td>%s</td>\n' % subject_kingdom
            complete_file += '    <td>%s</td>\n' % evalue
            complete_file += '    <td>%s</td>\n' % percent_identity
            complete_file += '    <td>%s</td>\n' % n_identical
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


        f = open(query_accession + '.html', 'w')
        complete_file = header_file + complete_file + '</table>\n'
        f.write(complete_file % (query_accession, query_accession))
        f.close()








if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input tabulated file")


    args = parser.parse_args()

    blast_tab2htm(args.input_tab)


