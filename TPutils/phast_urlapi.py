#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import time
import urllib2
import os
import urllib

def get_phast(accession):
    import genbank2refseq
    import re
    try:
        ready = False
        while not ready:
            if not os.path.exists(accession):
                os.mkdir(accession)
            f = urllib2.urlopen('http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?acc=%s' % accession)
            pattern = 'Result is ready for'
            page = f.read()
            if re.search(pattern, page):
                ready = True
            else:
                print ("not ready")
                time.sleep(2)

    except urllib2.HTTPError:
        ready = False
        while not ready:
            if not os.path.exists(accession):
                os.mkdir(accession)
            accession = genbank2refseq.genbank2refseq(accession)
            f = urllib2.urlopen('http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?acc=%s' % accession)
            pattern = 'Result is ready for'
            page = f.read()
            if re.search(pattern, page):
                ready = True
            else:
                print ("not ready")
                time.sleep(2)

    file_1 = 'http://phast.wishartlab.com/tmp/%s/summary.txt' % accession
    file_2 = 'http://phast.wishartlab.com/tmp/%s/detail.txt' % accession
    file_3 = 'http://phast.wishartlab.com/tmp/%s/image.png' % accession


    with open(os.path.join(accession, "summary.txt"), 'w') as f:
        f.write(urllib2.urlopen(file_1).read())
    with open(os.path.join(accession, "detail.txt"), 'w') as f:
        f.write(urllib2.urlopen(file_2).read())

    urllib.urlretrieve(file_3, os.path.join(accession, "image.png"))


def parse_phast_summary(phast_summary_file):
    import re
    with open(phast_summary_file, 'r') as f:
        rows = [row.strip() for row in f.readlines()]
        tag = False
        i = 0
        result =  {}
        for row in rows:
            if len(row) == 0:
                continue
            if "gi|" in row:
                genome = re.split('\s{2,}',row)
                accession = genome[0].split('|')[-2]
                genome_name = genome[-1]
                result[accession] = {}
            if "REGION" in row:
                headers = re.split('\s+', row)
                tag = True
                continue
            if tag:
                if i == 0:
                    i+=1
                    continue
                data = re.split('\s+', row)
                result[accession][int(data[0])] = {}
                result[accession][int(data[0])]["genome"] = genome_name
                for head, dat in zip(headers, data):
                    result[accession][int(data[0])][head] = dat
        return result

def format_phast_summary(dico):
    header = 'genome\t' \
             'accession\t' \
             'n\t'\
             'PHAGE+HYPO_PROTEIN_PERCENTAGE\t' \
             'BACTERIAL_PROTEIN_NUM\t' \
             'PHAGE_SPECIES_NUM\t' \
             'MOST_COMMON_PHAGE_NUM\t' \
             'COMPLETENESS(score)\t' \
             'REGION\t' \
             'PHAGE_HIT_PROTEIN_NUM\t' \
             'REGION_POSITION\t' \
             'REGION_LENGTH\t' \
             'TRNA_NUM\t' \
             'GC_PERCENTAGE\t' \
             'MOST_COMMON_PHAGE_PERCENTAGE\t' \
             'MOST_COMMON_PHAGE_NAME\t' \
             'HYPOTHETICAL_PROTEIN_NUM\t' \
             'SPECIFIC_KEYWORD\t' \
             'TOTAL_PROTEIN_NUM\t' \
             'ATT_SITE_SHOWUP\t'
    print (header)
    for accession in dico:
        print (accession)
        for n_phage in dico[accession]:
            #print n_phage
            phage_data = dico[accession][n_phage]
            #for i in phage_data:
            #    print i
            #print phage_data
            try:
                data = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (phage_data['genome'],
                accession,
                n_phage,
                phage_data['PHAGE+HYPO_PROTEIN_PERCENTAGE'],
                phage_data['BACTERIAL_PROTEIN_NUM'],
                phage_data['PHAGE_SPECIES_NUM'],
                phage_data['MOST_COMMON_PHAGE_NUM'],
                phage_data['COMPLETENESS(score)'],
                phage_data['REGION'],
                phage_data['PHAGE_HIT_PROTEIN_NUM'],
                phage_data['REGION_POSITION'],
                phage_data['REGION_LENGTH'],
                phage_data['TRNA_NUM'],
                phage_data['GC_PERCENTAGE'],
                phage_data['MOST_COMMON_PHAGE_PERCENTAGE'],
                phage_data['MOST_COMMON_PHAGE_NAME'],
                phage_data['HYPOTHETICAL_PROTEIN_NUM'],
                phage_data['SPECIFIC_KEYWORD'],
                phage_data['TOTAL_PROTEIN_NUM'],
                phage_data['ATT_SITE_SHOWUP'])
                print (data)
            except:
                #print "fai!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                data = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (phage_data['genome'],
                accession,
                n_phage,
                phage_data['PHAGE+HYPO_PROTEIN_PERCENTAGE'],
                phage_data['BACTERIAL_PROTEIN_NUM'],
                phage_data['PHAGE_SPECIES_NUM'],
                '-',
                phage_data['COMPLETENESS(score)'],
                phage_data['REGION'],
                phage_data['PHAGE_HIT_PROTEIN_NUM'],
                phage_data['REGION_POSITION'],
                phage_data['REGION_LENGTH'],
                phage_data['TRNA_NUM'],
                phage_data['GC_PERCENTAGE'],
                phage_data['FIRST_MOST_COMMON_PHAGE_PERCENTAGE'],
                phage_data['FIRST_MOST_COMMON_PHAGE_NUM'],
                phage_data['HYPOTHETICAL_PROTEIN_NUM'],
                phage_data['SPECIFIC_KEYWORD'],
                phage_data['TOTAL_PROTEIN_NUM'],
                phage_data['ATT_SITE_SHOWUP'])
                print (data)


'''
PHAGE+HYPO_PROTEIN_PERCENTAGE
MOST_COMMON_PHAGE_NAME(hit_genes_count)
BACTERIAL_PROTEIN_NUM
PHAGE_SPECIES_NUM
COMPLETENESS(score)
REGION
PHAGE_HIT_PROTEIN_NUM
REGION_POSITION
REGION_LENGTH
TRNA_NUM
GC_PERCENTAGE
genome
FIRST_MOST_COMMON_PHAGE_PERCENTAGE
FIRST_MOST_COMMON_PHAGE_NUM
HYPOTHETICAL_PROTEIN_NUM
SPECIFIC_KEYWORD
TOTAL_PROTEIN_NUM
ATT_SITE_SHOWUP

'''



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--ncbi_accessions', type=str, default=False, help="input ncbi accession(s)", nargs='+')
    parser.add_argument("-p", '--parse_summary_files', type=str, default=False, help="parse phase summary file(s)", nargs='+')

    args = parser.parse_args()

    if args.ncbi_accessions:
        for i in args.ncbi_accessions:
            get_phast(i)

    if args.parse_summary_files:
        all_summaries = {}
        for i in args.parse_summary_files:
            all_summaries.update(parse_phast_summary(i))
        format_phast_summary(all_summaries)
