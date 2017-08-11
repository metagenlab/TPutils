#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

def gbk2summary(records):
    from Bio.SeqUtils import GC

    for record in records:
        gc_content = round(GC(record.seq),2)
        accession = record.id
        #print dir(record)
        #print record.dbxrefs
        try:
            assembly = [i for i in record.dbxrefs if 'Assembly' in i][0].split(':')[1]
        except:
            try:
                assembly = [i for i in record.dbxrefs if 'Assembly' in i][0]
            except:
                assembly = '-'
        description = record.description
        length= len(record)

        try:
            date = record.annotations['date']
        except:
            date = '-'
        #try:
        #    annot_date = record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Software revision']
        #    annot_version = record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Date']
        #except:
        #    annot_data = '-'
        #    annot_version = '-'
        #['structured_comment']['Genome-Annotation-Data']
        refseq_annotation_data = []
        taxo_data = '\t'.join(record.annotations['taxonomy'])
        n_features = len([i for i in record.features if i.type=='CDS'])
        try:
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (assembly, accession, description, gc_content,length, n_features, date, taxo_data)
        except:
            print '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (assembly, accession, description, gc_content,length, n_features, date)

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", nargs='+')
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    all_records = []
    for genbank in args.input_gbk:
        input_handle = open(genbank, "rU")
        seq_records = list(SeqIO.parse(input_handle, "genbank"))
        all_records+=seq_records



    gbk2summary(all_records)
