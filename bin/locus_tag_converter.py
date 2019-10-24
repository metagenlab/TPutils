#!/usr/bin/python

def table_locus_converter(locus2locus, table_file):

    with open(table_file, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            for n,column in enumerate(data):
                try:
                    data[n] = locus2locus[column]
                except KeyError:
                    pass
            print '\t'.join(data)

def compare_gbk(file_a, file_b, table_file):

    from Bio import SeqIO
    handle_a = open(file_a, 'r')
    handle_b = open(file_b, 'r')

    records_a = [i for i in SeqIO.parse(handle_a, 'genbank')]
    records_b = [i for i in SeqIO.parse(handle_b, 'genbank')]
    locus2locus = {}
    for record in records_a:
        new_record = records_b[0]
        for ref_feature in record.features:

            if 'pseudo' in ref_feature.qualifiers:
               continue
            if ref_feature.type=='CDS':
                ref_locus_tag = ref_feature.qualifiers['locus_tag'][0]
                #print ref_locus_tag
            else:
                continue

            identical_CDS = 0
            match = False
            # count number of identical features (exact same location)
            for n, new_feature in enumerate(records_b[0].features):
                if new_feature.type == 'CDS':
                    if ref_feature.qualifiers['translation'] == new_feature.qualifiers['translation']:
                        print "%s\t%s\t%s" % (record.id,ref_locus_tag,new_feature.qualifiers['locus_tag'][0])
                        locus2locus[ref_locus_tag] = new_feature.qualifiers['locus_tag'][0]
                        identical_CDS +=1
                        new_record.features[n].qualifiers['locus_tag'] = ref_locus_tag
                        match = True
                        break
            if not match:
                #print '%s\t%s\t-' % (record.id, ref_locus_tag)
                locus2locus[ref_locus_tag] = '-'
    table_locus_converter(locus2locus, table_file)
    #with open('corresp.gbk', 'w') as tt:
    #    SeqIO.write(new_record,tt,'genbank')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", '--gbk1', type=str, help="gbk 1")
    parser.add_argument("-b", '--gbk2', type=str, help="gbk 2")
    parser.add_argument("-t", '--table', type=str, help="table")

    args = parser.parse_args()
    compare_gbk(args.gbk1, args.gbk2, args.table)

