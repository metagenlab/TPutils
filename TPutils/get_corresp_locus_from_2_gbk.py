#!/usr/bin/env python

def compare_gbk(file_a, file_b):

    from Bio import SeqIO
    handle_a = open(file_a, 'r')
    handle_b = open(file_b, 'r')

    records_a = [i for i in SeqIO.parse(handle_a, 'genbank')]
    records_b = [i for i in SeqIO.parse(handle_b, 'genbank')]
    new_record = records_b
    for record in records_a:
        new_record = records_b[0]

        for ref_feature in record.features:
            match_count = 0
            if 'pseudo' in ref_feature.qualifiers:
                continue
            if 'pseudogene' in ref_feature.qualifiers:
                continue
            if ref_feature.type=='CDS':
                ref_locus_tag = ref_feature.qualifiers['locus_tag'][0]
                #print ref_locus_tag
            else:
                continue

            identical_CDS = 0
            match = False
            # count number of identical features (exact same location)
            for n_rec, record_b in enumerate(records_b):
                for n, new_feature in enumerate(record_b.features):
                    if 'pseudo' in new_feature.qualifiers:
                        continue
                    if 'pseudogene' in new_feature.qualifiers:
                        continue
                    if new_feature.type == 'CDS':
                        #print(ref_feature)
                        #print(new_feature)
                        if ref_feature.qualifiers['translation'] == new_feature.qualifiers['translation']:
                            print ("%s\t%s\t%s\t%s" % (record.id,ref_locus_tag,new_feature.qualifiers['locus_tag'][0], match_count))
                            identical_CDS +=1
                            match_count+=1
                            #new_record.features[n].qualifiers['locus_tag'] = ref_locus_tag
                            match = True

                            break
            if not match:
                print ('%s\t%s\tNOMATCH' % (record.id, ref_locus_tag))
                        
    with open('corresp.gbk', 'w') as tt:
        SeqIO.write(new_record,tt,'genbank')

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()

	parser.add_argument("-a", '--gbk1', type=str, help="gbk 1")
	parser.add_argument("-b", '--gbk2', type=str, help="gbk 2")

	args = parser.parse_args()	
	compare_gbk(args.gbk1, args.gbk2)
