#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  

def gbk2faa(seq_records, outname):
    output_handle = open(outname, "w")

    for record in seq_records:
        for seq_feature in record.features :
            if seq_feature.type=="CDS" :
                #print seq_feature
                assert len(seq_feature.qualifiers['translation'])==1
                # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
                try:
                    output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                            seq_feature.qualifiers["db_xref"][0].split(":")[1],
                            seq_feature.qualifiers["protein_id"][0],
                            seq_feature.qualifiers["product"][0],
                            record.description,
                            seq_feature.qualifiers['translation'][0]))
                except:
                    output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                            seq_feature.qualifiers["db_xref"][0].split(":")[1],
                            seq_feature.qualifiers["protein_id"][0],
                            #seq_feature.qualifiers["note"][0],
                            record.description,
                            seq_feature.qualifiers['translation'][0]))



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk',type=str,help="input gbk file")
    parser.add_argument("-o", '--outname',type=str,help="putput_name", default = False)


    args = parser.parse_args()

    if not args.outname:
        outname = args.input_gbk.split(".")[0]+".faa"
    else:
        outname = args.outname

    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    gbk2faa(seq_records, outname)
