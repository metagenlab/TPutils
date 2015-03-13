#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to ffn
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------



def gbk2ffn(seq_records, outname):
    output_handle = open(outname, "w")

    for record in seq_records:
        for seq_feature in record.features:
            if seq_feature.type == "CDS":
                #print seq_feature
                try:
                    len(seq_feature.qualifiers['translation'])
                except:
                    print seq_feature
                    print "pseudogene?"
                    continue
                #assert len(seq_feature.qualifiers['translation'])==1
                # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
                try:
                    output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                            seq_feature.qualifiers["db_xref"][0].split(":")[1],
                            seq_feature.qualifiers["protein_id"][0],
                            seq_feature.qualifiers["product"][0],
                            record.description,
                            seq_feature.extract(record.seq)))
                except:
                    try:

                        output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                            seq_feature.qualifiers["db_xref"][0].split(":")[1],
                            seq_feature.qualifiers["protein_id"][0],
                            #seq_feature.qualifiers["note"][0],
                            record.description,
                            seq_feature.extract(record.seq)))
                    except:
                        output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                            seq_feature.qualifiers["locus_tag"][0],
                            seq_feature.qualifiers["locus_tag"][0],
                            record.description,
                            seq_feature.extract(record.seq)))

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()

    if not args.outname:
        outname = args.input_gbk.split(".")[0]+".faa"
    else:
        outname = args.outname

    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    gbk2ffn(seq_records, outname)
