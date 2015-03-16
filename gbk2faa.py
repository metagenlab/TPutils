#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def gbk2faa(seq_records, outname, format = False):
    import re
    output_handle = open(outname, "w")

    for record in seq_records:
        description = record.description
        description = re.sub(", complete genome\.", "", description)
        description = re.sub(", complete sequence\.", "", description)
        description = re.sub("strain ", "", description)
        description = re.sub("str\. ", "", description)
        description = re.sub(" complete genome sequence\.", "", description)
        description = re.sub(" complete genome\.", "", description)
        description = re.sub(" chromosome", "", description)
        description = re.sub(" DNA", "S.", description)

        for seq_feature in record.features:
            if seq_feature.type == "CDS":
                # check presence of a protein sequence
                try:
                    len(seq_feature.qualifiers['translation'])
                except:
                    print seq_feature
                    #sys.stderr.write(seq_feature.location.start)
                    #sys.stderr.write("pseudogene?")
                    continue
                #assert len(seq_feature.qualifiers['translation'])==1
                # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
                if format:
                    # output only locus tag and species name
                    try:
                        output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers["locus_tag"][0],
                                                description,
                                                seq_feature.qualifiers['translation'][0]))
                    except:
                        print seq_feature
                else:

                    try:
                        output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                seq_feature.qualifiers["product"][0],
                                description,
                                seq_feature.qualifiers['translation'][0]))
                    except:
                        try:

                            output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                seq_feature.qualifiers["protein_id"][0],
                                #seq_feature.qualifiers["note"][0],
                                description,
                                seq_feature.qualifiers['translation'][0]))
                        except:
                            output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                seq_feature.qualifiers["locus_tag"][0],
                                seq_feature.qualifiers["locus_tag"][0],
                                #seq_feature.qualifiers["note"][0],
                                description,
                                seq_feature.qualifiers['translation'][0]))


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-f", '--format', action='store_true', help="format header: >locus description", default=False)

    args = parser.parse_args()

    if not args.outname:
        outname = args.input_gbk.split(".")[0]+".faa"
    else:
        outname = args.outname

    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    gbk2faa(seq_records, outname, args.format)
