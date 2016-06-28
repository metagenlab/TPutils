#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

def gbk2fna(seq_records, outname):
    output_handle = open(outname, "w")

    if len(seq_records) > 1:
        for record in seq_records:
            output_handle.write(">%s\n%s\n" % (record.name, record.seq))
    else:
        try:
            gi = seq_records[0].annotations["gi"]
            source = seq_records[0].annotations["source"]
            output_handle.write(">%s gi|%s|%s\n%s\n" % (seq_records[0].name, gi, source, seq_records[0].seq))
        except:
            output_handle.write(">%s\n%s\n" % (seq_records[0].name, seq_records[0].seq))

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)


    args = parser.parse_args()
    input_handle = open(args.input_gbk, "rU")
    seq_records = list(SeqIO.parse(input_handle, "genbank"))
    if not args.outname:
        #outname = args.input_gbk.split(".")[0]+".fna"
        outname = seq_records[0].id.split('.')[0] + ".fna"
    else:
        outname = args.outname


    gbk2fna(seq_records, outname)
