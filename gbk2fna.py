#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


def gbk2fna(seq_records, outname):
    output_handle = open(outname, "w")

    for record in seq_records:
        try:
            gi = record.annotations["gi"]
            source = record.annotations["source"]
            output_handle.write(">%s gi|%s|%s\n%s\n" % (record.name, gi, source, record.seq))

        except:
            output_handle.write(">%s\n%s\n" % (record.id, record.seq))
 

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
    gbk2fna(seq_records, outname)
