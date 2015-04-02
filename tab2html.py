#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def tab2htm(input_file):
    with open(input_file, "r") as f:
        i = 0
        print '<table cellspacing="0" border="0">'

        for line in f:
            if i == 0:
                header = line.rstrip().split("\t")
                print '<tr>'
                print '    <th>%s</th>' % header[0]
                print '    <th>%s</th>' % header[1]
                print '    <th>%s</th>' % header[3]
                print '    <th>%s</th>' % header[4]
                print '    <th>%s</th>' % header[5]
                print '    <th>%s</th>' % header[6]
                print '    <th colspan="4">Download</th>'
                print '    </tr>'

                i+=1
                continue
            data = line.rstrip().split("\t")
            print '<tr>'
            print '    <td><a href="http://www.ncbi.nlm.nih.gov/nuccore/%s">%s<a></td>' % (data[2], data[0])
            print '    <td>%s</td>' % data[1]
            print '    <td>%s</td>' % data[3]
            print '    <td>%s</td>' % data[4]
            print '    <td>%s</td>' % data[5]
            print '    <td>%s</td>' % data[6]

            #print '<td colspan="6"><table width="800" border=0  class=table_genomes>'

            gbk = '<td style="vertical-align:top" nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/gbk/%s.gbk">' \
            '<span style="background-color:#64FE2E;color:black;padding:1px;">GBK</span></a>' \
            '</td>' % (data[0].split(".")[0])

            faa = '<td snowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/faa/%s.faa">' \
            '<span style="background-color:#64FE2E;color:black;padding:1px;">FAA</span></a>' \
            '</td>' % (data[0].split(".")[0])

            fna = '<td nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/fna/%s.fna">' \
            '<span style="background-color:#64FE2E;color:black;padding:1px;">FNA</span></a>' \
            '</td>' % (data[0].split(".")[0])

            ffn = '<td nowrap>' \
            '<a  style="text-decoration:none;" href="/assets/chlamdb/ffn/%s.ffn">' \
            '<span style="background-color:#64FE2E;color:black;padding:1px;">FFN</span></a>' \
            '</td>' % (data[0].split(".")[0])
            #print '<tr>'
            print gbk
            print faa
            print fna
            print ffn
            #print '</tr>'
            #print '</table></td>'
            #print '</tr>'

    print '</table>'

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input tabulated file")


    args = parser.parse_args()

    tab2htm(args.input_tab)


