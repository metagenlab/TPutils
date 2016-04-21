#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

def parse_barrnap_output(barnap_out):

    #print barnap_out
    rrna_list = []
    for line in barnap_out.split("\n"):
        s = line.rstrip().split("\t")
   
        if len(s)<9:
            pass
        else:
            liste = [s[0], s[3], s[4], s[6], s[8]]
            rrna_list.append(liste)
    return rrna_list
  
#fabrique une liste avec les 16Srrna
def search_16S_rrna(processed_list):
   
    rrna_16S = []
    rrna_23S = []
    rrna_5S = []
    for line in processed_list:
        if re.match('Name=16S_rRNA;product=16S ribosomal RNA', line[-1]) is not None:
            rrna_16S.append(line)
            
        elif re.match('Name=23S_rRNA;product=23S ribosomal RNA', line[-1]) is not None:
            rrna_23S.append(line)
            
        elif re.match('Name=5S_rRNA;product=5S ribosomal RNA', line[-1]) is not None:
            rrna_5S.append(line)

        else:
            print "unknown sequence", line[-1]

    return (rrna_16S, rrna_23S, rrna_5S)


#trouve le plus long 16S
def find_longest_16S(rrna_16S):
    #print 'rrna_list', rrna_16S
    try:
        longest_rrna = rrna_16S[0]
    except:
        return None
    for i in range(1,len(rrna_16S)):
        length_of_16S = int(rrna_16S[i][2]) - int(rrna_16S[i][1])
        length_longest = int(longest_rrna[2]) - int(longest_rrna[1])
        if length_of_16S > length_longest:
            longest_rrna = rrna_16S[i]
    return longest_rrna



def extract_seq(fasta_file, fasta_header_name, start, stop, header = ""):
    for one_fasta_entry in SeqIO.parse(fasta_file, "fasta"):
        if fasta_header_name == one_fasta_entry.name:
            sequence_def = one_fasta_entry[int(start):int(stop)]
            print ">" + sequence_def.name + header + "\n" + sequence_def.seq




if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input', type=str,help="input fasta")
    parser.add_argument("-1", '--only_16s', action="store_true",help="get only 16s")
    parser.add_argument("-2", '--only_23s', action="store_true",help="get only 23s")


    args = parser.parse_args()


    stdout_str, stderr_str, runcode = shell_command.shell_command("barrnap %s" % args.input)

    if 'not found' in stderr_str:
        raise(Exception('Barrnap was not found on path, please install it'))
        import sys
        sys.exit()


    rrna_16S, rrna_23S, rrna_5S  = search_16S_rrna(parse_barrnap_output(stdout_str))

    longest16 = find_longest_16S(rrna_16S)
    longest23 = find_longest_16S(rrna_23S)


    
    
    if args.only_16s:
        if longest16:
            extract_seq(args.input, longest16[0], longest16[1], longest16[2], "_16S")
    elif args.only_23s:
        if longest23:
            extract_seq(args.input, longest23[0], longest23[1], longest23[2], "_23S")
    else:
        if longest16:
            extract_seq(args.input, longest16[0], longest16[1], longest16[2], "_16S")
        if longest23:
            extract_seq(args.input, longest23[0], longest23[1], longest23[2], "_23S")


