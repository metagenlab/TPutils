#!/usr/bin/env python
# -*- coding: utf-8 -*-


def rename_hmm(hmm_list):
    import re



    for hmm in hmm_list:
        out_name = hmm.split('.')[0] + '_rename.hmm'
        o = open(out_name, 'w')
        with open(hmm) as f:
            for line in f:
                name_r = re.compile('^NAME')
                if name_r.match(line):
                    try:
                        name = line.rstrip().split(' ')[2]
                        print 'name', name
                        new = 'NAME  %s\n' % hmm.split('.')[0]
                    except:
                        print 'ECHEC---------------', line
                    o.write(new)
                else:
                    o.write(line)

def rename_hmm2(hmm_list):
    import re



    for hmm in hmm_list:
        data = [line for line in open(hmm, 'r')]
        name_r = re.compile('^NAME')
        for line in data:
            if name_r.match(line):
                out_name = line.rstrip().split('|')[1] + '.hmm'
        o = open(out_name, 'w')
        for line in data:
            if name_r.match(line):
                name = line.rstrip().split('|')[1]
                print 'name', name
                new = 'NAME  %s\n' % name #hmm.split('.')[0]
                o.write(new)
            else:
                o.write(line)
                    
if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", '--hmm_list', help="input hmm files", nargs='+')

    args = parser.parse_args()
    rename_hmm2(args.hmm_list)
