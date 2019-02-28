#!/usr/bin/env python

# Ublast
# Author: Trestan Pillonel (claire.bertelli[@]gmail.com)
# Date: 04.2015
# ---------------------------------------------------------------------------

import argparse
from Bio import Entrez
from shell_command import shell_command
import subprocess

Entrez.email = "claire.bertelli@chuv.ch"
def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def run_ublast(nb_hit, input, udb_list, evalue, output, out_q):
    import shell_command
    cmd='usearch8.0.1623_i86linux32 -maxhits %s -ublast %s -db %s ' \
        '-evalue %s -blast6out %s.tab -alnout %s_msa.fa ' \
        '-dbmatched %s_matched.fa'

    for udb in udb_list:
        print udb
        one_cmd = cmd % (nb_hit, input, udb, evalue, output, output, output)
        print one_cmd
        out, err, code = shell_command.shell_command(one_cmd)
        print "ok!"
        print out
        out_q.put(out)


def run_multiple_ublast(input, output, udb_list, evalue, nb_hit):
    import numpy as np
    import numpy
    from multiprocessing import Process, Queue
    import time
    from multiprocessing import cpu_count
    import os
    import MySQLdb



    '''
    :parameters for BLAST
    :return: blast results with taxonomic information
    scientific names and kingdom will only be retrieved if the taxid database has been installed locally ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    '''



    out_q = Queue()

    n_cpu = 2
    n_poc_per_list = int(numpy.ceil(len(udb_list)/float(n_cpu)))
    query_lists = _chunks(udb_list, n_poc_per_list)
    #print query_lists
    procs = []
    for one_udb_list in query_lists:
        proc = Process(target=run_ublast, args=(nb_hit, input, one_udb_list, evalue, output, out_q))
        procs.append(proc)
        proc.start()

    # Collect all results into a single result dict. We know how many dicts
    # with results to expect.
    ublast_result = []
    for i in range(n_cpu):
        ublast_result+=out_q.get()
    print "join proc"
    time.sleep(5)


    # Wait for all worker processes to finish
    for proc in procs:
        proc.join()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input', type=str, help="Fasta or multifasta sequence")
    parser.add_argument("-o", '--output', type=str, help="Prefix for the output file")
    parser.add_argument("-d", '--udbs', type=str, help="BLAST database to use", nargs='+')
    parser.add_argument("-e", '--evalue', default=0.00001, type=str, help="E-value cutoff to report results (default=10-5)")
    parser.add_argument("-n", '--nb_hit', default=10, type=int, help="Number of hits to keep (default=10)")


    args = parser.parse_args()

    print "Performing ublast..."
    run_multiple_ublast(args.input, args.output, args.udbs, args.evalue, args.nb_hit)
    print "UBLAST result now available"