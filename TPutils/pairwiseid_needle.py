#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

#Produce identity matrix from fasta with multiple (non aligned) sequences
#Alignment using the global aligment software needle from the EMBOSS package.
#DEPENDANCY: needle (EMBOSS) http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
#Build all pairwise alignments (which can be long). The multiprocessing module if configured to run 8 alignments in parallel.

#TODO refactoring code, multiple useless parts inherited from previous parwiseid.py
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 20014



from Bio import pairwise2
#from Bio.SubsMat import MatrixInfo as matlist
import re
#from Bio.SubsMat import read_text_matrix
#from Bio import SubsMat
from Bio.pairwise2 import format_alignment
import scipy
import numpy.ma as ma
from multiprocessing import Process, Queue
from multiprocessing import cpu_count
import numpy as np
from Bio import SeqIO
import re
from Bio.Seq import Seq
import sys
import math
from Bio import AlignIO
from io import StringIO
from shell_command import shell_command
from tempfile import NamedTemporaryFile
from Bio.Emboss.Applications import NeedleCommandline


#
# This matrix was created by Todd Lowe   12/10/92
#
# Uses ambiguous nucleotide codes, probabilities rounded to
#  nearest integer
#
# Lowest score = -4, Highest score = 5
#
#    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
#A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2
#T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2
#G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2
#C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2
#S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1
#W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1
#R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1
#Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1
#K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1
#M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1
#B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1
#V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1
#H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  
#D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1
#N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1


def check_alphabet(seq):
    for base in seq:
        if base not in ["A", "T", "G", "C", "a", "t", "g", "c"]:
            print ("strange base:", base )



            

DNA_matrix = {('B', 'N'): -1.0, ('S', 'W'): -4.0, ('G', 'G'): 5.0, ('S', 'N'): -1.0, ('K', 'G'): 1.0, ('H', 'W'): -1.0, ('W', 'R'): -2.0, ('C', 'M'): 1.0, ('V', 'M'): -1.0, ('Y', 'W'): -2.0, ('W', 'C'): -4.0, ('A', 'N'): -2.0, ('N', 'S'): -1.0, ('A', 'Y'): -4.0, ('Y', 'D'): -3.0, ('M', 'Y'): -2.0, ('D', 'G'): -1.0, ('V', 'N'): -1.0, ('Y', 'V'): -3.0, ('G', 'W'): -4.0, ('D', 'W'): -1.0, ('M', 'K'): -4.0, ('R', 'K'): -2.0, ('K', 'W'): -2.0, ('M', 'B'): -3.0, ('W', 'B'): -3.0, ('D', 'A'): -1.0, ('S', 'C'): 1.0, ('Y', 'G'): -4.0, ('H', 'B'): -2.0, ('C', 'K'): -4.0, ('H', 'K'): -3.0, ('Y', 'Y'): -1.0, ('G', 'V'): -1.0, ('K', 'V'): -3.0, ('C', 'Y'): 1.0, ('V', 'Y'): -3.0, ('W', 'A'): 1.0, ('G', 'D'): -1.0, ('K', 'D'): -1.0, ('T', 'N'): -2.0, ('N', 'C'): -2.0, ('W', 'W'): -1.0, ('T', 'W'): 1.0, ('A', 'R'): 1.0, ('M', 'R'): -2.0, ('V', 'H'): -2.0, ('R', 'M'): -2.0, ('S', 'S'): -1.0, ('D', 'H'): -2.0, ('H', 'R'): -3.0, ('D', 'C'): -4.0, ('K', 'C'): -4.0, ('S', 'A'): -4.0, ('W', 'V'): -3.0, ('B', 'A'): -4.0, ('N', 'W'): -1.0, ('G', 'T'): -4.0, ('R', 'N'): -1.0, ('K', 'T'): 1.0, ('C', 'W'): -4.0, ('A', 'K'): -4.0, ('W', 'G'): -4.0, ('T', 'G'): -4.0, ('D', 'K'): -1.0, ('A', 'B'): -4.0, ('K', 'B'): -1.0, ('Y', 'H'): -1.0, ('N', 'N'): -1.0, ('A', 'T'): -4.0, ('B', 'B'): -1.0, ('C', 'H'): -1.0, ('G', 'K'): 1.0, ('D', 'S'): -3.0, ('M', 'G'): -4.0, ('K', 'S'): -2.0, ('C', 'V'): -1.0, ('V', 'T'): -4.0, ('S', 'H'): -3.0, ('K', 'A'): -4.0, ('H', 'Y'): -1.0, ('Y', 'K'): -2.0, ('W', 'T'): 1.0, ('B', 'C'): -1.0, ('C', 'G'): -4.0, ('V', 'K'): -3.0, ('K', 'R'): -2.0, ('A', 'M'): 1.0, ('A', 'D'): -1.0, ('B', 'R'): -3.0, ('T', 'S'): -4.0, ('M', 'W'): -2.0, ('A', 'V'): -1.0, ('B', 'D'): -2.0, ('M', 'N'): -1.0, ('V', 'D'): -2.0, ('N', 'Y'): -1.0, ('M', 'A'): 1.0, ('N', 'K'): -1.0, ('C', 'T'): -4.0, ('V', 'V'): -1.0, ('W', 'D'): -1.0, ('H', 'V'): -2.0, ('B', 'S'): -1.0, ('D', 'N'): -1.0, ('Y', 'M'): -2.0, ('H', 'D'): -2.0, ('D', 'M'): -3.0, ('H', 'M'): -1.0, ('G', 'H'): -4.0, ('R', 'R'): -1.0, ('C', 'S'): 1.0, ('N', 'A'): -2.0, ('V', 'W'): -3.0, ('T', 'C'): -4.0, ('B', 'T'): -1.0, ('K', 'N'): -1.0, ('T', 'H'): -1.0, ('G', 'Y'): -4.0, ('R', 'A'): 1.0, ('C', 'D'): -4.0, ('D', 'V'): -2.0, ('M', 'H'): -1.0, ('S', 'V'): -1.0, ('M', 'C'): 1.0, ('R', 'S'): -2.0, ('C', 'R'): -4.0, ('S', 'M'): -2.0, ('H', 'T'): -1.0, ('S', 'D'): -3.0, ('K', 'M'): -4.0, ('D', 'D'): -1.0, ('R', 'B'): -3.0, ('B', 'G'): -1.0, ('C', 'C'): 5.0, ('V', 'G'): -1.0, ('W', 'K'): -2.0, ('G', 'N'): -2.0, ('R', 'T'): -4.0, ('A', 'A'): 5.0, ('W', 'Y'): -2.0, ('T', 'A'): -4.0, ('B', 'V'): -2.0, ('T', 'V'): -4.0, ('A', 'S'): -4.0, ('Y', 'N'): -1.0, ('M', 'S'): -2.0, ('R', 'C'): -4.0, ('B', 'H'): -2.0, ('C', 'B'): -1.0, ('D', 'T'): -1.0, ('G', 'M'): -4.0, ('S', 'T'): -4.0, ('D', 'B'): -2.0, ('S', 'K'): -2.0, ('V', 'R'): -1.0, ('S', 'B'): -1.0, ('B', 'W'): -3.0, ('K', 'K'): -1.0, ('H', 'C'): -1.0, ('N', 'T'): -2.0, ('H', 'H'): -1.0, ('R', 'D'): -1.0, ('C', 'A'): -4.0, ('V', 'A'): -1.0, ('A', 'H'): -1.0, ('R', 'V'): -1.0, ('A', 'C'): -4.0, ('V', 'S'): -1.0, ('T', 'T'): 5.0, ('M', 'M'): -1.0, ('D', 'R'): -1.0, ('M', 'D'): -3.0, ('V', 'B'): -2.0, ('W', 'H'): -1.0, ('G', 'C'): -4.0, ('S', 'R'): -2.0, ('R', 'W'): -2.0, ('H', 'S'): -3.0, ('Y', 'A'): -4.0, ('B', 'Y'): -1.0, ('H', 'A'): -1.0, ('N', 'D'): -1.0, ('Y', 'S'): -2.0, ('H', 'N'): -1.0, ('B', 'K'): -1.0, ('N', 'G'): -2.0, ('V', 'C'): -1.0, ('G', 'B'): -1.0, ('T', 'D'): -1.0, ('N', 'B'): -1.0, ('T', 'M'): -4.0, ('K', 'H'): -3.0, ('T', 'R'): -4.0, ('M', 'T'): -4.0, ('A', 'W'): 1.0, ('Y', 'R'): -4.0, ('G', 'S'): 1.0, ('R', 'G'): 1.0, ('S', 'Y'): -2.0, ('N', 'H'): -1.0, ('W', 'N'): -1.0, ('G', 'A'): -4.0, ('D', 'Y'): -3.0, ('R', 'Y'): -4.0, ('K', 'Y'): -2.0, ('S', 'G'): 1.0, ('N', 'R'): -1.0, ('Y', 'C'): 1.0, ('H', 'G'): -4.0, ('N', 'V'): -1.0, ('G', 'R'): 1.0, ('R', 'H'): -3.0, ('B', 'M'): -3.0, ('W', 'M'): -2.0, ('T', 'B'): -1.0, ('A', 'G'): -4.0, ('Y', 'B'): -1.0, ('W', 'S'): -4.0, ('T', 'K'): 1.0, ('C', 'N'): -2.0, ('M', 'V'): -1.0, ('N', 'M'): -1.0, ('Y', 'T'): 1.0, ('T', 'Y'): 1.0}






def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]



def _read_text_matrix(data_file):
    "Initially used to parse the EDAN substitution matrix"
    matrix = {}
    tmp = data_file.read().split("\n") 
    table=[] 
    for i in tmp: 
        table.append(i.split()) 
    # remove records beginning with ``#'' 
    for rec in table[:]: 
        if (rec.count('#') > 0):
            table.remove(rec) 
    
    # remove null lists
    while (table.count([]) > 0): 
        table.remove([])
    # build a dictionary 
    alphabet = table[0] 
    j = 0 
    for rec in table[1:]:
        row = alphabet[j]
        if re.compile('[A-z\*]').match(rec[0]):
             first_col = 1
        else: 
            first_col = 0
        i = 0 
        for field in rec[first_col:]:
            col = alphabet[i] 
            matrix[(row, col)] = float(field) 
            i += 1
        j += 1
    # delete entries with an asteris
    for i in matrix: 
        if '*' in i: 
             del(matrix[i]) 
    return matrix

def pairewise_identity(seq1, seq2):
    "return identity calculated as: n identical sites/n aligned sites (gaps non included)"
    A=list(seq1)
    B=list(seq2)
    identical_sites = 0
    aligned_sites = 0
    gaps=0
    for n in range(0, len(A)):
        if str(A[n]) != "-" and str(B[n]) != "-":
            aligned_sites+=1
        else:
            continue
        if str(A[n])==str(B[n]):
            identical_sites+=1
    try:
        identity = 100*(identical_sites/float(aligned_sites))
    except:
        identity = 0
    #print "aligned", float(aligned_sites)
    #print "identical", identical_sites
    if aligned_sites == 0:
        # return false if align length = 0
        return False
    else:
        return 100*(identical_sites/float(aligned_sites))

    
def global_align(seq_record1, seq_record2):
    """Global alignment using the Bio.pairwise2 package. 
    Check if sequences are nucleotide or amino acids using the _verify_alphabet function from the Bio.Alphabet module.
    """

    #from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    #from Bio.Alphabet import _verify_alphabet
    
    #gap_open = -10
    #gap_extend = -0.5
    #print type(seq_record1)
    #print type(seq_record2)
    
    seq_record1.seq = seq_record1.seq.upper()
    seq_record2.seq = seq_record2.seq.upper()
   
    '''
        temp_ref = NamedTemporaryFile(delete=False, mode='w')
        fastastr = StringIO()

        SeqIO.write(self.ref_locus_seqrecord, fastastr, 'fasta')

        temp_ref.write(fastastr.getvalue())
        temp_ref.flush()
    
    '''
   
   
    seq1_file = NamedTemporaryFile()
    fastastr = StringIO()
    SeqIO.write(seq_record1, fastastr, "fasta")
    seq1_file.write(fastastr.getvalue())
    seq1_file.flush()        
    seq2_file = NamedTemporaryFile()
    fastastr2 = StringIO()
    SeqIO.write(seq_record2, fastastr2, "fasta")
    seq2_file.write(fastastr2.getvalue())
    seq2_file.flush()

    #print seq_record1.seq.alphabet
    #print seq_record1.seq.alphabet
    
    
    #seq_record1.seq.alphabet = IUPAC.ambiguous_dna
    #seq_record2.seq.alphabet = IUPAC.ambiguous_dna


    #print "1", _verify_alphabet(seq_record1.seq)
    #print "2",_verify_alphabet(seq_record2.seq)
    
    if _verify_alphabet(seq_record1.seq) and _verify_alphabet(seq_record2.seq):
        #print "DNA!"
    #    alns = pairwise2.align.globalds(seq1, seq2, DNA_matrix, gap_open, gap_extend)
    #    print ">"+noms[id_seq1]
    #    print alns[0][0]
    #    print ">"+noms[id_seq2]
    #    print alns[0][1]
    #    return  alns[0]
        needle_cline=NeedleCommandline(asequence=seq1_file.name,
                                       bsequence=seq2_file.name,
                                       stdout=True,
                                       gapopen=10,
                                       gapextend=0.5,
                                       auto=True,
                                       aformat="srspair")
        stdout,stderr=needle_cline()
        #print stdout
        print (needle_cline)
        align = AlignIO.read(StringIO.StringIO(stdout), "emboss")
        return align

 
 
        

    #seq_record1.seq.alphabet = IUPAC.extended_protein
    #seq_record2.seq.alphabet = IUPAC.extended_protein
    #print seq1
    #print _verify_alphabet(seq1)

    if _verify_alphabet(seq_record1.seq) and _verify_alphabet(seq_record2.seq):
        #print "AA!"
    #    alns = pairwise2.align.globalds(seq1, seq2, matlist.blosum62, gap_open, gap_extend)
    #    return  alns[0]

        needle_cline=NeedleCommandline(asequence=seq1_file.name, bsequence=seq2_file.name, stdout=True, gapopen=10, gapextend=0.5, auto=True, aformat="srspair")
        stdout,stderr=needle_cline()
        align = AlignIO.read(StringIO.StringIO(stdout), "emboss")
        return align

    else:
        print (seq_record1.seq)
        check_alphabet(seq_record1.seq)
        print (seq_record2.seq)
        check_alphabet(seq_record2.seq)
        raise "unkown alphabet!"




    






    



        

def get_identity_from_2_seqrecord(record, seq1, seq2):
    """
    record = record containg all input fasta sequences
    seq1 and seq2: indexes of the to sequences to align and get identity

    """
    align = global_align(record[seq1], record[seq2])
    identity = pairewise_identity(align[0].seq, align[1].seq)
    return identity


def get_identity_from_combination_list(record, combinations, out_q):

    """
    record: record containg all input fasta sequences
    combinations: list of tuples containing all pairwise combinations of the record sequences: [(1,1), (1,2), (1,3),...]
    """
    
    #print "compb list", combinations

    result = []
    for one_comb in combinations:
        #print "one_comb", one_comb
        identity = get_identity_from_2_seqrecord(record, one_comb[0], one_comb[1])
        result.append(one_comb + (identity,))

    out_q.put(result)

def get_identity_matrix_from_multifasta(record):
    combinations = []
    
    # get a list of all combination of sequences to align
    # allows then to split the job into 8 process
    for i in range(0, len(record)):
        for y in range(i, len(record)):
            combinations.append((i,y))
    #print combinations

    out_q = Queue()
    n_cpu = 1
    
    n_poc_per_list = math.ceil(len(combinations)/float(n_cpu))
    #print "n_poc_per_list", n_poc_per_list
    #print range(0, len(combinations))
    # split the jobs into 8 lists
    query_lists = chunks(range(0, len(combinations)), int(n_poc_per_list))
    # start the 8 processes
    procs = []
    for one_list in query_lists:
        comb_list = [combinations[i] for i in one_list]
        #print len(comb_list)
        proc = Process(target=get_identity_from_combination_list, args=(record, comb_list, out_q))
        procs.append(proc)
        proc.start()

    complete_result = []
    for i in range(len(procs)):
        complete_result += out_q.get()
    import time
    time.sleep(5)

    for proc in procs:
        proc.join()

    #print "result length:", len(complete_result)

    identity_matrix = np.empty([len(record), len(record)])
    for i in complete_result:
        #print i
        identity_matrix[i[0],i[1]] = i[2]
        identity_matrix[i[1],i[0]] = i[2]

    return identity_matrix



def write_id_table(record, array, out_name):
    f = open(out_name, "w")
    ids = [i.name for i in record]
    f.write("\t" + "\t".join(ids) + "\n")
    for i in range(0,len(array)):
        name = ids[i]
        row = array[i,:]
        row = [str(round(i,2)) for i in row]
        f.write(name + "\t" + "\t".join(row) + "\n")





if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", '--seq1',type=str,help="seq1")
    parser.add_argument("-b", '--seq2',type=str,help="seq2")
    parser.add_argument("-m", '--multifasta',type=str,help="input multi fasta, align all against all")
    parser.add_argument("-o", '--out_name',type=str,help="output name, default = input_name_identity.tab")

    args = parser.parse_args()

    if args.seq1 and args.seq2:

        sys.exit()
        #alns = pairwise2.align.globalds(args.seq1, args.seq2, DNA_matrix, -10, -0.5)
        #print format_alignment(*alns[0])

    
    if not args.out_name:
        out_name = re.sub("\..*", "_identity.tab", args.multifasta)
    else:
        out_name = args.out_name

    if args.multifasta: 
        #handle = open(args.multifasta, "rU")
        multifasta = [record for record in SeqIO.parse(args.multifasta, "fasta")]
        id_matrix = get_identity_matrix_from_multifasta(multifasta)
        write_id_table(multifasta, id_matrix, out_name)

    #print SubsMat.make_log_odds_matrix(mat)
    #mat.print_full_mat(alphabet=["A", "T", "G", "C", "S", "W", "R","Y", "K", "M", "B", "V", "H", "D", "N"])
    #print matlist.blosum62
    
