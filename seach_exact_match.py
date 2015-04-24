#!/usr/bin/env python







if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input',type=str,help="input fasta")
    parser.add_argument("-k",'--kmer',type=str,help="exact match to search")
        
    args = parser.parse_args()
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    
    target_seq = Seq(args.kmer, generic_dna)

    
    from Bio import SeqIO
    handle = open(args.input, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        seq = record.seq
        #import time
        #time.sleep(3)
        for i in range(0,len(seq)):
            #if i%1000000 == 0:
            #    print i
            #print seq[i:i+len(args.kmer)]
            subseq = seq[i:i+len(target_seq)]
            #print subseq
            #print args.kmer
            if str(subseq) == str(target_seq):
                print "MATCH:", i, i + len(target_seq), str(target_seq)
            elif str(subseq) ==  str(target_seq.reverse_complement()):
                print "MATCH:", i, i + len(target_seq), str(target_seq), "(rev. complement)"
            elif str(subseq) ==  str(target_seq.complement()):
                print "MATCH:", i, i + len(target_seq), str(target_seq), "(complement)"
            else:
                pass
        
    handle.close()
