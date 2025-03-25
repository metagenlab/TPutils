#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


'''
def multi_delete(index_list, list_):
    from datetime import datetime
    indexes = sorted(index_list, reverse=True)
    for i, index in enumerate(indexes):
        if i%10000 == 0:
            print ("%s%%, %s" % (round((float(i)/len(indexes))*100,2), str(datetime.now())))
        del list_[i]
    return list_
'''
def multi_delete(indexes, list_):
    from datetime import datetime
    indexes = sorted(indexes, reverse=True)
    for i, index in enumerate(indexes):
        if i%10000 == 0:
            print ("%s%%, %s" % (round((float(i)/len(indexes))*100,2), str(datetime.now())))
        del list_[index]
    return list_

def compare_sequences(seq_list, threshold=0):
    align_length = len(seq_list[0])
    for seq in seq_list:
        if len(seq) != align_length:
            raise IOError('Input sequences should have the same length!')

    remove_list = []
    for i in range(0, align_length):
        if i%10000 == 0:
            print ("%s %%" % round((float(i)/align_length)*100,2) )
        missing_count = 0
        for seq in seq_list:
            if seq[i].lower() == 'n':
                missing_count+=1
            elif seq[i].lower() == '-':
                missing_count+=1
            else:
                continue
        if missing_count > threshold:
           remove_list.append(i) 
    return remove_list


def clean_align(input_records):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC    
    sequences = [list(i.seq) for i in input_records]
    columns_to_remove = compare_sequences(sequences, 3)
    print ('removing %s columns from alignment of length %s (%s%%)' % (len(columns_to_remove),
                                                                      len(sequences[0]),
                                                                      round(float(len(columns_to_remove))/len(sequences[0])*100, 2)))
    for i, record in enumerate(input_records):
        print (i)
        new_seq = Seq(''.join(multi_delete(columns_to_remove, sequences[i])), IUPAC.ambiguous_dna)
        input_records[i].seq = new_seq
    output_handle = open("test.fa", "w")    
    SeqIO.write(input_records, output_handle, "fasta")
    




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--multifasta', type=str, help="input multifasta")
    args = parser.parse_args()
    from Bio import SeqIO
    records = [i for i in SeqIO.parse(args.multifasta, "fasta")]
    clean_align(records)
    

    
