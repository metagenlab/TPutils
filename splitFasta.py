#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# split big multi fasta in multiple files with definied number of sequences (default 1000)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# -----------------------------------------------------


def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch


def split_fasta(fasta_file, number_per_file=1000):
    from Bio import SeqIO
    record_iter = SeqIO.parse(open(fasta_file), "fasta")
    for i, batch in enumerate(batch_iterator(record_iter, number_per_file)) :
         filename = "group_%i.fasta" % (i+1)
         handle = open(filename, "w")
         count = SeqIO.write(batch, handle, "fasta")
         handle.close()
         print "Wrote %i records to %s" % (count, filename)
    return i+1


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta',type=str,help="input fasta file")
    parser.add_argument("-n", '--number',type=int,help="number of sequences per file", default = 1000)

                              
    args = parser.parse_args()
    split_fasta(args.input_fasta, args.number)

