#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def GC_skew(seq, window=100):
    """Calculates GC skew (G-C)/(G+C) for multiple windows along the sequence.
    Returns a list of ratios (floats), controlled by the length of the sequence 
    and the size of the window.
    Does NOT look at any ambiguous nucleotides.
    """
    values = []
    for i in range(0, len(seq), window):
        s = seq[i: i + window] 
        g = s.count('G') + s.count('g') 
        c = s.count('C') + s.count('c')
        #print g-c
        #print g+c
        #print s
        try:
            skew = (g-c)/float(g+c)
        except:
            skew = 0
        values.append(skew) 
    return values



def plot_contig_len(handle, out_name):
    pp = PdfPages(out_name)
    sizes = [len(rec) for rec in SeqIO.parse(handle, "fasta")]
    pylab.hist(sizes, bins=20)
    pylab.title("%i contig sequences\nLengths %i to %i" \
            % (len(sizes),min(sizes),max(sizes)))
    pylab.xlabel("Sequence length (bp)")
    pylab.ylabel("Count")
    pylab.savefig(pp, format='pdf')
    pylab.close()
    pp.close()


def plot_GC_length(handle, out_name, xlim = False):
    pp = PdfPages(out_name)
    parsed_handle = [record for record in SeqIO.parse(handle, "fasta")]
    gc_values = [GC(rec.seq)  for rec in parsed_handle]
    
    #for i in parsed_handle:
    #    print i.id, GC(i.seq), len(i.seq)
    #print len(gc_values)
    len_values = [len(rec.seq) for rec in parsed_handle]
    #print len(len_values)
    #print gc_values

    #pylab.plot(gc_values)
    #pylab.title("%i orchid sequences\nGC%% %0.1f to %0.1f" \
    #            % (len(gc_values),min(gc_values),max(gc_values)))
    #pylab.xlabel("Genes")
    #pylab.ylabel("GC%")
    #pylab.show()
    pylab.plot(len_values,gc_values, 'ro')
    pylab.ylim(20.0,80.0)
    if xlim:
        pylab.xlim(0, int(xlim))
    pylab.savefig(pp, format='pdf')
    pylab.close()
    pp.close()


def gc_values(handle):
    parsed_handle = [record for record in SeqIO.parse(handle, "fasta")]
    gc_values = [GC(rec.seq)  for rec in parsed_handle]
    seq_ids = [rec.id  for rec in parsed_handle]
    for i in range(0, len(seq_ids)):
        print seq_ids[i] +"\t"+ str(gc_values[i])

def whole_gc(handle):
    seq = ""
    parsed_handle = [record for record in SeqIO.parse(handle, "fasta")]
    for record in parsed_handle:
        seq+= record.seq
    return GC(seq)



def plot_GC(handle, outname):
    pp = PdfPages(outname)
    parsed_handle = [record for record in SeqIO.parse(handle, "fasta")]

    for i in range(0, len(parsed_handle)):

        gc_skew = GC_skew(parsed_handle[i].seq)

        cumulated_skew = np.cumsum(gc_skew)
        pylab.plot(cumulated_skew)
        pylab.plot(gc_skew)
        pylab.title("contig %s"  % (str(parsed_handle[i].id)))
        pylab.savefig(pp, format='pdf')
        pylab.close()

    all_seq = ''
    for record in parsed_handle:
        all_seq+=record.seq
    gc_skew = GC_skew(all_seq)
    cumulated_skew = np.cumsum(gc_skew)
    pylab.plot(cumulated_skew)
    pylab.plot(gc_skew)
    pylab.title("concatenated contigs")
    pylab.savefig(pp, format='pdf')
    pp.close()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fasta',type=str,help="input fasta file ")
    parser.add_argument("-i",'--ind',action='store_true',help="get individual sequence GC percent ")
    parser.add_argument("-p",'--plot_gc',action='store_true',help="plot GC content")
    parser.add_argument("-l",'--len',action='store_true',help="plot sequence length versus GC")
    parser.add_argument("-d",'--len_dist',action='store_true',help="plot length distribution")
    parser.add_argument("-c",'--conc_gc',action='store_true',help="concatenated sequences GC")
    parser.add_argument("-o",'--outname',type=str,help="output pdf file name (for plots only)")
    parser.add_argument("-x",'--xlim',type=str,help="xlim for sequence size (optional)")


    args = parser.parse_args()
    handle = open(args.fasta)

    if args.plot_gc:
        plot_GC(handle, args.outname)
    if args.len:
        plot_GC_length(handle, args.outname, args.xlim)
    if args.len_dist:
        plot_contig_len(handle, args.outname)
    if  args.ind:
        gc_values(handle)

    if args.conc_gc:
        print whole_gc(handle)
