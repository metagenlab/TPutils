#!/usr/bin/env python

# Perform various GC calculations: GC skew, gc variations from the average, plots of length vs gc, circos input files,...
# TODO separate circos stuff from the rest
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

from Bio import SeqIO
from Bio.SeqUtils import GC, GC_skew
import pylab
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages



def circos_gc_var(record, windows=1000, shift=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%

    '''
    circos_string = ''
    from Bio.SeqFeature import FeatureLocation
    average_gc = GC(record.seq)
    gap_locations = []
    for feature in record.features:
        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))
    else:
        gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))
    if len(gap_locations) > 1:
        gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))


        for i in range(0, len(gap_locations)):
            if i == 0:
                seq = record.seq[0:gap_locations[i].start]
                chr_start = 0
            else:
                seq = record.seq[gap_locations[i-1].end:gap_locations[i].start]
                chr_start = gap_locations[i-1].end
            contig_name = record.name + "_%s" % (i +1)
            for i in range(0, len(seq), windows):
                start = i
                stop = i + windows
                #gc = ((GC(record.seq[start:stop])/average_gc) - 1)*100
                gc = GC(record.seq[start:stop]) - average_gc
                if stop > len(seq):
                    stop = len(seq)
                    if stop - start < 200:
                        break
                section_start = chr_start + start
                section_end = chr_start + stop
                circos_string += "%s %s %s %s\n" % (contig_name, section_start, section_end, gc)
    else:
        seq = record.seq
        contig_name = record.id #.split('.')[0]
        for i in range(0, len(seq), windows):
            start = i
            stop = i + windows
            #gc = ((GC(record.seq[start:stop])/average_gc) - 1)*100
            gc = GC(record.seq[start:stop]) - average_gc
            if stop > len(seq):
                stop = len(seq)
                if stop - start < 200:
                    break
            circos_string += "%s %s %s %s\n" % (contig_name, start+shift, stop+shift, gc)
    return circos_string


def circos_gc_skew(record, windows=1000, shift=0):
    '''
    :param record:
    :return: circos string with difference as compared to the average GC
    ex: average = 32
        GC(seq[3000:4000]) = 44
        diff = 44 - 32 = 12%

    '''
    from Bio.SeqFeature import FeatureLocation
    circos_string = ''
    #print "GENOME SIZE:", len(record.seq)

    gap_locations = []
    for feature in record.features:

        if feature.type == "assembly_gap":
            gap_locations.append(feature.location)
    if len(gap_locations) == 0:
        gap_locations.append(FeatureLocation(0, len(record.seq)))

    else:
        #gap_locations.append(FeatureLocation(gap_locations[-1].end + 1, len(record.seq)))
        gap_locations.append(FeatureLocation(len(record.seq), len(record.seq)))
    #print 'gap locations', gap_locations
    if len(gap_locations) > 1:
        for i in range(0, len(gap_locations)):
            if i == 0:
                seq = record.seq[0:gap_locations[i].start]
                chr_start = 0
            else:
                seq = record.seq[gap_locations[i-1].end:gap_locations[i].start]
                chr_start = gap_locations[i-1].end
            #print i, "seq", gap_locations[i-1].end, gap_locations[i].start, gap_locations[i].start - gap_locations[i-1].end

            try:
                values = GC_skew(seq, windows)
            except:
                print len(seq), seq
            contig_name = record.name + "_%s" % (i + 1)

            for i in range(0, len(values)):
                start = i *windows
                stop = start + windows
                #gc = ((GC(record.seq[start:stop])/average_gc) - 1)*100
                section_start = chr_start + start
                section_end = chr_start + stop
                circos_string += "%s %s %s %s\n" % (contig_name, section_start+shift, section_end+shift, values[i])
    else:
        #print 'no gaps!'
        try:
            values = GC_skew(record.seq, windows)
        except:
            values = GC_skew(record.seq, 2000)
        for i in range(0, len(values)):
            start = i *windows
            stop = start + windows

            circos_string += "%s %s %s %s\n" % (record.id, start+shift, stop+shift, values[i])
    return circos_string

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
    pylab.xlabel("Sequence length (bp)")
    pylab.ylabel("GC (%)")

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

def whole_gc(records):
    seq = ""
    for record in records:
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
    parser.add_argument("-z",'--circos', action='store_true', help="circos output")


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
        records = [record for record in SeqIO.parse(handle, "fasta")]
        print whole_gc(records)
    
    if args.circos:
        records = [record for record in SeqIO.parse(handle, "fasta")]
        print circos_gc_skew(records)
