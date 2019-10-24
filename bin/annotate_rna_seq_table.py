#!/usr/bin/env python

import get_chipseqpeak_neigbor_genes

def parse_rnaseq_file(rnaseq_file):

    accession2peak2position = {}

    with open(rnaseq_file) as f:
        for i, row in enumerate(f):
            if i == 0:
                continue
            else:
                data = row.rstrip().split('\t')
                if data[0] not in accession2peak2position:
                    accession2peak2position[data[0]] = {}
                    accession2peak2position[data[0]]['peak_%s' % i] = int(data[3])
                else:
                    accession2peak2position[data[0]]['peak_%s' % i] = int(data[3])
    return accession2peak2position

def annotate_chipseq(genbank_file, rna_seq_file):

    '''

    iterer CDS
        iterer peaks
        if max_peak within the window of the start codon:
            keep

    :param genbank_file:
    :param peak_file:
    :return:
    '''

    from Bio import SeqIO

    peak_center_positions = parse_peak_file(peak_file)
    with open(genbank_file) as g:
        records = [i for i in SeqIO.parse(g, 'genbank')]
        for record in records:

            record_peaks = peak_center_positions[record.id]
            peak2scores_dico = peak2scores(peak_file)
            locus2COG, locus2KEGG = get_annotation(record.name, 'chlamydia_04_16')
            print "peak_id\tpeak_max\tlocus_tag\tstart\tend\tstrand\tproduct\tpeak_length\tpeak_shape_score\tpeak_p_value" \
                  "\tCOG\tCOG_description\tCOG_categories\tko\tko_description\tkegg_pathways"
            for peak in record_peaks:
                #print peak, record_peaks[peak]
                for feature in record.features:
                    if feature.type == 'CDS':
                        #print feature.location.start
                        #print feature.location.end
                        if feature.location.strand == 1:
                            amont_position = int(feature.location.start) - max_amont
                            aval_position = int(feature.location.start) + max_aval
                        else:
                            amont_position = int(feature.location.end) + max_amont
                            aval_position = int(feature.location.end) - max_aval
                        #print peak, max([amont_position, aval_position]), min([amont_position, aval_position])
                        if record_peaks[peak] <= max([amont_position, aval_position]) and record_peaks[peak] >= min([amont_position, aval_position]):
                            #print 'peak!!!!!'
                            try:
                                cogg_data = locus2COG[feature.qualifiers['locus_tag'][0]]
                                cog_name = cogg_data[0][0]
                                cog_description = cogg_data[0][2]
                                cog_categories = ','.join([i[1] for i in cogg_data])
                            except:
                                cog_name = '-'
                                cog_description = '-'
                                cog_categories = '-'
                            try:
                                ko_data = locus2KEGG[feature.qualifiers['locus_tag'][0]]
                                #print ko_data
                                #print [i[2] for i in ko_data]
                                #print [i[3] for i in ko_data]
                                ko_name = ko_data[0][0]
                                ko_description = ko_data[0][1]
                                ko_pathways = [i[2] for i in ko_data]
                                ko_pathways_description = [i[3] for i in ko_data]

                                pathway_annot = ''
                                for i, path in enumerate(ko_pathways):
                                    pathway_annot+="%s(%s) " % (path, ko_pathways_description[i])


                            except:
                                ko_name = '-'
                                pathway_annot = '-'
                                ko_description = '-'

                            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (peak,
                                                                                record_peaks[peak],
                                                                                feature.qualifiers['locus_tag'][0],
                                                                                feature.location.start,
                                                                                feature.location.end,
                                                                                feature.location.strand,
                                                                                feature.qualifiers['product'][0],
                                                                                peak2scores_dico[peak][0],
                                                                                peak2scores_dico[peak][1],
                                                                                peak2scores_dico[peak][2],
                                                                                      cog_name,
                                                                                      cog_description,
                                                                                      cog_categories,
                                                                                                      ko_name,
                                                                                                      ko_description,
                                                                                                      pathway_annot)
                            #print feature
                            #import sys
                            #sys.exit()


if __name__ == '__main__':
    import argparse
    import generate_bsub_file
    import rerun_LFS_jobs

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_chipseq_table', type=str, help="input chipseq table file")
    parser.add_argument("-g", '--genbank', help="genbank file")
    parser.add_argument("-av", '--dist_amont', type=int, help="max distance en amont du start codon (default=400)")
    parser.add_argument("-am", '--dist_aval', type=int, help="max distance en aval du start codon (default=100)")

    args = parser.parse_args()


    get_peak_associated_genes(args.genbank, args.input_chipseq_table)
