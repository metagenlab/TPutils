#!/usr/bin/env python


import vcf
from Bio import SeqIO

'''
## FILTER=<ID=PASS,Description="All filters passed">
## ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
## FILTER=<ID=LowQual,Description="Low quality">
## FILTER=<ID=freqalt,Description="Set if true: GT='alt' & (FORMAT/AD[0:1])/(FORMAT/DP) < 0.75">
## FILTER=<ID=freqref,Description="Set if true: GT='ref' & (FORMAT/AD[0:0])/(FORMAT/DP) < 0.75">
## FILTER=<ID=lowcov,Description="Set if true: (FORMAT/DP)<10">
## FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
## FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
## FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
## FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
## FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
## FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
## FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
'''

def parse_gbk(gbk_file):

    record_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, 'genbank'))
    return record_dict


def get_neiboring_orf(position, feature_list):

    if int(position) > int(feature_list[-1].location.end):
        try:
            gene = feature_list[-1].qualifiers["gene"][0]
            locus_tag = feature_list[-1].qualifiers["locus_tag"][0]
        except KeyError:
            gene = '-'
            locus_tag = feature_list[-1].qualifiers["locus_tag"][0]
        return ['%s (%s)' % (locus_tag, gene), '-']
    first_feature = False
    for n, feature in enumerate(feature_list):
        if not first_feature:
            if int(position) < int(feature.location.start):

                try:
                    gene = feature.qualifiers["gene"][0]
                    locus_tag = feature.qualifiers["locus_tag"][0]
                except KeyError:
                    gene = '-'
                    locus_tag = feature.qualifiers["locus_tag"][0]
                return ["-", "%s (%s)" % (locus_tag, gene)]
        if feature.type != 'source':
            first_feature = True


        if int(position) > int(feature.location.end) and int(position) < int(feature_list[n + 1].location.start):
            try:
                gene1 = feature.qualifiers["gene"][0]
                locus_tag1 = feature.qualifiers["locus_tag"][0]
            except KeyError:
                gene1 = '-'
                locus_tag1 = feature.qualifiers["locus_tag"][0]
            try:
                gene2 = feature_list[n + 1].qualifiers["gene"][0]
                locus_tag2 = feature_list[n + 1].qualifiers["locus_tag"][0]
            except KeyError:
                gene2 = '-'
                locus_tag2 = feature_list[n + 1].qualifiers["locus_tag"][0]
            return ["%s (%s)" % (locus_tag1, gene1), "%s (%s)" % (locus_tag2, gene2)]
    return ["no CDS", "no CDS"]



def parse_vcf(vcf_file, gbk_file):
    from copy import copy
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import MutableSeq
    from Bio.Alphabet import generic_dna

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    gbk_dico = parse_gbk(gbk_file)

    intergenic = 0
    synonymous = 0
    non_synonymous = 0
    indel = 0

    header = ["contig", "contig_length", "position", "AC", "REF", "ALT", "MUTATION_location", "MUTATION_type", "ORF", "gene", "orf_before", "orf_after"] + list(vcf_reader.filters.keys())
    print ('\t'.join(header))
    for n, record in enumerate(vcf_reader):

        record_alt = copy(gbk_dico[record.CHROM])
        record_alt.seq = MutableSeq(str(record_alt.seq), generic_dna)
        # default to intergenic (modigied if feature match found)
        mut_location = 'Intergenic'
        mut_type = '-'
        orf_name = '-'
        gene = '-'
        for feature in record_alt.features:
            if int(record.POS) in feature and feature.type != "source":
                mut_location = feature.type
                orf_name = feature.qualifiers["locus_tag"][0]
                try:
                    gene = feature.qualifiers["gene"][0]
                except:
                    gene = '-'
                #print("SNP in feature--------------------------------------------:", feature.type, feature.location)
                if feature.type == 'CDS':

                    if len(record.ALT[0]) > 1:
                        mut_type = 'INDEL'
                        indel+=1
                        continue
                    else:
                        aa_seq_ref = str(feature.extract(record_alt.seq).translate())
                        # mutate reference sequence
                        record_alt.seq[int(record.POS)-1] = str(record.ALT[0])
                        # check if synonymous or not
                        aa_seq_alt = str(feature.extract(record_alt.seq).translate())

                        if str(aa_seq_ref) == str(aa_seq_alt):
                            mut_type = 'SYNONYMOUS'
                            synonymous += 1
                        else:
                            mut_type = 'NON SYNONYMOUS'
                            if len(record.FILTER) == 0:
                                non_synonymous+=1
                        #print (str(aa_seq_ref.translate()) == str(record_alt.translate()))
        if mut_location == 'Intergenic':
            orf_before, orf_after = get_neiboring_orf(int(record.POS), record_alt.features)
        else:
            orf_before, orf_after = ['-', '-']
        #print (filter_dico["PASS"].id)
        #print (filter_dico["PASS"].desc)
        contig = record.CHROM
        ac = record.INFO["AC"][0]
        position = record.POS
        ref = record.REF
        alt = record.ALT[0]
        filter_status = []

        # if any of the test failed, set PASS as failed
        if len(record.FILTER) != 0:
            record.FILTER.append('PASS')

        for filter_name in vcf_reader.filters:
            if filter_name in record.FILTER:
                if filter_name == 'PASS':
                    filter_status.append('NO')
                else:
                    filter_status.append('FAIL')
            else:
                if filter_name == 'PASS':
                    filter_status.append('YES')
                else:
                    filter_status.append('SUCCESS')

        print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig,
                                                   len(record_alt),
                                                   position,
                                                   ac,
                                                    ref,
                                                    alt,
                                                    mut_location,
                                                    mut_type,
                                                    orf_name,
                                                    gene,
                                                    orf_before,
                                                    orf_after,
                                                    '\t'.join(list(filter_status))))


    #print ('synonymous', synonymous)
    #print ('non synonymous', non_synonymous)
    #print ('indel', indel)


f = '/home/tpillone/work/projets/2018_diag_analyses/2018_11_19_TATRas/typing/gatk_gvcfs/full_genome_TATRas-control_assembled_genome/bwa/filtering/freq_cov_decomposed_normalized.vcf'
f = '/home/tpillone/work/projets/2018_diag_analyses/2018_11_19_TATRas/samples/TATRas-mutant-A/snps/gatk_gvcfs/TATRas-control_assembled_genome/bwa/freq.vcf'
g = '/home/tpillone/work/projets/2018_diag_analyses/2018_11_19_TATRas/samples/TATRas-control/annotation/TATRas-control.gbk'

parse_vcf(f, g)
