#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pairwiseid_needle
from blast_utils import run_prodigal

def check_blast_colocalization(record1, record2, seq_range = 0.15):
    start1 = int(record1.alignments[0].hsps[0].sbjct_start)
    end1 = int(record1.alignments[0].hsps[0].sbjct_end)

    start_range = sorted([start1 - start1*seq_range, start1 + start1*seq_range])
    end_range = sorted([end1 - end1*seq_range, end1 + end1*seq_range])

    end2 = record2.alignments[0].hsps[0].sbjct_start
    start2 = record2.alignments[0].hsps[0].sbjct_end

    if start2 > start_range[0] and start2 < start_range[1]:
        start = True
    else:
        print "start", start2, start_range
        start = False

    if end2 > end_range[0] and end2 < end_range[1]:
        end = True
    else:
        end = False
        print "end", end2, end_range

    if start == True and end == True:
        return True
    else:
        return False


def hmm_and_extract(query_fasta_prot, db_prot, seq_header = "query_seq"):
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    from Bio import SearchIO
    import os

    #print 'db_prot', db_prot
    #print 'querxy', query_fasta_prot
    db_prot_data = SeqIO.to_dict(SeqIO.parse(db_prot, "fasta"))
    # -E 0.000001
    # 0.000001
    # -T 150
    print
    qcode = os.path.basename(query_fasta_prot).split('.')[0]
    rcode = os.path.basename(db_prot).split('.')[0]
    # --max
    # -T <x>     : report sequences >= this score threshold in output
    hmmer_cline = 'hmmsearch -T 150 --tblout temp_hmm_%s_%s.tab %s %s' % (qcode, rcode, query_fasta_prot, db_prot)
    print hmmer_cline
    stdout, stderr, code = shell_command.shell_command(hmmer_cline)

    result_handle = open("temp_hmm_%s_%s.tab" % (qcode, rcode), 'r')
    hmmer_records = [i for i in SearchIO.parse(result_handle, 'hmmer3-tab')]

    try:
        best_hit_id = hmmer_records[0].hits[0].id
        #print dir(hmmer_records[0].hits[0])

        if float(hmmer_records[0].hits[0].bias) > 0:
            #print '################### bias!'
            for one_hit in hmmer_records[0].hits:
                score_diff = float(hmmer_records[0].hits[0].bitscore) - float(one_hit.bitscore)
                print 'score 1', 'score 2'
                print "score_diff", score_diff,  0.15*float(hmmer_records[0].hits[0].bitscore)
                if (float(one_hit.bias) < float(hmmer_records[0].hits[0].bias)) \
                        and (score_diff < (0.15*float(hmmer_records[0].hits[0].bitscore))):
                    best_hit_id = one_hit.id
                    break
        print 'best hit', best_hit_id
        best = db_prot_data[best_hit_id]
        best.name = best.name.split('_')[0]
        split_name = best.id.split('_')
        if len(split_name)>2:
            best.id = best_hit_id + '|' + split_name[0] + '_' + split_name[1]
        else:
            best.id = best_hit_id #+ '|' + split_name[0]

        return db_prot_data[best_hit_id]
    except:
        return False

def blastp_and_extract(query_fasta_prot, db_prot, seq_header = "query_seq"):
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio import SeqIO

    print 'db_prot', db_prot
    print 'querxy', query_fasta_prot
    db_prot_data = SeqIO.to_dict(SeqIO.parse(db_prot, "fasta"))
    
    blastp_cline = NcbiblastpCommandline(query=query_fasta_prot, db=db_prot, evalue=0.001, outfmt=5, out="temp.xml")
    print blastp_cline
    stdout, stderr = blastp_cline()

    result_handle = open("temp.xml")
    blast_records = NCBIXML.parse(result_handle)
    blast_records = list(blast_records)
    #best_hit = blast_records[0]
    #best_hit_length = len(best_hit.alignments[0].hsps[0].sbjct.replace("-", ""))
    best_hit = False
    for i in range(1, len(blast_records)):
        try:
            blast_record = blast_records[i]
            #print blast_record.query, "vs", blast_record.alignments[0].hit_def
            #print blast_record.alignments[0].hsps[0].sbjct.replace("-", "")
            if not best_hit:
                best_hit = blast_record.alignments[0].title.split(' ')[1]
            else:
                if blast_record.alignments[0].title.split(' ')[1] != best_hit:
                    print "Different hit!!!"
        except:
            print "NO hit!!!!!!!!!!!!!!"
    #print best_record.seq
    #best_record.seq = best_record.seq.remove("*")
    if not best_hit:
        return False
    else:

        best = db_prot_data[best_hit]
        best.name = best.name.split('_')[0]
        split_name = best.id.split('_')
        if len(split_name)>2:
            best.id = split_name[0] + '_' + split_name[1]
        else:
            best.id = split_name[0]



        return db_prot_data[best_hit]



def tblastn_and_extract(query_fasta_prot, db_nucl, seq_header = "query_seq"):
    '''return longest seq among the best hits (consider only the best hit of each query_fasta_prot)
    1. check if all best hits are localized in the same region
    2. return the longest (without gaps)

    '''
    from Bio.Blast.Applications import NcbitblastnCommandline
    from Bio.Blast import NCBIXML
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    tblastn_cline = NcbitblastnCommandline(query=query_fasta_prot, db=db_nucl, evalue=0.001, outfmt=5, out="temp.xml")
    stdout, stderr = tblastn_cline()

    result_handle = open("temp.xml")
    blast_records = NCBIXML.parse(result_handle)
    blast_records = list(blast_records)
    best_hit = blast_records[0]
    best_hit_length = len(best_hit.alignments[0].hsps[0].sbjct.replace("-", ""))
    for i in range(1, len(blast_records)):
        blast_record = blast_records[i]
        #print blast_record.query, "vs", blast_record.alignments[0].hit_def
        #print blast_record.alignments[0].hsps[0].sbjct.replace("-", "")

        if "*" in blast_record.alignments[0].hsps[0].sbjct:
            print "Achtung, stopcodon in the align"
            continue
            
        
        elif not check_blast_colocalization(best_hit, blast_record):
            print "problem with colocalization"
            continue
        else:
            hit_len = len(blast_record.alignments[0].hsps[0].sbjct.replace("-", ""))
            if hit_len > best_hit_length:
                best_hit_length = hit_len
                best_hit = blast_record


    descript = "(" + best_hit.query + " vs " + best_hit.alignments[0].hit_def + ")"
    print descript
    seq = best_hit.alignments[0].hsps[0].sbjct.replace("-", "")
    biorecord = SeqRecord(seq=Seq(seq, IUPAC.protein),
                          id=seq_header,
                          name=seq_header,
                          description=descript,
                          dbxrefs=[])
    return biorecord



def format_out(id_matrixes, genes_list, ids):
    out_str = "gene\t"
    #ids = [i.name for i in merged_record]
    out_str += "\t".join(ids) + "\n"
    for gene, matrix in zip(genes_list, id_matrixes):
        if matrix is None:
            continue
        identities = [str(round(y, 2)) for y in matrix[-1]]
        out_str += gene + "\t" + "\t".join(identities[:-1]) + "\n"
    return out_str



def main(protein_multi_fasta, fasta_files, seq_header, out_name, blast_p = False, hmmer = False, reannotate=False):
    import os
    import re
    genes_list = []
    all_mat = []
    target_id_list = []



    if reannotate:
        aa_fasta = []
        if blast_p or hmmer:
            for one_fasta in fasta_files:
                out = one_fasta.split('.')[0]+'_prodig.fa'
                aa_fasta.append(out)
                run_prodigal(one_fasta,output_name=out)
                if blast_p:
                    #shell_command.shell_command("formatdb -i %s -p T" % out)
                    pass
    else:
        aa_fasta = fasta_files
    report_handle = open('classification_report.txt', 'w')

    protein2genome2presence = {}

    marker2genome2best_hit = {}

    for one_protein in protein_multi_fasta:

        protein_id = os.path.basename(one_protein).split('.')[0]
        protein2genome2presence[protein_id] = {}
        marker2genome2best_hit[protein_id] = {}

        new_records = []
        # for each aa file, make blast
        for one_target_fasta in aa_fasta:

            target_id = re.sub('_prodig',
                               '',
                               os.path.basename(one_target_fasta).split('.')[0])

            #genes_list.append(os.path.basename(one_protein).split(".")[0])
            if not blast_p and not hmmer:
                new_record = [tblastn_and_extract(one_protein, one_target_fasta, seq_header)]
            if hmmer:
                new_record = hmm_and_extract(one_protein, one_target_fasta, seq_header)
                if new_record:
                    new_records.append(new_record)
                    protein2genome2presence[protein_id][target_id] = 1
                    marker2genome2best_hit[protein_id][target_id] = new_record.id
                else:
                    print 'no hits', one_protein, one_target_fasta
                    report_handle.write('No hit:\t%s\t%s\n' % (protein_id,
                                                               target_id))
                    protein2genome2presence[protein_id][target_id] = 0
                    marker2genome2best_hit[protein_id][target_id] = '-'
            if blast_p:
                new_record = blastp_and_extract(one_protein, one_target_fasta, seq_header)
                if new_record:
                    new_records.append(new_record)
            #handle = open(one_protein, "rU")
            #initial_records = [record for record in SeqIO.parse(handle, "fasta")]

        if len(new_records)>0:
            f = open("new_strain_" + os.path.basename(one_protein).split('.')[0] + '.faa', 'w')
            SeqIO.write(new_records, f, 'fasta')
            #id_matrix = pairwiseid_needle.get_identity_matrix_from_multifasta(new_records)
            #all_mat.append(id_matrix)
        else:
            #all_mat.append(None)
            print "no hit for %s" % one_protein

        #out_name = one_protein.split('.')[0] + '_new.fa'
        #f = open(out_name, "w")
        #f.write(format_out(all_mat, genes_list, target_id_list))
        #f.close()
    report_handle.close()
    print protein2genome2presence

    with open("presence_absence_matrix.tab", 'w') as m:
        protein_list = protein2genome2presence.keys()
        m.write('\t' + '\t'.join(protein_list) + '\n')
        for genome in protein2genome2presence[protein_list[0]].keys():
            one_list = [str(protein2genome2presence[i][genome]) for i in protein_list]
            m.write(genome + '\t' + '\t'.join(one_list) + '\n')
    with open("locus_table.tab", 'w') as n:
        protein_list = marker2genome2best_hit.keys()
        n.write('\t' + '\t'.join(protein_list) + '\n')
        for genome in marker2genome2best_hit[protein_list[0]].keys():
            one_list = [str(marker2genome2best_hit[i][genome]) for i in protein_list]
            n.write(genome + '\t' + '\t'.join(one_list) + '\n')

    return marker2genome2best_hit



if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--fasta_list', type=str, help="input fasta file: nucl (reannotation/6 frame translation) or protein (no reannotation)", nargs='+')
    parser.add_argument("-q", '--prot_fasta', type=str, help="input protein fasta/hmm profiles", nargs='+')
    parser.add_argument("-n", '--seq_header', type=str, help="seq header name")
    parser.add_argument("-o", '--out_name', type=str, help="output_identity_matrix", default = "id_matrix.txt")
    parser.add_argument("-p", '--blast_p', action="store_true", help="perform ORFing with prodigal and blastP search")
    parser.add_argument("-m", '--hmmer', action="store_true", help="perform ORFing with prodigal and hmm search search")
    parser.add_argument("-r", '--reanotate', action="store_true", help="reannotate with prodigual (default=False, aa input)")
    parser.add_argument("-s", '--six_trame_translation', default=False, help="biodb_name for six frame tranlsation complete ORFing (min size of 30aa)")


    args = parser.parse_args()


    if args.blast_p and args.hmmer:
        raise('use either blastp or hmm searches, not both!')
    elif args.blast_p:
        main(args.prot_fasta, args.fasta_list, args.seq_header, args.out_name, True)
    elif args.hmmer and not args.six_trame_translation:
        import biosql_own_sql_tables
        print 'hmm!---'
        marker2genome2best_hit = main(args.prot_fasta,
                                      args.fasta_list,
                                      args.seq_header,
                                      args.out_name,
                                      False,
                                      True,
                                      args.reanotate)

        dico_seq = {}
        for genome in marker2genome2best_hit[marker2genome2best_hit.keys()[0]]:
            genome_file = '%s.ffn' % genome
            print genome_file
            dico_seq.update(SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta')))

        for marker in marker2genome2best_hit:
            locus_list = []
            for genome in marker2genome2best_hit[marker]:
                if marker2genome2best_hit[marker][genome] != '-':
                    locus_list.append(marker2genome2best_hit[marker][genome])

            # récuperation des séquences avec mysql
            ''''
            records = biosql_own_sql_tables.locus_list2nucleotide_fasta("chlamydia_04_16", locus_list)

            with open('%s_best_hits.ffn' % marker, 'w') as m:
                #SeqIO.write(records, m, 'fasta')

                for i in records:
                    m.write(">%s|%s\n%s\n" % (i.id, i.description, str(i.seq)))
                #with open('%s_best_hits.ffn' % marker, 'w') as m:
                #    SeqIO.write(records, m, 'fasta')
            '''
            # a adapter: possibilité de fournir les ffn correspondant
            with open('%s_best_hits.ffn' % marker, 'w') as m:
                record_l = []
                for locus in locus_list:
                    record_l.append(dico_seq[locus])
                SeqIO.write(record_l, m, 'fasta')

    elif args.reanotate and not args.six_trame_translation:
        shell_command.shell_command("formatdb -i %s -p F" % args.fasta_db)
        main(args.prot_fasta, args.fasta_db, args.seq_header, args.out_name, args)
    elif args.six_trame_translation and args.hmmer:
        import universal_markers_hmm
        fasta_genome_list = universal_markers_hmm.biodatabase2six_frame_translation_all(args.six_trame_translation)
        print 'six frame trans!!!'
        print fasta_genome_list
        marker2genome2best_hit = main(args.prot_fasta,
                                      fasta_genome_list,
                                      args.seq_header,
                                      args.out_name,
                                      False,
                                      True,
                                      args.reanotate)

        dico_seq = {}
        for genome in marker2genome2best_hit[marker2genome2best_hit.keys()[0]]:
            genome_file = '%s.ffn' % genome
            print genome_file
            dico_seq.update(SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta')))

        for marker in marker2genome2best_hit:
            locus_list = []
            for genome in marker2genome2best_hit[marker]:
                if marker2genome2best_hit[marker][genome] != '-':
                    locus_list.append(marker2genome2best_hit[marker][genome])

            # récuperation des séquences avec mysql
            ''''
            records = biosql_own_sql_tables.locus_list2nucleotide_fasta("chlamydia_04_16", locus_list)

            with open('%s_best_hits.ffn' % marker, 'w') as m:
                #SeqIO.write(records, m, 'fasta')

                for i in records:
                    m.write(">%s|%s\n%s\n" % (i.id, i.description, str(i.seq)))
                #with open('%s_best_hits.ffn' % marker, 'w') as m:
                #    SeqIO.write(records, m, 'fasta')
            '''
            # a adapter: possibilité de fournir les ffn correspondant
            with open('%s_best_hits.ffn' % marker, 'w') as m:
                record_l = []
                for locus in locus_list:
                    record_l.append(dico_seq[locus])
                SeqIO.write(record_l, m, 'fasta')

    else:
        print 'wrong combination of args'



