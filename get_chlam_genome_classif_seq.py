#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pairwiseid_needle

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


def run_prodigal(fasta_seq):
    from Bio import SeqIO
    import shell_command
    import StringIO
    from tempfile import NamedTemporaryFile
    cmd = "prodigal -q -a temp.faa -i %s" % fasta_seq
    print cmd
    sdt_out, sdt_err, err = shell_command.shell_command(cmd)
    print sdt_out
    print sdt_err
    shell_command.shell_command('sed -i "s/*//g" temp.faa')
    #print sdt_out, sdt_err, err
    #shell_command.shell_command("seqret -sequence %s -feature -fformat gff -fopenfile temp.gff -osformat genbank -auto -outseq temp.gbk" % fasta_seq)
    #print sdt_out
    #fasta_file = NamedTemporaryFile()
    #fasta = open("temp.faa", 'w')
    #fasta.write(sdt_out)
    
    #print "genbank", genbank
    #for i in genbank:
    #    print "record", i

    #test = open("test.gbk", 'w')
    #test.write(sdt_out)
    
    #for i in genbank:
    #    print i
    #records = [i for i in genbank]
    #print records
    #SeqIO.write(genbank, fasta, "fasta")
    #fasta.close()
    #fasta_file.flush()
    return "temp.faa"


def blastp_and_extract(query_fasta_prot, db_prot, seq_header = "query_seq"):
    from Bio.Blast.Applications import NcbiblastpCommandline
    from Bio.Blast import NCBIXML
    from Bio import SeqIO

    print 'db_prot', db_prot
    print 'querxy', query_fasta_prot
    db_prot_data = SeqIO.to_dict(SeqIO.parse(db_prot, "fasta"))
    
    blastp_cline = NcbiblastpCommandline(query=query_fasta_prot, db=db_prot, evalue=0.001, outfmt=5, out="temp.xml")
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



def main(protein_multi_fasta, fasta_db, seq_header, out_name, blast_p = False):
    genes_list = []
    all_mat = []
    for one_fasta in protein_multi_fasta:
        print one_fasta
        genes_list.append(os.path.basename(one_fasta).split(".")[0])
        if not blast_p:
            new_record = [tblastn_and_extract(one_fasta, fasta_db, seq_header)]
        else:
            new_record = [blastp_and_extract(one_fasta, fasta_db, seq_header)]
        handle = open(one_fasta, "rU")
        initial_records = [record for record in SeqIO.parse(handle, "fasta")]
        print "new_record", new_record
        if new_record[0]:
            merged_record = initial_records + new_record
            f = open("new_strain_" + os.path.basename(one_fasta), 'w')
            SeqIO.write(merged_record, f, 'fasta')
            id_matrix = pairwiseid_needle.get_identity_matrix_from_multifasta(merged_record)
            all_mat.append(id_matrix)
        else:
            all_mat.append(None)
            print "no hit for %s" % one_fasta
    seq_ids = [i.name for i in initial_records]

    f = open(out_name, "w")
    f.write(format_out(all_mat, genes_list, seq_ids))
    f.close()



if __name__ == '__main__':

    import argparse
    import shell_command
    from Bio import SeqIO
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--fasta_db', type=str, help="input fasta file")
    parser.add_argument("-q", '--prot_fasta', type=str, help="input protein fasta", nargs='+')
    parser.add_argument("-n", '--seq_header', type=str, help="seq header name")
    parser.add_argument("-o", '--out_name', type=str, help="output_identity_matrix", default = "id_matrix.txt")
    parser.add_argument("-p", '--blast_p', action="store_true", help="perform ORFing with prodigal and blastP search")


    args = parser.parse_args()

    if args.blast_p:
        blast_db = run_prodigal(args.fasta_db)
        shell_command.shell_command("formatdb -i %s -p T" % blast_db)
        main(args.prot_fasta, blast_db, args.seq_header, args.out_name, True)
    else:
        shell_command.shell_command("formatdb -i %s -p F" % args.fasta_db)
        main(args.prot_fasta, args.fasta_db, args.seq_header, args.out_name)


