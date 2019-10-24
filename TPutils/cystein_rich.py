#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input',type=str,help="input genbank")
    parser.add_argument("-l",'--locus',type=str,help="locus_tag_prefix")      
    args = parser.parse_args()

    target_aa = ['U', 'C', 'u', 'c']

    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')

    sql = 'select locus_tag, SP, TM from orthology_detail_chlamydia_03_15'

    sql2 = 'select protein_id, locus_tag from orthology_detail_chlamydia_03_15'

    protein_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    #print protein_id2locus_tag
    data = server.adaptor.execute_and_fetchall(sql,)
    locus_tag2SP_TM = {}
    for i in data:
        locus_tag2SP_TM[i[0]] = [i[1], i[2]]


    from Bio import SeqIO
    handle = open(args.input, "rU")
    print 'protein_id\tlocus\ttransmembrane_domains\tsignal_peptide\tcystein(%)\tn_C_U\tprotein_length\tdescription'
    for record in SeqIO.parse(handle, "fasta") :
        target_n = 0
        protein_length = len(record.seq)
        for aa in record.seq:
            if aa in target_aa:
                target_n +=1
        # gi|297376554|ref|ADI38384.1|
        try:
            protein_id = record.id.split('|')[3]
        except:
            protein_id = record.id
        try:

            locus_tag = protein_id2locus_tag[protein_id]
        except:
            locus_tag = protein_id
        SP = locus_tag2SP_TM[locus_tag][0]
        TM = locus_tag2SP_TM[locus_tag][1]


        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (protein_id, locus_tag, TM, SP, 100*float(target_n)/protein_length, target_n, protein_length, record.description)
        
