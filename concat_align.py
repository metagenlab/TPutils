#! /usr/bin/env python

all_locus = [
"wcw_1683",
"wcw_1636",
"wcw_1751",
"wcw_0277",
"wcw_1006",
"wcw_1572",
"wcw_0128",
"wcw_1710",
"wcw_1223",
"wcw_1671"
]



if __name__ == '__main__':
    import argparse
    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Align import MultipleSeqAlignment
    import shell_command
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",'--fasta',type=str,help="input fasta files ", nargs="+")


    keep1 =[
        60,
        137,
        67,
        68,
        64,
        55,
        50,
        292,
        291,
        125,
        135,
        284]

    
    keep2 =[
        136,
        137,
        287,
        285,
        286,
        139,
        138,
        241,
        290,
        62,
        66,
        68 ,
        67,
        64,
        55,
        48,
        283,
        124,
        122,
        123,
        47,
        46,
        45,
        52 ,
        49,
        50,
        291,
        292,
        125,
        135,
        289,
        288,
        60,
        291,
        284,
        135,
        125,
        292]

    args = parser.parse_args()

    taxons = []
    all_seq_data = {}
    for one_fasta in args.fasta:
        all_seq_data[one_fasta] = {}
        alignment = AlignIO.read(open(one_fasta), "fasta")
        for record in alignment:
            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record
    print taxons
    print all_seq_data


    concat_data = {}

    for one_fasta in args.fasta:
        for taxon in taxons:
            if taxon not in all_seq_data[one_fasta]:
                print taxon, "not in", one_fasta
                seq = Seq("-"*len(all_seq_data[one_fasta][all_seq_data[one_fasta].keys()[0]]))
                all_seq_data[one_fasta][taxon] = SeqRecord(seq, id=taxon)
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    for i in concat_data:
        print len(concat_data[i])#, concat_data[i]
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])

    AlignIO.write(MSA, "msa.fa", "fasta")

    '''

    import manipulate_biosqldb
    server, db =manipulate_biosqldb.load_db("chlamydia_02_15")
    for locus in all_locus:

        orthogroup = manipulate_biosqldb.locus_tag2orthogroup_id(server, locus, "chlamydia_02_15")
        print locus, orthogroup
        a,b,c =shell_command.shell_command("cp /home/trestan/Dropbox/dev/django/test_1/assets/chlamydia_02_15_fasta_by_taxons/%s.txt /home/trestan/Dropbox/dev/concat_align_test" % orthogroup)
        print a,b,c



    for i in args.fasta:
        from Bio import SeqIO
        fasta = SeqIO.parse(i, "fasta")
        new_fasta = []
        for record in fasta:
            if int(record.id) not in keep2:
                continue
            else:
                new_fasta.append(record)
        out_name = i.split(".")[0] + "_subset.txt"
        with open(out_name, "w") as f:
            SeqIO.write(new_fasta,f, "fasta")

    '''
