


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