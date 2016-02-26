


def run_prodigal(fasta_seq, output_name='temp.faa'):
    from Bio import SeqIO
    import shell_command
    import StringIO
    from tempfile import NamedTemporaryFile
    cmd = "prodigal -q -a %s -i %s" % (output_name, fasta_seq)
    print cmd
    sdt_out, sdt_err, err = shell_command.shell_command(cmd)
    print sdt_out
    print sdt_err
    shell_command.shell_command('sed -i "s/*//g" %s' % output_name)
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
    return output_name


class Blast():

    def __init__(self, query, database, protein=False):
        import NCBIXML
        import shell_command
        self.query = query
        self.database = database
        self.protein = protein
        self.blast_path_var= '$BLASTDB:/temp/blastdb.temp'


    def id_generator(self, size=6, chars=False): # + string.digits
        import random
        import string

        if not chars:
            chars = string.ascii_lowercase
        return ''.join(random.choice(chars) for _ in range(size))

    def format_database(self):
        import shell_command

        new_database = self.id_generator(8)

        if self.protein:
            shell_command.shell_command('formatdb -i %s -t /tmp/%s.temp -o T -p F -n /tmp/%s.temp -b T' % (self.database, new_database, new_database))
        else:
            shell_command.shell_command('formatdb -i %s -t /tmp/%s.temp -o T -p T -n /tmp/%s.temp -b T' % (self.database, new_database, new_database))

        self.database = new_database


    def run_tblastn(self):



        tblastn_cline = NcbitblastnCommandline(query='dnaa.temp',
                                             db=contig_file,
                                             evalue=0.001,
                                             outfmt=5,
                                             out="dnaa_blast2.xml")

        stdout, stderr = tblastn_cline()

        result_handle = open("dnaa_blast2.xml", 'r')
        blast_records = [i for i in NCBIXML.parse(result_handle)]
