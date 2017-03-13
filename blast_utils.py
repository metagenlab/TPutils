


def run_prodigal(fasta_seq, output_name='temp.faa'):
    from Bio import SeqIO
    import shell_command
    import StringIO
    from tempfile import NamedTemporaryFile
    # -q quiet
    # -a Write protein translations to the selected file
    # -i Specify input file
    # -c:  Closed ends.  Do not allow genes to run off edges. # not activated
    cmd = "prodigal -q -a %s -i %s" % (output_name, fasta_seq)

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


class Hmm():

    def __init__(self, hmm_profiles, database, output_dir=False, call_genes=False, score_cutoff=100, cut_tc=False, filter_bitscore='best'):
        from tempfile import NamedTemporaryFile

        assert isinstance(hmm_profiles, list)

        if isinstance(database, list) and len(database)>1:
            # if multiple databases, compare the bitscore obtained for each profile against each protein
            # use either the best biscrore or the median bitscrore to filter results
            self.bitscore_litering = filter_bitscore
            self.multiple_databases = True
            self.profile2scores = {}

        self.hmm_profiles = hmm_profiles
        self.database = database
        # -T <x>     : report sequences >= this score threshold in output
        # out, query and db
        #self.hmmer_cmd = 'hmmsearch -T %s -E 1e-10 -o %s %s %s'
        if not cut_tc:
            self.hmmer_cmd = 'hmmsearch -T %s -o' % score_cutoff + ' %s %s %s'
        else:
            self.hmmer_cmd = 'hmmsearch --cut_tc -o %s %s %s'
        self.hmmer_score_cutoff = score_cutoff
        self.hmmer_output_list = []

        if call_genes:
            #if not output_dir:

            temp_prodigal_file = NamedTemporaryFile()
            self.database = temp_prodigal_file.name
            # add content to temporary file
            #temp_file.write(str(self))

            run_prodigal(database, output_name=self.database)

    def run_hmmer(self, profiles=False):
        from tempfile import NamedTemporaryFile
        import shell_command


        if not profiles:
            profiles = self.hmm_profiles

        header = ["profile_id",
                "profile_length",
                "best_hit_id",
                "bias",
                "bitscore",
                "evalue",
                "query_start",
                "query_end",
                "query_coverage",
                "hit_start",
                "hit_end"]
        results = []#[header]
        for profile in profiles:
            temp_file = NamedTemporaryFile()
            self.hmmer_output_list.append(temp_file.name)
            if not isinstance(self.database, list):
                cmd = self.hmmer_cmd % (temp_file.name, profile, self.database)
                print cmd
                stout, sterr, code = shell_command.shell_command(cmd) # self.hmmer_score_cutoff,
                if code != 0:
                    import sys
                    sys.stdout.write("\n%s\n%s\n" % (stout, sterr))
                    sys.exit()

                parsed_data = self._parse_hmmsearch(temp_file.name)

                if len(parsed_data) == 0:
                    print 'No domains respecting score threshold for %s, continue...' % profile
                    continue

                if not isinstance(parsed_data[0], dict):
                    results.append(['%s' % parsed_data[0], '-', '-', '-', '-', '-', '-', '-', '-', '-'])
                else:
                    hsp_list = parsed_data
                    for x in range(0,len(hsp_list)):
                        #results += '\t'.join([str(hsp_list[x][i]) for i in header])
                        #results += '\n'
                        results.append([str(hsp_list[x][i]) for i in header])
            else:
                # multiple databases: performing bitscore filtering
                self.biodb2best_hits = {}
                for database in self.database:
                    stout, sterr, code = shell_command.shell_command(self.hmmer_cmd % (self.hmmer_score_cutoff, temp_file.name, profile, self.database))
                    if code != 0:
                        import sys
                        sys.stdout.write("\n%s\n%s\n" % (stout, sterr))
                        sys.exit()

                    parsed_data = self._parse_hmmsearch(temp_file.name)


                    '''
                    if not isinstance(parsed_data[0], dict):
                        pass
                    else:
                        # all hsp have the same bitscore, only use the first hsp
                        if parsed_data[0]['profile_id'] not in self.profile2scores:
                            self.profile2scores[parsed_data[0]['profile_id']] = [parsed_data[0]['bitscore']]
                        else:
                            self.profile2scores[parsed_data[0]['profile_id']].append(parsed_data[0]['bitscore'])
                        hsp_list = parsed_data
                        for x in range(0,len(hsp_list)):
                            results += '\t'.join([str(hsp_list[x][i]) for i in header])
                            results += '\n'
                    '''


        return results

    def _parse_hmmsearch(self, hmmsearch_result):
        from Bio import SearchIO

        result_handle = open(hmmsearch_result, 'r')
        hmmer_records = [i for i in SearchIO.parse(result_handle, 'hmmer3-text')]

        try:
            best_hit_id = hmmer_records[0].hits[0].id
        except IndexError:
            return [hmmer_records[0].id]
        else:

            hsp_list = []
            #print dir(hmmer_records[0])
            profile_id =  hmmer_records[0].id
            '''
            'append', 'bias', 'bitscore', 'description', 'description_all', 'domain_exp_num', 'domain_obs_num',
            'evalue', 'filter', 'fragments', 'hsps', 'id',
            'id_all', 'index', 'is_included', 'map', 'pop', 'query_description', 'query_id', 'sort'

            '''
            for hsp in hmmer_records[0].hits[0].hsps:
                result = {
                    "profile_id" : hmmer_records[0].id,
                    "profile_length" : hmmer_records[0].seq_len,
                    "best_hit_id" : hmmer_records[0].hits[0].id,
                    "bias" : hmmer_records[0].hits[0].bias,
                    "bitscore" : hmmer_records[0].hits[0].bitscore,
                    "hit_start" : hsp.hit_start,
                    "hit_end" : hsp.hit_end,
                    "query_start" : hsp.query_start,
                    "query_end" : hsp.query_end,
                    "evalue" : hmmer_records[0].hits[0].evalue,
                    "query_coverage" : round((hsp.query_end-hsp.query_start)/float(hmmer_records[0].seq_len),2)
                    }
                hsp_list.append(result)
            return hsp_list

    def filter_hmmer_results(self):
        pass



class Blast():

    def __init__(self, query, database, protein=False, formatdb=False, best_hit_only=True):
        import os

        self.query = query
        self.database = database
        self.protein = protein
        self.best_hit_only = best_hit_only
        self.formatdb = formatdb
        self.working_dir = os.getcwd()
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
            shell_command.shell_command('formatdb -i %s -t /tmp/%s.temp -o T -p T -n /tmp/%s.temp -b T' % (self.database, new_database, new_database))
        else:
            shell_command.shell_command('formatdb -i %s -t /tmp/%s.temp -o T -p F -n /tmp/%s.temp -b T' % (self.database, new_database, new_database))

        self.database_path = '/tmp/%s.temp' % new_database

    def run_blastp(self):
        from Bio.Blast.Applications import NcbiblastpCommandline
        import os

        outpath = os.path.join(self.working_dir, 'blast_result.tab')
        blastp_cline = NcbiblastpCommandline(query= self.query,
                                            db=self.database,
                                            evalue=0.005,
                                            outfmt=6,
                                            out=outpath)
        stdout, stderr = blastp_cline()
        print stderr

        with open(outpath, 'r') as result_handle:

            self.best_hit_list = []
            for line in result_handle:
                if line.split('\t')[0] in self.best_hit_list:
                    continue
                else:
                    self.best_hit_list.append(line.rstrip().split('\t'))



    def run_tblastn(self):


        '''
        tblastn_cline = NcbitblastnCommandline(query='dnaa.temp',
                                             db=contig_file,
                                             evalue=0.001,
                                             outfmt=5,
                                             out="dnaa_blast2.xml")

        stdout, stderr = tblastn_cline()

        result_handle = open("dnaa_blast2.xml", 'r')
        blast_records = [i for i in NCBIXML.parse(result_handle)]
        '''
        pass