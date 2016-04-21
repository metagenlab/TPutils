#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Edit embl file from RAST to work with embl submission system. Example usage: ' \
                                                 'RASTembl_format.py -i 2718.4.embl -l locux_prefixXYZ -o "Cardiobacterium hominis" -t 573 -p ' \
                                                 ' PRJEB13157 -d "Draft genome sequence of Cardiobacterium hominis"')
    parser.add_argument("-i",'--input',type=str,help="input embl", required=True)
    parser.add_argument("-l",'--locus',type=str,help="locus_tag_prefix", required=True)
    parser.add_argument("-o",'--organism',type=str,help="organism name", required=True)
    parser.add_argument("-t",'--taxon_id',type=str,help="taxonomic id", required=True)
    parser.add_argument("-p",'--project',type=str,help="project ID", required=True)
    parser.add_argument("-d",'--description',type=str,help="description", required=True)

    args = parser.parse_args()

    from Bio import SeqIO
    from Bio.SeqFeature import Reference
    from Bio.SeqFeature import FeatureLocation
    import datetime
    today = datetime.date.today()
    handle = open(args.input, "rU")
    gbk_records = [i for i in SeqIO.parse(handle, "embl")]

    print 'n records', len(gbk_records)

    type_list = []
    updated_records = []
    contig_list = []
    x = 1

    for record in gbk_records:

        '''
        ID   XXX; XXX; linear; genomic DNA; STD; PROK; 1445 BP.
        XX
        AC   XXX;
        XX
        AC * _gnlCHUVkpeu_contig000001
        XX
        PR   Project:PRJEB12345; ok
        XX
        DE   XXX; ok
        XX
        RN   [1] ok
        RA   Submitter, A.; ok
        RT   "Bacullis sp. strain XYZ genome annotated using Prokka."; ok
        RL   Submitted (18-Apr-2016) to the INSDC. ok
        XX
        '''

        record.id="XXX"
        record.name = 'XXX'
        contig_name = record.description.split('Contig ')[1].split(' ')[0]
        contig_list.append(contig_name)
        record.description = args.description
        record.dbxrefs.append("Project:%s" % args.project)
        record.annotations['accessions'] = ['XXX', 'contig']
        record.annotations["data_file_division"] = 'XXX'
        record.annotations["references"] = [Reference()]
        record.annotations["references"][0].authors='XXX'
        record.annotations["references"][0].location = [FeatureLocation(0, len(record))]
        record.annotations["references"][0].title = ''
        record.annotations["references"][0].journal='Submitted (%s) to the INSDC.' % today.strftime('%d-%b-%Y')
        new_features = []
        for i in range(0, len(record.features)):
            type_list.append(record.features[i].type)
            if record.features[i].type == 'source':
                del record.features[i].qualifiers['project']
                del record.features[i].qualifiers['genome_md5']
                del record.features[i].qualifiers['genome_id']
                record.features[i].qualifiers['db_xref'] = ['taxon:%s' % args.taxon_id]
                record.features[i].qualifiers['organism'] = [args.organism]
            elif record.features[i].type in ["CDS", "tRNA","rRNA"]:
                record.features[i].qualifiers["locus_tag"]= ["%s_%04d" % (args.locus, x)]#"CHUV0807_%s" % x
                x+=1
            else:
                raise IOError('Unkown feature type: %s, consider to modify this script to take it into account' % record.features[i].type)
            if record.features[i].type == "CDS":
                if 'x' in record.features[i].qualifiers['translation'][0]:
                    record.features[i].qualifiers['translation'][0] = record.features[i].qualifiers['translation'][0].replace('x','')
                if len(record.features[i].location)/3 != len(record.features[i].qualifiers['translation'][0])+1:
                    print 'Discrepancy between featurelocation length and translation length, removing feature:'
                    print record.features[i]
                    x-=1
                    continue
            new_features.append(record.features[i])
        record.features = new_features
        updated_records.append(record)
    print 'n locus: %s' % x
    out_name = args.input.split(".")[0] + "_locus.tmp"
    handle2 = open(out_name, "w")
    SeqIO.write(updated_records, handle2, "embl")
    handle2.close()
    out_handle =  open(args.input.split(".")[0] + "_locus.embl", "w")
    for type in set(type_list):
        if type not in ["CDS", "tRNA","rRNA", 'source']:
            raise IOError('Unexpected feature type: %s' % type)
    #import sys
    i = 0
    with open(args.input.split(".")[0] + "_locus.tmp") as f:
        match = False
        for line in f:
            if match:
                line += 'AC * _Contig_%s_%s\nXX\n' % (i+1, contig_list[i])
                i+=1
                match=False
            if 'ID   ' in line:
                line = 'ID   XXX; XXX; linear; XXX; XXX; XXX; XXX.\n'
            if 'AC   XXX;' in line:
                line = 'AC   ;\n'
                match = True
            if 'RA   XXX;' in line:
                line +='RT   ;\n'
            if 'OS   .' in line or 'OC   .' in line:
                continue

            out_handle.write(line)
    import shell_command
    shell_command.shell_command('rm %s' % (args.input.split(".")[0] + "_locus.tmp"))