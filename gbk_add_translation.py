#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", '--gbk',type=str,help="gbk file")
    parser.add_argument("-f", '--faa',type=str,help="faa file")
    parser.add_argument("-n", '--fna',type=str,help="fna file")
    
    args = parser.parse_args()
    genbank = [i for i in SeqIO.parse(args.gbk, "genbank")]

    if args.gbk is not None and args.faa is None and args.fna is None:
	for record in genbank:
             for feature in record.features: 
                 if feature.type == 'CDS':
                    try:
                        len(feature.qualifiers['translation'])
                    except:
                        import re
                        try:
                            print "adding translation to %s" % feature.qualifiers['locus_tag']
                            seq = str(feature.extract(record.seq).translate())
                            if seq[-1] == '*':
                                seq = seq[0:-1]
			    feature.qualifiers['translation'] = [seq]
                        except KeyError:
                            print '----------------No locus tag for:'
                            print feature
        out_name = args.gbk.split('.')[0] + '_translations.gbk'
        output_handle = open(out_name, 'w')            
        SeqIO.write(genbank,output_handle, "genbank")                   

    if args.faa is not None and args.gbk is not None:
        
        faa = SeqIO.parse(args.faa, "fasta")
        
        locus2seq = {}
        for record in faa:
            locus2seq[record.name] = record.seq
        #print locus2seq

        
        for record in genbank:
            i = 1
            for feature in record.features:
                i+=1
                if feature.type == 'CDS':
                    print i
                    try:
                        feature.qualifiers['translation'] = [str(locus2seq[feature.qualifiers['protein_id'][0]])]
                    except:
                        try:
                            feature.qualifiers['translation'] = [str(locus2seq[feature.qualifiers['locus_tag'][0]])]
                        except:
                            
                            print feature
                        
        output_handle = open("test.gbk", 'w')            
        SeqIO.write(genbank,output_handle, "genbank")
    elif args.faa is None and args.gbk is not None and args.fna is not None:
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna, generic_protein
        fna = [i for i in SeqIO.parse(args.fna, "fasta")]
        print fna[0]
        for record_gbk, record_fna in zip(genbank, fna):
            
            record_gbk.seq = Seq(str(record_fna.seq), generic_dna)
            #print record_gbk.seq
            for feature in record_gbk.features:
                if feature.type == 'CDS':
                    seq = str(feature.extract(record_fna.seq).translate())
		    if seq[-1] == '*':
                        seq = seq[0:-1]                               
                    feature.qualifiers['translation'] = [seq]
                    #print feature
        output_handle = open("test.gbk", 'w')            
        SeqIO.write(genbank,output_handle, "genbank")

    else:
        print 'Either provide gbk file fith corresponding faa or fna file'
