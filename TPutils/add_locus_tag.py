#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input',type=str,help="input genbank")
    parser.add_argument("-l",'--locus',type=str,help="locus_tag_prefix")
    parser.add_argument("-c",'--check_locus',action='store_true', help="check if a genbank file contain locus tags")
    parser.add_argument("-g",'--replace_gene',action='store_true', help="Replace missing locus_tag by gene name")


    args = parser.parse_args()

    from Bio import SeqIO
    handle = open(args.input, "rU")
    gbk_records = [i for i in SeqIO.parse(handle, "genbank")]


    if not args.check_locus and not args.replace_gene:
        for record in  gbk_records:
            x = 1
            for i in range(0, len(record.features)):
                print ("x", x)

                if record.features[i].type == "gene":
                    print (record.features[i].type, record.features[i+1].type)
                    #xprint record.features[i].qualifiers
                    record.features[i].qualifiers["locus_tag"]= ["%s_%s" % (args.locus,x)]
                    y = 1
                    try:
                        while record.features[i+y].type != "gene":
                            print ('y', y)
                            #if record.features[i+1].type == "CDS" or "RNA" in record.features[i+1].type:
                            record.features[i+y].qualifiers["locus_tag"] = ["%s_%s" % (args.locus,x)]
                            x+=1
                            y+=1
                    except:
                        print ('n locus', x)

                    y = 1
            if x == 1:
                y = 0
                try:
                    while record.features[y].type != "gene":
                        print ('y', y)
                        # if record.features[i+1].type == "CDS" or "RNA" in record.features[i+1].type:
                        record.features[y].qualifiers["locus_tag"] = ["%s_%s" % (args.locus,x)]
                        x+=1
                        y+=1
                except:
                    pass
            out_name = args.input.split(".")[0] + "_locus.gbk"
            handle2 = open(out_name, "w")
            SeqIO.write(record, handle2, "genbank")

    elif args.replace_gene:
        import sys
        count_CDS = 0
        count_no_locus = 0
        from Bio import SeqIO
        handle = open(args.input, "rU")
        for record in gbk_records :
            for feature in record.features:
                if feature.type=='CDS':
                    count_CDS += 1
                    try:
                        test = feature.qualifiers['locus_tag']
                        # print 'locus ok'
                    except:
                        print ('############## no locus:', feature.qualifiers['gene'])
                        feature.qualifiers["locus_tag"] = feature.qualifiers['gene']
                        count_no_locus += 1
                    #sys.exit()

            print ('Features with updated locus:', count_no_locus)
            print ('Total number of CDS:', count_CDS)
            out_name = args.input.split(".")[0] + "_locus.gbk"
            handle2 = open(out_name, "w")
            SeqIO.write(record, handle2, "genbank")

    else:
        import gbk_check
        for record in gbk_records:
            count_no_locus, count_CDS = gbk_check.count_missing_locus_tags(record)
            print ('Features without locus:', count_no_locus)
            print ('Total number of CDS:', count_CDS)
