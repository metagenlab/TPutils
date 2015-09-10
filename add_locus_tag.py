#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--input',type=str,help="input genbank")
    parser.add_argument("-l",'--locus',type=str,help="locus_tag_prefix")
    parser.add_argument("-c",'--check_locus',action='store_true', help="check if a genbank file contain locus tags")


    args = parser.parse_args()
    if not args.check_locus:
        from Bio import SeqIO
        handle = open(args.input, "rU")
        for record in SeqIO.parse(handle, "genbank") :

            x = 1
            for i in range(0, len(record.features)):
                print "x", x

                if record.features[i].type == "gene":
                    print record.features[i].type, record.features[i+1].type
                    #xprint record.features[i].qualifiers
                    record.features[i].qualifiers["locus_tag"]= ["%s_%s" % (args.locus,x)]
                    y = 1
                    try:
                        while record.features[i+y].type != "gene":
                            print 'y', y
                            #if record.features[i+1].type == "CDS" or "RNA" in record.features[i+1].type:
                            record.features[i+y].qualifiers["locus_tag"] = ["%s_%s" % (args.locus,x)]
                            x+=1
                            y+=1
                    except:
                        print 'n locus', x

                    y = 1
            if x == 1:
                y = 0
                try:
                    while record.features[y].type != "gene":
                        print 'y', y
                        #if record.features[i+1].type == "CDS" or "RNA" in record.features[i+1].type:
                        record.features[y].qualifiers["locus_tag"] = ["%s_%s" % (args.locus,x)]
                        x+=1
                        y+=1
                except:
                    pass

                #print record.features[i+1]

            out_name = args.input.split(".")[0] + "_locus.gbk"
            handle2 = open(out_name, "w")
            SeqIO.write(record, handle2, "genbank")
    else:
        import sys
        from Bio import SeqIO
        handle = open(args.input, "rU")
        for record in SeqIO.parse(handle, "genbank") :
            for feature in record.features:
                if feature.type=='CDS':
                    try:
                        test = feature.qualifiers['locus_tag']
                        #print 'locus ok'
                    except:
                        print '############## no locus ######################'
                        print feature
                    #sys.exit()
            print 'locus ok'