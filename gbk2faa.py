#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------


def gbk2faa(seq_records, pformat=False, lformat=False, remove_redundancy=False, get_translation=False, outname=False):
    import re, sys

    all_locus_ids = []
    all_prot_ids = []
    for record in seq_records:
        input_handle = open(record, "rU")
        seq_records = list(SeqIO.parse(input_handle, "genbank"))

        length_records = [len(i.seq) for i in seq_records]
        longest_record = length_records.index(max(length_records))

        if not outname:
            outname = seq_records[longest_record].id.split('.')[0] + ".faa"

        output_handle = open(outname, "w")

        if pformat and lformat:
            raise('You have to chose either protein id or locus tag as header, not both!')

        for record in seq_records:
            description = record.description
            description = re.sub(", complete genome\.", "", description)
            description = re.sub(", complete sequence\.", "", description)
            description = re.sub("strain ", "", description)
            description = re.sub("str\. ", "", description)
            description = re.sub(" complete genome sequence\.", "", description)
            description = re.sub(" complete genome\.", "", description)
            description = re.sub(" chromosome", "", description)
            description = re.sub(" DNA", "S.", description)




            for seq_feature in record.features:
                if 'pseudo' in seq_feature.qualifiers:
                    continue
                if seq_feature.type == "CDS":
                    try:
                        if seq_feature.qualifiers["protein_id"][0] in all_prot_ids:
                            print '%s (%s), Duplicated protein id: %s' % (record.id,
                                                                     record.description,
                                                                     seq_feature.qualifiers["protein_id"][0])
                            if remove_redundancy:
                                print 'skipping...'
                                continue
                        else:
                            all_prot_ids.append(seq_feature.qualifiers["protein_id"][0])
                    except KeyError:
                        # no protein ids
                        pass
                    try:
                        if seq_feature.qualifiers["locus_tag"][0] in all_locus_ids:
                            print '%s (%s), Duplicated locus id: %s, skipping' % (record.id,
                                                                     record.description,
                                                                     seq_feature.qualifiers["locus_tag"][0])
                            # if remove redundancy, skip writing
                            if remove_redundancy:
                                print 'skipping...'
                                continue
                        else:
                            all_locus_ids.append(seq_feature.qualifiers["locus_tag"][0])
                    except KeyError:
                        # no protein ids
                        pass
                    # check presence of a protein sequence
                    try:
                        len(seq_feature.qualifiers['translation'])
                    except:
                        #print seq_feature
                        #sys.stderr.write(seq_feature.location.start)
		        sys.stderr.write("%s (%s) - %s: no sequence, is it a pseudogene, or genbank without translated CDS?\n" % (record.id,
		                                                             record.description,
		                                                             seq_feature.qualifiers["locus_tag"][0]))
			if get_translation:
                                import re
                                seq = str(seq_feature.extract(record.seq).translate())
				if seq[-1] == '*':
                                    seq = seq[0:-1]
				seq_feature.qualifiers['translation'] = [seq]
				print seq
			else:
				continue


                        
                    #assert len(seq_feature.qualifiers['translation'])==1
                    # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]

                    if pformat:
                        # output only protein_id and species name
                        try:
                            output_handle.write(">%s %s\n%s\n" % (
                                    seq_feature.qualifiers["protein_id"][0],
                                                    description,
                                                    seq_feature.qualifiers['translation'][0]))
                        except:
                            try:

                                output_handle.write(">%s %s\n%s\n" % (
                                    seq_feature.qualifiers["protein_id"][0],
                                                    description,
                                                    seq_feature.qualifiers['translation'][0]))
                            except:
                                #print seq_feature
                                pass
                    elif lformat:
                        # output only locus tag and species name
                        try:
                            output_handle.write(">%s %s\n%s\n" % (
                                    seq_feature.qualifiers["locus_tag"][0],
                                                    description,
                                                    seq_feature.qualifiers['translation'][0]))
                        except:
                            try:

                                output_handle.write(">%s %s\n%s\n" % (
                                    seq_feature.qualifiers["protein_id"][0],
                                                    description,
                                                    seq_feature.qualifiers['translation'][0]))
                            except:
                                #print seq_feature
                                pass
                    else:

                        try:
                            output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                                    seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                    seq_feature.qualifiers["protein_id"][0],
                                    seq_feature.qualifiers["product"][0],
                                    description,
                                    seq_feature.qualifiers['translation'][0]))
                        except:
                            try:

                                output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                    seq_feature.qualifiers["db_xref"][0].split(":")[1],
                                    seq_feature.qualifiers["protein_id"][0],
                                    #seq_feature.qualifiers["note"][0],
                                    description,
                                    seq_feature.qualifiers['translation'][0]))
                            except:
                                try:
                                    output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                                        seq_feature.qualifiers["locus_tag"][0],
                                        seq_feature.qualifiers["locus_tag"][0],
                                        #seq_feature.qualifiers["note"][0],
                                        description,
                                        seq_feature.qualifiers['translation'][0]))
                                except:
                                    if 'gene' in seq_feature.qualifiers:
                                        output_handle.write(">%s %s %s \n%s\n" % (
                                        seq_feature.qualifiers["gene"][0],
                                        seq_feature.qualifiers["product"][0],
                                        description,
                                        seq_feature.qualifiers['translation'][0]))
                                    else:
                                        output_handle.write(">%s %s\n%s\n" % (
                                        seq_feature.qualifiers["product"][0],
                                        description,
                                        seq_feature.qualifiers['translation'][0]))                                        



if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_gbk', type=str, help="input gbk file", nargs='+')
    #parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-f", '--lformat', action='store_true', help="format header: >locus description", default=False)
    parser.add_argument("-p", '--pformat', action='store_true', help="format header: >protein_id description", default=False)
    parser.add_argument("-o", '--outname', help="outname (optinal)", default=False)
    parser.add_argument("-r", '--remove', action='store_true', help="remove redundancy (protein id or locus tags persent mor than once)", default=False)
    parser.add_argument("-t", '--translate', action='store_true', help="translate from DNA if translation not available (not for pseudo tagged features)", default=False)

    args = parser.parse_args()






    if args.lformat:
        gbk2faa(args.input_gbk, lformat=args.lformat, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
    elif args.pformat:
        gbk2faa(args.input_gbk, pformat=args.pformat, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
    else:
        gbk2faa(args.input_gbk, False, False, remove_redundancy=args.remove, get_translation=args.translate, outname=args.outname)
