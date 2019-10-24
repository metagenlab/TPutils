#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def parse_blast(blast_result_file_list):
    blast2data = {}
    for one_blast_file in blast_result_file_list:
        with open(one_blast_file, 'r') as f:
            for line in f:
                line = line.split('\t')
                # if new record
                if one_blast_file not in blast2data:
                    blast2data[one_blast_file] = {}
                    # keep identity, start, end, name
                    blast2data[one_blast_file][line[0]] = [float(line[4]), int(line[11]), int(line[12]), line[0]]
                else:
                     blast2data[one_blast_file][line[0]] = [float(line[4]), int(line[11]), int(line[12]), line[0]]
    return blast2data

def remove_record_redundancy(blast2data):
    # check if multiple proteins in the database have hits again the same region
    # keep only the best hit

    for one_blast in blast2data.keys():
        for ref_gene in blast2data[one_blast].keys():

            for query_gene in blast2data[one_blast].keys():
                overlap = False
                if ref_gene == query_gene:
                    continue
                # check if position is overlapping
                try:
                    if are_blast_hits_overlaping(blast2data[one_blast][ref_gene], blast2data[one_blast][query_gene]):


                        if blast2data[one_blast][ref_gene][0] > blast2data[one_blast][query_gene][0]:
                            del blast2data[one_blast][query_gene]

                        else:
                            del blast2data[one_blast][ref_gene]

                            break
                except KeyError:
                    pass
    return blast2data


def get_colocating_blast(blast2data):

    result_list = []

    colocalisation_done = []

    #print blast2data.keys()

    for i, one_blast_record in enumerate(blast2data):
        # copare with all other blast files
        #print i, one_blast_record
	
        for n, first_blast in enumerate(blast2data[one_blast_record]):

      
	    if blast2data[one_blast_record][first_blast] in colocalisation_done:
                continue

            temp_res = ['-']*len(blast2data)
            #print 'first blast %s & %s/%s' % (first_blast, n,  len(blast2data[one_blast_record]))

            for y, second_blast_record in enumerate(blast2data):

                # if blast was already seen, skip it (may be replaced by a remove?)
                if blast2data[one_blast_record][first_blast] not in colocalisation_done:
                    # keep the name of the blast hit

                    for second_blast in blast2data[second_blast_record]:
		        #if blast2data[second_blast_record][second_blast][3] == 'ClpC' or blast2data[one_blast_record][first_blast][3] == 'ClpC':
		        #    #print blast2data[one_blast_record][first_blast], blast2data[second_blast_record][second_blast], are_blast_hits_overlaping(blast2data[one_blast_record][first_blast], blast2data[second_blast_record][second_blast])
                        if second_blast == first_blast:
                            if float(blast2data[second_blast_record][second_blast][0])>=80:
                                temp_res[y] = "%s (%s)" % (blast2data[second_blast_record][second_blast][3],blast2data[second_blast_record][second_blast][0])
				#colocalisation_done.append(blast2data[second_blast_record][second_blast])
                            continue

                        try:
                            coloc = are_blast_hits_overlaping(blast2data[one_blast_record][first_blast], blast2data[second_blast_record][second_blast])


                            if coloc:
                                #print one_blast_record, second_blast_record, blast2data[one_blast_record][first_blast][3], blast2data[second_blast_record][second_blast][3], coloc
                                # keep description
                                if float(blast2data[second_blast_record][second_blast][0])>=80:
                                    temp_res[y] = "%s (%s)" % (blast2data[second_blast_record][second_blast][3], blast2data[second_blast_record][second_blast][0])
				    # keep second hit in memory (and not not consider it next time)
                                    colocalisation_done.append(blast2data[second_blast_record][second_blast])
                                    # break the loop (we can only have one hit colocating, as hits matching the same region were removed)

                        except:
                            print "one_blast_record", one_blast_record
                            print "first_blast", first_blast
                            print "second_blast_record", second_blast_record
                            print "second_blast", second_blast
            #print temp_res
            if temp_res not in result_list:
                result_list.append(temp_res)
    for i in result_list:
        if i.count('-') == len(i):
            pass
        else:
            print '\t'.join(i)






def are_blast_hits_overlaping(blast1, blast2):

    # test si start de blast1 = start blast2
    # ou si start de blast1 se trouve entre start2 et stop2
    sorted_coordinates1 = sorted(blast1[1:3])
    sorted_coordinates2 = sorted(blast2[1:3])
    if sorted_coordinates2[1] <= sorted_coordinates1[1] and sorted_coordinates2[1]>= sorted_coordinates1[0]:
	length_1 = sorted_coordinates2[1]-sorted_coordinates2[0] 
	length_2 = sorted_coordinates1[1]-sorted_coordinates1[0]
	#print min(length_1, length_2)/float(max(length_1, length_2))
	if min(length_1, length_2)/float(max(length_1, length_2)) >= 0.5:
            return True
	else:
	    return False
    # test du cas inverse
    if sorted_coordinates1[1] <= sorted_coordinates2[1] and sorted_coordinates1[1]>= sorted_coordinates2[0]:
	length_1 = sorted_coordinates2[1]-sorted_coordinates2[0] 
	length_2 = sorted_coordinates1[1]-sorted_coordinates1[0]
	#print min(length_1, length_2)/float(max(length_1, length_2))
	if min(length_1, length_2)/float(max(length_1, length_2)) >= 0.5:
            return True
	else:
	    return False

        return True
    else:
        return False


def main(blast_files):

    blast2data_raw = parse_blast(blast_files)

    blast2data = remove_record_redundancy(blast2data_raw)

    get_colocating_blast(blast2data)





def locus2description(fasta_files):
    from Bio import SeqIO

    locus2d = {}
    for fasta in fasta_files:
        with open(fasta) as f:
            for i in SeqIO.parse(f, 'fasta'):
                locus2d[i.name] = ' '.join(i.description.split(' ')[1:])

    return locus2d




if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    import gbk2accessiontodefinition
    parser = argparse.ArgumentParser()
    parser.add_argument("-b",'--blast_files', type=str, help="blast files", nargs="+")
    parser.add_argument("-f",'--fasta_files', type=str, help="fasta files", nargs="+")

    args = parser.parse_args()

    if args.fasta_files:
        l2d = locus2description(args.fasta_files)
        with open('comp_table.tab', 'r') as f:
            for i in f:

                line = i.rstrip().split('\t')
                for y, value in enumerate(line):
                    if 'SaCp_00021' in value:
                        print line, '#################################'
                    if value == '-':
                        pass
                    else:
                        try:
                            line[y] = l2d[value.split(' ')[0]]
                        except:
                            pass
                print line
                print '\t'.join(line)


    if args.blast_files:

        main(args.blast_files)

