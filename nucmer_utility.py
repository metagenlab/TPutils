#!/usr/bin/python



def execute_promer(fasta1,
                   fasta2,
                   algo="nucmer",
                   coords=True,
                   minimum_align_length=100,
                   minimum_identity=30):

    '''



    :param fasta1: reference fata
    :param fasta2: query fasta(s)
    :param algo: nucmer/promer
    :param coords: execute show-coords
    :param minimum_align_length: default:100
    :param minimum_identity: default:30
    :return: coords file list, delta file list
    '''

    import shell_command
    import os

    #cmd1 = 'promer --mum -l 5 %s %s' % (fasta1, fasta2)

    delta_files = []
    coord_files = []
    for one_fasta in fasta2:
        if algo == 'nucmer':
            # 03.2017 => changed from -mum to -mumreference. 01.06.17 removed -mumreference
            cmd1 = 'nucmer -b 200 -c 65 -g 90 -l 20 -p %s %s %s' % (os.path.basename(one_fasta).split('.')[0],
                                                                         fasta1,
                                                                         one_fasta)
            print (cmd1)
            a, b, c = shell_command.shell_command(cmd1)
            if c != 0:
                raise Exception("%s" % b)

            if coords:
                cmd2 = 'show-coords -T -r -c -L %s -I %s %s.delta > %s.coords' % (minimum_align_length,
                                                                                  minimum_identity,
                                                                                  os.path.basename(one_fasta).split('.')[0],
                                                                                  os.path.basename(one_fasta).split('.')[0])
                a, b, c = shell_command.shell_command(cmd2)
                if c != 0:
                    raise Exception("%s" % b)


        elif algo == 'promer':
            # promer --mum -l 5
            cmd1 = 'promer -l 5 -p %s %s %s' % (os.path.basename(one_fasta).split('.')[0],
                                                fasta1,
                                                one_fasta)

            a, b, c = shell_command.shell_command(cmd1)
            if c != 0:
                raise Exception("%s" % b)

            if coords:
                # show-coords -T -r -c -L 100 -I 30 out.delta
                cmd2 = 'show-coords -T -r -c -L %s -I %s %s.delta > %s.coords' % (minimum_align_length,
                                                                                  minimum_identity,
                                                                                  os.path.basename(one_fasta).split('.')[0],
                                                                                  os.path.basename(one_fasta).split('.')[0])

                a, b, c = shell_command.shell_command(cmd2)
                if c != 0:
                    raise Exception("%s" % b)



        delta_files.append('%s.delta' % os.path.basename(one_fasta).split('.')[0])
        if coords:
            coord_files.append('%s.coords' % os.path.basename(one_fasta).split('.')[0])

    return coord_files, delta_files

def delta_file2start_stop_list(delta_input,
                               contigs_add,
                               algo='nucmer',
                               minimum_identity=85):

    '''

    :param coords_input:
    :param contig_add:
    :param algo:
    :return:
    '''


    import re
    import shell_command

    print ('contig add', contigs_add)

    out_name = delta_input.split('.')[0]

    cmd2 = 'show-coords -T -r -c -L 100 -I %s %s > gaps_%s.coords' % (minimum_identity,
                                                                      delta_input,
                                                                      out_name)
    a, b, c = shell_command.shell_command(cmd2)
    if c != 0:
        print ('nucmer error!!!')
        print (a,b,c)
        print (cmd2)
        import sys
        sys.exit()

    with open("gaps_%s.coords" % out_name, 'rU') as infile:
        rawLinks = justLinks(infile.readlines());

    if algo == 'promer':
        shift = 4
    elif algo == 'nucmer':
        shift = 0

    contig2start_stop_list = {}

    i = 1
    for n, row in enumerate(rawLinks):

        l = re.split(r'\t+', row.rstrip('\n'))

        start = int(l[0]) + contigs_add[l[9+shift]][0]
        stop = int(l[1]) + contigs_add[l[9+shift]][0]
        contig_ref = l[9+shift]


        if contig_ref not in contig2start_stop_list:



            contig2start_stop_list[contig_ref] = {}
            contig2start_stop_list[contig_ref]["start"] = [start]
            contig2start_stop_list[contig_ref]["stop"] = [stop]
        else:
            contig2start_stop_list[contig_ref]["start"].append(start)
            contig2start_stop_list[contig_ref]["stop"].append(stop)


    return contig2start_stop_list



def justLinks(coords_input):
    ##Find first row after header
    import re
    i = 0
    for row in coords_input:
        if len(re.split(r'\t+', row)) > 3:
            headerRow = i + 1
            break
        else:
            i += 1

    return coords_input[headerRow:None]


def id_generator(size=6, chars=False ): # + string.digits
    import string
    import random
    if not chars:
        chars = string.ascii_lowercase
    return ''.join(random.choice(chars) for _ in range(size))

def get_contigs_coords(ref_fasta):
    '''

    given a fragmented fasta
    return coordinates of each 'contig' if contantenanted in the order provided

    :param self:
    :param ref_fasta:
    :return:
    '''

    from Bio import SeqIO
    import re
    contigs_add = {}
    fasta_data1 = [i for i in SeqIO.parse(open(ref_fasta), "fasta")]
    contig_end = 0
    for record in fasta_data1:
        name = re.sub("\|", "", record.name)
        contig_start = contig_end+1
        contig_end = contig_start+len(record)
        contigs_add[name] = [contig_start, contig_end]
    return contigs_add


def coord_file2circos_heat_file(coords_input, contigs_add, algo="nucmer"):
    '''

    :param coords_input:
    :param contig_add:
    :param algo:
    :return:
    '''

    import re

    with open(coords_input, 'rU') as infile:
        rawLinks = justLinks(infile.readlines());

    if algo == 'promer':
        shift = 4
    elif algo == 'nucmer':
        shift = 0

    link_file = coords_input.split('.')[0] + '.heat'
    with open(link_file, 'w') as f:
        i = 1
        hit_list = []
        query_list = []

        for n, row in enumerate(rawLinks):
            l = re.split(r'\t+', row.rstrip('\n'))
            if n == 0:
                contig = re.sub("\|", "", l[9+shift])
                start = int(l[0]) + contigs_add[contig][0]
                end = int(l[0])+1 + contigs_add[contig][0]
                f.write("%s\t%s\t%s\t30\tz=0\n" % (contig, start, end))
                f.write("%s\t%s\t%s\t100\tz=0\n" % (contig, start, end))

            contig = re.sub("\|", "", l[9+shift])
            start = int(l[0]) + contigs_add[contig][0]
            end = int(l[1]) + contigs_add[contig][0]
            identity = l[6]
            f.write('%s\t%s\t%s\t%s\tz=%s\n' % (contig, start, end, identity, identity))
            i += 1
            hit_list.append(re.sub("\|", "", l[9+shift]))
            query_list.append(re.sub("\|", "", l[10+shift]))
    return (hit_list, query_list)

def get_gaps_from_start_stop_lists(contig2start_stop_lists, contigs_add, min_gap_size=1000):
    import pandas as pd
    import numpy as np


    gap_data = {}
    for contig in contig2start_stop_lists:
        print (contig)
        #print contig2start_stop_lists[contig]["start"]
        data = pd.DataFrame({'start': contig2start_stop_lists[contig]["start"],
                'stop': contig2start_stop_lists[contig]["stop"] })
        data.start = data.start.astype(np.int64)
        data.stop = data.stop.astype(np.int64)
        #print 'add', data['start'][0]+1
        data_sort = data.sort_values(by=["start"])
        print (data_sort)
        index_start = 0
        comparison_index = 1
        stop = False

        while index_start < len(data_sort['start']):
            # deal with unaligned begining of contig
            if index_start == 0 and int(data_sort['start'][0]) != int(contigs_add[contig][0]+1):
                print ('gap at the begining of the contig!', [int(contigs_add[contig][0]), int(data_sort['start'][0])-1])
                if (int(data_sort['start'][0]) -1) - int(contigs_add[contig][0]) > min_gap_size:
                    gap_data[contig] = [[int(contigs_add[contig][0]), int(data_sort['start'][0])-1]]
            # deal with unaligned end of contigs

            '''
            if index_start == (len(data_sort['start'])-2) and int(data_sort['stop'][index_start+1]) != int(contigs_add[contig][1]):
                print 'gap at the end of a contig!'
                if int(contigs_add[contig][1]) - int(data_sort['stop'][0]) > min_gap_size:
                    if contig not in gap_data:
                        gap_data[contig] = [[data_sort['stop'][index_start], int(contigs_add[contig][1])]]
                    else:
                        gap_data[contig].append([data_sort['stop'][index_start], int(contigs_add[contig][1])])
            '''
            #print 'index', index_start
            #print 'compt', comparison_index
            if index_start+comparison_index >= len(data_sort['start']):
                break
            #print 'index', index_start
            #print 'comp index', comparison_index
            # cas si fragment suivant est compris a l'interieur du 1er
            # comparer avec le fragment suivant
            if int(data_sort['stop'][index_start+comparison_index])<= int(data_sort['stop'][index_start]):
                #print 'compris interne', index_start, comparison_index
                #print '1', data_sort['start'][index_start], data_sort['stop'][index_start]
                #print '2', data_sort['start'][index_start+comparison_index], data_sort['stop'][index_start+comparison_index]

                comparison_index+=1
                continue
            # si le start du 2eme se trouve dans le premier et est plus long
            # skipper le 1er fragment
            elif int(data_sort['start'][index_start+comparison_index])<=int(data_sort['stop'][index_start]) and  \
                int(data_sort['stop'][index_start+comparison_index]) > int(data_sort['stop'][index_start]):
                index_start+=comparison_index
                comparison_index=1
                continue
            elif int(data_sort['start'][index_start+comparison_index])-int(data_sort['stop'][index_start]) > min_gap_size:
                #print '############# Potential gap!!! ', contig, data_sort['stop'][index_start], "-" , data_sort['start'][index_start+comparison_index]
                if contig not in gap_data:
                    gap_data[contig] = [[data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]]]
                else:
                    gap_data[contig].append([data_sort['stop'][index_start], data_sort['start'][index_start+comparison_index]])


                #print index_start, comparison_index
                #print '1', data_sort['start'][index_start], data_sort['stop'][index_start]
                #print '2', data_sort['start'][index_start+comparison_index], data_sort['stop'][index_start+comparison_index]
                #print data_sort

                index_start+=comparison_index
                comparison_index=1
            else:
                #print 'add'
                index_start+=comparison_index
                comparison_index=1
        if max(data_sort['stop']) < contigs_add[contig][1] and (contigs_add[contig][1]- max(data_sort['stop'])) > min_gap_size:
            print ('end gap???--------------------------')
            if contig not in gap_data:
                gap_data[contig] = [[max(data_sort['stop']), contigs_add[contig][1]]]
            else:
                gap_data[contig].append([max(data_sort['stop']), contigs_add[contig][1]])
        return gap_data





def get_circos_gap_file_from_gap_start_stop(contig2start_stop_lists,
                                            contigs_add,
                                            out_labels='circos_gaps_labels.txt',
                                            out_highlight="circos_gaps_highlight.txt",
                                            min_gap_size=1000,
                                            gap_merge_distance=0):

    import numpy as np
    import pandas as pd
    h = open(out_highlight, 'w')
    g = open(out_labels, 'w')

    gap_data = get_gaps_from_start_stop_lists(contig2start_stop_lists, contigs_add, min_gap_size=min_gap_size)

    for contig in gap_data:
        print (gap_data[contig])
        if len(gap_data[contig])>1:
            start_list = []
            stop_list = []
            for gap in gap_data[contig]:
                start_list.append(gap[0])
                stop_list.append(gap[1])
            data = pd.DataFrame({'start': start_list,
            'stop': stop_list })
            data.start = data.start.astype(np.int64)
            data.stop = data.stop.astype(np.int64)
            #print 'add', data['start'][0]+1
            data_sort = data.sort(columns=["start"])
            print (data_sort)
            data_sort.to_csv("%s_alignments.csv" % contig)
            start = data_sort['start'][0]
            stop = data_sort['stop'][0]

            for i in range(0, len(data_sort['start'])-1):
                print (i, len(data_sort['start']))
                # if less than 10kb between 2 gaps, merge
                if data_sort['start'][i+1]-data_sort['stop'][i] < gap_merge_distance:
                    stop = data_sort['stop'][i+1]
                    # if the last gap considered, write it
                    if (i+2) == len(data_sort['start']):
                        # gap highlight
                        h.write("%s\t%s\t%s\n" % (contig,start, stop))
                        # gap lables (start and stop)
                        g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                        g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                # more than 10kb between the 2 gaps
                else:
                    if i+2 == len(data_sort['start']):
                        # last gaps compared. as not merged, write both
                        h.write("%s\t%s\t%s\n" % (contig,start, stop))
                        # write labels
                        g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                        g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                        h.write("%s\t%s\t%s\n" % (contig,data_sort['start'][i+1], data_sort['stop'][i+1]))
                        # write labels
                        g.write("%s\t%s\t%s\t%s\n"% (contig,data_sort['start'][i+1], data_sort['start'][i+1]+5, data_sort['start'][i+1]))
                        g.write("%s\t%s\t%s\t%s\n"% (contig,data_sort['stop'][i+1], data_sort['stop'][i+1]+5, data_sort['stop'][i+1]))
                    else:
                        # not merged, write gap and reinitialize start and stop of the next gap
                        h.write("%s\t%s\t%s\n" % (contig,start, stop))
                        # write labels
                        g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
                        g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))
                        start = data_sort['start'][i+1]
                        stop = data_sort['stop'][i+1]

        else:
            start = gap_data[contig][0][0]
            stop = gap_data[contig][0][1]
            h.write("%s\t%s\t%s\n" % (contig, start, stop))
            # write labels
            g.write("%s\t%s\t%s\t%s\n"% (contig,start, start+5, start))
            g.write("%s\t%s\t%s\t%s\n"% (contig,stop, stop+5, stop))

    h.close()
    g.close()



