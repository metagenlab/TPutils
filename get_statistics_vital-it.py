#!/usr/bin/env python

def parse_output_files(output_file):

    '''
    CPU time :               604.10 sec.
    Max Memory :             554.28 MB
    Average Memory :         374.82 MB
    Total Requested Memory : 6000.00 MB
    Delta Memory :           5445.72 MB
    (Delta: the difference between total requested memory and actual max usage.)
    Max Swap :               4117 MB

    :param output_file:
    :return:
    '''

    import re

    cpu_time = re.compile(".*CPU time :(.*)")
    max_memory = re.compile(".*Max Memory :(.*)")

    with open(output_file, 'r') as f:
        for line in f:
                #print(line)
                if re.match(cpu_time, line):
                    n_secondes = line.split()[-2]
                if re.match(max_memory, line):
                    n_mega = line.split()[-2]
    return n_secondes, n_mega

def plot_scores_distribution(file_list, outname):


    import pairwiseid_plots
    import numpy

    max_mem_list = []
    for one_file in file_list:
        max_mem_list.append(float(parse_output_files(one_file)[1]))

    #pairwiseid_plots.density_plot([max_mem_list], label_list=['mem'], output_path=outname)
    print ('%s\t%s\t%s\t%s' % (outname.split('.')[0],
                               min(max_mem_list),
                               max(max_mem_list),
                               numpy.median(max_mem_list)))

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--output_files', type=str, help="vital-it output files", default=False, nargs='+')
    parser.add_argument("-o", '--output_name', type=str, help="output name", default="out.svg")

    args = parser.parse_args()

    plot_scores_distribution(args.output_files, args.output_name)




