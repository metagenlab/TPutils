#!/usr/bin/env python

import sys 
import re 

m = re.compile("(.*)_(S[0-9]+)_[A-Z][0-9]+_(R[1|2])_")

prefix_and_sample_id2new_file_prefix = {}
file_prefix2count = {}
prefix2files = {}
for row in sys.stdin:
    row = row.rstrip()
    s = re.search(m, row)
    if s is None:
        raise IOError("No match for this row:\n%s" % row)
    prefix = s.groups()[0]
    sample_id = s.groups()[1]
    R = s.groups()[2]     

    file_prefix = prefix.split("/")[-1]
    if R == 'R1':
        if file_prefix in file_prefix2count:
            file_prefix2count[file_prefix] +=1
            prefix_and_sample_id2new_file_prefix["%s_%s" % (prefix, sample_id)] = "%s-%s" % (file_prefix, file_prefix2count[file_prefix])
        else:
            file_prefix2count[file_prefix] = 0
            prefix_and_sample_id2new_file_prefix["%s_%s" % (prefix, sample_id)] = file_prefix
    
    if prefix not in prefix2files:
        prefix2files[prefix] = {}
    if sample_id not in prefix2files[prefix]:
        prefix2files[prefix][sample_id] = {}
    if R in prefix2files[prefix][sample_id]:
        raise IOError("Two files with identical prefixes?\n- %s\n- %s" % (row, prefix2files[prefix][sample_id][R]))
    prefix2files[prefix][sample_id][R] = row 

for prefix in  prefix2files:
    for sample_id in prefix2files[prefix]:

        if 'R1' in  prefix2files[prefix][sample_id] and 'R2' in prefix2files[prefix][sample_id]:
            sys.stdout.write("%s\t%s\t%s\t%s\t%s\n" % (prefix,
                                                       sample_id,
                                                       prefix2files[prefix][sample_id]["R1"],
                                                       prefix2files[prefix][sample_id]["R2"],
                                                       prefix_and_sample_id2new_file_prefix["%s_%s" % (prefix, sample_id)]))
                                            
        elif 'R1' in  prefix2files[prefix][sample_id] and not 'R2' in prefix2files[prefix][sample_id]:
            sys.stdout.write("%s\t%s\t%s\t-\t%s\t-\n" % (prefix,
                                                             sample_id,
                                                             prefix2files[prefix][sample_id]["R1"],
                                                             prefix_and_sample_id2new_file_prefix["%s_%s" % (prefix, sample_id)]))        
        elif 'R2' in  prefix2files[prefix][sample_id] and not 'R1' in prefix2files[prefix][sample_id]:
            raise IOError("Missing R1 fastq for %s\n%s" % (prefix, prefix2files[prefix]))
        else:
            raise IOError("Unexpected file naming for %s\n%s" % (prefix, prefix2files[prefix]))
    
