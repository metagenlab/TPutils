#!/usr/bin/env python

import sys 
import re 

m = re.compile("(.*)_[A-Z][0-9]_(R[1|2])_")

prefix2files = {}
for row in sys.stdin:
    row = row.rstrip()
    s = re.search(m, row)
    prefix = s.groups()[0]
    R = s.groups()[1]
    if prefix not in prefix2files:
        prefix2files[prefix] = {}
    if R in prefix2files[prefix]:
        raise IOError("Two files with identical prfixes?\n- %s\n- %s" % (row, prefix2files[prefix][R]))
    prefix2files[prefix][R] = row 

for prefix in  prefix2files:
    if 'R1' in  prefix2files[prefix] and 'R2' in prefix2files[prefix]:
        sys.stdout.write("%s\t%s\t%s\n" % (prefix,
                              prefix2files[prefix]["R1"],
                              prefix2files[prefix]["R2"]))
    elif 'R1' in  prefix2files[prefix] and not 'R2' in prefix2files[prefix]:
        sys.stdout.write("%s\t%s\t-\n" % (prefix,
                                          prefix2files[prefix]["R1"]))        
    elif 'R2' in  prefix2files[prefix] and not 'R1' in prefix2files[prefix]:
        raise IOError("Missing R1 fastq for %s\n%s" % (prefix, prefix2files[prefix]))
    else:
        raise IOError("Unexpected file naming for %s\n%s" % (prefix, prefix2files[prefix]))
    