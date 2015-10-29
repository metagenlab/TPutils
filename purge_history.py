#!/usr/bin/env python
import sys
import re

m = re.compile("#[0-9]{10}")

p = open("purged_bash_history.txt", "w")

with open(".trest_bash_history", "r") as f:
    new_lines = {}
    all_lines = [i for i in f]
    new_lines_ordered = []
    for i in range(0,len(all_lines)):
        #print all_lines[i]
        #print all_lines[i+1]
        if m.match(all_lines[i]):
            continue
        if i%10000 == 0:
            print i, len(new_lines_ordered)
        if all_lines[i] not in new_lines:
            new_lines[all_lines[i]] = all_lines[i-1]
            new_lines_ordered.append([all_lines[i-1], all_lines[i]])
    #print len(new_lines_ordered)
        
    for i in new_lines_ordered:
        p.write("%s%s" % (i[0], i[1]))
        
