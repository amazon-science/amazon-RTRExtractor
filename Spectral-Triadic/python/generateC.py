import os
import sys
import graph_tools



fname = sys.argv[1]
mat = sys.argv[2]
f_input = open(fname,'r')

list_lines = f_input.readlines()

row = 0
for each_line in list_lines:
    tokens = (each_line.strip()).split('&')
    col = 0
    for entry in tokens:
        entry = entry.strip()
        if len(entry) == 0:
            entry = 0;
        to_print = mat+'['+str(row)+']['+str(col)+'] = '+str(entry)+';'
        print to_print
        col = col+1
    row = row + 1
    f_input.close()
