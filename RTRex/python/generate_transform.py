import os
import graph_tools
import triangle_counters
import four_vertex_counters
import five_vertex_counters
import sys
import numpy as np

transform = np.zeros(shape=(11,11))

for i in range(1,12):
    fname = sys.argv[1]+"pattern"+str(i)+".edges"
    G = graph_tools.graph()
    G.Read_edges(fname)
    all_counts = four_vertex_counters.four_vertex_count(G,'','',False)

    for j in range(0,11):
        transform[j][i-1] = all_counts[j]


final = ''
for i in range(0,11):
    row = ''
    for j in range(0,11):
        row = row + str(int(transform[i][j])) + ' '
    final = final + row + ';\n'

print(final)

# inverse = np.linalg.inv(transform)


