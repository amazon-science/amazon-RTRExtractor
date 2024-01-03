import os
import graph_tools
import triangle_counters
import four_vertex_counters
import five_vertex_counters
import sys
import numpy as np

transform_5 = np.zeros(shape=(34,34))

transform_4 = np.matrix('1 1 1 1 1 1 1 1 1 1 1 ;'
                        '0 1 2 2 3 3 3 4 4 5 6 ;'
                        '0 0 1 0 0 0 1 1 2 2 3 ;'
                        '0 0 0 1 3 3 2 5 4 8 12 ;'
                        '0 0 0 0 1 0 0 1 0 2 4 ;'
                        '0 0 0 0 0 1 0 1 0 2 4 ;'
                        '0 0 0 0 0 0 1 2 4 6 12 ;'
                        '0 0 0 0 0 0 0 1 0 4 12 ;'
                        '0 0 0 0 0 0 0 0 1 1 3 ;'
                        '0 0 0 0 0 0 0 0 0 1 6 ;'
                        '0 0 0 0 0 0 0 0 0 0 1')

for j in range(0,11):
# first 11 X 11 block is same as transform for 4
    for i in range(0,11):
        transform_5[i,j] = transform_4[i,j]
# the jth pattern doesn't contain any other subgraphs
    for i in range(11,34):
        transform_5[i,j] = 0

# dealing with the 12th pattern, wedge+edge

transform_5[0,11] = 1
transform_5[1,11] = 3
transform_5[2,11] = 2
transform_5[3,11] = 1
transform_5[4,11] = 0

for i in range(5,11):
    transform_5[i,11] = 0

transform_5[11,11] = 1
transform_5[12,11] = 0

for i in range(13,34):
    transform_5[i,11] = 0

# dealing with the 13th pattern, triangle+edge

transform_5[0,12] = 1
transform_5[1,12] = 4
transform_5[2,12] = 3
transform_5[3,12] = 3
transform_5[4,12] = 1

for i in range(5,12):
    transform_5[i,12] = 0

transform_5[11,12] = 3
transform_5[12,12] = 1

for i in range(13,34):
    transform_5[i,12] = 0



for j in range(13,34):
    num = j-12
    fname = "../graphs/pattern"+str(num)+".edges"
    G = graph_tools.graph()
    G.Read_edges(fname)
    all_disc_counts = five_vertex_counters.all_five_disconnected(G,'','')
    all_con_counts = five_vertex_counters.all_five_vertex(G,'','',False)

    for i in range(0,13):
        transform_5[i,j] = all_disc_counts[i]
    
    for i in range(13,34):
        transform_5[i,j] = all_con_counts[i-13]

for i in range(0,34):
    row = '\''
    for j in range(0,34):
        row = row + str(int(transform_5[i,j])) + ' '
    row = row + ';\''
    print(row)

trans_inv = np.linalg.inv(transform_5)
print('\n\n')

for i in range(0,34):
    row = '\''
    for j in range(0,34):
        row = row + str(trans_inv[i,j]) + ' '
    row = row + ';\''
    print(row)


