import os
import graph_tools
import triangle_counters
import four_vertex_counters
import five_vertex_counters
import sys
import numpy as np



f_mat = open('matrices_latex.txt','w')
f_c = open('matrix_c.txt','w')

transform = np.zeros(shape=(21,21))

for i in range(1,22):
    fname = sys.argv[1]+"pattern"+str(i)+".edges"
    G = graph_tools.graph()
    G.Read_edges(fname)
    all_counts = five_vertex_counters.all_five_vertex(G,'','',False)

    for j in range(0,21):
        transform[j][i-1] = all_counts[j]

print(transform)
inverse = np.linalg.inv(transform)


for i in range(0,21):
    latex_code = ''
    for j in range(0,21):
       if int(transform[i][j]) == 0:
           entry = ' '
       else:
           entry = str(int(transform[i][j]))
       latex_code = latex_code + entry  + ' & '
       f_c.write(sys.argv[2]+'['+str(i)+']['+str(j)+'] = '+entry+';\n')
    latex_code = latex_code[:-2] + '\\\\ \n'
    f_mat.write(latex_code)

f_mat.write('\n\n\n\n\n')
f_c.write('\n\n\n\n\n')

for i in range(0,21):
    latex_code = ''
    for j in range(0,21):
       if int(inverse[i][j]) == 0:
            entry = ' '
       else:
            entry = str(int(inverse[i][j]))
       latex_code = latex_code + entry + ' & '
       f_c.write(sys.argv[3]+'['+str(i)+']['+str(j)+'] = '+entry+';\n')
    latex_code = latex_code[:-2] + '\\\\ \n'
    f_mat.write(latex_code)



f_mat.close()
f_c.close()


