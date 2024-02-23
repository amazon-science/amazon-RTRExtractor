import graph_tools
import triangle_counters
import four_vertex_counters
import five_vertex_counters
import sys

# Usage:
# python tester_five.py <GRAPH NAME>

name = sys.argv[1]

G = graph_tools.graph()
G.Read_edges(name)
four_counts = four_vertex_counters.four_vertex_count(G,'','',False)
all_counts = five_vertex_counters.all_five_vertex(G,'','',False)
all_disc_counts = five_vertex_counters.all_five_disconnected(G,'','')


print "------------------------"
print long(four_counts[0]), "3-stars"
print long(four_counts[1]), "3-paths"
print long(four_counts[2]), "tailed-triangles"
print long(four_counts[3]), "4-cycles"
print long(four_counts[4]), "diamonds"
print long(four_counts[5]), "4-cliques"

print "---------------"
print long(all_disc_counts[0]), "Ind set"
print long(all_disc_counts[1]), "Only edge"
print long(all_disc_counts[2]), "Matching"
print long(all_disc_counts[3]), "Only wedge"
print long(all_disc_counts[4]), "Only triangle"
print long(all_disc_counts[5]), "Only 3-star"
print long(all_disc_counts[6]), "Only 3-path"
print long(all_disc_counts[7]), "Only Tailed tri"
print long(all_disc_counts[8]), "Only 4-cycle"
print long(all_disc_counts[9]), "Only Diamond"
print long(all_disc_counts[10]), "Only 4-clique"
print long(all_disc_counts[11]), "Wedge+edge"
print long(all_disc_counts[12]), "Triangle+edge"



print "----------------"
print all_counts[0], "stars"
print all_counts[1], "prongs"
print all_counts[2], "four paths"
print all_counts[3], "fork-tailed triangle"
print all_counts[4], "long-tailed triangle"
print all_counts[5], "double-tailed triangle"
print all_counts[6], "tailed 4-cycles"
print all_counts[7], "five cycles"
print all_counts[8], "hourglasses"
print all_counts[9], "cobras"
print all_counts[10], "stingrays"
print all_counts[11], "hatted 4-cycles"
print all_counts[12], "three wedge collisions"
print all_counts[13], "three triangles collisions"
print all_counts[14], "tailed 4-cliques"
print all_counts[15], "triangle strips"
print all_counts[16], "diamond-wedge collisions"
print all_counts[17], "wheels"
print all_counts[18], "hatted 4-cliques"
print all_counts[19], "almost cliques"
print all_counts[20], "five cliques"

f = open("nums_test.txt","w")

f.write(str(int(four_counts[0]))+"\n")
f.write(str(int(four_counts[1]))+"\n")
f.write(str(int(four_counts[2]))+"\n")
f.write(str(int(four_counts[3]))+"\n")
f.write(str(int(four_counts[4]))+"\n")
f.write(str(int(four_counts[5]))+"\n")
f.write(str(int(all_counts[0]))+"\n")
f.write(str(int(all_counts[1]))+"\n")
f.write(str(int(all_counts[2]))+"\n")
f.write(str(int(all_counts[3]))+"\n")
f.write(str(int(all_counts[4]))+"\n")
f.write(str(int(all_counts[5]))+"\n")
f.write(str(int(all_counts[6]))+"\n")
f.write(str(int(all_counts[7]))+"\n")
f.write(str(int(all_counts[8]))+"\n")
f.write(str(int(all_counts[9]))+"\n")
f.write(str(int(all_counts[10]))+"\n")
f.write(str(int(all_counts[11]))+"\n")
f.write(str(int(all_counts[12]))+"\n")
f.write(str(int(all_counts[13]))+"\n")
f.write(str(int(all_counts[14]))+"\n")
f.write(str(int(all_counts[15]))+"\n")
f.write(str(int(all_counts[16]))+"\n")
f.write(str(int(all_counts[17]))+"\n")
f.write(str(int(all_counts[18]))+"\n")
f.write(str(int(all_counts[19]))+"\n")
f.write(str(int(all_counts[20]))+"\n")
