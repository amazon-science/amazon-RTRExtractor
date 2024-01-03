import itertools
import sys
import triangle_counters
import four_vertex_counters
import numpy as np


def all_five_vertex(G,fname='',gname='',to_list=False):
    five_counts = four_path_based(G,'','',False)
    prong_count = prong(G,'','',False)
    wedge_stuff = wedge_collision(G,'','',False)
    stars = fourstar(G,'','',False)

    all_counts = []
    all_counts.append(stars)
    all_counts.append(prong_count)
    all_counts.append(five_counts[4])
    all_counts.append(five_counts[7])
    all_counts.append(five_counts[6])
    all_counts.append(five_counts[5])
    all_counts.append(five_counts[0])
    all_counts.append(five_counts[11])
    all_counts.append(five_counts[8])
    all_counts.append(five_counts[9])
    all_counts.append(five_counts[10])
    all_counts.append(five_counts[1])
    all_counts.append(wedge_stuff[0])
    all_counts.append(wedge_stuff[1])
    all_counts.append(five_counts[2])
    all_counts.append(five_counts[12])
    all_counts.append(wedge_stuff[3])
    all_counts.append(wedge_stuff[2])
    all_counts.append(five_counts[3])
    all_counts.append(wedge_stuff[4])
    all_counts.append(wedge_stuff[5])

    return all_counts

def all_five_disconnected(G,fname='',gname=''):
    all_counts = []
    size_info = G.Size()
    n = size_info[0]
    m = size_info[1]/2
    w = size_info[2]
    
    order = G.DegenOrdering()
    DG = G.Orient(order)
    
    tri_info = triangle_counters.triangle_info(DG)
    t = sum(tri_info[0].values())/3
    
    four_counts = four_vertex_counters.four_vertex_count(G,'','',False)
    
    all_counts.append((n*(n-1)*(n-2)*(n-3)*(n-4))/120)
    all_counts.append(m*((n-2)*(n-3)*(n-4))/6)
    all_counts.append(((m*(m-1))/2 - w)*(n-4))
    all_counts.append(w*((n-3)*(n-4))/2)
    all_counts.append(t*((n-3)*(n-4))/2)
    all_counts.append(four_counts[0]*(n-4));   
    all_counts.append(four_counts[1]*(n-4));   
    all_counts.append(four_counts[2]*(n-4));   
    all_counts.append(four_counts[3]*(n-4));   
    all_counts.append(four_counts[4]*(n-4));   
    all_counts.append(four_counts[5]*(n-4));   
    all_counts.append(w*(m-2) - 3*t - 3*four_counts[0] - 2*four_counts[1])
    all_counts.append(t*(m-3) - four_counts[2])

    return all_counts

def fourstar(G,fname='',gname='',to_list=False):
    stars = 0
    for node1 in G.vertices:
        deg = len(G.adj_list[node1])
        stars = stars + (deg*(deg-1)*(deg-2)*(deg-3))/24
    return stars

def wedge_collision(G,fname='',gname='',to_list=False):
    three_wed = 0
    diamond_wed = 0
    three_tri = 0
    wheel = 0
    five_clique = 0
    almost_clique = 0

    wedges = dict()
    for node1 in G.vertices:
        for (nbr1, nbr2) in itertools.combinations(G.adj_list[node1],2):
            if nbr2 > nbr1:
                swap = nbr1
                nbr1 = nbr2
                nbr2 = swap

#             print 'outer', node1, nbr1, nbr2
            if (nbr1,nbr2) not in wedges:
                wedges[(nbr1,nbr2)] = set({node1})
            else:
                wedges[(nbr1,nbr2)].add(node1)
#                 print 'adding', nbr1, nbr2, wedges[(nbr1,nbr2)]

    for pair in wedges.keys():
        num = len(wedges[pair])
        three_wed = three_wed + num*(num-1)*(num-2)/6

        if to_list:
            print(pair, wedges[pair])

        if pair[0] in G.adj_list[pair[1]]:
            connected = True
        else:
            connected = False

        for (nbr1, nbr2, nbr3) in itertools.combinations(wedges[pair],3):
            num_edges = 0
            if to_list:
                print(pair, nbr1, nbr2, nbr3)
            if nbr1 in G.adj_list[nbr2]:
                num_edges += 1
            if nbr1 in G.adj_list[nbr3]:
                num_edges += 1
            if nbr2 in G.adj_list[nbr3]:
                num_edges += 1



            diamond_wed += num_edges
            wheel = wheel+ (num_edges*(num_edges-1))/2
            almost_clique += (num_edges*(num_edges-1)*(num_edges-2))/6
            
            if connected:
                three_tri += 1
                almost_clique += (num_edges*(num_edges-1))/2
                five_clique += (num_edges*(num_edges-1)*(num_edges-2))/6
#                 if num_edges == 3:
#                     dummy = [int(pair[0]), int(pair[1]), int(nbr1), int(nbr2), int(nbr3)]
#                     dummy = sorted(dummy)
#                     print "found from",pair[0],pair[1]
#                     print "clique", dummy[0], dummy[1], dummy[2], dummy[3], dummy[4]

#             print connected, num_edges

           
    wheel = wheel/2
    five_clique = five_clique/10
    almost_clique = almost_clique/4  

    return (three_wed, three_tri, wheel, diamond_wed, almost_clique, five_clique)

def prong(G,fname='',gname='',to_list=False):
    prong = 0

    for node1 in G.vertices:
        for node2 in G.adj_list[node1]:
            for (nbr1, nbr2) in itertools.combinations(G.adj_list[node2],2):
                if nbr1 == node1 or nbr2 == node1:
                    continue
                for node3 in G.adj_list[node1]:
                    if node3 == node2 or node3 == nbr1 or node3 == nbr2:
                        continue
                    prong = prong+1
    return prong
  


def four_path_based(G,fname='',gname='',to_list=False):
    
    tailed_cycle = 0
    hatted_cycle = 0

    tailed_clique = 0
    hatted_clique = 0
    four_paths = 0
    double_tailed_tri = 0
    long_tailed_tri = 0
    fork_tailed_tri = 0
    hourglass = 0
    cobra = 0
    stingray = 0
    fivecycle = 0
    triangle_strip = 0

    for node1 in G.vertices:
        for node2 in G.adj_list[node1]:
            for node3 in G.adj_list[node2]:
                if node3 == node1:
                    continue
                for node4 in G.adj_list[node3]:
                    if node4 == node2 or node4 == node1:
                        continue;
                    for node5 in G.adj_list[node4]:
                        if node5 == node3 or node5 == node1:
                            continue
                        if node5 == node2:
                            for node6 in G.adj_list[node2]:
                                if node6 == node1 or node6 == node3 or node6 == node4:
                                    continue
                                fork_tailed_tri = fork_tailed_tri+1
                            continue
                        four_paths = four_paths+1
                        
                        if G.isEdge(node2,node4):
                            double_tailed_tri = double_tailed_tri+1

                        if G.isEdge(node3,node5):
                            long_tailed_tri = long_tailed_tri+1

                        if G.isEdge(node1,node3) and G.isEdge(node5,node3):
                            hourglass = hourglass + 1

                        if G.isEdge(node2,node5) and G.isEdge(node3,node5):
                            cobra = cobra + 1

                        if G.isEdge(node2,node5) and G.isEdge(node2,node4):
                            stingray = stingray+1
                        
                        if G.isEdge(node5,node1):
                            fivecycle = fivecycle+1
                        
                        if G.isEdge(node1,node5) and G.isEdge(node2,node5) and G.isEdge(node2,node4):
                            triangle_strip = triangle_strip+1
                        
                        if G.isEdge(node2,node5):
                            tailed_cycle = tailed_cycle+1
                            if to_list:
                                print("Tailed cycle:",node1,node2,node3,node4,node5)
                            
                            if G.isEdge(node3,node5) and G.isEdge(node2,node4):
                                tailed_clique = tailed_clique+1
                                if to_list:
                                    print("Tailed clique:",node1,node2,node3,node4,node5)

                            if G.isEdge(node1,node3):
                                hatted_cycle = hatted_cycle+1
                                if to_list:
                                    print("Hatted cycle:",node1,node2,node3,node4,node5)

                                if G.isEdge(node3,node5) and G.isEdge(node2,node4):
                                    hatted_clique = hatted_clique+1
                                    if to_list:
                                        print("Hatted clique:",node1,node2,node3,node4,node5)

    return (tailed_cycle/2,hatted_cycle/2,tailed_clique/6,hatted_clique/4,four_paths/2,double_tailed_tri/2,long_tailed_tri/2,fork_tailed_tri/4,hourglass/8,cobra/2,stingray/2,fivecycle/10,triangle_strip/2)


