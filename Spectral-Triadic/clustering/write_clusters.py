## A wrapper code that writes the clusters according to the original name.
## USAGE:
##        python3 write_clusters.py <PATH TO MAP FILE> <PATH TO DECOMPOSITION WITH NEW INTEGER LABELS>
##
## The assumption is that:
##       <PATH TO MAP FILE>: This is a map file where each line has the old (original) name and the new (integer) name, sorted by the new name. The second entry in each line is simply the line number, so it is redundant.
##       <PATH TO DECOMPOSITION>: This file has, in each line, a different cluster according to the new integer name.
##
## The procedure is quite simple. It simply stores the old names in a list actual_names, where the ith entry is the old map that was mapped to i. 
## On reading a cluster, it just maps the integer index (i) to the old name (ith entry in the list actual_names)
##
## C. Seshadhri, Aug 2023

import sys
import os

if(len(sys.argv) < 3): # incorrect arguments, so print usage statement and exit
    print('USAGE: python3 write_clusters.py <PATH TO MAP FILE> <PATH TO DECOMPOSITION>')
    sys.exit(-1)

mapping = sys.argv[1] # the first argument is the mapping file
decomposition = sys.argv[2] # the second argument is the decomposition, with the new integer labels

actual_names = list() # this is the list of actual names

with open(mapping,'r') as f_map:
    for line in f_map.readlines(): # reading each line of the map file
        tokens = line.split() # tokenize each line by whitespace
        actual_names.append(tokens[0]) # the first token is the actual name. Note that we simply put the first token of ith line in the ith position of actual_names.
                                       # We are assuming that the map file is sorted by the old (integer) label

with open('real-'+decomposition,'w') as f_realclus: # open another file to write the real clusters in
    with open(decomposition,'r') as f_decomp: # open the input clusters file
        for line in f_decomp.readlines(): # each line has a different cluster
            tokens = line.strip().split() # each token is a different vertex of that cluster
            for vert in tokens: # loop over all the tokens, which are the clusters
                if (int(vert) >= len(actual_names)): # some error here. If the vertex ID is outside the list of actual names, print it out and continue
                    print("Error, out of bounds: "+vert+"\n") # printing the error message
                else:
                    f_realclus.write(actual_names[int(vert)]+" ") # we can pick out the original/actual name from the list, so print it with a space after
            f_realclus.write("\n") # put in a newline to end the cluster




