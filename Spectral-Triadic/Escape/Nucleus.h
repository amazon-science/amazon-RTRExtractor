#ifndef ESCAPE_TRIADIC_H_
#define ESCAPE_TRIADIC_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include <set>
#include <inttypes.h>
using namespace std;

using namespace Escape;

// Structure that stores triangle information of graph

struct TriangleInfo
{
    Count total; // total number of triangles
    Count *perVertex;  // array storing number of triangles for each vertex
    Count *perEdge;   // arry storing number of triangles for each edge

};

struct CEdges{
    EdgeInfo *edges;
};

//use with a Graph or CGraph as argument
template <class T>
TriangleInfo newTriangleInfo(const T* graph)
{
  TriangleInfo ret;
  ret.total = 0;
  ret.perVertex = new Count[graph->nVertices];
  ret.perEdge = new Count[graph->nEdges];
  return ret;
}


void delTriangleInfo(TriangleInfo& info)
{
  delete[] info.perVertex;
  delete[] info.perEdge;
}


// Structure for storing all the triangles of graph
// Interpretation of the structure requires knowledge of the exact CGraph cg used to construct it
// Every index pos in the nbors list of cg corresponds to some edge (i,j). Note that each undirected edge appears as (i,j) and (j,i)
// Then trioffsets[pos] contains the starting index (in *triangles) of all the triangles that (i,j) is incident to.
// Basically, the portion of the array trioffsets[pos] to trioffsets[pos+1] in triangles is a list of vertices v1, v2,... such 
//that each vi forms a triangles with (i,j)
//
// We only store triangles for one copy of edge (i,j), where i < j in the degree ordering.
struct TriangleList
{
    VertexIdx total;
    VertexIdx *triangles;
    EdgeIdx *trioffsets;
};

struct Cluster
{
    EdgeIdx id;
    EdgeIdx c;
    VertexIdx n_vertices;
    EdgeIdx n_edges;
    EdgeIdx parent;
    std::set<VertexIdx> vertices;
};


// This wedge enumeration algorithm produces all triangles, and is more
// efficient than the previous version. We use binary search to find edges, and also pass
// the original graph to cut down on queries.
// Input: a pointer gout to a CGraph labeled according to degree
// Output: a TriangleInfo for g. The ordering of edges in perEdge (of TriangleInfo) is that same as g.

TriangleInfo betterWedgeEnumerator(CGraph *gout)
{
   TriangleInfo ret;   // output 
   ret.total = 0;      // initialize outout
   ret.perVertex = new EdgeIdx[gout->nVertices+1];
   ret.perEdge = new EdgeIdx[gout->nEdges+1]; 

   for (VertexIdx i=0; i < gout->nVertices; ++i)
       ret.perVertex[i] = 0;

   for (EdgeIdx j=0; j < gout->nEdges; ++j)
       ret.perEdge[j] = 0; 

   VertexIdx end1, end2;
   EdgeIdx loc;

   for (VertexIdx i=0; i < gout->nVertices; ++i) // loop over vertices
       for (EdgeIdx j = gout->offsets[i]; j < gout->offsets[i+1]; ++j)   // loop over neighbor of i
           for (EdgeIdx k = j+1; k < gout->offsets[i+1]; ++k)         // loop over another neighbor of i
           {
               end1 = gout->nbors[j];     // we are now looking at wedge (i, end1, end2), centered at i
               end2 = gout->nbors[k];

               // note that end1 < end2 because of the labeled ordering

               loc = gout->getEdgeBinary(end1,end2);
               if (loc != -1)        // (end1, end2) is present
               {
                   ret.total++;       // found a triangle! So update total.

                   ret.perVertex[i]++; // update all per vertex counts
                   ret.perVertex[end1]++;
                   ret.perVertex[end2]++;

                   ret.perEdge[j]++;  // update all per edge counts. Note that location used is same as position in g->nbors
                   ret.perEdge[k]++;
                   ret.perEdge[loc]++;
               }
           }

   return ret;
}


// This stores all triangles in a TriangleList structure, corresponding to CGraph g. 
// The number of triangles numtri is required for initialization

TriangleList storeAllTriangles(CGraph *g, EdgeIdx numtri)
{
    TriangleList ret;

    ret.total = 3*numtri;
    ret.triangles = new VertexIdx[3*numtri+1];
    ret.trioffsets = new EdgeIdx[g->nEdges+1];

    EdgeIdx posj = 0;

    EdgeIdx current = 0;
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        VertexIdx degi = g->offsets[i+1] - g->offsets[i];
        for (posj = g->offsets[i]; posj < g->offsets[i+1]; posj++)
        {
            VertexIdx j = g->nbors[posj];
//             if (i == g->nVertices-1)
//                 printf("j is %" PRId64 "\n",j);
            VertexIdx degj = g->offsets[j+1] - g->offsets[j];
            ret.trioffsets[posj] = current;
            if (degj < degi || (degj == degi && j <= i))
                continue;
            
            for (EdgeIdx ptr = g->offsets[i]; ptr < g->offsets[i+1]; ptr++)
            {
                VertexIdx nbr = g->nbors[ptr];
                if (g->getEdgeBinary(nbr,j) != -1)
                {
//                     if (i == g->nVertices-1)
//                         printf("%" PRId64 " at ptr %" PRId64 ": nbr is %" PRId64 "\n",current,ptr,nbr);
                    ret.triangles[current] = nbr;
                    current++;
                }
            }
        }
    }
    
    ret.trioffsets[posj] = current; // posj is the index after the last edge in g. we set this final offset to the current position

    return ret;
}


//Auxiliary function that given an edge creates a new cluster in the hierarchy, it will look through the neighbors of the edge, if it finds a non-clustered edge it will ad it as a child
//if it finds a cluster edge it will find the representative of that cluster and add it as a child

void findConnectedComponent(
    CGraph *g,
    EdgeIdx e, 
    EdgeIdx *neighbors, 
    EdgeIdx *neighbors_offset, 
    EdgeIdx *deletion_time, 
    EdgeIdx *parent, 
    EdgeIdx *compressed_parent, 
    EdgeIdx *children, 
    EdgeIdx *children_offset, 
    EdgeIdx &children_counter,
    bool *in_cluster)
{
    // start creating the stack add the neighbors of e
    EdgeIdx *stack = new EdgeIdx[g->nEdges+1]; // store list of edges to be explored
    EdgeIdx stack_index = -1; // largest index in stack

    for (EdgeIdx idxi = neighbors_offset[deletion_time[e]]; idxi < neighbors_offset[deletion_time[e]+1]; idxi++)
        {
            EdgeIdx neighbor = neighbors[idxi];
            stack[stack_index+1] = neighbor;
            stack_index++;
        }   

    in_cluster[e]=true;

    //while the stack is not empty
    while(stack_index >= 0){
        EdgeIdx edge = stack[stack_index];
        stack_index--;
        //check if the edge was already visited by e (its compressed parent is e)
        if(compressed_parent[edge]==e){
            continue;
        }
        //check if belongs already to a cluster
        if(in_cluster[edge]) {
            //need to find the cluster it belongs to
            EdgeIdx index = edge;
            while(compressed_parent[index] != index) {
                EdgeIdx new_index = compressed_parent[index];
                compressed_parent[index] = e;
                index = new_index;
            }
            //once in the top of the cluster, check if the cluster is already e
            if(index != e){
                //otherwise we add it as child of e
                parent[index] = e;
                //printf("set parent\n");
                compressed_parent[index] = e;
                children[children_counter] = index;
                children_counter++;
            }
            continue;
        }
        //else the edge doesn't belong to any cluster
        //add edge to the cluster
        parent[edge] = e;
        compressed_parent[edge] = e;
        children[children_counter] = edge;
        children_counter++;
        in_cluster[edge] = true;
        //add the neighbours to the stack
        for (EdgeIdx idxi = neighbors_offset[deletion_time[edge]]; idxi < neighbors_offset[deletion_time[edge]+1]; idxi++) //loop over all the children of the edge
        {
            EdgeIdx neighbor = neighbors[idxi];
            stack[stack_index+1] = neighbor;
            stack_index++;
        }    
    }
    free(stack);
}


//This function creates the final cluster rooted at root using the hierarchy
void completeCluster(
        Cluster *new_cluster, 
        EdgeIdx *root_to_cluster, 
        Cluster *clusters, 
        EdgeIdx root,
        EdgeIdx c, 
        EdgeIdx *kappa, 
        CGraph *g, 
        EdgeIdx *sources, 
        EdgeIdx *children, 
        EdgeIdx *children_offset,
        EdgeIdx *processing_time){

    // start creating the stack and add the root of the cluster
    EdgeIdx *stack = new EdgeIdx[g->nEdges+1]; // store list of edges to be explored
    stack[0] = root;
    EdgeIdx stack_index = 0; // largest index in stack

    //while the stack is not empty
    while(stack_index >= 0){
        EdgeIdx edge = stack[stack_index];
        stack_index--;
        //Check if the edge is from a higher level
        if (kappa[edge] > c) {
            //then we need to add this cluster to the current
            Cluster cluster_to_combine = clusters[root_to_cluster[edge]];
            new_cluster->n_edges += cluster_to_combine.n_edges;
            clusters[cluster_to_combine.id].parent = new_cluster->id;
            for (set<VertexIdx>::iterator itr = cluster_to_combine.vertices.begin(); itr != cluster_to_combine.vertices.end(); itr++)
            {
                if(new_cluster->vertices.count(*itr) == 0) {
                    new_cluster->vertices.insert(*itr);
                    new_cluster->n_vertices++;
                }
            }
        } else {
            //we just add the edge to the cluster and the children to the stack
            new_cluster->n_edges++;
            VertexIdx src = sources[edge];
            if(new_cluster->vertices.count(src) == 0) {
                new_cluster->vertices.insert(src);
                new_cluster->n_vertices++;
            }
            VertexIdx dest = g->nbors[edge];
            if(new_cluster->vertices.count(dest) == 0) {
                new_cluster->vertices.insert(dest);
                new_cluster->n_vertices++;
            }

            for (EdgeIdx idxi = children_offset[processing_time[edge]]; idxi < children_offset[processing_time[edge]+1]; idxi++) {
                stack[stack_index+1] = children[idxi];
                stack_index++;
            }
        }
    }
    free(stack);
}

//Main Function to generate the nuclei
//The output is a list of clusters of type Cluster
Cluster* nuclei(CGraph *g, TriangleList* tlist, bool noMerge, EdgeIdx &nClusters)
{
    EdgeIdx *triCount = new EdgeIdx[g->nEdges+1];  // store triangle count of each edge
    EdgeIdx *toDelete = new EdgeIdx[g->nEdges+1]; // store list of edges to be deleted

    EdgeIdx ind_toDelete = 0; // largest index in toDelete

    bool *deleted = new bool[g->nEdges+1];  // store flags for deleted edges
    VertexIdx *kappa = new VertexIdx[g->nEdges+1]; //store the kappa value (max c such that e belong to the c-truss) of each edge
    EdgeIdx *deletion_time = new EdgeIdx[g->nEdges+1]; //store the deletion time of every edge (that is the deletion ordering)

    EdgeIdx deletionCount = 0; //count the number of deleted edges so far

    //store the edge indexes of the children of every edge, the children_offset array uses the deletion_time[index] instead of the index directly
    EdgeIdx *neighbors = new EdgeIdx[tlist->total+1];
    EdgeIdx *neighbors_offset = new EdgeIdx[g->nEdges+1];
    
    neighbors_offset[0]=0; //the first edge starts at 0
    EdgeIdx neighbors_index = 0; //position of the next element to be added in the children array

    //store only the edges that are deleted for every value of c
    EdgeIdx *deleted_edges = new EdgeIdx[g->nEdges+1];
    EdgeIdx *deleted_offsets = new EdgeIdx[g->nEdges+1];
    EdgeIdx deleted_index = 0;
    deleted_offsets[0]=0;

    for (EdgeIdx i=0; i < g->nEdges+1; i++) // initialize all edges as not deleted and ancestors to -1
    {
        deleted[i] = false;
    }

    //Before starting we will precompute the source of every edge based on the id of the edge, this will help with other computations later
    EdgeIdx *sources = new EdgeIdx[g->nEdges+1];
    for (EdgeIdx i=0; i < g->nVertices; i++) {
        for (EdgeIdx j=g->offsets[i]; j < g->offsets[i+1]; j++) {
            sources[j]=i;
        }
    }

    EdgeIdx nEdges = g->nEdges;
    EdgeIdx c = 1; //set threshold to 1
    EdgeIdx count;
    EdgeIdx totalCount = 0;

    for (VertexIdx i=0; i < g->nVertices; i++)
            for (EdgeIdx posj=g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over all edges
            {
                VertexIdx j = g->nbors[posj]; // edge (i,j)
                if (i > j)   // only look at pairs ordered properly
                    continue;
                triCount[posj] = tlist->trioffsets[posj+1] - tlist->trioffsets[posj];  // store the number of triangles that edge (i,j) participates in. Information is exactly stored in tlist
            }
    printf("nEdges = %" PRId64 "\n",nEdges/2);
    while(totalCount < nEdges/2)
    {
        printf("c = %" PRId64 "\n",c-1);
        count = 0;
        for (VertexIdx i=0; i < g->nVertices; i++)
            for (EdgeIdx posj=g->offsets[i]; posj < g->offsets[i+1]; posj++) //looping over all edges
            {
                VertexIdx j = g->nbors[posj]; // edge (i,j)
                if (i > j || deleted[posj])   // only look at pairs ordered properly and not deleted
                    continue;
                if (triCount[posj] < c) // edge should be deleted
                {
                    EdgeIdx current;
                    current = posj;
                    toDelete[ind_toDelete] = current;   // storing current edge in toDelete
                    ind_toDelete++;  // update index
                    //printf("added edge (%ld,%ld) to deletion\n",i,j);
                }
            }

        while(ind_toDelete >= 1) // while toDelete is non-empty
        {
            EdgeIdx toRemove = toDelete[ind_toDelete-1]; // get edge to remove
            ind_toDelete--; // update last index in toDelete
            if(deleted[toRemove]){ //ignore if already deleted
                continue;
            }
            VertexIdx i = sources[toRemove];
            VertexIdx j = g->nbors[toRemove];
            //printf("deleting edge (%ld,%ld)\n",i,j);
            deleted[toRemove] = true;
            kappa[toRemove] = c-1;
            deletion_time[toRemove] = deletionCount;
            deletionCount++;

            deleted_edges[deleted_index] = toRemove;
            deleted_index++;

            count++;
            totalCount++;
        
    //         printf("%" PRId64 " %" PRId64 "\n",i,j); 
            for (EdgeIdx indk = tlist->trioffsets[toRemove]; indk < tlist->trioffsets[toRemove+1]; indk++) // looping over triangles that edge participates in
            {
                VertexIdx k = tlist->triangles[indk]; // (i,j,k) forms triangle
                //printf("k = %ld \n",k);
                EdgeIdx locik, locjk;
                if (i < k) {// get position of ordered edge (i,k)
                    locik = g->getEdgeBinary(i,k);
                } else {
                    locik = g->getEdgeBinary(k,i);
                }
                if (locik == -1) // something is wrong. (i,k) has to be edge in g
                {
                    printf("Error: edge i,k is not present in graph, but triangle (i,j,k) is stored in tlist\n");
                    exit(1);
                }
                if (!deleted[locik]){
                    neighbors[neighbors_index]=locik;
                    neighbors_index++;
                }

                if (j < k) {// get position of ordered edge (j,k)
                    locjk = g->getEdgeBinary(j,k);
                    //printf("added (%" PRId64 ",%" PRId64 ") as an ancestor of (%" PRId64 ",%" PRId64 ")\n",i,j,j,k);
                } else {
                    locjk = g->getEdgeBinary(k,j);
                    //printf("added (%" PRId64 ",%" PRId64 ") as an ancestor of (%" PRId64 ",%" PRId64 ")\n",i,j,k,j);
                }
                if (locjk == -1) // something is wrong. (j,k) must be edge in g
                {
                    printf("Error: edge j,k not present in graph, but triangle (i,j,k) is stored in tlist\n");
                    exit(1);
                }
                if (!deleted[locjk]){
                    neighbors[neighbors_index]=locjk;
                    neighbors_index++;
                }

                if (deleted[locik])  // if (i,k) has been already deleted, then this triangle has been deleted, so continue
                    continue;

                if (deleted[locjk]) // if (j,k) has been deleted, triangle (i,j,k) has been deleted, so continue
                    continue;


                // we delete triangle (i,j,k), so decrement triangle counts appropriately
                triCount[locik]--; // decrement count for (i,k)
                if (triCount[locik] < c) // (i,k) now participates in less than c triangles
                {
                    toDelete[ind_toDelete] = locik;
                    ind_toDelete++;
                    //printf("Added edge (%ld,%ld) to deletion\n",i,k);
                }
                triCount[locjk]--; // decrement count for (j,k)
                if (triCount[locjk] < c) // (j,k) now participates in less than c triangles
                {
                    toDelete[ind_toDelete] = locjk;
                    ind_toDelete++;
                    //printf("Added edge (%ld,%ld) to deletion\n",j,k);
                }
            }
            //update the offset array
            neighbors_offset[deletionCount]=neighbors_index;
        }
        deleted_offsets[c]=deleted_index;
        c++;
        printf("neighbors_index: %" PRId64 "\n", neighbors_index);
        printf("deleted: %" PRId64 "\n", count);
        printf("total deleted: %" PRId64 "\n", totalCount);
    }
    printf("Construct all the clusters\n");
    //We construct the clusters here, we start from the latest time and we list all the cluster at that level, the clusters are form by edges
    //but we will get the list of nodes that participates in such edges to make the output smaller.

    EdgeIdx *parent = new EdgeIdx[g->nEdges+1]; 
    EdgeIdx *compressed_parent = new EdgeIdx[g->nEdges+1]; 
    bool *in_cluster = new bool[g->nEdges+1];
    EdgeIdx *processing_time = new EdgeIdx[g->nEdges+1]; 
    EdgeIdx processing_time_counter = 0;
    EdgeIdx *children = new EdgeIdx[g->nEdges+1];
    EdgeIdx *children_offset = new EdgeIdx[g-> nEdges+1]; //The offset use processing_time[edge]
    EdgeIdx children_counter =  0;
    for (EdgeIdx i=0; i < g->nEdges+1; i++) // initialize all nodes parents to themselves
    {
        parent[i] = i;
        compressed_parent[i] = i;
        in_cluster[i] = false;
    }

    EdgeIdx *root_by_level = new EdgeIdx[g->nEdges+1];
    EdgeIdx *root_offsets = new EdgeIdx[g->nEdges+1]; //the roots for level c start at root_offsets[maxC - c]
    EdgeIdx rootIndex = 0;
    

    //get the highest c
    c = c-2;
    EdgeIdx maxC = c;
    while (c > 0)
    {
        for(EdgeIdx idxi = deleted_offsets[c]; idxi < deleted_offsets[c+1]; idxi++) //Iterate over each edge removed in the layer
        {
            //Get edge
            EdgeIdx e = deleted_edges[idxi]; //e is the i-th edge
            processing_time[e] = processing_time_counter;
            children_offset[processing_time_counter] = children_counter; 
            processing_time_counter++;
            
            //check if the edge has been already added to the tree
            if(in_cluster[e]) {
                continue;
            }

            //Otherwise we need to find all the descendants edges and clusters.
            findConnectedComponent(g, e, neighbors, neighbors_offset, deletion_time, parent, compressed_parent, children, children_offset, children_counter, in_cluster);

        }
        //get the roots for the level c
        root_offsets[maxC - c] = rootIndex;
        for(EdgeIdx idxi = deleted_offsets[c]; idxi < deleted_offsets[c+1]; idxi++) {//Iterate over each edge removed in the layer
            EdgeIdx e = deleted_edges[idxi];
            if(parent[e]==e || kappa[parent[e]]<c){
                root_by_level[rootIndex] = e;
                rootIndex++;
            }
        }

        printf("Generated Tree for c = %" PRId64 ":\n",c);
        printf("There are %" PRId64 " clusters at this level \n", (rootIndex - root_offsets[maxC - c]));
        c--;
    }
    root_offsets[maxC] = rootIndex;
    for(EdgeIdx idxi = deleted_offsets[0]; idxi < deleted_offsets[1]; idxi++){
        EdgeIdx e = deleted_edges[idxi];
        processing_time[e] = processing_time_counter;
        children_offset[processing_time_counter]= children_counter;
        processing_time_counter++;
    }
    children_offset[processing_time_counter] = children_counter; 
    //At this stage the tree contains the hierarchy of clusters, now need to go back and use the tree to compute the clusters

    //Generate the clusters
    nClusters = rootIndex;
    Cluster *clusters = new Cluster[nClusters+1];
    EdgeIdx *root_to_cluster = new EdgeIdx[nEdges+1];
    EdgeIdx clusterCounter = 0;
    c = maxC;
    while (c>0){
        printf("getting clusters for c = %" PRId64 "\n",c);
        for (EdgeIdx idxi = root_offsets[maxC-c]; idxi < root_offsets[maxC-c+1]; idxi++){
            //get the root of the cluster
            EdgeIdx root = root_by_level[idxi];
            Cluster *new_cluster = new Cluster();
            new_cluster->id = clusterCounter;
            root_to_cluster[root] = clusterCounter;
            new_cluster->c = c;
            new_cluster->n_edges = 0;
            new_cluster->n_vertices = 0;
            new_cluster->parent=-1;
            completeCluster(new_cluster, root_to_cluster, clusters, root, c, kappa, g, sources, children, children_offset, processing_time);
            clusters[clusterCounter] = *new_cluster;
            clusterCounter++;
        }
        c--;
    }
    return clusters;
}


#endif



