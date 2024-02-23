#ifndef ESCAPE_DECOMP_H_
#define ESCAPE_DECOMP_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/ClusterStructures.h"
#include <stack>
#include <list>
#include <queue>
#include <algorithm>
#include <chrono>
#include <vector>
#include <map>
#include <iostream>
#include <inttypes.h>


// all relevant structures are in ClusterStructures.h

using namespace Escape;
using namespace std;


// The common neighbor algorithm that produces the triangle weight
// Input: a pointer g to a CGraph that is sorted by ID (so binary search is possible).
// Output: a WeightedTriangleInfo for g. The ordering of edges in perEdge (of TriangleInfo) is that same as g.
// 
//   The algorithm does the common neighbor procedure for finding triangles. So it loops over every edge (u,v)
// and checks for common neighbors. The main trick is to search for the entries of the smaller list/degreee in
// the neighbor list of the other vertex. 

//   To ensure that every triangle is discovered exactly once, we use some tiebreaking rules. We process edge (u,v) only if
// u < v (in terms of id). Furthermore, we only consider the wedge (u,v,w) if w has the largest id. So triangle (u,v,w)
// is discovered exactly once: when we start with u, process edge (u,v), and intersect neighbor lists.
//
//   The weight of a triangle u, v, w is 1/{d(u)*d(v)*d(w)}, where d(u) is the degree of u. This is 
// based on defining the edge weight according to the normalized adjacency matrix (D^{-1/2}MD^{-1/2}).
//  


WeightedTriangleInfo commonNbr(CGraph *g)
{
   printf("Computing all triangle weighted info\n");

   WeightedTriangleInfo ret;   // output 
   // initialize output
   ret.total = 0;      
   ret.total_wgt = 0;  
   ret.perVertex = new Weight[g->nVertices+1];
   ret.perEdge = new Weight[g->nEdges+1]; 

   // initial all per vertex values to zero
   for (VertexIdx i=0; i < g->nVertices; ++i)
       ret.perVertex[i] = 0;

   // initial all per edge values to zero
   for (EdgeIdx j=0; j < g->nEdges; ++j)
       ret.perEdge[j] = 0; 

   // getting the partner map, so that we can update triangle weights more efficiently
   g->getPartnerMap(); // calling getPartnerMap, which returns an array

   // vertex variables
   VertexIdx v, w, lower, higher;
   EdgeIdx loc;
   Count degu, degv, degw; 
   Weight wgt;

   for (VertexIdx u=0; u < g->nVertices; ++u) // loop over vertices
       for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; ++i)   // loop over neighbor of u
       {
           v = g->nbors[i]; // get the neighbor v, present at location i in the nbors list

           degu = g->offsets[u+1] - g->offsets[u]; // get the degree of u
           degv = g->offsets[v+1] - g->offsets[v]; // get the degree of v

           // we will not process (u,v) if u > v. 
           if(u > v)
               continue;

           // we will now determine the neighbor with lower degree
           lower = (degu < degv) ? u : v; // checking if degu < degv
           higher = (degu < degv)? v : u; // getting the other vertex

           // loop of neighbors of lower, and look up these neighbors in the adj list for the other vertex
           for (EdgeIdx j = g->offsets[lower]; j < g->offsets[lower+1]; ++j)
           {
               w = g->nbors[j]; // get the neighbor w, preset at location j in the nbors list
               // we only consider this wedge if w has the largest id
               if (w <= u || w <= v)
                   continue;
               loc = g->getEdgeBinary(higher,w); // check if w is present in higher's nbors list
               if (loc != -1) // edge (upper, w) is present
               {
                   degw = g->offsets[w+1] - g->offsets[w]; // get the degree of w
                   wgt = 1/double(degu*degv*degw); // the weight of the triangle (u,v,w)
                   ret.total++; // increment the total number of triangles
                   ret.total_wgt += wgt; // add the weight to the total weight

                   ret.perEdge[i] += wgt; // add to the triangle weight of the edge (u,v), at location i
                   ret.perEdge[g->partnerMap[i]] += wgt; // edge (v,u), at the partner location for i
                   
                   ret.perEdge[j] += wgt; // (lower,w), at location j
                   ret.perEdge[g->partnerMap[j]] += wgt; // edge (w,lower), at the partner location for j

                   ret.perEdge[loc] += wgt; // (higher,w) at location loc
                   ret.perEdge[g->partnerMap[loc]] += wgt; // edge (w,higher), at the partner location for loc
                   
                   ret.perVertex[u] += wgt; // add to the triangle weight of vertex u
                   ret.perVertex[v] += wgt; // add to the triangle weight of vertex v
                   ret.perVertex[w] += wgt; // add to the triangle weight of vertex w
               }
           }
       }
    return ret;
}


/* deleteAndClean procedure: This is the real workhorse. It get a list of edges to delete (either from extraction or those that have already been detected
to be below the cleaning threshold). It deletes these edges, removes the corresponding triangles, and cleans all *adjacent* edges. Note that these are the
only edges whose triangle weights are affected. This process is propagated until convergence. Meaning, adjacent edges below the cleaning threshold
are deleted, which may lead to other adjacent edges being deleted, so on and so forth.

Input:
    g: A pointer to a CGraph
    edgeStatus: A pointer to a char array that stores the current status of the edges (Y = in the graph, N = not in the graph, S = in the stack). So the current graph has all 'N' edges already deleted.
    triInfo: A pointer to weighted triangle info structure, where weight are of the *current* graph. (So only triangles with non 'N' edges contribute to the weight.)
    toDelete: A pointer to a stack of edges to be deleted.
    eps: The cleaning threshold

Input assumption:
    - edgeStatus and toDelete are consistent. This means that edgeStatus[i] = 'S' iff the corresponding edge is in toDelete.
    - edgeStatus and triInfo are consistent. So triInfo only contains the weight of triangles whose edgeStatus is not 'N'

Output: void BUT the following inputs get changed
    The input assumptions will remain valid at the end of the procedure.
    - edgeStatus: edges will get deleted and be labeled 'N'. All edges in toDelete will get labeled 'N'.
    - triInfo: As edges get deleted, triInfo will be updated accordingly.
    - toDelete: This will become empty

The algorithm is as follows. We go over every edge in the stack. We list all the triangles that edge participates in. All the triangles
(with non 'N' edges) need to be deleted, so we reduce the triangle weight on all edges of the triangle. If any of these edges has
a triangle weight that goes below the cleaning threshold, we push it onto the slack toDelete. We keep continuing until the stack is empty.
*/

void deleteAndClean(CGraph* g, char* edgeStatus, WeightedTriangleInfo* triInfo, stack<Pair>* toDelete, double eps)
{
    // keep looping until the stack is empty
    while(!toDelete->empty())
    {
        Pair edge = toDelete->top(); // get the next edge to be deleted
        toDelete->pop(); // remove it from stack

        VertexIdx u = edge.first;  // get the endpoints of the edge
        VertexIdx v = edge.second; 
        EdgeIdx loc_uv = g->getEdgeBinary(u,v); // get the location of this edge
        if (loc_uv == -1) // edge was not found, so there is some error
        {
            printf("Error in deleteAndClean: edge (u,v) in stack, but v is not u's neighbor\n");
            printf("u = %" PRId64 ", v = %" PRId64 "\n",u,v);
        }

        EdgeIdx partloc_uv = g->partnerMap[loc_uv]; // get the location of the partner (v,u)

        if (edgeStatus[loc_uv] != 'S' || edgeStatus[partloc_uv] != 'S') // the edge status of either the edge or the partner is NOT 'S'. There's a problem, so throw an error
        {
            printf("Error in deleteAndClean: edge (u,v) in stack toDelete, but edgeStatus does not show that\n");
            printf("u = %" PRId64 ", v = %" PRId64 ", loc_uv = %" PRId64 ", partloc_uv = %" PRId64 ", edgeStatus[loc_uv] = %c, edgeStatus[partloc_uv] = %c\n",u,v,loc_uv,partloc_uv,edgeStatus[loc_uv],edgeStatus[partloc_uv]);
        }

        // at this point, we have pick out the edge at the top of the toDelete stack and found all the relevant information
        // it's time to get the triangles incident to (u,v) and remove them from the graph
       
        // first, we get the degree of u and v 
        Count degu = g->offsets[u+1] - g->offsets[u]; // get the degree of u
        Count degv = g->offsets[v+1] - g->offsets[v]; // get the degree of v

        // we will now determine the vertex with lower degree
        VertexIdx lower = (degu < degv) ? u : v; // checking if degu < degv
        VertexIdx higher = (degu < degv)? v : u; // getting the other vertex

        // just to have the degrees clearly
        Count deglower = (degu < degv) ? degu : degv; // get the lower degree
        Count deghigher = (degu < degv) ? degv : degu; // get the lower degree

        // loop of neighbors of lower, and look up these neighbors in the adj list for the other vertex
        for (EdgeIdx j = g->offsets[lower]; j < g->offsets[lower+1]; ++j)
        {
            if (edgeStatus[j] == 'N') // the current edge is deleted, so continue to the next edge
                continue;

            VertexIdx w = g->nbors[j]; // get the neighbor w, present at location j in the nbors list
            if (w == u || w == v) // seeing the vertex again, so continue
                continue;
            
            EdgeIdx loc_hiw = g->getEdgeBinary(higher,w); // check if w is present in higher's nbors list
            if (loc_hiw == -1) // edge was not found, so we just continue
                continue;

            if (edgeStatus[loc_hiw] == 'N') // edge is already deleted, so continue
                continue;

            // at this stage, the triangle (u,v,w) is present (none of the edges are deleted)
            // The edges are at indices loc_uv, j, and loc_hiw. Specifically, (lower,w) is at j and (higher,w) is at loc_hiw

            // let's get the indices of the partner edges
            EdgeIdx partj = g->partnerMap[j];
            EdgeIdx partloc_hiw = g->partnerMap[loc_hiw];

            // let us get the weight of the triangle (u,v,w)
            Count degw = g->offsets[w+1] - g->offsets[w]; // get the degree of w
            double wgt = 1/double(degu*degv*degw); // the weight of the triangle (u,v,w)

            // we decrement this weight from all the relevant triInfo content
            // first, the edges
            triInfo->perEdge[loc_uv] -= wgt;
            triInfo->perEdge[j] -= wgt;
            triInfo->perEdge[loc_hiw] -= wgt;
            // now for the partners
            triInfo->perEdge[partloc_uv] -= wgt;
            triInfo->perEdge[partj] -= wgt;
            triInfo->perEdge[partloc_hiw] -= wgt;

            // next, the vertices
            triInfo->perVertex[u] -= wgt;
            triInfo->perVertex[v] -= wgt;
            triInfo->perVertex[w] -= wgt;

            // finally, update totals
//             printf("Decrementing total\n");
            triInfo->total--;
            triInfo->total_wgt -= wgt;

            // now, we check if any of the edges (u,w) or (v,w) have gone below the cleaning threshold. If so, we push them onto the stack
            // It doesn't matter which copy of the edge we push, but we need to mark the edgeStatus properly

            // dealing with the edge (lower,w). We check that the edge is not on the stack and has weight below the cleaning threshold
            double edge_wgt = 1/double(deglower*degw);
            if(edgeStatus[j] != 'S' && triInfo->perEdge[j] < eps*edge_wgt) // (lower,w) has triangle weight below the cleaning threshold
            {
                Pair next; // set up a pair to push onto the stack
                next.first = lower;
                next.second = w;
                toDelete->push(next); // pushed the pair onto the stack

                // now update the edgeStatus for the edge and it's partner
                edgeStatus[j] = 'S';
                edgeStatus[partj] = 'S';
            }
            
            // we repeat the above process with the edge (higher,w). Maybe I should have written a function...? But I would have to pass all the data structures, so it feels like it would just add clunkiness
            edge_wgt = 1/double(deghigher*degw);
            if(edgeStatus[loc_hiw] != 'S' && triInfo->perEdge[loc_hiw] < eps*edge_wgt) // (higher,w) has triangle weight below the cleaning threshold
            {
                Pair next; // set up a pair to push onto the stack
                next.first = higher;
                next.second = w;
                toDelete->push(next); // pushed pair onto the stack

                // now update the edgeStatus for the edge and it's partner
                edgeStatus[loc_hiw] = 'S';
                edgeStatus[partloc_hiw] = 'S';
            }

            // at this point, we have processed the triangle (u,v,w) by removing it's weight from all the incident edges. We have checked the two candidate edges that may need to be removed
            // Processing is finished, so the loop moves to the next iteration (the next w)
        }
        // at this point, we have processed all the neighbors of (the lower endpoint of) edge (u,v)

        // as a final check, we make sure that triangle weight on edge (u,v) is actually zero. If not, there is some error
        // this is a double, so we check if it's more than 1/n
        if(abs(triInfo->perEdge[loc_uv]) > 1/double(g->nVertices) || abs(triInfo->perEdge[partloc_uv]) > 1/double(g->nVertices))
        {
            printf("Error in deleteAndClean. The triangle weight on an edge about to be deleted is non-zero.\n");
            printf("u = %" PRId64 ", v = %" PRId64 ", loc_uv = %" PRId64 ", partloc_uv = %" PRId64 ", triangle weight on (u,v) = %f, triangle weight on (v,u) = %f\n",u,v,loc_uv,partloc_uv,triInfo->perEdge[loc_uv], triInfo->perEdge[partloc_uv]);
        }

        // the actual edge removal
        edgeStatus[loc_uv] = 'N';
        edgeStatus[partloc_uv] = 'N';
    }
}

/* initialClean procedure: it applies the repeated Jaccard based cleaning procedure. Every edge
with Jaccard similarity less than a parameter epsilon is removed from the graph. This process
is done iteratively until convergence. When the process ends, all remaining edges have
Jaccard value about epsilon.

Input:
    g: a pointer to a CGraph
    trInfo: a pointer to the corresponding WeightedTriangleInfo object
    eps: the cleaning threshold

Output:
    A char/byte array (with the same index as g->nbors) indicating which edges remain.
    A value of Y means edge is still present (not deleted), and N means edges was deleted during cleaning.
    We use a char array, since we will also store other flags. This is used in subsequent functions to mark intermediate steps of cleaning.

*/

char* initialClean(CGraph* g, WeightedTriangleInfo* triInfo, double eps)
{
    printf("Starting the initial clean\n");
    // set up the return value, which is a boolean array
    char* edgeStatus = new char[g->nEdges+1];
    // no edge is deleted, so initialize all entries as Y
    for(EdgeIdx i = 0; i <= g->nEdges; i++)
        edgeStatus[i] = 'Y';


    stack<Pair> toDelete; // setting up a stack to store the edges that need to be deleted
    Count degu, degv;

    // loop over all the edges and find the ones whose Jaccard similarity is below the threshold eps
    // all these edges will be pushed onto the slack toDelete, and marked as "on the stack" ('S') in the array edgeStatus

    for(VertexIdx u=0; u < g->nVertices; u++)
       for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; ++i)   // loop over neighbor of u
       {
           VertexIdx v = g->nbors[i]; // get the neighbor v, present at location i in the nbors list

           degu = g->offsets[u+1] - g->offsets[u]; // get the degree of u
           degv = g->offsets[v+1] - g->offsets[v]; // get the degree of v

           // we will not process (u,v) if u > v. 
           if(u > v)
               continue;

           // create a pair storing the edge
           Pair edge;
           edge.first = u; // first end is u
           edge.second = v; // second endpoint is v

           double edgeWt = 1/double(degu*degv); // calculate the weight of the edge
//            printf("The edge data: u = %" PRId64 ", v = %" PRId64 ", i = %" PRId64 ", edgeWt = %f, triangle weight = %f\n",u,v,i,edgeWt,triInfo->perEdge[i]);
           if(edgeStatus[i] != 'S' && triInfo->perEdge[i] < eps*edgeWt) // if the edge is not already on stack and triangle weight is less than eps*edgeWt, this edge should be cleaned
           {
//                printf("The initial pushes: u = %" PRId64 ", v = %" PRId64 ", i = %" PRId64 ", edgeWt = %f, triangle weight = %f\n",u,v,i,edgeWt,triInfo->perEdge[i]);
               toDelete.push(edge); // we push the edge onto the stack
               edgeStatus[i] = 'S'; // the edge is marked S, to denote it is in the stack
               edgeStatus[g->partnerMap[i]] = 'S'; // also mark the partner as 'S'
           }
               
       }

    // we have populated the set toDelete with the edges below the threshold
    // we will now call the real workhorse to do the actual edge deletions and cleaning

    deleteAndClean(g, edgeStatus, triInfo, &toDelete, eps); // call the function, and set ret to be the array indicating which edges remain

    // print out total triangle count
    printf("Total triangle count, post cleaning = %" PRId64 "\n",triInfo->total);


    // edgeStatus is updated by the previous function, so return it as the output
    return edgeStatus;
}

/* populateStats procedure: given a cluster object, this function just populates the properties of the cluster

Input:
    g: Pointer to a CGraph
    clus: Pointer to a cluster object

Output:
    void, but it changes the clus object

It simply goes over all edges incident to clus, and checks which of them are contained inside the cluster
*/
void populateStats(CGraph* g, Cluster* clus)
{
    clus->nVertices = (clus->vertices).size(); // call the size function to get the number of vertices
    clus->nEdges = 0; // initialize the value to zero
    clus->cut = 0; // initialize the value to zero
    for (auto iter = begin(clus->vertices); iter != end(clus->vertices); iter++) // iterate over the set of vertices given by clus->vertices
    {
        VertexIdx u = *iter; // dereference to get the vertex
        for (EdgeIdx i = g->offsets[u]; i < g->offsets[u+1]; i++) // loop over all the neighbors of the vertex u
        {
            VertexIdx v = g->nbors[i]; // get the vertex v
            if ((clus->vertices).find(v) != (clus->vertices).end()) // cluster has the neighbor v
                clus->nEdges++; // increment the number of edges
            else // this is a cut edge, since the other endpoint lies outside the cluster
                clus->cut++;
        }
    }
    clus->nEdges = clus->nEdges/2; // divide by two, because each edge is counted twice
}





/* disjointExtract procedure: this follows the spectral triadic decomposition procedure.
First, it computes the weighted triangle info, and then runs the initial clean operation. At this stage,
every edge has a triangle weight at least eps times the edge weight.
Then, the extract happens. It process vertices in increasing order of degree, and extract
clusters centered at the vertices. This is done repeatedly until the graph is empty.

Input:
    g: Pointer to a CGgraph
    eps: the decomposition parameter, which is the cleaning threshold
    cluster_map: As reference, keeps a mapping from vertex id to cluster

Output:
    A set of Cluster objects, that partition all vertices in the graph.
    There may be singleton clusters for vertices that do not produce anything non-trivial.

The algorithm works as follows. When a vertex u is processed, we look at all non-deleted
neighbors v. If the degree of v is at most (the degree of u)/eps, then we add it to the cluster. 
Next, we iterate over all *current* triangles involving two neighbors v, w of u. (If any triangle
has a removed edge, we ignore it.) Each such triangle involves a
vertex x that is two hops away from u. For all x, we compute the total triangle weight of such triangles.
This is done by iterating of all triangles (v,w,x) and updating an unordered map keyed by x.

Finally, we have a map that gives, for all (relevant x), the triangle weight incident on x by
two vertices in the neighborhood of u. (Note that there are precisely the triangles that are not
in the current cluster.) We sort the x's in decreasing order of this triangle weight.

We add these vertices x to the cluster one by one, keeping track of (i) the triangle weight incident on x
and (ii) the total triangle weight of the current cluster. We maximize the ratio: 
    
    (total triangle weight inside neighborhood of u and triangle weight contributed by two-hop vertices)/(total triangle weight of all vertices in the cluster)

We simply process ther vertices x in decreasing order (as said before), and find the maximum. So it's just a linear
loop over all the x vertices. 

At this point, we have created a new cluster. We put all edges inside the cluster into a stack, and send this
to deleteAndClean. This function will remove all these edges and clean to convergence. We repeat this loop
to get the next cluster, so on and so forth. At the end, all edges will have been deleted.

C. Seshadhri, Aug 2023
*/

vector<Cluster> disjointExtract(CGraph* g, double eps, std::map<VertexIdx, VertexIdx>& cluster_map)
{

    auto start = std::chrono::high_resolution_clock::now();

    vector<Cluster> decomposition; // the return value, which is a list of clusters

    for (VertexIdx i=0; i < g->nVertices; i++) { //initialize the cluster map to -1 (unclustered)
        cluster_map[i] = -1; 
    }
    
    //The partner map is already generated in the commonNbr() step, so it is redundant (commenting this saves 30min in large graphs)
    //printf("Getting partner map\n");
    //g->getPartnerMap();

    WeightedTriangleInfo triInfo = commonNbr(g); // get the weighted triangle info

    auto stop = std::chrono::high_resolution_clock::now();
    cout << "Time for getting the weights: " << chrono::duration<double>(stop - start).count() << "s \n";
    start = std::chrono::high_resolution_clock::now();


    char* edgeStatus = initialClean(g, &triInfo, eps); // applying the initial clean, to get the current edgeStatus. This array has a 'Y' for every existing edge, and an 'N' for deleted edge. No other symbol should be present

    stop = std::chrono::high_resolution_clock::now();
    cout << "Time for the initial clean: " << chrono::duration<double>(stop - start).count() << "s \n";
    start = std::chrono::high_resolution_clock::now();

    // our next step is to sort all the vertices in increase ordering of degree, to set up the order of extraction

    Pair *deg_info = new Pair[g->nVertices]; // allocate memory for the vertices and their degrees

    // Construct array of pairs, storing vertex label and degree
    for (VertexIdx i=0; i < g->nVertices; i++) // looping over all vertices
    {
        deg_info[i].first = i; // storing the vertex id as the first item in pair
        deg_info[i].second = g->offsets[i+1] - g->offsets[i]; // storing the degree as the second item in pair
    }
    
    // sort the pairs by degree, using comparator (if degree is same, sort by old vertex label)
    std::sort(deg_info,deg_info + g->nVertices,pairCompareSecond); 
       
    // at this stage, the first entries of deg_info are the vertices sorted by increasing degree

    double* triwgt = new double[g->nVertices]; // allocating memory to store triangle weight of two-hop neighborhood
    for (VertexIdx i=0; i < g->nVertices; i++)
        triwgt[i] = 0; // initializing to zero

    // we create a bitmap to store vertex neighborhoods. Whenever a vertex is processed, we create this bitmap. At the end, we clear it out.
    bool* nbdBitMap = new bool[g->nVertices]; // a bitmap to store the neighborhood of a vertex
    bool* clustered = new bool[g->nVertices]; // a bitmap to mark the set of clustered vertices
    for (VertexIdx i=0; i < g->nVertices; i++)
    {
        nbdBitMap[i] = false;
        clustered[i] = false;
    }

    int n_clusters = 0; //keeps track of the total number of clusters found so far

    for (VertexIdx i=0; i < g->nVertices; i++) // looping of entries of deg_info
    {
        VertexIdx u = deg_info[i].first; // get the next vertex
        Count degu = deg_info[i].second; // get the degree of u, as the second entry of pair

//         printf("Seed degree %" PRId64 "\n",degu);

        if (clustered[u]) // vertex u has been clustered, so continue
            continue;

        Cluster next_clus; // setting up the next cluster
        next_clus.vertices.insert(u); // insert u into the cluster
        cluster_map[u] = n_clusters; //update the value of u in the map
        clustered[u] = true; // mark u as clustered

        set<VertexIdx> two_hop; // creating set with two hop neighborhood
        double internalWgt = 0; // this tracks the internal triangle weight of the cluster


        // an initial loop over the neighbors of u to set up the neighborhood bitmap, and populate the cluster
        for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) // looping over neighbors of u
        {
            VertexIdx v = g->nbors[j]; // storing the neighbor of u
            
            if (edgeStatus[j] == 'N') // edge is removed, so proceed to next edge
                continue;

            next_clus.vertices.insert(v); // insert v into the cluster
            nbdBitMap[v] = true; // insert into the bitmap
            cluster_map[v] = n_clusters; //update the value of v in the map
            clustered[v] = true; // mark v as clustered
        }

        // now for the real loop over the neighbors
        for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) // looping over neighbors of u
        {
            VertexIdx v = g->nbors[j]; // storing the neighbor of u
            
            if (edgeStatus[j] == 'N') // edge is removed, so proceed to next edge
                continue;
            
            Count degv = g->offsets[v+1] - g->offsets[v]; // get the degree of v
            if ((double)degv > (double)degu/eps) // v has much higher degree, then continue
                continue;

            // we now try to find the 2-hop neighbors that are also in the cluster. 
            // We find all the edges (v,w) in the neighborhood of u. From (v,w), we enumerate all triangles (v,w,x) to get the potential candidates for the cluster
            // To find w, we do an enumeration of another neighbor of u, so the index starts from j. So each edge (v,w) is generated exactly once.
            // We check if (v,w) is an 'Y' edge, and if so, enumerate triangles incident on (v,w)
            for (EdgeIdx k=j; k < g->offsets[u+1]; ++k) // looping over neighbors of v
            {
                VertexIdx w = g->nbors[k]; // storing the next neighbor 

                if (edgeStatus[k] == 'N') // edge is removed, so proceed to next edge
                    continue;
                Count degw = g->offsets[w+1] - g->offsets[w]; // get the degree of w

                if ((double)degw > (double)degu/eps) // v has much higher degree, then continue
                    continue;

                EdgeIdx loc = g->getEdgeBinary(u,w); // look for w in the neighbors of u
                if (loc == -1 || edgeStatus[loc] == 'N') // if the edge is not present, or has been deleted, continue
                    continue;

                // at this point, (u,v,w) is a triangle where all edges are labeled 'Y'

                internalWgt += 1/(degu*degv*degw); // add the triangle weight to the internal weight

                // We will enumerate over the triangles that (v,w) participates in. We start by finding the lower degree vertex.
                VertexIdx lower = (degv < degw) ? v : w; // the lower degree vertex is v
                VertexIdx higher = (degv < degw) ? w : v; // the higher degree vertex is w

                // enumerate over neighbors of lower, and check if they are neighbors of higher
                for (EdgeIdx ell = g->offsets[lower]; ell < g->offsets[lower+1]; ell++)
                {
                    if(edgeStatus[ell] == 'N') // this edge has been removed, continue
                        continue;
                    VertexIdx x = g->nbors[ell]; // getting the actual vertex
                    Count degx = g->offsets[x+1] - g->offsets[x]; // get the degree of x

                    if(x == lower || x == higher || x == u) // x is just a copy, or x is u, we can ignore
                        continue;

                    EdgeIdx loc_hix = g->getEdgeBinary(higher,x); //
                    if(loc == -1 || edgeStatus[loc_hix] == 'N') // if edge is not present, or has been deleted
                        continue;
                    
                    // at this stage (v,w,x) is a triangle. x is not u. x is either a neighbor of u or a two-hop neighbor 
                    // if x is a neighbor, then this is triangle that we should add to the internal triangle weight. 
                    // To ensure that we do not overcount the triangles, we only perform this addition if x is the highest label vertex

                    double wgt = 1/(degx*degv*degw); // the weight of the triangle

                    if (nbdBitMap[x] == true) // so x is a neighbor
                    {
                        if (x > v && x > w) // so x is the highest label vertex
                            internalWgt += wgt; // add this triangle to the internal weight
                        continue; // so we do not process this triangle further, and more on to the next triangle
                    }

                    // at this point, (v,w,x) is a triangle and x is a two-hop neighbor
                    two_hop.insert(x); // insert x into the potential two-hop set
                    triwgt[x] += wgt; // increment the triangle weight at x by this triangles weight
                    
                    // loop ends. The triangle (v,w,x) has been processed.
                }
                
                // at this point, we have processed all triangles incident to (v,w), finding all the two-hop neighbors that triangles incident to (v,w) lead to.
            }
            // at this point, we have processed all choice of w, for a given v. Meaning, we have processed all triangles incident to (u,v) and following up on other triangles they lead to.
        }
        // at this point, we have processed all neighbors of u. We have determine all the two-hop neighbors relevant to the cluster.
        // The data structures two_hop and triwgt have been populated. We now have to determine which vertices in two_hop will enter the cluster

        // The first step is to put the information from two_hop and triwgt into an array of pairs. We sort in decreasing order of contributions to triangle weight

        wgtPair* candidates = new wgtPair[two_hop.size()]; // an array of candidates to add into the cluster. We allocate the right amount of memory, since we know exactly how many vertices we're dealing with
        VertexIdx index = 0; //setting up an index to populate candidates
        double potentialWgt = 0; // this track the potential weight we can gain by the two-hop vertices

        for (auto iter = begin(two_hop); iter != end(two_hop); ++iter) // loop over the set two_hop using an iterator
        {
            VertexIdx curr = *iter; // get the current vertex
            wgtPair curr_pair; // the current pair
            curr_pair.vertex = curr; // store the current vertex
            curr_pair.wgt = triwgt[curr]; // store the current triangle weight
            // now, we store this pair and update the index
            candidates[index] = curr_pair;
            index++;
         
            potentialWgt += triwgt[curr]; // add the potential weight that can gained.
            triwgt[curr] = 0; // VERY IMPORTANT! We need to reset the triwgt array, since it is reused in every cluster iteration
        }
        // we have populated the candidates array. We just need to sort it and process accordingly.
        sort(candidates,candidates+index,wgtCompareDecreasing); // sort the candidates in decreasing order of their weight

        // we now set up the processing of adding vertices to the cluster. The array candidates has all the two-hop vertices sorted in decreasing order of their potential contributions
        // Suppose we just take u and its neighborhood. We will gain internalWgt, and lose potentialWgt. 
        // So our initial objective is the ratio internalWgt/(internalWgt + potentialWgt).
        //
        // When we insert vertex x (from candidates) into our cluster, we gain triwgt[x]. Technically, this has been zeroed out. 
        // This is stored as the second entry in the pair of candidates. But the denominator now includes triInfo.perVertex[x].
        // Suppose we include a subset X from candidates. Our objective is the ratio:
        //
        //              internalWgt + \sum_{x \in X} (tri-wgt of x into nbd of u)
        //              _________________________________________________________
        //              internalWgt + potentialWgt + \sum_{x \in X} (total triangle weight incident to x)
        //
        // Technically, there is an overcount. Some triangles incident to x can also be counted in the potentialWgt. Triangles incident to two x's are counted twice. But for algorithmic simplicity, we will ignore this.
        // To optimize over all subset X, we will do a simple heuristic. We have all candidates in reverse sorted order of their tri-wgt values, and we simply optimize over all prefixes of this order.
        //

        double ratio = internalWgt/(internalWgt + potentialWgt); // the current ratio. We are trying to maximize
        int maxind = -1; // setting up the index of the maximum. -1 means that we don't take any vertices from candidates

        double num_sum = 0, den_sum = 0; // setting up the numerator and denominator sums

        VertexIdx size = two_hop.size();
        for (index = 0; index < size; index++) // loop over the entries in candidates
        {
            num_sum += candidates[index].wgt; // adding the candidate weight
            den_sum += triInfo.perVertex[candidates[index].vertex]; // get the total triangle weight incident to that vertex

            double new_ratio = (internalWgt + num_sum)/(internalWgt + potentialWgt + den_sum); // the next ratio on including the vertex
            if (new_ratio > ratio) // so we have increase our objective
            {
                ratio = new_ratio; // update the value of the ratio to the improved version
                maxind = index; // record the current index
            }
        }

        // at this point, maxind will be the index of the prefix among candidates that optimizes (maximizes) the ratio.
        // So we loop over the candidates up to maxind and add of the those

        delete[] candidates; // free the memory 

        for (index=0; index <= maxind; index++) // loop over candidates, just going up maxind
        {
            VertexIdx next_vert = candidates[index].vertex; // get the vertex to be added
            next_clus.vertices.insert(next_vert); // insert the vertex into the cluster
            cluster_map[next_vert] = n_clusters; //update the value of next_very in the map
            clustered[next_vert] = true; // mark the vertex as clustered
        }
        
        // finally, a loop over the neighbors of u to clear up the neighborhood bitmap
        for (EdgeIdx j=g->offsets[u]; j < g->offsets[u+1]; ++j) // looping over neighbors of u
        {
            VertexIdx v = g->nbors[j]; // storing the neighbor of u
            nbdBitMap[v] = false; // reset the bitmap
        }

        // we have constructed the cluster, in the set next_clus.
        // We now write out the basic stats of the cluster, done in a separate function. Note that the function modifies the stats in the next_clus object.
        populateStats(g,&next_clus); 

        // We now add it the cluster to the decomposition.
        decomposition.push_back(next_clus); // replaced by push_back as it is a vector
        n_clusters++; //increase the counter of clusters
//         printf("Cluster size %" PRId64 "\n",next_clus.vertices.size());

        // we now need to remove all the edges incident to the decomposition, and perform the corresponding cleaning
        stack<Pair> toDelete; // setting up a stack of pairs to delete
        for (auto iterc = begin(next_clus.vertices); iterc != end(next_clus.vertices); iterc++) // loop over vertices the cluster
        {
            VertexIdx z = *iterc; // get the vertex under consideration
            for (EdgeIdx j = g->offsets[z]; j < g->offsets[z+1]; j++) // loop over the neighbors of z
            {
                if (edgeStatus[j] != 'Y') // this edge is either already removed or on the stack
                    continue; // move on the next edges

                Pair next; // set up a pair to put on the stack
                next.first = z;
                next.second = g->nbors[j]; // putting the other vertex in the pair

                toDelete.push(next); // put this pair on the stack
                edgeStatus[j] = 'S'; // mark this edge on the stack
                edgeStatus[g->partnerMap[j]] = 'S'; // mark the partner edge as being on the stack
            }
        }

        // we have set up the stack toDelete, and updated edgeStatus accordingly. 
        // So we call the deleteAndClean function to update the graph
        deleteAndClean(g, edgeStatus, &triInfo, &toDelete, eps);
    }

    stop = std::chrono::high_resolution_clock::now();
    cout << "Time for getting the clusters: " << chrono::duration<double>(stop - start).count() << "s \n";

    // we have looped over all the vertices and decomposed the graph. We return the decomposition.
    return decomposition;
}


/* expandClusters procedure: This procedure expands the existing clusters adding isolated vertices that are strongly connected to a cluster.

Input:
    g: Pointer to a CGgraph
    decomposition: The decomposition result from the triadic algorithm (or any decomposotion) of type vector<Cluster>
    cluster_map: A map from vertex to cluster id showing the cluster membership

Output:
    A set of Cluster objects, that partition all vertices in the graph.

For each non-clustered vertex the algorithm will find the cluster id of all the neighbours and add that vertex to the cluster with more neighbours
if there are at least threshold neighbours.

Daniel Paul Pena, Sep 2023
*/
vector<Cluster> expandClusters(CGraph* g, vector<Cluster> decomposition, std::map<VertexIdx, VertexIdx> cluster_map, int threshold, std::map<VertexIdx, double>& scores)
{
    VertexIdx size = decomposition.size();
    for (VertexIdx i = 0; i < size; i++) { //loop over all vertices

        if(decomposition[i].nVertices > 1){
            continue;
        }
        VertexIdx v  = *decomposition[i].vertices.cbegin(); //v is the only vertex of a cluster of size 1
        std::map<VertexIdx, VertexIdx> counts;

        VertexIdx best = -1;
        VertexIdx maxScore = threshold; //We require at least threshold neighbours in the cluster
        //the vertex is not clustered
        EdgeIdx number_edges = g->offsets[v+1] - g->offsets[v];
        for (EdgeIdx j=g->offsets[v]; j < g->offsets[v+1]; ++j) // looping over neighbors of v
        {
            VertexIdx u = g->nbors[j];
            if (cluster_map[u] == -1){
                continue;
            }
            VertexIdx cluster_id = cluster_map[u]; //neighbour u belongs to cluster cluster_id
            if (counts.count(cluster_id) == 0) {
                counts[cluster_id] = 0;
            }
            counts[cluster_id]++;

            //Update the best neighboring cluster if neccesary
            if (counts[cluster_id] > maxScore){
                best = cluster_id;
                maxScore = counts[cluster_id];
            }
        }
        //best is the closest neighbour for vertex v, check it is not -1
        if (best == -1){
            //In this case no cluster went above the threshold, we keep the vertex in its own cluster
            continue;
        }
        //we need to add v to the cluster best and remove it from the original cluster
        decomposition[best].vertices.insert(v);
        decomposition[i].vertices.erase(v);
        scores[v] = (double) maxScore/number_edges;
    }
    //to finish we recompute the stats of all the clusters
    for (VertexIdx i = 0; i < size; i++) {
        populateStats(g, &decomposition[i]);
    }

    return decomposition;

}

/* connectedComp procedure: This finds the connected components of the graph, and stores them as a list of clusters.

Input:
    g: Pointer to a CGgraph

Output:
    A set of Cluster objects, that partition all vertices in the graph.

The procedure is a standard bfs.

C. Seshadhri, Aug 2023
*/

vector<Cluster> connectedComp(CGraph* g)
{
    vector<Cluster> decomposition;

    bool* visited = new bool[g->nVertices]; // a bitmap storing the visited status
    for (VertexIdx u=0; u < g->nVertices; u++) // initializing the visited array
        visited[u] = false;

    queue<VertexIdx> q; // the initially empty queue of vertices

    for (VertexIdx u=0; u < g->nVertices; u++) // loop over the vertices
    {
        if (visited[u]) // u is already visited
            continue; // move on to the next vertex

        Cluster next_clus; // start the next cluster
        next_clus.vertices.insert(u); // insert u into the cluster
        visited[u] = true; // mark u as visited
        q.push(u); // push u into the queue to initialize

        while(!q.empty()) // while the queue is non-empty
        {
            VertexIdx next = q.front(); // getting the next vertex of the q
            q.pop(); // also remove the element
            for (EdgeIdx j = g->offsets[next]; j < g->offsets[next+1]; j++) // loop over the neighbors of j
            {
                VertexIdx nbr = g->nbors[j]; // get the neighbor under question
                if (!visited[nbr]) // the neighbor has not been visited
                {
                    q.push(nbr); // add the neighbor to the q
                    visited[nbr] = true; // we have visited the neighbor
                    next_clus.vertices.insert(nbr); // get the neighbor into the cluster
                }
            }
        }
        // we have constructed the cluster, in the set next_clus.
        // We now write out the basic stats of the cluster, done in a separate function. Note that the function modifies the stats in the next_clus object.
        populateStats(g,&next_clus); 

        // We now add it the cluster to the decomposition.
        decomposition.push_back(next_clus); // replaced by push_back as it is a vector
    }

    return decomposition;

}















#endif


