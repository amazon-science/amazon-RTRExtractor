#ifndef ESCAPE_CLUSTERSTRUCT_H
#define ESCAPE_CLUSTERSTRUCT_H

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include <set>
#include <unordered_set>
#include <list>

using namespace Escape;
using namespace std;

// Structure that stores weighted triangle information of graph

struct WeightedTriangleInfo
{

    Count total; // total number of triangles
    double total_wgt; // total weight of triangles
    double *perVertex;  // array storing total weight of triangles on each vertex
    double *perEdge;   // array storing weight of triangles on each edge

};

// Structure that stores a cluster

struct Cluster
{
    unordered_set<VertexIdx> vertices; // set storing vertices of a cluster
    Count nVertices; // the size of the cluster in vertices
    Count nEdges; // the number of edges inside the cluster
    Count cut; // the number of edges leaving the cluster
};

// Structure storing pair that has a vertex and a triangle weight

struct wgtPair
{
    VertexIdx vertex; // the vertex being stored
    double wgt; // the corresponding weight
};

// comparator that compares the weight in wgtPair objects. It will sort in DECREASING order
bool wgtCompareDecreasing(wgtPair firstPair, wgtPair nextPair)
{
    return firstPair.wgt > nextPair.wgt;
}

// comparator that compares the number of vertices in Cluster objects. It will sort in DECREASING order
bool clusterCompareDecreasing(Cluster first, Cluster second)
{
    return first.nVertices > second.nVertices;
}



#endif



