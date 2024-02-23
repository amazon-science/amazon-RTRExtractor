#include "Escape/TriangleProgram.h"


using namespace Escape;


//unfortunate that C++14 does not yet permit this to be declared inside
//countTriangles.
template <bool doPerVertex, bool doPerEdge>
struct CountFunctor
{
  Count total;
  Count *perVertex;
  Count *perEdge;

  void operator ()(VertexIdx u, VertexIdx v, VertexIdx w
    , EdgeIdx vu, EdgeIdx vw, EdgeIdx uw)
  {
    if (doPerVertex)
    {
      ++perVertex[u];
      ++perVertex[v];
      ++perVertex[w];
    }

    if (doPerEdge)
    {
      ++perEdge[vu];
      ++perEdge[vw];
      ++perEdge[uw];

      //We need to also increment counts for edges
      //going the other way. Looking up edges in here is bad.
      //What is the correct solution?  It's particularly bad for
      //undirected graphs.
    }

    ++total;
  }
};


//The per-edge count here is wrong.  Do not use until fixed.

//this may be useful to expose in its own right.
template <bool degreeOrdered, bool directed, bool doPerVertex, bool doPerEdge>
static Count countTriangles(const CGraph* gOut
  , const CGraph *gIn
  , Count *perVertex
  , Count *perEdge)
{
  auto cf = CountFunctor<doPerVertex, doPerEdge>{0, perVertex, perEdge};
  triangleProgram<degreeOrdered, true, !directed>(gOut, cf, gIn);
  return cf.total;
}


template <bool degreeOrdered, bool directed, bool doPerVertex>
static Count countTriangles1(const CGraph* gOut
  , const CGraph *gIn
  , Count *perVertex
  , Count *perEdge)
{
  if (perEdge)
    return countTriangles<degreeOrdered, directed, doPerVertex, true>(gOut, gIn, perVertex, perEdge);
  else
    return countTriangles<degreeOrdered, directed, doPerVertex, false>(gOut, gIn, perVertex, perEdge);
}


template <bool degreeOrdered, bool directed>
static Count countTriangles2(const CGraph* gOut, const CGraph *gIn
  , Count *perVertex, Count *perEdge)
{
  if (perVertex)
    return countTriangles1<degreeOrdered, directed, true>(gOut, gIn, perVertex, perEdge);
  else
    return countTriangles1<degreeOrdered, directed, false>(gOut, gIn, perVertex, perEdge);
}


template <bool degreeOrdered>
static Count countTriangles3(bool directed, const CGraph* gOut, const CGraph *gIn
  , Count *perVertex, Count *perEdge)
{
  if (directed)
    return countTriangles2<degreeOrdered, true>(gOut, gIn, perVertex, perEdge);
  else
    return countTriangles2<degreeOrdered, false>(gOut, gIn, perVertex, perEdge);
}


//This is just a wrapper that selects the right template parameters using the
//helper functions above.
Count Escape::countTriangles(const CGraph* gOut
  , bool degreeOrdered
  , bool directed
  , const CGraph *gIn
  , Count *perVertex
  , Count *perEdge)
{
  if (degreeOrdered)
    directed = true;

  if (!directed)
    gIn = gOut;
  else if (!gIn)
  {
    //compute transpose, todo
    printf("please implement this\n");
    exit(1);
  }

  if (perVertex)
    std::fill(perVertex, perVertex + gOut->nVertices, 0);

  if (perEdge)
    std::fill(perEdge, perEdge + gOut->nEdges, 0);

  if (degreeOrdered)
    return countTriangles3<true>(directed, gOut, gIn, perVertex, perEdge);
  else
    return countTriangles3<false>(directed, gOut, gIn, perVertex, perEdge);
}
