#include "Escape/GraphIO.h"
#include <inttypes.h>

using namespace Escape;


static bool isBlankLine( const char* line )
{
  while( *line )
  {
    if( !isspace( *line ) )
      return false;
    ++line;
  }
  return true;
}


static ErrorCode loadGraph_Escape(const char *path, Graph& graph, int undirected)
{
  FILE* f = fopen(path, "r");
  if (!f)
  {
    fprintf(stderr, "could not open file %s\n", path);
    return ecInvalidInput;
  }

  graph.nVertices = 0;
  EdgeIdx iEdge = 0;
  char line[1024];
  while (fgets(line, sizeof(line), f))
  {
    //Ignore comment lines.
    if (line[0] == '#')
      continue;

    if (isBlankLine(line))
      continue;
      
    int64_t i1, i2;
    sscanf(line, "%" PRId64 "%" PRId64 "", &i1, &i2);
    if (graph.nVertices == 0)
    {
      graph.nVertices = i1;
      if (undirected)
          graph.nEdges = 2*i2;
      else
          graph.nEdges = i2;
      graph.srcs = new VertexIdx[graph.nEdges];
      graph.dsts = new VertexIdx[graph.nEdges];
    }
    else
    {
      graph.srcs[iEdge] = i1;
      graph.dsts[iEdge] = i2;
      ++iEdge;
      if (undirected)
      {
          graph.srcs[iEdge] = i2;
          graph.dsts[iEdge] = i1;
          ++iEdge;
      }
    }
  }
  fclose(f);

  if (iEdge < graph.nEdges)
  {
    fprintf(stderr, "expected %" PRId64 " edges, only got %" PRId64 "\n", graph.nEdges, iEdge);
    return ecIOError;
  }
  
  return ecNone;
}



ErrorCode Escape::loadGraph(const char *path, Graph& graph, int undirected, IOFormat fmt)
{
  switch (fmt)
  {
    case IOFormat::escape:
      return loadGraph_Escape(path, graph, undirected);
    
    default:
      return ecUnsupportedFormat;
  }
}
