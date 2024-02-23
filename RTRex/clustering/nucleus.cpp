#include "Escape/GraphIO.h"
#include "Escape/Digraph.h"
#include "Escape/Nucleus.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <inttypes.h>
using namespace std::chrono;

using namespace Escape;

int main(int argc, char *argv[])
{
  Graph g;
  auto start = high_resolution_clock::now();
  if (loadGraph(argv[1], g, 1, IOFormat::escape))
    exit(1);

  printf("Loaded graph\n");
  auto stop = high_resolution_clock::now();
  cout << "time: " << chrono::duration<double>(stop - start).count() << "s \n";

  CGraph cg = makeCSR(g);
  cg.sortById();
  printf("Converted to CSR\n");
  stop = high_resolution_clock::now();
  cout << "time: " << chrono::duration<double>(stop - start).count() << "s \n";

  printf("Relabeling graph\n");
  VertexIdx *inverse = new VertexIdx[cg.nVertices];
  CGraph cg_relabel = cg.renameByDegreeOrderMapping(inverse);
  cg_relabel.sortById();

  printf("Creating DAG\n");
  CDAG dag = degreeOrdered(&cg_relabel);

  (dag.outlist).sortById();
  (dag.inlist).sortById();

  printf("Original data\n");
  printf("----------------------\n");
  printf("%" PRId64 " edges\n",cg_relabel.nEdges);
  printf("----------------------\n");

  printf("Counting triangles\n");
  TriangleInfo tri_info = betterWedgeEnumerator(&(dag.outlist));
  printf("Found %" PRId64 " triangles\n",tri_info.total);
  stop = high_resolution_clock::now();
  cout << "time: " << chrono::duration<double>(stop - start).count() << "s \n";

  printf("Getting all triangles\n");
  TriangleList allTris = storeAllTriangles(&cg_relabel,tri_info.total);
  stop = high_resolution_clock::now();
  cout << "time: " << chrono::duration<double>(stop - start).count() << "s \n";

  printf("Compute nuclei\n");
  EdgeIdx nClusters = 0;
  Cluster* clusters = nuclei(&cg_relabel,&allTris, true, nClusters); 
  stop = high_resolution_clock::now();
  cout << "time: " << chrono::duration<double>(stop - start).count() << "s \n";

  printf("Writting clusters\n");
  ofstream myfile (argv[2]);
  if (myfile.is_open())
  {
    for (EdgeIdx idxj=0; idxj<nClusters; idxj++) {
      Cluster cluster = clusters[idxj];
      myfile << "Cluster id: "<< cluster.id << ", C: " << cluster.c << ", n_vertices: " << cluster.n_vertices << ", n_edges: " << cluster.n_edges << ", parent: "<< cluster.parent <<"\n";
      std::set<VertexIdx>::iterator itr;
      for (itr = cluster.vertices.begin(); itr != cluster.vertices.end(); itr++){
          myfile << (inverse[*itr]) << " ";
      }
      myfile << "\n";
    }
    myfile.close();
  } else {
    printf("error writting output\n");
  }


}

