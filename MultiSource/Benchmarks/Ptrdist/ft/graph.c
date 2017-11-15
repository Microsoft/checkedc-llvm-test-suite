/*
 * graph.c
 *
 * The author of this software is Alain K\"{a}gi.
 *
 * Copyright (c) 1993 by Alain K\"{a}gi and the University of Wisconsin
 * Board of Trustees.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 *
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR THE UNIVERSITY OF
 * WISCONSIN MAKE ANY REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING
 * THE MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR
 * PURPOSE.
 *
 * ------------------------------------------------------------------------
 *
 * The graph etc.
 *
 * ------------------------------------------------------------------------
 *
 * $Id$
 *
 */

#include <assert.h>
#include <stdio_checked.h>
#include <stdlib_checked.h>
#include "graph.h"

#pragma BOUNDS_CHECKED ON

#define TRUE 1
#define FALSE 0

#ifdef __MINGW32__
#define random() rand()
#endif

#define GET_WEIGHT	(random() + 1) % MAX_WEIGHT

/*
 * Local variables.
 */
static int generatedEdges;

/*
 * Local functions.
 */
_Ptr<Vertices>  GenTree(int vertex);
_Ptr<Vertices>  AddEdges(_Ptr<Vertices>  graph, int nVertex, int nEdge);
_Ptr<Vertices>  PickVertex(_Ptr<Vertices>  graph, int whichVertex);
void      Connect(_Ptr<Vertices>  vertex1, _Ptr<Vertices>  vertex2);
int       Duplicate(_Ptr<Vertices>  vertex1, _Ptr<Vertices>  vertex2);
_Ptr<Vertices>  NewVertex(void);
_Ptr<Edges>  NewEdge(void);
void      PrintNeighbors(_Ptr<Vertices>  vertex);

/*
 * Local variables.
 */
static id = 1;

/*
 * Generates random connected undirected graphs.  Unfortunately this current
 * implementation does not generate the graphs with a uniform distribution.
 *
 * Apparently a good reference is Tinhofer G., ,
 * C. Hanser, Verlag, M\"{u}nchen 1980.
 */
_Ptr<Vertices> 
GenGraph(int nVertex, int nEdge)
{
  _Ptr<Vertices>  graph = 0;

  assert(nEdge + 1 >= nVertex);
  assert(nEdge <= nVertex * (nVertex - 1) / 2);

  generatedEdges = 0;

  graph = GenTree(nVertex);
  graph = AddEdges(graph, nVertex, nEdge - nVertex + 1);
  return(graph);
}

_Ptr<Vertices> 
GenTree(int nVertex)
{
  int       i;
  int       weight;
  _Ptr<Vertices>  vertex = 0;
  _Ptr<Vertices>  graph = 0;
  _Ptr<Edges>  edge = 0;

  graph = NewVertex();
  NEXT_VERTEX(graph) = graph;

  for(i = 1; i < nVertex; i++)
  {
    vertex = NewVertex();
    edge = NewEdge();

    /*
     * The newly created vertex has one edge ...
     */
    EDGES(vertex) = edge;

    /*
     * ... which is connected to the graph so far generated.  The connection
     * point in the graph is picked at random.
     */
    VERTEX(edge) = PickVertex(graph, random() % i);
    weight = GET_WEIGHT;
    WEIGHT(edge) = weight;
    SOURCE(edge) = vertex;

    /*
     * Link the new vertex into the graph.
     */
    NEXT_VERTEX(vertex) = NEXT_VERTEX(graph);
    NEXT_VERTEX(graph) = vertex;

    /*
     * Add an edge to the vertex randomly picked as the connection point.
     */
    edge = NewEdge();
    WEIGHT(edge) = weight;
    SOURCE(edge) = VERTEX(EDGES(vertex));
    VERTEX(edge) = vertex;
    NEXT_EDGE(edge) = EDGES(VERTEX(EDGES(vertex)));
    EDGES(VERTEX(EDGES(vertex))) = edge;
   }

  return(graph);
}

_Ptr<Vertices> 
AddEdges(_Ptr<Vertices>  graph, int nVertex, int nEdge)
{
  int       i;
  _Ptr<Vertices>  vertex1 = 0;
  _Ptr<Vertices>  vertex2 = 0;

  assert(graph != NULL_VERTEX);
  assert(nEdge >= 0);

  for(i = 0; i < nEdge; i++)
  {
    do
    {
      vertex1 = PickVertex(graph, random() % nVertex);
      vertex2 = PickVertex(NEXT_VERTEX(vertex1), (random() % nVertex) - 1);
      assert(vertex1 != vertex2);
    }
    while(Duplicate(vertex1, vertex2));

    Connect(vertex1, vertex2);
  }
  
  return(graph);
}

_Ptr<Vertices> 
PickVertex(_Ptr<Vertices>  graph, int whichVertex)
{
  int       i;

  for(i = 0; i < whichVertex; i++)
  {
    graph = NEXT_VERTEX(graph);
  }

  return(graph);
}

void
Connect(_Ptr<Vertices>  vertex1, _Ptr<Vertices>  vertex2)
{
  int    weight;
  _Ptr<Edges>  edge = 0;

  weight = GET_WEIGHT;

  edge = NewEdge();
  WEIGHT(edge) = weight;
  SOURCE(edge) = vertex1;
  VERTEX(edge) = vertex2;
  NEXT_EDGE(edge) = EDGES(vertex1);
  EDGES(vertex1) = edge;
  
  edge = NewEdge();
  WEIGHT(edge) = weight;
  SOURCE(edge) = vertex2;
  VERTEX(edge) = vertex1;
  NEXT_EDGE(edge) = EDGES(vertex2);
  EDGES(vertex2) = edge;
}

int
Duplicate(_Ptr<Vertices>  vertex1, _Ptr<Vertices>  vertex2)
{
  _Ptr<Edges>  edge = 0;

  edge = EDGES(vertex1);

  while(edge != NULL_EDGE)
  {
    if(VERTEX(edge) == vertex2)
    {
      return(TRUE);
    }

    edge = NEXT_EDGE(edge);
  }

  return(FALSE);
}

_Ptr<Vertices> 
NewVertex(void)
{
  _Ptr<Vertices>  vertex = 0;

  vertex = calloc(1, sizeof(Vertices));

  if(vertex == NULL)
  {
    fprintf(stderr, "Could not malloc\n");
    exit(1);
  }

  ID(vertex) = id++;
  EDGES(vertex) = NULL;
  NEXT_VERTEX(vertex) = NULL;

  return(vertex);
}

_Ptr<Edges> 
NewEdge(void)
{
  _Ptr<Edges>  edge = 0;

  edge = calloc(1, sizeof(Edges));

  if(edge == NULL)
  {
    fprintf(stderr, "Could not malloc\n");
    exit(1);
  }

  WEIGHT(edge) = 0;
  VERTEX(edge) = NULL;
  NEXT_EDGE(edge) = NULL;
  
  return(edge);
}

void
PrintGraph(_Ptr<Vertices>  graph)
{
  _Ptr<Vertices>  vertex = 0;

  assert(graph != NULL);

  vertex = graph;
  do
  {
    printf("Vertex %d is connected to:", ID(vertex));
    PrintNeighbors(vertex);
    printf("\n");
    vertex = NEXT_VERTEX(vertex);
  }
  while(vertex != graph);
}

void
PrintNeighbors(_Ptr<Vertices>  vertex)
{
  _Ptr<Edges>  edge = 0;

  edge = EDGES(vertex);
  while(edge != NULL)
  {
    printf(" %d(%d)[%d]", ID(VERTEX(edge)), WEIGHT(edge), ID(SOURCE(edge)));
    edge = NEXT_EDGE(edge);
  }
}

