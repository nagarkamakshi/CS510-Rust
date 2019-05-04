#include <bits/stdc++.h>


struct Edge
{
	int src, dest, weight;
};

struct Graph
{
	int V, E;
	// graph is represented as an array of edges.
	struct Edge* edge;
};

// Creates a graph with V vertices and E edges
struct Graph* createGraph(int V, int E)
{
	struct Graph* graph = new Graph;
	graph->V = V;
	graph->E = E;
	graph->edge = new Edge[E];
	return graph;
}

int main(){

  int v=5, e=10;
  int graph[v][v];
  for (int i=0;i<v;i++){
    for(int j=i;j<v;j++){
      if(i==j){
        graph[i][j]=0;
      }else{
      graph[i][j]= rand()%((v*(v-1))/2);
       graph[j][i] = graph[i][j];
     }
    }
  }
  for (int i=0;i<v;i++){
    for(int j=0;j<v;j++){
        printf("%d   ", graph[i][j]);
      }
      printf("\n");
    }


  return 0;
}
