
// Shortest Path Calculation using Dijkstra,
// BellmanFord and FloydWarshall Algorithm.
// The impemetation is extended from CLRS book,
// websites such as wikipedia and geeksforgeeks.

#include <iostream>
#include <list>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>
#include <vector>
#include <climits>
#include <fstream>
#include<limits.h>
#include <iomanip>

#define INF INT_MAX
using namespace std;

int V_max = 100;
int weight_max = 10;
int vnum[] = {50,100,500,1000};
int vnum_size = 4;
std::vector<int> myvector;

// Graph class represents a directed graph
// using adjacency list representation
class Graph
{
public:
    int V;    // No. of vertices
    int E;    // No. of edges
    // Pointer to an array containing
    // adjacency lists
    list<pair<int, int> > *adj;

    Graph(int V);   // Constructor

    // function to add an edge to graph
    void addEdge(int v, int w, int weight);
    void printEdge();
    int numVertex();
};

// This function returns number of vertex in
// Graph
int Graph::numVertex()
{
    return this->V;
}

// This is constructor for Graph
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<pair<int, int> >[V];
}

// This function add an edge along with the
// weight of the edge.
void Graph::addEdge(int v, int w, int weight)
{
    cout << "construct edge: (" << v << ", " << w << ", " << weight << ")" << endl;
    adj[v].push_back(make_pair(w, weight)); // Add w to v’s list.
}

// Helper function to print an edge
void Graph::printEdge()
{
    for (int k = 0; k < V; k++){
        for (list<pair<int, int> >::const_iterator i = adj[k].begin(); i != adj[k].end(); i++){
            cout << k << "," << i->first << "," << i->second;
            cout << endl;
        }
    }
}

// generate sample of edges randomly
vector<int> randSample(int E, int E_max)
{
    vector<int> edgeVector;
    for (int i = 0; i < E_max; i++){
        edgeVector.push_back(i);
    }
    random_shuffle(edgeVector.begin(), edgeVector.end());
    vector<int> randVector(edgeVector.begin(), edgeVector.begin() + E);
    return randVector;
}

// generate graph of connecting random edges
Graph genGraph(int V, int E, int E_max)
{
    srand(time(NULL));
    cout << "V = " << V << endl;
    cout << "E = " << E << endl;
    //generate graph
    Graph g(V);
    g.E = E;
    //random edges in the vector
    vector<int> Edges = randSample(E, E_max);
    //construct edges
    for (int i = 0; i < Edges.size(); i++){
        int v = Edges[i]/(V-1);
        int w = Edges[i]%(V-1);
        if (w >= v) w++;
        int weight = rand() % weight_max + 1;
        g.addEdge(v, w, weight);
    }
    return g;
}

// Helper function for Dijstra function
int minDistance(int dist[], bool sptSet[], int V)
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
    {
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

// To print the final solution
int printSolution(int dist[], int V)
{
    cout<<"\n Vertex   Distance from Source"<<endl;
    for (int i = 0; i < V; i++)
        cout<<i<<"\t "<<dist[i]<<endl;
        //printf("%d \t\t %d\n", i, dist[i]);
}

// The main function that finds shortest distances from src to
// all other vertices using Dijkstra's algorithm.
float dijkstra(Graph g, int src)
{
    clock_t t1, t2;
    t1 = clock();
    int V = g.numVertex();
    int dist[V];     // The output array.  dist[i] will hold the shortest
                     // distance from src to i

    bool sptSet[V]; // sptSet[i] will be true if vertex i is included in shortest
                    // path tree or shortest distance from src to i is finalized

    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < V-1; count++)
    {
        // Pick the minimum distance vertex from the set of vertices not
        // yet processed. u is always equal to src in first iteration.
        int u = minDistance(dist, sptSet, V);

        // Mark the picked vertex as processed
        sptSet[u] = true;

        if (dist[u] == INT_MAX) break;
        // Update dist value of the adjacent vertices of the picked vertex.
        for (list<pair<int, int> >::const_iterator i = g.adj[u].begin(); i != g.adj[u].end(); i++)
        {
            int v = i->first;
            // Update dist[v] only if is not in sptSet, there is an edge from
            // u to v, and total weight of path from src to  v through u is
            // smaller than current value of dist[v]
            if (!sptSet[v] && dist[u]+ i->second < dist[v])
                dist[v] = dist[u] + i->second;
        }
    }
    t2 = clock();

    // print the constructed distance array
    printSolution(dist, V);
    for (int i=0;i<V;i++){
        myvector.push_back(dist[i]);
    }

    float diff = ((float)t2-(float)t1);
    return diff / CLOCKS_PER_SEC;
}

// function for dijkstra all pair shorted path
float DijstraAPSP(Graph g)
{
  clock_t t1, t2;
  t1 = clock();
  int V = g.numVertex();
  int graph1[V][V];
  for(int i = 0; i<V;i++){
          dijkstra(g,i);
  }
  t2 = clock();
  float diff = ((float)t2-(float)t1);
  /*for (int i =0;i<int(myvector.size());i++){
          printf("my vector has  %d \n", myvector[i]);
  } */
  int k=0;
  for(int i =0;i<V;i++){
    for (int j = 0; j<V;j++)
    {
        graph1[i][j] = myvector[k];
        k= k+1;
          //printf("my vector has is %d \n", myvector[i*j]);
      }
  }
  //Print the result in matrix format
  for (int i = 0; i < V; i++)
  {
    for (int j = 0; j < V; j++)
    {
      if (graph1[i][j] == INF)
            printf("%7s", "INF");
      else
            printf ("%7d", graph1[i][j]);
    }
      printf("\n");
  }
  return diff / CLOCKS_PER_SEC;
}
// The main function that finds shortest distances from src to
// all other vertices using Bellman-Ford algorithm.  The function
// also detects negative weight cycle
float BellmanFord(Graph g, int src)
{
    clock_t t1, t2;
    t1 = clock();
    int V = g.numVertex();
    int E = g.E;
    int dist[V];

    // Step 1: Initialize distances from src to all other vertices
    // as INFINITE
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX;
    dist[src] = 0;

    // Step 2: Relax all edges |V| - 1 times. A simple shortest
    // path from src to any other vertex can have at-most |V| - 1
    // edges
    for (int i = 1; i <= V-1; i++)
    {
        int change = 0;
        for (int u = 0; u < V; u++)
        {
            if (dist[u] == INT_MAX) continue;
            for (list<pair<int, int> >::const_iterator i = g.adj[u].begin(); i != g.adj[u].end(); i++)
            {
                int v = i->first;
                int weight = i->second;
                if (dist[u] + weight < dist[v]){
                    dist[v] = dist[u] + weight;
                    change = 1;
                }
            }
        }
        if (!change) break;
    }

    // Step 3: check for negative-weight cycles.  The above step
    // guarantees shortest distances if graph doesn't contain
    // negative weight cycle.  If we get a shorter path, then there
    // is a cycle.
    /*for (int i = 0; i < E; i++)
    {
        int u = graph->edge[i].src;
        int v = graph->edge[i].dest;
        int weight = graph->edge[i].weight;
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
            printf("Graph contains negative weight cycle");
    }*/

    t2 = clock();
    printSolution(dist, V);
    float diff = ((float)t2-(float)t1);
    return diff / CLOCKS_PER_SEC;
}

// The main function that finds shortest distances from each vertex to
// all other vertices using Floyd Warshall's algorithm.
float floydWarshall (Graph g)
{
    /* dist[][] will be the output matrix that will finally have the shortest
      distances between every pair of vertices */
    clock_t t1, t2;
    t1 = clock();
    int V = g.numVertex();
    int dist1[V][V];
      for(int i =0;i<V;i++){
        for (int j = 0; j<V;j++){
          if(i==j){
            dist1[i][j] = 0;
        }
        else{
        dist1[i][j] = INF;
        }
      }
    }
    /* Initialize the solution matrix same as input graph matrix. Or
       we can say the initial values of shortest distances are based
       on shortest paths considering no intermediate vertex. */
    for(int i=0;i<V;i++){
        for (list<pair<int, int> >::const_iterator k = g.adj[i].begin(); k != g.adj[i].end(); k++){
              //cout<<"i = "<<i<<"\t" << "k value = "<<k->first <<"\t" << "weight= "<< k->second<<"\t"<<'\n' ;
              dist1[i][k->first]= k->second;
            }
        }
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of an iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of an iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */
   for (int k = 0; k < V; k++)
    {
        // Pick all vertices as source one by one
        for (int i = 0; i < V; i++)
        {
            // Pick all vertices as destination for the
            // above picked source
            for (int j = 0; j < V; j++)
            {
                // If vertex k is on the shortest path from
                // i to j, then update the value of dist[i][j]
                if ( dist1[i][k] != INF
                    && dist1[k][j] != INF
                    && (dist1[i][k] + dist1[k][j] < dist1[i][j]))
                    {
                    dist1[i][j] = dist1[i][k] + dist1[k][j];
                    }
            }
        }
    }
        t2 = clock();
    // Print the shortest distance matrix
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
        {
            if (dist1[i][j] == INF)
                printf("%7s", "INF");
            else
                printf ("%7d", dist1[i][j]);
        }
        printf("\n");
    }


    float diff = ((float)t2-(float)t1);
    return diff / CLOCKS_PER_SEC;
  //  printSolution1(dist1);
}

int main()
{
    float seconds, seconds2,seconds3,sec;
    ofstream outputFile;
    outputFile.open("output.csv");
    outputFile << "V,E,Dijkstra,Bellman-Ford,FloydWarshall,DijstraAPSP"<<endl;
    for (int i = 0; i < vnum_size; i++){
        int V = vnum[i];
        int E_max = V * V - V;
      for (int E = V+10; E <= 50*V; E+=5*V){
            outputFile << V << "," << E << ",";
            Graph g = genGraph(V, E, E_max);
            //g.printEdge();
            seconds = dijkstra(g, 0);
            //std::cout << std::setprecision(9) << seconds << '\n';
            //cout<<"Dijkstra's Run time : "<<seconds*1000<<endl;
            seconds2 = BellmanFord(g, 0);
            //cout<<"BellmanFord Run time : "<<seconds2*1000<<endl;
            seconds3= floydWarshall(g);
            //cout<<"floydWarshall Run time : "<<seconds3*1000<<endl;
            sec = DijstraAPSP(g);
            outputFile << seconds << "," << seconds2 << ","<< seconds3 << ","<< sec << endl;
       }
    }
    outputFile.close();
    return 0;
}
