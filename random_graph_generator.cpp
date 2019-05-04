#include<iostream>
#include<stdlib.h>

using namespace std;

// A function to generate random graph.
void GenerateRandGraphs(int NOE, int NOV)
{
	int i, j, edge[NOE][2], count;
	i = 0;
	// Build a connection between two random vertex.
	while(i < NOE)
	{
		edge[i][0] = rand()%NOV+1;
		edge[i][1] = rand()%NOV+1;

		if(edge[i][0] == edge[i][1])
			continue;
		else
		{
			for(j = 0; j < i; j++)
			{
				if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) || (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0]))
					i--;
			}
		}
		i++;
	}

	// Print the random graph.
	cout<<"\nThe generated random random graph is: ";
	for(i = 0; i < NOV; i++)
	{
		count = 0;
		cout<<"\n\t"<<i+1<<"-> { ";
		for(j = 0; j < NOE; j++)
		{
			if(edge[j][0] == i+1)
			{
				cout<<edge[j][1]<<"   ";
				count++;
			}
			else if(edge[j][1] == i+1)
			{
				cout<<edge[j][0]<<"   ";
				count++;
			}
			else if(j == NOE-1 && count == 0)
				cout<<"Isolated Vertex!";
		}
		cout<<" }";
	}
}

int main()
{
	int n, i, e, v;

	cout<<"Random graph generation: ";

	// Randomly assign vertex and edges of the graph.
	v = 11+rand()%10;
	cout<<"\nThe graph has "<<v<<" vertexes.";
	e = rand()%((v*(v-1))/2);
	cout<<"\nThe graph has "<<e<<" edges.";

	// A function to generate a random undirected graph with e edges and v vertexes.
	GenerateRandGraphs(e, v);
}
