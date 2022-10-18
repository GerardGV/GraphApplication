#include "pch.h"
#include "Graph.h"
#include <vector>
#include <queue>
using namespace std;

// =============================================================================
// SpanningTreePrim ============================================================
// =============================================================================


struct funct
{
	bool operator()(CEdge* left, CEdge* right)
	{
		return left->m_Length > right->m_Length;
	}
};

/* 
bool compEdgesLenght(CEdge left, CEdge right)
{
	return left.m_Length < right.m_Length;
}
*/


//cola de prioridad y resultado

CSpanningTree SpanningTreePrim(CGraph& graph) {

	struct comparator {
		bool operator()(CEdge* pE1, CEdge* pE2) {
			return pE1->m_Length > pE2->m_Length;
		}
	};
	priority_queue<CEdge*, std::vector<CEdge*>, comparator> queue;

	CSpanningTree tree(&graph);

	//solo tenemos añadir las aristas del vertice que visitamos i solo cogemos las exploran nuevos vertices
	for (CEdge& e : graph.m_Edges) queue.push(&e);

	int label = 1;
	//for (CVertex& v : graph.m_Vertices) v.m_KruskalLabel = label++;

	//we get the edge with less lenght
	CEdge* pE = queue.top();

	//the origin has been visited
	pE->m_pOrigin->m_PrimInTree = true;
	int visited = 1;
	while (visited < graph.m_Vertices.size()){

		//reverse has not been used
		if (pE->m_reverseUsed != true)
		{
			//origin has to be visited, destination not visited
			if (pE->m_pOrigin->m_PrimInTree == true && pE->m_pDestination->m_PrimInTree != true) {
				tree.m_Edges.push_back(pE);
				//destination now its visited
				pE->m_pDestination->m_PrimInTree = true;
				visited++;
				//we indicated that we are not going to use the reverse
				pE->m_pReverseEdge->m_reverseUsed = true;
			}
		}


		//we control accesing to the queue when it is empty
		if (!queue.empty())
		{
			pE = queue.top();

			queue.pop();
		}
		else
			break;
		
	}
	return tree;
} 



 		/*
		//we start exploring the first vertex
		CVertex* vertex = &*(graph.m_Vertices.begin());

		//we add the edge to the graph because it is now visited
		solucio.m_pGraph->m_Vertices.push_back(*vertex);

		//list of reverse edges we don't want
		list<string> bannedEdges;

		//list of reverse edges we don't want
		while (solucio.m_pGraph->m_Vertices.size() != graph.m_Vertices.size())
		{
			//we select the first edge like the best
			//CEdge* bestEdge = (*(*vertex).m_Edges.begin());
			


			//fill pqueue_less
			if (vertex->m_visited != true)
			{
				for (auto it = vertex->m_Edges.begin(); it != vertex->m_Edges.end(); it++)
				{
					CEdge* tmp = (*it);
					//we look if the edge is a reverse
					//size_t found = tmp->m_Name.find("$");

					cout << " trying to ADD EDGE: " << (*(*it)).m_Name << " , LENGHT: " << (*(*it)).m_Length << endl;
					cout << " DESTINATION: " << (*(*it)).m_pDestination->m_Name << endl;
					//if (!solucio.m_pGraph->MemberP(tmp->m_pReverseEdge))
					//auto found = find(bannedEdges.begin(), bannedEdges.end(), tmp->m_pReverseEdge->m_Name.c_str());
					if (tmp->m_pDestination->m_visited==false)
					{
						cout << "ADDING EDGE: " << (*(*it)).m_Name << " , LENGHT: " << (*(*it)).m_Length << endl;

						pCandidates.push(tmp);
					}
					
					
											
				}
				vertex->m_visited = true;
			}
			

			//cout <<endl << "BEST EDGE: " << (*pCandidates).top()->m_Name << endl;
			//we add to the tree the best edge that is on top of the binary heap to the graph
			solucio.m_pGraph->m_Edges.push_back(*pCandidates.top());
			CEdge* tmp = pCandidates.top();
			solucio.Add(tmp);

			cout << "bestEdge: " << tmp->m_Name << ", LENGHT " << tmp->m_Length << endl;
			//we save the next vertex if it has not been explored
			if (pCandidates.top()->m_pDestination->m_visited != true)
			{
				vertex = pCandidates.top()->m_pDestination;
				//we add the edge to the graph because it is now visited
				solucio.m_pGraph->m_Vertices.push_back(*vertex);
			}

			//we take out the best edge because we have already added to the graph
			pCandidates.pop();


		}


		//loop to get the best edges
		
		return solucio;
		//return NULL;
	}

	return NULL;
}
*/
/*
CSpanningTree SpanningTreePrim(CGraph& graph)
{
	//if the graph has edges, we wiill able to apply the algorithm
	if (!graph.m_Edges.empty())
	{

		vector<CEdge*> candidates;
		
		CSpanningTree solucio(new CGraph(graph));
		solucio.m_pGraph->Clear();

		//we get all the posible edges
		for (CEdge& e : graph.m_Edges)
		{
			pCandidates.push(&e);
		}

		//origin vertex it's visited
		CEdge* tmp = pCandidates.top();
		tmp->m_pOrigin->m_visited = true;

		//we add the vertex to the solution
		solucio.m_pGraph->m_Vertices.push_back(*tmp->m_pOrigin);

		//now, we will get the BEST amoung ALL THE EDGES until we visit all the vertex
		while (solucio.m_pGraph->m_Vertices.size() != graph.m_Vertices.size())
		{
			cout << " trying to ADD EDGE: " << (*tmp).m_Name << " , LENGHT: " << (*tmp).m_Length << endl;
			cout << " DESTINATION: " << (*tmp).m_pDestination->m_Name << endl;
			if (tmp->m_eliminated != true)
			{

				if (tmp->m_pDestination->m_visited == false)
				{
					cout << "ADDING EDGE: " << (*tmp).m_Name << " , LENGHT: " << (*tmp).m_Length << endl;

					//we add the best edge
					solucio.m_pGraph->m_Edges.push_back(*tmp);
					solucio.Add(tmp);

					//we won't be able to use the reverse vertex
					tmp->m_pReverseEdge->m_eliminated = true;

					//we put the destination vertex visited
					tmp->m_pDestination->m_visited = true;

					//we add the vertex to the solution
					solucio.m_pGraph->m_Vertices.push_back(*tmp->m_pDestination);
				}
			}


			//we take out the best edge because we have already added to the graph
			pCandidates.pop();

			//we get the next best edge
			tmp = pCandidates.top();

		}//while

		/*
		//we start exploring the first vertex
		CVertex* vertex = &*(graph.m_Vertices.begin());

		//we add the edge to the graph because it is now visited
		solucio.m_pGraph->m_Vertices.push_back(*vertex);

		//list of reverse edges we don't want
		list<string> bannedEdges;

		//list of reverse edges we don't want
		while (solucio.m_pGraph->m_Vertices.size() != graph.m_Vertices.size())
		{
			//we select the first edge like the best
			//CEdge* bestEdge = (*(*vertex).m_Edges.begin());



			//fill pqueue_less
			if (vertex->m_visited != true)
			{
				for (auto it = vertex->m_Edges.begin(); it != vertex->m_Edges.end(); it++)
				{
					CEdge* tmp = (*it);
					//we look if the edge is a reverse
					//size_t found = tmp->m_Name.find("$");

					cout << " trying to ADD EDGE: " << (*(*it)).m_Name << " , LENGHT: " << (*(*it)).m_Length << endl;
					cout << " DESTINATION: " << (*(*it)).m_pDestination->m_Name << endl;
					//if (!solucio.m_pGraph->MemberP(tmp->m_pReverseEdge))
					//auto found = find(bannedEdges.begin(), bannedEdges.end(), tmp->m_pReverseEdge->m_Name.c_str());
					if (tmp->m_pDestination->m_visited==false)
					{
						cout << "ADDING EDGE: " << (*(*it)).m_Name << " , LENGHT: " << (*(*it)).m_Length << endl;

						pCandidates.push(tmp);
					}



				}
				vertex->m_visited = true;
			}


			//cout <<endl << "BEST EDGE: " << (*pCandidates).top()->m_Name << endl;
			//we add to the tree the best edge that is on top of the binary heap to the graph
			solucio.m_pGraph->m_Edges.push_back(*pCandidates.top());
			CEdge* tmp = pCandidates.top();
			solucio.Add(tmp);

			cout << "bestEdge: " << tmp->m_Name << ", LENGHT " << tmp->m_Length << endl;
			//we save the next vertex if it has not been explored
			if (pCandidates.top()->m_pDestination->m_visited != true)
			{
				vertex = pCandidates.top()->m_pDestination;
				//we add the edge to the graph because it is now visited
				solucio.m_pGraph->m_Vertices.push_back(*vertex);
			}

			//we take out the best edge because we have already added to the graph
			pCandidates.pop();


		}


		//loop to get the best edges
		
		return solucio;
		//return NULL;
	}

	return NULL;
}
*/