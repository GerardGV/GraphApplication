#include "pch.h"
#include "Graph.h"
#include <queue>

//#include <unordered_map>

// =============================================================================
// Dijkstra ====================================================================
// =============================================================================

void Dijkstra(CGraph& graph, CVertex* pStart)
{
	//we mark Start vertex as visited and we set the distance to 0
	pStart->m_DijkstraDistance = 0;
	pStart->m_DijkstraVisit = true;
	CVertex* pActual = pStart;


	int visitats = 1;


	//we calculate the dijkstra distance of the neighbours through the conection with the actual Vertex
	while (visitats != graph.m_Vertices.size()) {


		for (CEdge* e : pActual->m_Edges) {


			if (e->m_Length + pActual->m_DijkstraDistance < e->m_pDestination->m_DijkstraDistance) {

				(*e).m_pDestination->m_DijkstraDistance = (*e).m_Length + pActual->m_DijkstraDistance;
				
			}
		}

		//me mark the actual vertex as visited
		pActual->m_DijkstraVisit = true;
		visitats++;
		


		//we look for the next vertex, the nearest, to explore it
		if (visitats != graph.m_Vertices.size())
		{
			
			double min = numeric_limits<double>::max();
			
			for (auto it =graph.m_Vertices.begin(); it != graph.m_Vertices.end(); it++)
			{
				if (min > (*it).m_DijkstraDistance && (*it).m_DijkstraVisit == false)
				{
					min = (*it).m_DijkstraDistance;
					pActual = &(*it);
				}
			}
		}
		



	}


}

// =============================================================================
// DijkstraQueue ===============================================================
// =============================================================================

void DijkstraQueue(CGraph& graph, CVertex *pStart)
{
	
	//we define the condicion to sorted the priority queue
	struct comparator {
		bool operator()(CVertex* pE1, CVertex* pE2) {
			return pE1->m_DijkstraDistance > pE2->m_DijkstraDistance;
		}
	};

	priority_queue<CVertex*, std::vector<CVertex*>, comparator> nextNearVertex;

	for (CVertex& it : graph.m_Vertices)
	{
		it.m_DijkstraDistance = numeric_limits<double>::max();
		it.m_DijkstraVisit = false;
	}

	//we push the start node and we put its distane to 0
	nextNearVertex.push(pStart) ;
	pStart->m_DijkstraDistance = 0;
	

	//we iterate until visit all the vertex
	while (!nextNearVertex.empty())
	{
		//we go to the nearest vertex and we count +1 vertex visited
		pStart = nextNearVertex.top();
		//cout << "AUX Name: " << pStart->m_Name << endl;
		nextNearVertex.pop();

		//we only explore the vertex if it has not been wxplore before
		if (pStart->m_DijkstraVisit == false)
		{		
			//we calculate news distances Vertex from the actual vertex, pStart, to its neighbours, and we push the neighbours as candidates to be the next vertex that we are going to explore
			for (auto e = pStart->m_Edges.begin(); e != pStart->m_Edges.end();e++)
			{
				//we update the dijkstra lengh of destination vertex with the lenght until the last vertex + lenght of the edge if it is less than the actual lenght
				if ((*(*e)).m_pDestination->m_DijkstraDistance > (*(*e)).m_pOrigin->m_DijkstraDistance + (*(*e)).m_Length)
				{
					//we update the distance 
					(*(*e)).m_pDestination->m_DijkstraDistance = (*(*e)).m_pOrigin->m_DijkstraDistance + (*(*e)).m_Length;

					//we update the origin edge
					(*(*e)).m_pDestination->m_pDijkstraPrevious = (*e);
				}

				//we push neighbours vertexs to the candidates lis to be the next vertex that we are going to explore
				nextNearVertex.push((*(*e)).m_pDestination);


			}

		}

		//we mark the actual Vertex as visited
		pStart->m_DijkstraVisit = true;

		
	}

	
}
