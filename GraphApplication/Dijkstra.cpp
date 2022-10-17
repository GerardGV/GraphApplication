#include "pch.h"
#include "Graph.h"
#include <queue>

//#include <unordered_map>

// =============================================================================
// Dijkstra ====================================================================
// =============================================================================

void Dijkstra(CGraph& graph, CVertex* pStart)
{
	pStart->m_DijkstraDistance = 0;
	//vector<CVertex*> candidates;
	pStart->m_DijkstraVisit = true;
	//visited.push_back(*pStart);
	CVertex* pActual = pStart;

	//list<CVertex*> candidates;
	//candidates.push_back(pActual);
	//cout << "Il vertice di partenza è: " << pStart->m_Name << " e la sua distanza: " << pStart->m_DijkstraDistance << endl;
	//auto it = candidates.begin();
	//int indexMin = 0;
	int visitats = 1;

	while (visitats != graph.m_Vertices.size()) {

		//candidates.erase(next(candidates.begin(), indexMin));

		//we check if we find a best dijkstra distance to the neighbours vertex through actual vertex
		for (CEdge* e : pActual->m_Edges) {

			/*
			if(e->m_pDestination->m_DijkstraVisit==false)
				candidates.push_back(e->m_pDestination);
			*/
			//cout << "La distanza del vertice: " << e->m_pDestination->m_Name << " è: " << e->m_pDestination->m_DijkstraDistance << endl;

			if (e->m_Length + pActual->m_DijkstraDistance < e->m_pDestination->m_DijkstraDistance) {

				(*e).m_pDestination->m_DijkstraDistance = (*e).m_Length + pActual->m_DijkstraDistance;
				
				//cout << "Il vertice attuale è: " << pStart->m_Name << " e la sua distanza: " << dist << endl;
			}
		}

		//me mark the actual vertex as visited
		pActual->m_DijkstraVisit = true;
		visitats++;
		


		//we look for the next vertex, the nearest, to explore it
		if (visitats != graph.m_Vertices.size())
		{
			pActual =&graph.m_Vertices.front();
			double min = pActual->m_DijkstraDistance;
			auto it = graph.m_Vertices.begin();
			it++;
			for (;it != graph.m_Vertices.end(); it++)
			{
				if (min > (*it).m_DijkstraDistance && (*it).m_DijkstraVisit == false)
				{
					min = (*it).m_DijkstraDistance;
					pActual = &(*it);
				}
			}
		}
		


		/*
		if (!candidates.empty())
		{
			pActual = candidates[0];
			double min = candidates[0]->m_DijkstraDistance;

			for (int index = 1; index < candidates.size(); index++)
			{
				if (min > candidates[index]->m_DijkstraDistance && candidates[index]->m_DijkstraVisit == false)
				{
					min = candidates[index]->m_DijkstraDistance;
					pActual = candidates[index];
					indexMin = index;
				}
			}
		}
		*/

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
