#include "pch.h"
#include "Graph.h"
#include <set>


// =============================================================================
// SalesmanTrackBacktracking ===================================================
// =============================================================================


//we create these global variables that we don't modify in each call of the recursive function
CTrack CamiMesCurt2(NULL);
double LongitudCamiMesCurt2;

CVertex* pDesti2;
CTrack CamiActual2(NULL);
double LongitudCamiActual2;


list<CVertex*> VISITES2;

//we have created global variable to take advantage of the local memory in the functions, saving time accesing to the values
struct NodeCami2 {
	CEdge* m_pEdge;
	NodeCami2* m_pAnterior;
};


CTrack SalesmanTrackBacktrackingRec(NodeCami2* pAnterior, CVertex* pActual)
{
	if (pActual == pDesti2) {
		if (LongitudCamiActual2 < LongitudCamiMesCurt2) {
			
			bool correct = true;
			
			//we check that we passed all the vertex that we have to visit
			for (auto it = VISITES2.begin(); it != VISITES2.end(); it++)
			{
				//if we didn't pass through one vertex of visits, the way is no a possible answer
				if ((*it)->m_count == 0)
				{
					correct = false;
					break;
				}
			}
			
			//if the way is a possible answer
			if (correct)
			{
				//if the actual way has less cost, we will take it as shortest way
				if (LongitudCamiMesCurt2 > LongitudCamiActual2)
				{
					CamiMesCurt2.Clear();
					while (pAnterior) {
						CamiMesCurt2.m_Edges.push_front(pAnterior->m_pEdge);
						pAnterior = pAnterior->m_pAnterior;
					}
					LongitudCamiMesCurt2 = LongitudCamiActual2;
				}
				

			}
			

		}
	}
	else 
		//we are not going to explore a vertex if the current lengh is already longer that the shortest way
		if (LongitudCamiActual2 < LongitudCamiMesCurt2) {
		
			//pActual->m_JaHePassat = true;
		//we pass +1 time trhough the vertex
		pActual->m_count++;

		//node to  save actual way
		NodeCami2 node;
		node.m_pAnterior = pAnterior;
		for (CEdge* pE : pActual->m_Edges) {
			//if (!pE->m_pDestination->m_JaHePassat) {
				
			node.m_pEdge = pE;
			LongitudCamiActual2 += pE->m_Length;
			NodeCami2* tmp = &node;
			SalesmanTrackBacktrackingRec(&node, pE->m_pDestination);
			LongitudCamiActual2 -= pE->m_Length;
		//}
		}

		//we get to the previous state
		pActual->m_JaHePassat = false;
		pActual->m_count--;
	}

	return CamiMesCurt2;
}


CTrack SalesmanTrackBacktracking(CGraph &graph, CVisits &visits)
{
	if (!graph.m_Edges.empty())
	{
		pDesti2 = visits.m_Vertices.back();
		//we take out the inicial node and the destination
		VISITES2= visits.m_Vertices;
		VISITES2.pop_back();
		VISITES2.pop_front();

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();
		return SalesmanTrackBacktrackingRec(NULL, VISITES2.front());
	}


	return NULL;
}



// =============================================================================
// SalesmanTrackBacktrackingGreedy =============================================
// =============================================================================


CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits)
{
	return CTrack(&graph);
}
