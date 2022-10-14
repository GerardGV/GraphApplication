#include "pch.h"
#include "Graph.h"


// =============================================================================
// Dijkstra ====================================================================
// =============================================================================

void Dijkstra(CGraph& graph, CVertex *pStart)
{
	

	//we equal all the lenghts of the grapfh to infinite
	for (CEdge* e : pStart->m_Edges ) {
		
		if (e->m_Name != pStart->m_Name)
		{
			e->m_Length = numeric_limits<double>::max();
		}
		
	}

}

// =============================================================================
// DijkstraQueue ===============================================================
// =============================================================================

void DijkstraQueue(CGraph& graph, CVertex *pStart)
{
}
