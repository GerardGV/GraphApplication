#include "pch.h"
#include "Graph.h"
// SalesmanTrackGreedy =========================================================

CTrack SalesmanTrackGreedy(CGraph& graph, CVisits &visits)
{
	
	//result track
	CTrack res= CTrack(&graph);
	CTrack aux = CTrack(&graph);

	//creating v, the first visit
	//auto firstVisita = visits.m_Vertices.begin();
	//CVertex* v = *firstVisita;


	CVertex* v = visits.m_Vertices.front();

	//we create a list of candidates without the first and the last candidate
	list<CVertex*> candidats(visits.m_Vertices);
	candidats.pop_back();
	candidats.pop_front();

	
	


	while (!candidats.empty())
	{
		//select v1, it is which has less Dijkstra distance
		CVertex* v1 = candidats.front();

		//apply Dijkstra
		DijkstraQueue(graph, v);

		//select v1, it is which has less Dijkstra distance
		v1 =  candidats.front();
		auto it = candidats.begin();
		it++;
		for (;it != candidats.end();it++)
		{
			if(v1->m_DijkstraDistance > (*(*it)).m_DijkstraDistance)
			{
				v1 = (*it);
			}
		}

		//add the edges between v and v1 to make get the 
		aux.m_Edges.push_front(v1->m_pDijkstraPrevious);
		CVertex* vMiddle=v1;
		while (vMiddle->m_pDijkstraPrevious->m_pOrigin->m_Name != v->m_Name)
		{
			//we go to the next vertex between v and v1
			vMiddle = vMiddle->m_pDijkstraPrevious->m_pOrigin;
			aux.m_Edges.push_front(vMiddle->m_pDijkstraPrevious);
		}

		//we join the track unti v1 with the answer track
		res.Append(aux);

		//we have to empty the edges of the auxiliar because we have already add them
		aux.Clear();

		//take out candidate from Cnadidates
		for (auto it = candidats.begin(); it != candidats.end(); it++)
		{
			if ((*it)->m_Name == v1->m_Name)
			{
				candidats.erase(it);
				break;
			}
		}

		//v1->m_JaHePassat = true;
		
		v = v1;
	}

	//apply Dijkstra
	DijkstraQueue(graph, v);

	//add the edges between v and the final destination 
	aux.m_Edges.push_front(visits.m_Vertices.back()->m_pDijkstraPrevious);
	CVertex* vMiddle = visits.m_Vertices.back();
	while (vMiddle->m_pDijkstraPrevious->m_pOrigin->m_Name != v->m_Name)
	{
		//we go to the next vertex between v and v1
		vMiddle = vMiddle->m_pDijkstraPrevious->m_pOrigin;
		aux.m_Edges.push_front(vMiddle->m_pDijkstraPrevious);
	}

	//we join the track unti v1 with the answer track
	res.Append(aux);

	return res;

	//return CTrack(&graph);
}
