#include "pch.h"
#include "Graph.h"
#include <set>
#include <queue>


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


//versio de SALESMAN backtracking pur amb CTRACK per guarda cami actual

/*
void SalesmanTrackBacktrackingRec(CVertex* pActual)
{
	//we are not going to explore a vertex if the current lengh is already longer that the shortest way
	if (LongitudCamiActual2 < LongitudCamiMesCurt2)
	{
		if (pActual == pDesti2) {
		 

			bool correct(true);
			//int correct(0);

			//we check that we passed all the vertex that we have to visit
			for (CVertex* v1 : VISITES2)
			{
				//if we didn't pass through one vertex of visits, the way is no a possible answer
				
				if (v1->m_count == 0)
				{
					correct = false;
					break;
				}
				

				//solucio proposta per el profe, comprobar que surt una aresta dels vertex que hem de visitar
				//for (CEdge* edge:v1->m_Edges)
				//{
					//if (edge->m_Used) {
						
						//correct++;
						//break;
					//}
				//}
				

			}

			//if the way is a possible answer
			if (correct)
			{
				
				//CamiMesCurt2.Clear();
				//while (pAnterior) {
					//CamiMesCurt2.m_Edges.push_front(pAnterior->m_pEdge);
					//pAnterior = pAnterior->m_pAnterior;
				//}
				

				CamiMesCurt2 = CamiActual2;

				LongitudCamiMesCurt2 = LongitudCamiActual2;


			}
			else
			{
				//node to  save actual way
				//NodeCami2 node;
				//node.m_pAnterior = pAnterior;

				for (CEdge* pE : pActual->m_Edges) {
					//if (!pE->m_pDestination->m_JaHePassat) {

					if (!pE->m_Used)
					{
						pE->m_Used = true;
						//node.m_pEdge = pE;

						CamiActual2.m_Edges.push_back(pE);

						LongitudCamiActual2 += pE->m_Length;


						//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
						//cout << "Edge: " << pE->m_Name << endl;
						SalesmanTrackBacktrackingRec(pE->m_pDestination);
						//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;

						CamiActual2.m_Edges.pop_back();
						LongitudCamiActual2 -= pE->m_Length;
						pE->m_Used = false;
					}
				}

			}
		
		}
		else
		{

			//pActual->m_JaHePassat = true;
			
			//we pass +1 time trhough the vertex
			pActual->m_count++;

			//node to  save actual way
			//NodeCami2 node;
			//node.m_pAnterior = pAnterior;
			for (CEdge* pE : pActual->m_Edges) {

				if (!pE->m_Used)
				{
					pE->m_Used = true;
					//node.m_pEdge = pE;

					CamiActual2.m_Edges.push_back(pE);

					LongitudCamiActual2 += pE->m_Length;


					//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
					//cout << "Edge: " << pE->m_Name << endl;
					SalesmanTrackBacktrackingRec(pE->m_pDestination);
					//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;

					CamiActual2.m_Edges.pop_back();
					LongitudCamiActual2 -= pE->m_Length;
					pE->m_Used = false;
				}

			}

			//we get to the previous state
			pActual->m_count--;
			
		}
	}
	

	//return CamiMesCurt2;
}
*/

//versio per fer crida de recursiu amb Ctrack per guardar cami actual
/*
CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits)
{
	if (!graph.m_Edges.empty())
	{
		pDesti2 = visits.m_Vertices.back();
		//we take out the inicial node and the destination
		VISITES2 = visits.m_Vertices;
		VISITES2.pop_back();
		VISITES2.pop_front();

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();
		//NodeCami2 first;
		//we add the first
		//first.m_pEdge = visits.m_Vertices.front()->m_Edges.front();
		//we took the first neighbour that is the destination of the first edge of the beginner 

		//actual = CTrack(&graph);

		SalesmanTrackBacktrackingRec(visits.m_Vertices.front());

		return CamiMesCurt2;
	}


	return NULL;
}
*/

//versio versio de SALESMAN backtracking puramb amb llista enllaçada NodeCami2 per guarda cami actual
 
void SalesmanTrackBacktrackingRec(NodeCami2* pAnterior, CVertex* pActual)
{
	//we are not going to explore a vertex if the current lengh is already longer that the shortest way
	if (LongitudCamiActual2 < LongitudCamiMesCurt2)
	{

		if (pActual == pDesti2) {


			//this variable is to control that we have visited all the vertex that we must
			bool correct = true;

			//we check that we passed all the vertex that we have to visit
			for (CVertex* v1 : VISITES2)
			{
				//if we didn't pass through one vertex of visits, the way is no a possible answer
				if (v1->m_count == 0)
				{
					correct = false;
					break;
				}
			}

			//if the way is a possible answer
			if (correct)
			{
				//we update the optimal way with the way that we have found
				CamiMesCurt2.Clear();
				while (pAnterior) {
					CamiMesCurt2.m_Edges.push_front(pAnterior->m_pEdge);
					pAnterior = pAnterior->m_pAnterior;
				}
				LongitudCamiMesCurt2 = LongitudCamiActual2;



			}
			else
			{
				//node to  save actual way
				NodeCami2 node;
				node.m_pAnterior = pAnterior;
				for (CEdge* pE : pActual->m_Edges) {

					if (!pE->m_Used)
					{
						pE->m_Used = true;
						node.m_pEdge = pE;
						LongitudCamiActual2 += pE->m_Length;
						//cout << "Edge: " << pE->m_Name << endl;
						SalesmanTrackBacktrackingRec(&node, pE->m_pDestination);
						//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;
						LongitudCamiActual2 -= pE->m_Length;
						pE->m_Used = false;
					}
				}

			}
			
		}
		else
		{

			pActual->m_count++;

			//node to  save actual way
			NodeCami2 node;
			node.m_pAnterior = pAnterior;
			for (CEdge* pE : pActual->m_Edges) {


				if (!pE->m_Used)
				{
					pE->m_Used = true;
					node.m_pEdge = pE;
					LongitudCamiActual2 += pE->m_Length;
					//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
					//cout << "Edge: " << pE->m_Name << endl;
					SalesmanTrackBacktrackingRec(&node, pE->m_pDestination);
					//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;
					LongitudCamiActual2 -= pE->m_Length;
					pE->m_Used = false;
				}

			}

			//we get to the previous state
			pActual->m_count--;
			
		}
	}


}


CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits)
{
	if (!graph.m_Edges.empty())
	{
		pDesti2 = visits.m_Vertices.back();
		//we take out the inicial node and the destination
		VISITES2 = visits.m_Vertices;
		VISITES2.pop_back();
		VISITES2.pop_front();

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();

		SalesmanTrackBacktrackingRec(NULL, visits.m_Vertices.front());

		return CamiMesCurt2;
	}


	return NULL;
}



// =============================================================================
// SalesmanTrackBacktrackingGreedy =============================================
// =============================================================================


struct comparator {
	bool operator()(CVertex* pE1, CVertex* pE2) {
		return pE1->m_DijkstraDistance > pE2->m_DijkstraDistance;
	}
};


//we are going to use a priority queue because we are going to explore the vertex in order of the Dijkstra Distance. 
	//from the cheapest to the most expensive
vector<priority_queue<CVertex*, std::vector<CVertex*>, comparator>> dijkstraDistancesMatrix;



struct index{

	index(CVertex* indexFirst)
	{
		indexPrevi->first = indexFirst;
	}

	pair<CVertex*,CVertex*>* indexPrevi;

	index* m_pAnterior;
};


CVisits* VISITES;



void SalesmanTrackBacktrackingGreedyRec(index* indexAterior, CVertex* pActual) {

	//if the current lenght is bigger we are not to continue
	if (LongitudCamiActual2 < LongitudCamiMesCurt2) {
		
		if (pActual == pDesti2)
		{
			//we update the optimal way with the current
			LongitudCamiMesCurt2 = LongitudCamiActual2;

			//we update the optimal way with the way that we have found
			CamiMesCurt2.Clear();


			while (indexAterior) {

				//we apply Dijkstra to get all the edges of the section that are part othe optimal way
				DijkstraQueue(*VISITES->m_pGraph, indexAterior->indexPrevi->first);//this is because each time we had execute Dijkstra to create the matrix, we update the DijkstraPrecious edges

				auto it = VISITES->m_Vertices.begin();

				advance(it, indexAterior->indexPrevi->second->m_indexMatrix);

				//we compare that the origin of the edge is the vertex from where we had calculated all dijkstraDistances
				while ((*it)->m_pDijkstraPrevious->m_pOrigin != indexAterior->indexPrevi->first)
				{
					//we create the optimal way adding all its sections
					CamiMesCurt2.m_Edges.push_front((*it)->m_pDijkstraPrevious);
				}

				//we go to the next node in the linked list
				indexAterior = indexAterior->m_pAnterior;
			}


		}
		else {

			//we create a node of the linked list where we save section of the optimal way
			index indexActual(pActual);
			//indexActual.indexPrevi->first = pActual;

			pair<CVertex*, CVertex*> actualSection;
			while (!dijkstraDistancesMatrix[pActual->m_indexMatrix].empty())
			{
				//we add the weight of the optimal way section to the total weight of the way
				LongitudCamiActual2+= dijkstraDistancesMatrix[pActual->m_indexMatrix].top()->m_DijkstraDistance;

				//we update the destination vertex of the section
				indexActual.indexPrevi->second= dijkstraDistancesMatrix[pActual->m_indexMatrix].top();

				//we sent the nearest vertex to the current vertex to the recursive function to explore it
				SalesmanTrackBacktrackingGreedyRec(&indexActual, dijkstraDistancesMatrix[pActual->m_indexMatrix].top());//the nearest vertex is on the top of the priority queue
			
				//we take out the weight of the optimal way section to do a step back
				LongitudCamiActual2 -= dijkstraDistancesMatrix[pActual->m_indexMatrix].top()->m_DijkstraDistance;

				//we take it out and we try with the next near vertex
				dijkstraDistancesMatrix[pActual->m_indexMatrix].pop();

			}
			
		}
	}

	
}



CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits)
{
	//we control that the graf has edges
	if (!graph.m_Edges.empty())
	{
		//int index = 0;
		int row = 0;
		for (CVertex* v : visits.m_Vertices) {

			//we calculate all the minimum weight from the vertex v until all the vertex of the graph 
			DijkstraQueue(graph, v);

			//we create a priority queue to simulate a column in the matrix, each element of the vector Dijkstra distance
			//simulates a row in the matrix
			priority_queue<CVertex*, std::vector<CVertex*>, comparator> colomn;
			dijkstraDistancesMatrix.push_back(colomn);
			v->m_indexMatrix = row;
			//index++;
			

			//we sort  in each row with the priority queue the vertexs that we have to visited depending of the weight, each vertex is a diferent column
			for (CVertex* v: visits.m_Vertices) {
				dijkstraDistancesMatrix[row].push(v);

			}

			//we take out the first column because it is gonna refer to the same vertex of thw row because columns are priorities queue an the distance to itself is 0
			dijkstraDistancesMatrix[row].pop();


			//we go to the next vertex
			row++;
		}

		//we get the destination and we select all the visits that can be or not visitied with the current way
		pDesti2 = visits.m_Vertices.back();
		
		//we put the vertex that we need to visit as global because we need them in case we get to a possible solution
		VISITES = &visits;

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();

		//we pass the current graph to a global variable because we are not going to modify it
		//graphToRec = &graph;

		//we look for the best way to visit all the nodes
		SalesmanTrackBacktrackingGreedyRec(NULL, visits.m_Vertices.front());
	}

	return CamiMesCurt2;
}
