#include "pch.h"
#include "Graph.h"
#include <queue>
#include <iostream>
#include <iomanip> 
#include <random>
#include <chrono>

#define GRADIENT_ITERATION 10000 
#define SEED 1
// SalesmanTrackProbabilistic ==================================================

class Cell {
public:
	double dijkstraDistance;
	list<CEdge*> tram;
};



CTrack SalesmanTrackProbabilistic(CGraph& graph, CVisits& visits)
{

	if (!visits.m_Vertices.empty())
	{

		//MATRIX CREATION

		//we create ways matrix with n vextors = n-1 vistis because we dont want the dijkstra distance from final vertex to the others
		vector<vector<Cell>> matriu(visits.m_Vertices.size() - 1);
		int index = 0;

		for (auto vOriginDijks = visits.m_Vertices.begin(); vOriginDijks != next(visits.m_Vertices.begin(), visits.m_Vertices.size() - 1); vOriginDijks++)
		{
			DijkstraQueue(graph, *vOriginDijks);

			//we save in each vertex its index from the matrix
			(*vOriginDijks)->m_indexMatrix = index;

			for (CVertex* vDestinationDijks : visits.m_Vertices) {

				Cell matrixCell;
				//we have to save in each vell of the matrix the DIJKSTRA way and de distance of this way
				matrixCell.dijkstraDistance = vDestinationDijks->m_DijkstraDistance;

				//we dont have dijkstraPrevios in the vertex from we explore because is th origin of Dijkstra
				if (*vOriginDijks != vDestinationDijks)
				{
					while (vDestinationDijks != *vOriginDijks)
					{
						matrixCell.tram.push_front(vDestinationDijks->m_pDijkstraPrevious);
						vDestinationDijks = vDestinationDijks->m_pDijkstraPrevious->m_pOrigin;
					}
				}

				matriu[index].push_back(matrixCell);

			}


			index++;
		}

		//we put the index of the column that represents the destination vertex in the matrix
		visits.m_Vertices.back()->m_indexMatrix = index;

		//-------------------------------------------------------------------

		//we choose randome begining

		deque<int> range;//we need some place to save the vertes that still can be choose

		for (CVertex* v : visits.m_Vertices)
		{
			range.push_back(v->m_indexMatrix);
		}
		range.pop_front();
		range.pop_back();

		deque<int> copyRange(range);

		int newVertex;
		int indexVertex;
		double firstLenght = 0;

		vector<int> firstSolution;
		//we always start at INICI vertex
		firstSolution.push_back(0);
		while (!copyRange.empty())
		{
			//srand(SEED);//change the seed of randome generator
			index = (rand() % copyRange.size());
			newVertex = copyRange[index];

			//we take out the elemenet that has been chosen
			copyRange.erase(next(copyRange.begin(), index));

			firstLenght += matriu[firstSolution.back()][newVertex].dijkstraDistance;
			firstSolution.push_back(newVertex);
		}

		//adding the lenght from the las vertex in the way to the destination vertex
		firstLenght += matriu[firstSolution.back()][visits.m_Vertices.back()->m_indexMatrix].dijkstraDistance;
		firstSolution.push_back(visits.m_Vertices.back()->m_indexMatrix);


		//Descens de Gradient

		//same as before and at the end of the loop we are going to compare the best solution got it in Descens de Grdient with the first solution that we created after the creation of the matrix
		double GradientLenght = std::numeric_limits<double>::max();
		vector<int> GradientSolution;
		for (int i = 0; i < GRADIENT_ITERATION*visits.m_Vertices.size(); i++)
		{
		other:;
			//we get back all the possible vertex to visit to calculate a new possible solution
			copyRange = range;

			double currentLenght = 0;
			int newVertex;
			int indexGrad;
			vector<int> currentSolution;
			currentSolution.push_back(0);
			
			while (!copyRange.empty() && currentLenght < GradientLenght)
			{
				if (currentLenght > GradientLenght)
					goto other;
				//srand(SEED);//change the seed of randome generator
				indexGrad = (rand() % copyRange.size());
				newVertex = copyRange[indexGrad];

				//we take out the elemenet that has been chosen
				copyRange.erase(next(copyRange.begin(), indexGrad));

				currentLenght += matriu[currentSolution.back()][newVertex].dijkstraDistance;
				currentSolution.push_back(newVertex);
			}

			//adding the lenght from the las vertex in the way to the destination vertex
			currentLenght += matriu[currentSolution.back()][visits.m_Vertices.back()->m_indexMatrix].dijkstraDistance;
			currentSolution.push_back(visits.m_Vertices.back()->m_indexMatrix);

			//we save the best solution
			if (currentLenght < GradientLenght)
			{
				GradientLenght = currentLenght;
				GradientSolution = currentSolution;

			}

		}

		//we create the solution the optimal one
		if (GradientLenght < firstLenght) {

			//we create the CTrack with th first section already addes
			CTrack solution = CTrack(&graph, matriu[GradientSolution[0]][GradientSolution[1]].tram);

			//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
			for (int element = 1; element < GradientSolution.size() - 1; element++) {

				for (CEdge* e : matriu[GradientSolution[element]][GradientSolution[element + 1]].tram)
					solution.m_Edges.push_back(e);
				//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
			}

			return solution;
		}
		else {
			//we create the CTrack with th first section already addes
			CTrack solution = CTrack(&graph, matriu[firstSolution[0]][firstSolution[1]].tram);

			//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
			for (int element = 1; element < firstSolution.size() - 1; element++) {

				for (CEdge* e : matriu[firstSolution[element]][firstSolution[element + 1]].tram)
					solution.m_Edges.push_back(e);
				//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
			}

			return solution;
		}

		
	}


	return CTrack(&graph);
}
