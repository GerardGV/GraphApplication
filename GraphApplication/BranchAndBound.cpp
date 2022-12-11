#include "pch.h"
#include "Graph.h"
#include <queue>
#include <iostream>
#include <iomanip> 
#include <numeric>



// SalesmanTrackBranchAndBound1 ===================================================

class CBBNode1{
public:
	vector<int> m_indexes;
	vector<bool> m_visited;

	//B&B1
	double m_dijkstraLenght;

	//B&B2
	double m_cotaInferior;
	double m_cotaSuperior;


	//B&b1
	CBBNode1(double dijkstraDistance, vector<int> indexes, int nVertex) :m_dijkstraLenght(dijkstraDistance), m_indexes(indexes) {

		m_visited.resize(nVertex, false);

		//first one has been visited
		m_visited[0] = true;

	}

	//B&B2
	CBBNode1(double dijkstraDistance, vector<int> indexes, int nVertex, double cotaInferior, double cotaSuperior) :m_dijkstraLenght(dijkstraDistance), m_indexes(indexes), m_cotaInferior(cotaInferior), m_cotaSuperior(cotaSuperior) {

		m_visited.resize(nVertex, false);

		//first one has been visited
		m_visited[0] = true;

	}



	/*	
	CBBNode1(double dijkstraDistance, CEdge* pEdge, CVertex* pDestination) :m_dijkstraLenght(dijkstraDistance)
		, cotaInferior(m_dijkstraLenght + pDestination->m_Point.Distance(pEdge->m_pDestination->m_Point)) {}
	*/

	CBBNode1(const CBBNode1& node)
	{
		m_indexes = node.m_indexes;

		m_visited = node.m_visited;
	
	}

};


struct comparator1 {
	bool operator()(const CBBNode1* s1, const CBBNode1* s2) {
		return s1->m_dijkstraLenght > s2->m_dijkstraLenght;
	}
};


struct comparator2 {
	bool operator()(const CBBNode1* s1, const CBBNode1* s2) {
		return s1->m_cotaInferior > s2->m_cotaInferior;
	}
};

class Cell {
	public:
	double dijkstraDistance;
	list<CEdge*> tram;
};





CTrack SalesmanTrackBranchAndBound1(CGraph& graph, CVisits& visits)
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


	//------------------------------------------------------------------

	priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator1> allPossibleSolution;

	//we explore the first node and we add ways, choosign with the priority queue the next one that we are going to explore, it dependes of DijkstraDistance to each one
	//but with don't add the cell v1 column v1 because the dijkstraDistance is 0 and that generate a loop adn we are not going to explore the destination node so size-1
	for (int col = 1; col < visits.m_Vertices.size(); col++) {

		//we create te vector where we are going to seave the indexes of th tram
		vector<int> vIndex;
		vIndex.push_back(0);
		vIndex.push_back(col);
		allPossibleSolution.push(new CBBNode1(matriu[0][col].dijkstraDistance, vIndex, visits.m_Vertices.size()));

	}




	//we are going to iterate until find a solution,  the first one it's the optimal in B&B, or not get any solution
	while (!allPossibleSolution.empty())
	{
		//we take the tram with less weight(DijkstraDistance) and we compare if the las node in this is the destination
		CBBNode1* topNode = allPossibleSolution.top();
		allPossibleSolution.pop();

		//we mark the next vertex in the list of index as visited
		topNode->m_visited[topNode->m_indexes.back()] = true;

		//if the last index is equal to the index of the destination that means that we are in a solution
		if (topNode->m_indexes.back() == visits.m_Vertices.back()->m_indexMatrix)
		{
			//we check if the possible solution it's right
			for (bool visitVertex : topNode->m_visited)
			{
				if (!visitVertex)
				{
					goto noSolution;
				}
			}

			//we create the CTrack with th first section already addes
			CTrack solution = CTrack(&graph, matriu[topNode->m_indexes[0]][topNode->m_indexes[1]].tram);

			//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
			for (int element = 1; element < topNode->m_indexes.size() - 1; element++) {

				for (CEdge* e : matriu[topNode->m_indexes[element]][topNode->m_indexes[element + 1]].tram)
					solution.m_Edges.push_back(e);
				//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
			}


			//we release the memory
			delete topNode;

			return solution;

		}

	noSolution:

		//we unmark the next vertex in the list of index as visited
		//visited[topNode->m_indexes.back()] = false;

		//we are not going to explore the destination vertex
		if (topNode->m_indexes.back() != visits.m_Vertices.back()->m_indexMatrix)
		{//if we don't get a solution we push all the new possible section of the solution way plus the current lenght
			for (int col = 0; col < visits.m_Vertices.size(); col++)
			{
				//we control that we don't go to the same vertex, if we don't put this conditional, we will always in the same vertex because Dijkstra distance to it is 0
				if (col != topNode->m_indexes.back() && !topNode->m_visited[col])
				{
					//CBBNode(double* dijkstraDistance, vector<int> indexes)
					//CBBNode* newNode(topNode);

					CBBNode1* newNode = new CBBNode1(*topNode);


					newNode->m_indexes.push_back(col);

					newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

					allPossibleSolution.push(newNode);
				}

			}
		}


	}

	return CTrack(&graph);

}


// SalesmanTrackBranchAndBound2 ===================================================


CTrack SalesmanTrackBranchAndBound2(CGraph& graph, CVisits &visits)
{

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

	//if we only have orgin and destination with return Dijkstra way
	if (visits.m_Vertices.size() == 2) {

		//we create the CTrack with th first section already addes
		return CTrack(&graph, matriu[0][1].tram);


	}
	else {

		//we put the index of the column that represents the destination vertex in the matrix
		visits.m_Vertices.back()->m_indexMatrix = index;

		//-----------------------------------------

		priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator2> allPossibleSolution;


		//creation of cotas


		//each element of next vectors is the DujkstraDistance to one node that we have to visited.
		vector <double> min;
		vector <double> max;
		min.resize(visits.m_Vertices.size(), numeric_limits<double>::max());

		max.resize(visits.m_Vertices.size(), 0.0);

		for (int row = 0; row < visits.m_Vertices.size() - 1; row++)
		{
			//we check if we find another max or min DijkstraDistance in the new row to nodes
			for (int col = 0; col < visits.m_Vertices.size(); col++)
			{
				if (row != col) {
					if (matriu[row][col].dijkstraDistance < min[col])
					{
						min[col] = matriu[row][col].dijkstraDistance;
					}

					if (matriu[row][col].dijkstraDistance > max[col])
					{
						max[col] = matriu[row][col].dijkstraDistance;
					}
				}

			}
		}

		//summation of all the minimums and maximums ways to each vertex
		double cotaInferior = accumulate(min.begin(), min.end(), 0.0);
		double cotaSuperior = accumulate(max.begin(), max.end(), 0.0) + pow(10, -6);

		double globalCotaSuperior = cotaSuperior;

		//It is created the first node with the initial cotas
		vector<int> vIndex;
		vIndex.push_back(0);
		allPossibleSolution.push(new CBBNode1(0, vIndex, visits.m_Vertices.size(), cotaInferior, cotaSuperior));


		//we are going to iterate until find a solution,  the first one it's the optimal in B&B, or not get any solution
		while (!allPossibleSolution.empty())
		{
			//we take the tram with less weight(DijkstraDistance) and we compare if the las node in this is the destination
			CBBNode1* topNode = allPossibleSolution.top();
			allPossibleSolution.pop();

			//we mark the next vertex in the list of index as visited
			topNode->m_visited[topNode->m_indexes.back()] = true;

			//if the last index is equal to the index of the destination that means that we are in a solution
			if (topNode->m_indexes.back() == visits.m_Vertices.back()->m_indexMatrix)
			{
				//we check if the possible solution it's right
				for (bool visitVertex : topNode->m_visited)
				{
					if (!visitVertex)
					{
						goto noSolution;
					}
				}

				//we create the CTrack with th first section already addes
				CTrack solution = CTrack(&graph, matriu[topNode->m_indexes[0]][topNode->m_indexes[1]].tram);

				//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
				for (int element = 1; element < topNode->m_indexes.size() - 1; element++) {

					for (CEdge* e : matriu[topNode->m_indexes[element]][topNode->m_indexes[element + 1]].tram)
						solution.m_Edges.push_back(e);
					//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
				}


				//we release the memory
				delete topNode;

				return solution;

			}

		noSolution:


			//we are not going to explore the destination vertex
			if (topNode->m_indexes.back() != visits.m_Vertices.back()->m_indexMatrix)
			{//if we don't get a solution we push all the new possible section of the solution way plus the current lenght

	

				for (int col = 0; col < visits.m_Vertices.size(); col++)
				{

					//we control that we don't go to the same vertex, if we don't put this conditional, we will always in the same vertex because Dijkstra distance to it is 0
					//It is controled that we are not goint to go to visited vertex and we are not going to create a way with only the origin and the destination origin->destination
					if (col != topNode->m_indexes.back() && !topNode->m_visited[col] && (col != visits.m_Vertices.size()-1 || topNode->m_indexes.size() != 1))
					{


						//calculation of the new cotaSuperior to knwo if it is interesting to continous exploring this solutions
						cotaInferior = topNode->m_cotaInferior - min[col] + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

						if (globalCotaSuperior > cotaInferior)
						{


							//cota inferior calculation
							cotaSuperior = topNode->m_cotaSuperior - max[col] + matriu[topNode->m_indexes.back()][col].dijkstraDistance + pow(10, -6);

							//when it is found a smaller cotaSuperior we take it as the new global

							if (globalCotaSuperior > cotaSuperior)
							{
								globalCotaSuperior = cotaSuperior;
							}

						
							//Creation of the new solution or node
							CBBNode1* newNode = new CBBNode1(*topNode);
							newNode->m_cotaInferior = cotaInferior;
							newNode->m_cotaSuperior = cotaSuperior;

							newNode->m_indexes.push_back(col);

						
							newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

							allPossibleSolution.push(newNode);
						}


					}

				}
			}


		}

	}

	return CTrack(&graph);

}



// SalesmanTrackBranchAndBound3 ===================================================


CTrack SalesmanTrackBranchAndBound3(CGraph& graph, CVisits &visits)
{

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

	//if we only have orgin and destination with return Dijkstra way
	if (visits.m_Vertices.size() == 2) {

		//we create the CTrack with th first section already addes
		return CTrack(&graph, matriu[0][1].tram);


	}
	else {

		//we put the index of the column that represents the destination vertex in the matrix
		visits.m_Vertices.back()->m_indexMatrix = index;

		//-----------------------------------------

		priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator2> allPossibleSolution;


		//creation of cotas


		//each element of next vectors is the DujkstraDistance to one node that we have to visited.
		vector <double> min;
		vector <double> max;
		min.resize(visits.m_Vertices.size(), numeric_limits<double>::max());

		max.resize(visits.m_Vertices.size(), 0.0);

		for (int row = 0; row < visits.m_Vertices.size() - 1; row++)
		{
			//we check if we find another max or min DijkstraDistance in the new row to nodes
			for (int col = 0; col < visits.m_Vertices.size(); col++)
			{
				if (row != col) {
					if (matriu[row][col].dijkstraDistance < min[col])
					{
						min[col] = matriu[row][col].dijkstraDistance;
					}

					if (matriu[row][col].dijkstraDistance > max[col])
					{
						max[col] = matriu[row][col].dijkstraDistance;
					}
				}

			}
		}

		//summation of all the minimums and maximums ways to each vertex
		double cotaInferior = accumulate(min.begin(), min.end(), 0.0);
		double cotaSuperior = accumulate(max.begin(), max.end(), 0.0) + pow(10, -6);

		double globalCotaSuperior = cotaSuperior;

		//It is created the first node with the initial cotas
		vector<int> vIndex;
		vIndex.push_back(0);
		allPossibleSolution.push(new CBBNode1(0, vIndex, visits.m_Vertices.size(), cotaInferior, cotaSuperior));


		//we are going to iterate until find a solution,  the first one it's the optimal in B&B, or not get any solution
		while (!allPossibleSolution.empty())
		{

		noSolution:
			//we take the tram with less weight(DijkstraDistance) and we compare if the las node in this is the destination
			CBBNode1* topNode = allPossibleSolution.top();
			allPossibleSolution.pop();

			//we mark the next vertex in the list of index as visited
			topNode->m_visited[topNode->m_indexes.back()] = true;

			//if the last index is equal to the index of the destination that means that we are in a solution
			if (topNode->m_indexes.back() == visits.m_Vertices.back()->m_indexMatrix)
			{
				//we check if the possible solution it's right
				/*
				for (bool visitVertex : topNode->m_visited)
				{
					if (!visitVertex)
					{
						goto noSolution;
					}
				}
				*/

				//we get the unique vertex that we visited
				auto ip = std::unique(topNode->m_indexes.begin(), topNode->m_indexes.begin() + topNode->m_indexes.size());

				topNode->m_indexes.resize(distance(topNode->m_indexes.begin(), ip));

				//we check if we have the same number of unique values as the number of vertex we have to visited
				if (topNode->m_indexes.size() != visits.m_Vertices.size())
					goto noSolution;

				//we create the CTrack with th first section already addes
				CTrack solution = CTrack(&graph, matriu[topNode->m_indexes[0]][topNode->m_indexes[1]].tram);

				//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
				for (int element = 1; element < topNode->m_indexes.size() - 1; element++) {

					for (CEdge* e : matriu[topNode->m_indexes[element]][topNode->m_indexes[element + 1]].tram)
						solution.m_Edges.push_back(e);
					//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
				}


				//we release the memory
				delete topNode;

				return solution;

			}

	


			
			//if we don't get a solution we push all the new possible section of the solution way plus the current lenght

			for (int row = 0; row < visits.m_Vertices.size() - 1; row++)
			{
				//It's not consider the rows of visited vertex
				//bool a1 = topNode->m_visited[row] == false;
				//bool a2 = row == topNode->m_indexes.back();
				if (!topNode->m_visited[row]) {

					//we check if we find another max or min DijkstraDistance in the new row to nodes
					for (int col = 0; col < visits.m_Vertices.size(); col++)
					{
						//It is not considere the column of the last visited vertex

						if (row != col && !topNode->m_visited[col]) {
							if (matriu[row][col].dijkstraDistance < min[col])
							{
								min[col] = matriu[row][col].dijkstraDistance;
							}

							if (matriu[row][col].dijkstraDistance > max[col])
							{
								max[col] = matriu[row][col].dijkstraDistance;
							}
						}




					}
				}

			}


			//loop to check the row of the last vertex in the way, It is not considered going to the final vertex
			for (int col = 0; col < visits.m_Vertices.size() - 1; col++) {

				//It is not considered the diagonal of the matrix because it is 0, we don't consider column of visited vertex
				if (topNode->m_indexes.back() != col && !topNode->m_visited[col]) {

					if (matriu[topNode->m_indexes.back()][col].dijkstraDistance < min[col])
					{
						min[col] = matriu[topNode->m_indexes.back()][col].dijkstraDistance;
					}

					if (matriu[topNode->m_indexes.back()][col].dijkstraDistance > max[col])
					{
						max[col] = matriu[topNode->m_indexes.back()][col].dijkstraDistance;
					}
				}
			}

			for (int col = 0; col < visits.m_Vertices.size(); col++)
			{

				//we control that we don't go to the same vertex, if we don't put this conditional, we will always in the same vertex because Dijkstra distance to it is 0
				//It is controled that we are not goint to go to visited vertex and we are not going to create a way with only the origin and the destination origin->destination
				if (col != topNode->m_indexes.back() && !topNode->m_visited[col] && (col != visits.m_Vertices.size() - 1 || topNode->m_indexes.size() != 1))
				{


					//calculation of the new cotaSuperior to knwo if it is interesting to continous exploring this solutions
					cotaInferior = topNode->m_cotaInferior - min[col] + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

					if (globalCotaSuperior > cotaInferior)
					{


						//cota inferior calculation
						cotaSuperior = topNode->m_cotaSuperior - max[col] + matriu[topNode->m_indexes.back()][col].dijkstraDistance + pow(10, -6);

						//when it is found a smaller cotaSuperior we take it as the new global

						if (globalCotaSuperior > cotaSuperior)
						{
							globalCotaSuperior = cotaSuperior;
						}


						//Creation of the new solution or node
						CBBNode1* newNode = new CBBNode1(*topNode);
						newNode->m_cotaInferior = cotaInferior;
						newNode->m_cotaSuperior = cotaSuperior;

						newNode->m_indexes.push_back(col);


						newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

						allPossibleSolution.push(newNode);
					}


				}

			}
			


		}

	}

	return CTrack(&graph);
}
