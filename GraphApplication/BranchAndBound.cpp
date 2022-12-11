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

	CBBNode1(double dijkstraDistance, vector<int> indexes, int nVertex, int secondVertex):m_dijkstraLenght(dijkstraDistance), m_indexes(indexes) { 
		
		m_visited.resize(nVertex, false); 

		//first one has been visited
		m_visited[0] = true;

		//we have possible solution with 3 vertex always so we put the first and the second one as visited
		m_visited[secondVertex] = true;
	}

	CBBNode1(double dijkstraDistance, vector<int> indexes, int nVertex) :m_dijkstraLenght(dijkstraDistance), m_indexes(indexes) {

		m_visited.resize(nVertex, false);

		//first one has been visited
		m_visited[0] = true;

	}


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
		return (s1->m_dijkstraLenght + (pow(10, -6) + s1->m_cotaSuperior)) > (s2->m_dijkstraLenght + (pow(10, -6) + s2->m_cotaSuperior));
	}
};


class Cell {
	public:
	double dijkstraDistance;
	list<CEdge*> tram;
};


/*
vector<int> correctSolution = { 0,1,2,4,3,5,2,6, 7};

CTrack SalesmanTrackBranchAndBound1(CGraph& graph, CVisits& visits)
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

	//we put the index of the column that represents the destination vertex in the matrix
	visits.m_Vertices.back()->m_indexMatrix = index;




	priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator1> allPossibleSolution;

	//we explore the first node and we add ways, choosign with the priority queue the next one that we are going to explore, it dependes of DijkstraDistance to each one
	//but with don't add the cell v1 column v1 because the dijkstraDistance is 0 and that generate a loop adn we are not going to explore the destination node so size-1
	for (int vertex1 = 1; vertex1 < visits.m_Vertices.size()-1 ; vertex1++) {

		//we create te vector where we are going to seave the indexes of th tram
		vector<int> vIndex;
		vIndex.resize(3, 0);
		//vIndex.push_back(0);
		vIndex[1]=vertex1;

		//if we don't have a simple graph, we start with 3 elements vectors of posible solutions
		for (int vertex2 = 0; vertex2 < visits.m_Vertices.size(); vertex2++)
		{

			//to avoid adding v1->v2->v2
			if (vertex1 != vertex2) {

				vIndex[2]=vertex2;

				/*
				CBBNode1* node = new CBBNode1((matriu[0][vertex1].dijkstraDistance + matriu[vertex1][vertex2].dijkstraDistance), vIndex, visits.m_Vertices.size());
				node->m_visited[1] = true;
				allPossibleSolution.push(node);
				


				allPossibleSolution.push(new CBBNode1((matriu[0][vertex1].dijkstraDistance + matriu[vertex1][vertex2].dijkstraDistance), vIndex, visits.m_Vertices.size(), vertex2));

			}

		}
		

	}




	//we are going to iterate until find a solution,  the first one it's the optimal in B&B, or not get any solution
	while (!allPossibleSolution.empty())
	{
		//we take the tram with less weight(DijkstraDistance) and we compare if the las node in this is the destination
		CBBNode1* topNode = allPossibleSolution.top();
		allPossibleSolution.pop();

		//we mark the next vertex in the list of index as visited
		topNode->m_visited[topNode->m_indexes.back()] = true;

		
		if (topNode->m_indexes.size() == 9 && topNode->m_indexes == correctSolution)
		{
			cout << "We got it" << endl;
		}
		

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
				//we control that we are not adding a solution that goes from v1 -> v0 -> v1 -> v0 again, avoid loop solution
				if (col != topNode->m_indexes.back() && col != topNode->m_indexes[topNode->m_indexes.size()-2] && topNode->m_indexes[topNode->m_indexes.size()-1] != topNode->m_indexes[topNode->m_indexes.size() - 3])
				{
					
					CBBNode1* newNode = new CBBNode1(*topNode);


					newNode->m_indexes.push_back(col);

					
					if (newNode->m_indexes.size() == 9)
					{
						if(newNode->m_indexes == correctSolution)
							cout << "Solution created" << endl;
					}
					
					newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

					allPossibleSolution.push(newNode);
				}

			}
		}


	}

	return CTrack(&graph);
	
}

*/


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
				if (col != topNode->m_indexes.back() && col != topNode->m_indexes[0])
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

vector<int> correctSolution = { 0,1,2,5,3,4,2,6,7 };

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
	double cotaSuperior = accumulate(max.begin(), max.end(), 0.0);
	
	double globalCotaSuperior = cotaSuperior;

	//we explore the first node and we add ways, choosign with the priority queue the next one that we are going to explore, it dependes of DijkstraDistance to each one
	//but with don't add the cell v1 column v1 because the dijkstraDistance is 0 and that generate a loop adn we are not going to explore the destination node so size-1
	for (int col = 1; col < visits.m_Vertices.size()-1; col++) {

		//we create te vector where we are going to seave the indexes of th tram
		vector<int> vIndex;
		vIndex.push_back(0);
		vIndex.push_back(col);
		
		allPossibleSolution.push(new CBBNode1(matriu[0][col].dijkstraDistance, vIndex, visits.m_Vertices.size(), cotaInferior, cotaSuperior));

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

			//we acces to the las CVertex in the current solution
			auto lastNode= visits.m_Vertices.begin();
			advance(lastNode, topNode->m_indexes.back());

			
			for (int col = 0; col < visits.m_Vertices.size(); col++)
			{
				
				//we control that we don't go to the same vertex, if we don't put this conditional, we will always in the same vertex because Dijkstra distance to it is 0
				if (col != topNode->m_indexes.back() && col != topNode->m_indexes[0])
				{

					if (topNode->m_indexes.size() >= 3) {
						
						if (col == topNode->m_indexes[topNode->m_indexes.size() - 2] && topNode->m_indexes[topNode->m_indexes.size() - 1] == topNode->m_indexes[topNode->m_indexes.size() - 3])
							goto avoidLoop;
						
						
					}


					//Acces to the Cvertex we try to add
					auto it = visits.m_Vertices.begin();
					advance(it, col);

					//calculation of the new cotaSuperior to knwo if it is interesting to continous exploring this solutions
					cotaInferior = topNode->m_cotaInferior - min[(*it)->m_indexMatrix] + (*lastNode)->m_Point.Distance((*it)->m_Point);

					if (globalCotaSuperior > cotaInferior)
					{
						//CBBNode(double* dijkstraDistance, vector<int> indexes)
						//CBBNode* newNode(topNode);


						//cota inferior calculation
						cotaSuperior = topNode->m_cotaSuperior - max[(*it)->m_indexMatrix] + (*lastNode)->m_Point.Distance((*it)->m_Point);

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

						if (newNode->m_indexes.size() == 9)
						{
							if (newNode->m_indexes == correctSolution)
								cout << "Solution created" << endl;
						}

						newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][col].dijkstraDistance;

						allPossibleSolution.push(newNode);
					}
					
					
				}
				avoidLoop:;
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
	double cotaSuperior = accumulate(max.begin(), max.end(), 0.0);

	double globalCotaSuperior = cotaSuperior;

	//we explore the first node and we add ways, choosign with the priority queue the next one that we are going to explore, it dependes of DijkstraDistance to each one
	//but with don't add the cell v1 column v1 because the dijkstraDistance is 0 and that generate a loop adn we are not going to explore the destination node so size-1
	for (int col = 1; col < visits.m_Vertices.size() - 1; col++) {

		//we create te vector where we are going to seave the indexes of th tram
		vector<int> vIndex;
		vIndex.push_back(0);
		vIndex.push_back(col);

		allPossibleSolution.push(new CBBNode1(matriu[0][col].dijkstraDistance, vIndex, visits.m_Vertices.size(), cotaInferior, cotaSuperior));

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

			/*
			vector<int>::iterator ip;


			vector<int> tmp=topNode->m_indexes;

			// Using std::unique
			ip = std::unique(tmp.begin(), tmp.begin() + tmp.size());
			// Resizing the vector so as to remove the undefined terms
			tmp.resize(distance(tmp.begin(), ip));


			
			*/

			min.clear();
			max.clear();
			min.resize(visits.m_Vertices.size(), numeric_limits<double>::max());
			max.resize(visits.m_Vertices.size(), 0.0);
			

			for (int row = 0; row < visits.m_Vertices.size() - 1; row++)
			{
				//It's not consider the rows of visited vertex
				if (!topNode->m_visited[row]) {

					//we check if we find another max or min DijkstraDistance in the new row to nodes
					for (int col = 0; col < visits.m_Vertices.size(); col++)
					{
						//It is not considere the column of the last vertex with visited, It is not consider going from last vertex to the destination vertex
						if (col != topNode->m_indexes.back() && (row != topNode->m_indexes.back() || col != visits.m_Vertices.size()-1)) {


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
				}

			}
			//we acces to the las CVertex in the current solution
			auto lastNode = visits.m_Vertices.begin();
			advance(lastNode, topNode->m_indexes.back());


			for (int newVertex = 0; newVertex < visits.m_Vertices.size(); newVertex++)
			{



				//we control that we don't go to the same vertex, if we don't put this conditional, we will always in the same vertex because Dijkstra distance to it is 0
				if (newVertex != topNode->m_indexes.back() && newVertex != topNode->m_indexes[0])
				{

					if (topNode->m_indexes.size() >= 3) {

						if (newVertex == topNode->m_indexes[topNode->m_indexes.size() - 2] && topNode->m_indexes[topNode->m_indexes.size() - 1] == topNode->m_indexes[topNode->m_indexes.size() - 3])
							goto avoidLoop;


					}


					//Acces to the Cvertex we try to add
					auto it = visits.m_Vertices.begin();
					advance(it, newVertex);

					//calculation of the new cotaSuperior to knwo if it is interesting to continous exploring this solutions
					cotaInferior = topNode->m_cotaInferior - min[(*it)->m_indexMatrix] + (*lastNode)->m_Point.Distance((*it)->m_Point);

					if (globalCotaSuperior > cotaInferior)
					{
						//CBBNode(double* dijkstraDistance, vector<int> indexes)
						//CBBNode* newNode(topNode);


						//cota inferior calculation
						cotaSuperior = topNode->m_cotaSuperior - max[(*it)->m_indexMatrix] + (*lastNode)->m_Point.Distance((*it)->m_Point);

						//when it is found a smaller cotaSuperior we take it as the new global
						if (globalCotaSuperior > cotaSuperior)
						{
							globalCotaSuperior = cotaSuperior;
						}

						//Creation of the new solution or node
						CBBNode1* newNode = new CBBNode1(*topNode);
						newNode->m_cotaInferior = cotaInferior;
						newNode->m_cotaSuperior = cotaSuperior;

						newNode->m_indexes.push_back(newVertex);

						if (newNode->m_indexes.size() == 9)
						{
							if (newNode->m_indexes == correctSolution)
								cout << "Solution created" << endl;
						}

						newNode->m_dijkstraLenght = topNode->m_dijkstraLenght + matriu[topNode->m_indexes.back()][newVertex].dijkstraDistance;

						allPossibleSolution.push(newNode);
					}


				}
			avoidLoop:;
			}
		}


	}

	return CTrack(&graph);
}
