/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/*

This file contains a class for co events -> events that span several nodes of gene tree(s)

Created the: 25-01-2016
by: Wandrille Duchemin

Last modified the: 07-07-2016
by: Wandrille Duchemin

*/

#include "CoEvent.h"

/*
Takes:
	- id (int): an index in GeneList

Returns:
	(bool): true if there is at leats an adjacency with the gene
*/
bool CoEvent::HasNeighbors(int id)
{
	for(unsigned i = 0; i < AdjList.size(); i++)
	{
		if(AdjList[i].n1 == id)
			return true;
		else if(AdjList[i].n2 == id)
			return true; 
	}
	return false;
}


/*
Takes:
	- SourceId (int) : id of the source for the connex component

Returns:
	vector <int> list if the ids of the connex component that contains SourceId

*/
vector < int > CoEvent::exploreOneConnexComponent(int SourceId)
{
	//cout << "exploring a connex component from " << SourceId << endl;
	//1. initialize a bool array to keep track of the nodes we've seen already

	bool crible [GeneList.size()];
	for(unsigned i = 0 ; i < GeneList.size(); i++)
		crible[i] = false;

	//2. initialize exporation
	vector <int> ConnexComponent;
	vector <int> CurrentNodes;

	crible[SourceId] = false;
	ConnexComponent.push_back(SourceId);
	CurrentNodes.push_back(SourceId);

	while(CurrentNodes.size() > 0)
	{
		int currentId = CurrentNodes.back();
		CurrentNodes.pop_back();

		vector <int> neighbors = getGeneNeighbors(currentId);
		for(unsigned j=  0; j < neighbors.size(); j++)
		{
			if(!crible[neighbors[j]]) // first time we see that node
			{
				crible[neighbors[j]] = true;
				ConnexComponent.push_back(neighbors[j]);
				CurrentNodes.push_back(neighbors[j]);
			}
		}
		//cout << "current component: " << ConnexComponent.size() << " todo: " << CurrentNodes.size() << endl;
	}

	return ConnexComponent;
}

/*
Takes:
	- id (int): an index in GeneList

Returns:
	(int): number of gene that share an adjacency with the gene (actually, number of adjacency where the gene appears)
*/
int CoEvent::getGeneNumberOfNeighbors(int id)
{
	int c = 0;
	for(unsigned i = 0; i < AdjList.size(); i++)
	{
		if(AdjList[i].n1 == id)
			c += 1;
		else if(AdjList[i].n2 == id)
			c += 1;
	}
	return c;
}


/*
Returns:
	(int): number of genes composing that coevent
*/
int CoEvent::getNumberOfGene()
{
	return GeneList.size();
}

/*
Returns:
	(int): number of adjacencies composing that coevent
*/
int CoEvent::getNumberOfAdj()
{
	return AdjList.size();
}


/*
Takes:
	- i (int): an index in GeneList

Returns:
	(Gene): the gene at index i (thows error if the index is invalid)
*/
Gene CoEvent::getGene(int i)
{
	if((i < 0) || (i >= GeneList.size()))
		throw Exception("CoEvent::getGene : gave invalid index.");
	return GeneList[i];
}

/*
Takes:
	- i (int): an index in AdjList

Returns:
	(AdjNode): the adjacency at index i (thows error if the index is invalid)
*/
AdjNode CoEvent::getAdj(int i)
{
	if((i < 0) || (i >= AdjList.size()))
		throw Exception("CoEvent::getGene : gave invalid index.");
	return AdjList[i];	
}

/*
Takes:
	- gfam (int): Gene family index
	- nodeid (int): id the gene family's reconciled tree

Returns:
	(int): index of the Gene in GeneList. -1 if the gene is absent.
*/
int CoEvent::getGeneIndex(int gfam, int nodeid)
{
	for(unsigned i =0; i < GeneList.size(); i++)
	{
		if( GeneList[i].Gfam == gfam)
		{
			if(GeneList[i].NodeId == nodeid)
			{
				return i;
			}
		}
	}
	return -1;
}

/*
Takes:
	- eqclassid (int): equivalence class index
	- nodeid (int): id the equivalence class's adjacency tree

Returns:
	(int): index of the AdjNode in AdjList. -1 if the gene is absent.

*/
int CoEvent::getAdjIndex(int eqclassid, int nodeid)
{
	for(unsigned i =0; i < AdjList.size(); i++)
	{
		if( AdjList[i].EqClassId == eqclassid)
		{
			if(AdjList[i].NodeId == nodeid)
			{
				return i;
			}
		}
	}
	return -1;	
}

//adding Gene or Adj

/*
Adds a gene to the CoEvent,
also checks if the gene already exists.

Takes:
	- Gfam (int): id of the gene family of the gene
	- NodeId (int): id of the node in the reconciled tree

Returns:
	(int): the index of the Gene in GeneList
*/
int CoEvent::addGene(int Gfam, int NodeId)
{

	int ind = getGeneIndex(Gfam, NodeId);
	if( ind != -1)//the gene already exists -> no need to add it
		return ind;


	//adding element to gene list
	GeneList.push_back(Gene());

	GeneList.back().Gfam = Gfam;
	GeneList.back().NodeId = NodeId;
	GeneList.back().UTS = -1;
	GeneList.back().LTS = -1;

	//grObs.createNode(&GeneList.back()); // creating the node in the associated graph

	return GeneList.size() -1 ;
}


/*
Adds a gene to the CoEvent,
also checks if the gene already exists and sets the UTS/LTS.

Takes:
	- Gfam (int): id of the gene family of the gene
	- NodeId (int): id of the node in the reconciled tree
	- Rtree (ReconciledTree *): pointer to the reconciled tree

Returns:
	(int): the index of the Gene in GeneList
*/
int CoEvent::addGene(int Gfam, int NodeId, ReconciledTree * Rtree)
{
	int ind = addGene(Gfam,NodeId); //creating the Gene and the node -> getting its index in GeneList

	int uts;
	int lts;

	int TSstatus = Rtree->getTimeSliceStatus(); // 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	if(TSstatus == 0) //no time slice -> no need to update them
		return ind;
	else if(TSstatus == 1)// precise time slice -> uts == lts
	{
		uts = Rtree->getNodeUpperBoundaryTS(NodeId);
		lts = uts;
	}
	else //bounded time slices
	{
		uts = Rtree->getNodeUpperBoundaryTS(NodeId);
		lts = Rtree->getNodeLowerBoundaryTS(NodeId);
	}

	GeneList[ind].UTS = uts;
	GeneList[ind].LTS = lts;


	if( GeneList.size() == 1)//only one Gene (this one) -> directly set the TSs
	{
		UpperTS = uts;
		LowerTS = lts;
	}
	else //checking and updating boundaries
	{
		if(UpperTS > uts)
			UpperTS = uts;
		if(LowerTS < lts)
			LowerTS = lts;
	}

	return ind;
}

/*
also checks if the AdjNode already exists, as well as the Genes. 

Takes:
	- EqClassId (int): id of the Equivalence Class of the adjacency
	- NodeId (int): id of the node of the adjacency in the adjacency tree of the equivalence class
	- Atree (AdjTree *): pointer to the adjacency tree where the AdjNode is 
	- Gfam1 (int): Gene family id of the first Gene of the Adjacency
	- Gfam2 (int): Gene family id of the second Gene of the Adjacency

Returns:
	(int): index of the adjacency in AdjList
*/
int CoEvent::addAdj(int EqClassId, int NodeId, AdjTree * Atree, int Gfam1, int Gfam2)
{
	//cout << "CoEvent::addAdj " << EqClassId << " "<< NodeId << endl;
	//int ind = getAdjIndex(EqClassId,NodeId);
	//cout << ind << endl; 
	//if( ind != -1) // the Adjacency already exists in this CoEvent -> nothing to do
	//	return ind;
	int ind = -1;
	//1. add the genes
	pair <int,int> gnodes = Atree->getNodeNodeIds(NodeId);

	int n1 = addGene( Gfam1 , gnodes.first );
	int n2 = addGene( Gfam2 , gnodes.second );

	//2. add the adjacency
	AdjList.push_back(AdjNode());

	AdjList.back().EqClassId = EqClassId;
	AdjList.back().NodeId = NodeId;	
	AdjList.back().n1 = n1;
	AdjList.back().n2 = n2;

	ind = AdjList.size() - 1 ;

	//3. set the species and event if this is the first adj
	if( AdjList.size() == 1)
	{
		species = Atree->getNodeSpecies(NodeId);
		event = Atree->getNodeEvent(NodeId);
	}
	return ind;
}


/*
Adds an adjacency to the CoEvent.
also checks if the AdjNode already exists, as well as the Genes. Sets the UTS/LTS 

Takes:
	- EqClassId (int): id of the Equivalence Class of the adjacency
	- NodeId (int): id of the node of the adjacency in the adjacency tree of the equivalence class
	- Atree (AdjTree *): pointer to the adjacency tree where the AdjNode is 
	- Gfam1 (int): Gene family id of the first Gene of the Adjacency
	- Gfam2 (int): Gene family id of the second Gene of the Adjacency
	- Rtree1 (ReconciledTree *): pointer to the reconciled tree of the first half of the adjacency
	- Rtree2 (ReconciledTree *): pointer to the reconciled tree of the second half of the adjacency

Returns:
	(int): index of the adjacency in AdjList
*/
int CoEvent::addAdj(int EqClassId, int NodeId, AdjTree * Atree, int Gfam1, int Gfam2, ReconciledTree * Rtree1, ReconciledTree * Rtree2)
{
	//cout << "CoEvent::addAdj " << EqClassId << " "<< NodeId << endl;
	//int ind = getAdjIndex(EqClassId,NodeId);
	//cout << ind << endl; 
	//if( ind != -1) // the Adjacency already exists in this CoEvent -> nothing to do
	//	return ind;
	int ind = -1;


	//1. add the genes
	pair <int,int> gnodes = Atree->getNodeNodeIds(NodeId);

	int n1 = addGene(Gfam1,gnodes.first,Rtree1);
	int n2 = addGene(Gfam2,gnodes.second,Rtree2);

	//2. add the adjacency
	AdjList.push_back(AdjNode());

	AdjList.back().EqClassId = EqClassId;
	AdjList.back().NodeId = NodeId;	
	AdjList.back().n1 = n1;
	AdjList.back().n2 = n2;

	ind = AdjList.size() - 1 ;

	//grObs.link(&GeneList[n1], &GeneList[n2], &AdjList[ind]);//creating the edge corresponding to that adjacency

	//3. set the species and event if this is the first adj
	if( AdjList.size() == 1)
	{
		species = Atree->getNodeSpecies(NodeId);
		event = Atree->getNodeEvent(NodeId);
	}
	return ind;
}

//removing gene or Adj

/*
removes Gene at index i

Takes:
	- i (int): index of the node to delete in GeneList
*/
void CoEvent::removeGene(int i)
{
	if(( i< 0) || (i > GeneList.size()-1))//invalid id -> does nothing
		return;

	//grObs.deleteNode(&GeneList[i]);

	GeneList.erase(GeneList.begin()+i);	

}

//void removeGeneAndLinkedAdj(int i); // removes Gene at index i as well as all NodeAdj implicating this Gene; Also updates UTS and LTS of the CoEvent

/*
removes Adj at index i

Takes:
	- i (int): index of the node to delete in AdjList
*/
void CoEvent::removeAdj(int i)
{
	if(( i< 0) || (i > GeneList.size()-1))//invalid id -> does nothing
		return;

	int n1 = AdjList[i].n1;
	int n2 = AdjList[i].n2;

	//grObs.unlink(&GeneList[n1],&GeneList[n2]); // removing from the graph

	AdjList.erase(AdjList.begin()+i); // removing from the list
}

/*
Removes NodeAdj at index i as well as all Gene that are not participating to another edge in the CoEvent.
Also updates UTS and LTS of the CoEvent

Takes:
	- i (int): index of the node to delete in AdjList
*/
void CoEvent::removeAdjAndLinkedGene(int i)
{

	if(( i< 0) || (i > GeneList.size()-1))//invalid id -> does nothing
		return;

	int n1 = AdjList[i].n1;
	int n2 = AdjList[i].n2;

	//grObs.unlink(&GeneList[n1],&GeneList[n2]); // removing from the graph

	AdjList.erase(AdjList.begin()+i); // removing from the list

	if( n2 > n1) //making sure we delete the element with the biggest index first
	{
		int tmp = n1;
		n1 = n2;
		n2 = n1;
	}

	bool removedgene = false;
	if( !HasNeighbors(n1))
	{
		removeGene(n1);
		removedgene = true;
	}

	if( !HasNeighbors(n2) )
	{
		removeGene(n2);
		removedgene = true;
	}

	if(removedgene)
		updateTimeSlices();
}

/*
Updates the upper and lower time slices of the coevent
*/
void CoEvent::updateTimeSlices()
{
	if(UpperTS == -1) // this means that time slices are not active in that CoEvent so we ignore the command
		return;

	if(GeneList.size() == 0) // no gene...
		return;

	//reinitializing the UTS and LTS
	UpperTS = GeneList[0].UTS;
	LowerTS = GeneList[0].LTS;

	for(unsigned i = 1 ; i <  GeneList.size() ; i++)
	{
		if(GeneList[i].UTS < UpperTS)
			UpperTS = GeneList[i].UTS;
		if(GeneList[i].LTS < LowerTS)
			LowerTS = GeneList[i].LTS;
	}
}


bool CoEvent::isCompatible(int NodeId, AdjTree * Atree)
{

	if(GeneList.size() == 0)
		return true; // always compatible if the coevent is empty

	if(Atree->getNodeSpecies(NodeId) != species)
		return false;

	if(Atree->getNodeEvent(NodeId) != event)
		return false;

	return true;
}

bool CoEvent::isCompatible(int NodeId, AdjTree * Atree, int UTS, int LTS)
{
	if(! isCompatible(NodeId,Atree))
		return false;

	if(UpperTS != -1)
		if(LTS > UpperTS)
			return false;

	if(LowerTS != -1)
		if(UTS < LowerTS)
			return false;

	return true;
}

/*
Returns:
	(int) : number of connex components (actual number of coevents)
*/
int CoEvent::getNbConnexComponents()
{
	if(nbConnexComponents == -1)
	{
		//cout << "in CoEvent::getLinearScore -> looking for connex"<<endl;
		vector < vector<int> > connex = getConnexComponents();
		nbConnexComponents = connex.size();
	}
	return nbConnexComponents;
}

/*
Takes:
	- Dcost (double): cost of a single duplication
	- Tcost (double): cost of a single transfer
	- Lcost (double): cost of a single loss

Returns:
	(double): cost of the coevent (or rather cost to reimburse) equal to (1 - NumberOfGene) * EventCost
*/
double CoEvent::getLinearScore(double Dcost, double Tcost, double Lcost)
{
	//cout << "in CoEvent::getLinearScore"<<endl;

	//determining which score applies
	double s;

	if(event == D)
		s = Dcost;
	else if(event == L)
		s = Lcost;
	else if(event == R)
		s = Tcost;
	else
		return 0;

	if(nbConnexComponents == -1)
	{
		//cout << "in CoEvent::getLinearScore -> looking for connex"<<endl;
		vector < vector<int> > connex = getConnexComponents();
		nbConnexComponents = connex.size();
	}
	//cout << "coev lin score:" <<  GeneList.size() << "-"<< nbConnexComponents << endl;

	s *= - (double)(GeneList.size() - nbConnexComponents);

	return s;
}


/*
Takes:
	- Gid (int): index of a gene in GeneList

Returns:
	( vector<int> ): list of the neighbors of Gid's id
*/
vector <int> CoEvent::getGeneNeighbors(int Gid)
{
	if(( Gid < 0 ) || ( Gid >= GeneList.size() ))
		throw Exception("CoEvent::getGeneDegree : invalid id given.");

	vector<int> neighbors;
	for(unsigned i = 0; i <  AdjList.size(); i++)
	{
		if( AdjList[i].n1 == Gid )
			neighbors.push_back(AdjList[i].n2);
		else if( AdjList[i].n2 == Gid )
			neighbors.push_back(AdjList[i].n1);
	}
	return neighbors;
}

/*
Returns;
	( vector < vector <int> > ): each vector<int> represents a connex component in the graph formed by the co-event genes and adjacencies. The <int> are Gene index in GeneList
*/
vector < vector <int> > CoEvent::getConnexComponents()
{
	//cout << "in CoEvent::getConnexComponents" << endl;
	//1. initialize a bool array to keep track of the nodes we've seen already
	int nbGene = GeneList.size();
	bool crible [nbGene];

	for(unsigned i = 0 ; i < nbGene; i++)
		crible[i] = false;

	//cout << "crible initiated" << endl;


	vector < vector <int> > ConnexComponents;

	//2. explore the graph
	int i = 0;
	while( i < nbGene )
	{
		if(!crible[i]) // gene unknown
		{
			vector <int> connex = exploreOneConnexComponent(i);

			ConnexComponents.push_back(connex);

			for(unsigned j = 0; j < connex.size(); j++)
				crible[ connex[j] ] = true;
		}
		i++;
		//cout << "explored one conex component" << endl;
	}

	return ConnexComponents;
}

/*
Removes all Genes and adjs that sports GfamId.
May leaves some nodes hanging

Takes:
 - GfamId (int): the id of a GeneFamily

Returns:
	(bool): true if at least a gene / adj has been removed

*/
bool CoEvent::removeGeneFam(int GfamId)
{

	vector <int> removedId;

	size_t i = 0;
	while(i < GeneList.size())
	{
		if(GeneList[i].Gfam == GfamId)
			removedId.push_back(i);
		
		i++;
	}


	if(removedId.size() == 0)
		return false; //no gene with GfamId -> nothing to change here

	for(i=0 ; i < removedId.size() ; i++)
		removeGene(removedId[i]);

	i = 0; /// SO inefficient
	while(i < AdjList.size())
	{
		bool remove = false;
		for(size_t j = 0 ; j < removedId.size() ; j++)
		{	
			if( ( AdjList[i].n1 == removedId[j] ) || ( AdjList[i].n2 == removedId[j] ) )
			{
				remove = true;
				break;
			}
		}
		if(remove)
			removeAdj(i);
		else
			i++;
	}

	return true;
}

bool CoEvent::removeECF(int ECFId)
{
	bool changed = false;
	unsigned i = 0 ;
	while(i < AdjList.size())
	{
		if(AdjList[i].EqClassId == ECFId)
		{
			removeAdjAndLinkedGene(i);
			changed = true;
		}
		else
			i++;
	}
	return changed;
}