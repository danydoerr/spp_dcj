#ifndef RECONCILED_TREE_H_
#define RECONCILED_TREE_H_

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

This file contains a class for reconciled trees

Created the: 26-10-2015
by: Wandrille Duchemin

Last modified the: 08-09-2016
by: Wandrille Duchemin

*/

#include <string>
#include <vector>
#include <map>

#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Exceptions.h>
#include <fstream>

#include "MySpeciesTree.h"
#include "MyGeneTree.h"
#include "CladeReconciliation.h" //object representing the reconciliation of a clade split as outputted by TERA
#include "ReconciliationEvent.h" //object representing and event of reconciliation as outputted by TERA

#include "XMLUtils.h"

using namespace bpp;

//Variables used to get and set some node properties
const string spe = "S";  //name of the node property containing their species
const string ev = "Ev";  //name of the node property containing their event
const string clnum = "Cnum";//name of the node property containing their clade id
const string porder = "Po";//name of the node property containing their post order
//const string alpha = "Al";//name of the node property containing their alpha status (whether they are in a dead lineage or no)
const string ts = "TS";//name of the node property containing their time slice
const string uts = "UTS";//upper time slice boundary
const string lts = "LTS";//lower time slice boundary


//Code of the reconciliation events

const int C		= 0; //Current (leaf)
const int S		= 1; //Speciation (in the species tree)
const int L		= 2; //Loss (in the species tree)
const int D		= 3; //Duplication (in the species tree)
const int Sout	= 4; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
const int R		= 5; //Transfer reception
const int N		= 6; //no event (to account for time slices)
const int Bout	= 7; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)


class ReconciledTree: public TreeTemplate<Node>
{
protected:

	map<int,vector<int> > CladeIdToNodeIds;
	//map<int, Node * > NodeIdToNodeP;


	int nbNodes;
	int MyNextId()
	{
		return nbNodes;
	}

	int TimeSliceStatus; // 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)


	vector <string> cladeToLeafListAssociationAux(int currentNode, map <int, vector <string> > & cladeToLeafList );

	MyGeneNode * getMyGeneTreeTopologyAux(Node * currentNode);

	int countAnyEvent( int evttype);


	/*
	Adds the reconciliation to the given NodeId associated with idU.
	This function additionally creates and annotates a loss node when necessary (dependent on the event).

	Takes:
	 - RE (ReconciliationEvent) : a reconciliation event specifying the species it occurs in and the type of event
	 - NodeId (int) : node id
	 - VERBOSE (bool) (default = false)

	Returns:
	 (int) node id to branch the next nodes to
	*/
	int addReconciliationEvent(ReconciliationEvent RE, int NodeId, Node * currentNode, bool VERBOSE = false);


	/*
	This function look for T and TFD events in order to add the correct Reception event to the tree.
	If warranted, it annotates NodeId as a Reception and create its son.
	
	Takes:
	 - PRE (ReconciliationEvent) : the previous reconciliation event undergone to arrive to NodeId
	 - NodeId (int) : node id
	 - VERBOSE (bool) (default = false)
	
	Returns:
	 (int) node id to branch the next nodes to
	
	*/
	int accountForPreviousSplit(ReconciliationEvent PRE, int NodeId,bool VERBOSE=false);

	/*
	Creates or annotate a leaf node at NodeId
	
	Checks if he parent is a transfer (Sout or Bout). If it is, the function will annotate NodeId as a reception node and create a leaf node.
	Otherwise just annotate the node as a leaf.
	
	Takes:
	 - NodeId (int): id of the node to annotate with species and clade already set
	
	 */
	void CreateLeafNode(int NodeId);


	/*
	add the reconciliation and return idUl and idUr or -1,-1 if idU is a leaf
	
	Takes:
		- CR (CladeReconciliation) : a clade reconciliation of id idU
		- VERBOSE (bool) (default = false)

	Returns:
		pair<int,int> : idUl and idUr or -1,-1 if idU is a leaf

	*/
	pair<int,int> addCladeReconciliation(CladeReconciliation CR, bool VERBOSE=false);


	/*
	Takes:
		- Nodeid (int): a node id
		- evtLine (string): recPhyloXML line describing a reconciliation event
	*/
	void setNodeDetailsFromXMLLine(int Nodeid, string evtLine, map< string , int > & speciesNameToNodeId  , bool VERBOSE);

	//same as previous; uses the interpreted xml line instead
	void setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value, map< string , int > & speciesNameToNodeId  , bool VERBOSE);

	/*
	read and add a clade from a recPhyloXML stream
	
	Takes:
		- fileIN (ifstream): streaam to a recPhyloXML file
		- VERBOSE (bool) (default = false)

	*/
	void addXMLClade(ifstream& fileIN, int Nodeid, map< string , int > & speciesNameToNodeId  , bool VERBOSE = false); // will need a species tree


	void readPhyloXMLFile(ifstream& fileIN, MySpeciesTree * Stree, map< string , int > & speciesNameToNodeId , bool VERBOSE );

	void removeNoEventNodes();

	void PreComputeBoundedTimeSlices(map <int,int> idXToLTS);
	/*
	Recursively sets the upper and lower boundaries of a subtree.
	Takes:
		- NodeId (int) : id of the root of the subtree
		- idXToUTS (map <int,int>) : map of species id to maximum associated time slice
		- idXToLTS (map <int,int>) : map of species id to minimum associated time slice
	*/
	void setSubtreeBoundaryTS(int NodeId, 	map <int,int> idXToUTS, map <int,int> idXToLTS);

	// used to fix the holes than can be left after all null nodes are deleted in BTS mode
	void addIntermediaryNullBTS(int NodeId);


	string NodeString(int nodeid, bool hideLosses = false);

public:
	ReconciledTree(): TreeTemplate<Node>(), nbNodes(0)
	{} // empty tree

	ReconciledTree(Node * root, int TSstatus);

	/**
	 * Constructor that takes a vector of CladeReconciliation
	 */
	ReconciledTree(vector<CladeReconciliation> * CladeRecs,bool VERBOSE = false);


	//Constructor that uses a stream from a recPhyloXML file
	ReconciledTree(ifstream& fileIN, MySpeciesTree * Stree, map< string , int > & speciesNameToNodeId , bool VERBOSE = false);

	//Constructor that uses a recPhyloXML file name
	ReconciledTree(string phyloxmlFileName, MySpeciesTree * Stree , map< string , int > & speciesNameToNodeId, bool VERBOSE = false);


	~ReconciledTree()
	{
	CladeIdToNodeIds.clear();
	}//fairly simple destruction



	//access to NodeId with CladeId
	vector<int> getNodeIdsFromCladeId(int CladeId);
	void addNodeIdToCladeIdMap(int CladeId, int NodeId);


	int getTimeSliceStatus();

	//getter of node properties
	int getNodeSpecies(int nodeid); // retrieve the specified node species. 
	int getNodePostOrder(int nodeid);
	int getNodeTimeSlice(int nodeid);
	int getNodeUpperBoundaryTS(int nodeid);
	int getNodeLowerBoundaryTS(int nodeid);
	int getNodeEvent(int nodeid); // retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
	int getNodeCladeNum(int nodeid); // retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance

	//setter of node properties
	void setNodeSpecies(int nodeid, int speciesid); // set the specified node species.
	void setNodePostOrder(int nodeid, int postorder);
	void setNodeTimeSlice(int nodeid, int timeslice);
	void setNodeUpperBoundaryTS(int nodeid, int timeslice);
	void setNodeLowerBoundaryTS(int nodeid, int timeslice);
	void setNodeEvent(int nodeid, int eventid); // assigns an event to the node according to the eventid
	void setNodeCladeNum(int nodeid, int cladenum); // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance

	void setNodeSpecies(Node * NP, int speciesid); // set the specified node species.
	void setNodePostOrder(Node * NP, int postorder);
	void setNodeTimeSlice(Node * NP, int timeslice);
	void setNodeUpperBoundaryTS(Node * NP, int timeslice);
	void setNodeLowerBoundaryTS(Node * NP, int timeslice);
	void setNodeEvent(Node * NP, int eventid); // assigns an event to the node according to the eventid
	void setNodeCladeNum(Node * NP, int cladenum); // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance



	void resetNodeCladeNum(int nodeid);//unset the clade id of the node. Always use this function because it also updates the CladeIdToNodeIds map

	void printNode(int nodeid);//prints node information
	void printMe();//prints all nodes
	string NewickString(bool hideLosses = false);


	void computeBoundedTimeSlices(MySpeciesTree * SpeciesTree);
	void setToNonTimeSliced();
	void SubdivideTree();//expects that node already have an exact set time slice

	int getNumberOfDuplication();
	int getNumberOfLoss();
	int getNumberOfTransfer();


	bool haveSameSpecies(Node  * n1, Node * n2);// just check species
	bool areTSCompatible(Node  * n1, Node * n2);// just check TS
	bool areCompatible(Node  * n1, Node * n2); //check if the node are compatible (ie. same species and comparable timeslice). Presumes that both nodes have the same timeslice status as the tree
	bool areCompatible(int id1, int id2); //check if the node are compatible (ie. same species and comparable timeslice). Presumes that both nodes have the same timeslice status as the tree


	bool isAncestor(int AncId, int SonId);//Checks if AncId is the ancestor of SonId
	vector <Node *> pathFromAncestorToSon(int AncId, int SonId); // Get the path from ancestor to son.
	vector <int> idPathFromAncestorToSon(int AncId, int SonId); // Get the path from ancestor to son.

	int getIdWithName(const string &name);


	bool isRealLeaf(int id);
	bool isExtant(int evtcode);
	bool isSpeciation(int evtcode);
	bool isLoss(int evtcode);
	bool isDup(int evtcode);
	bool isSout(int evtcode);
	bool isRec(int evtcode);
	bool isNull(int evtcode);
	bool isBout(int evtcode);

	ReconciledTree * cloneSubtree(int newRootId);

	vector <string> getRealLeavesNames();

	MyGeneTree getMyGeneTreeTopology();

	map <int, vector <string> > cladeToLeafListAssociation();

	vector <int> getLeafOrSpeciationDescendants(int nodeId);

	void makeSingleLeafTree( string name, int speciesId );

};
//the post order attribute is not really used. Why not discard it?
//overload my node accesser?



#endif /*RECONCILED_TREE_H_*/