#ifndef ADJ_TREE_H_
#define ADJ_TREE_H_

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

This file contains a class for an adjacency history as a tree

Created the: 14-12-2015
by: Wandrille Duchemin

Last modified the: 16-06-2016
by: Wandrille Duchemin

*/

#include <string>

#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Exceptions.h>
//#include <Bpp/BppBoolean.h>


#include "ReconciledTree.h" // for the common event code


#include "XMLUtils.h" // from AdjRecphyloXML reading

using namespace bpp;
using namespace std;

//Variables used to get and set some node properties
//const string spe = "S";  //name of the node property containing their species
const string nodeid1 = "nid1"; // name of the node property contnaining the node id in Gfamily1
const string nodeid2 = "nid2"; // name of the node property contnaining the node id in Gfamily2
//const string ev = "Ev";  //name of the node property containing their event
const string coev = "coev";//name of the node property containing the coevolution status of the event


//Code of the reconciliation events

/*const int C		= 0; //Current (leaf)
const int S		= 1; //Speciation (in the species tree)
const int L		= 2; //Loss (in the species tree)
const int D		= 3; //Duplication (in the species tree)
const int Sout	= 4; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
const int R		= 5; //Transfer reception
//const int N		= 6; //no event (to account for time slices)
const int Bout	= 7; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
*/
const int Abreak = 8;//adjacency break




class AdjTree: public TreeTemplate<Node>
{
protected:
	int Gfamily1;
	int Gfamily2;
	bool RootIsGain;

	pair <string, string> leafNamesfromAname(string Aname);

	vector <string> getC1LeafListAux( int currentId, int gfamIndex, vector < vector <string> > & LeafList );
	

	int countNbBreakAux(Node * n);
	int countNbGainAux(Node * n);

	void printNode(Node * n);
	void printNodeRec(Node * n,int level);

	string NodeString(Node * n, bool hideLosses = false);

	void addABreakInDeadNode(int Nodeid);


	/*
	Takes:
		- Nodeid (int): a node id
		- evtLine (string): recPhyloXML line describing a reconciliation event
	*/
	void setNodeDetailsFromXMLLine(int Nodeid, string evtLine, bool VERBOSE);

	//same as previous; uses the interpreted xml line instead
	void setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value, bool VERBOSE);

	/*
	read and add a clade from a recPhyloXML stream
	
	Takes:
		- fileIN (ifstream): streaam to a recPhyloXML file
		- VERBOSE (bool) (default = false)

	*/
	void addXMLClade(ifstream& fileIN, int Nodeid, bool VERBOSE = false); // will need a species tree

	void readPhyloXMLFile(ifstream& fileIN, bool VERBOSE );


	void setNodeNamesAux(int sens1 , int sens2, Node * current);

public:

	AdjTree(){}
	AdjTree(int gfam1, int gfam2): Gfamily1(gfam1) , Gfamily2(gfam2)
	{}

	AdjTree(int gfam1, int gfam2, bool rootgain): Gfamily1(gfam1) , Gfamily2(gfam2), RootIsGain(rootgain)
	{}

	AdjTree(bool rootgain): RootIsGain(rootgain)
	{}

	AdjTree(Node * root, int rootgain, int gfam1, int gfam2);

	AdjTree(ifstream& fileIN,bool rootgain, int gfam1, int gfam2, bool VERBOSE = false);


	~AdjTree(){}


	int getGfamily1();
	int getGfamily2();
	bool isRootGain();

	void setGfamily1(int newval);
	void setGfamily2(int newval);
	void setRootGain(bool newval);
	

	int getNodeSpecies(int id);
	pair <int,int> getNodeNodeIds(int id);
	bool getNodeCoEventStatus(int id);
	int getNodeEvent(int id);

	int getNodeSpecies(Node * n);
	pair <int,int> getNodeNodeIds(Node * n);
	bool getNodeCoEventStatus(Node * n);
	int getNodeEvent(Node * n);


	void setNodeSpecies(int id, int sp);
	void setNodeNodeIds(int id, pair <int,int> nodeids);
	void setNodeCoEventStatus(int id, bool co);
	void setNodeEvent(int id, int e);

	Node * createNode(pair <int,int> nodeids, int sp, bool co, int e);
	Node * createBreakNode();
	
	void refine(ReconciledTree * Rtree1, ReconciledTree * Rtree2, int sens1 = 0, int sens2 = 0);
	void setNodeNames(int sens1 , int sens2);

	void printMe();
	string NewickString(bool hideLosses = false);


	bool isExtant(int evtcode);
	bool isSpeciation(int evtcode);
	bool isLoss(int evtcode);
	bool isDup(int evtcode);
	bool isSout(int evtcode);
	bool isRec(int evtcode);
	bool isNull(int evtcode);
	bool isBout(int evtcode);
	bool isAdjBreak(int evtcode);

	bool isLivingAdjBreak(int Nodeid, int evtcode = -1);


	int countNbBreak();
	int countNbGain();

	AdjTree * cloneSubtree(int newRootId);

	void getC1LeafList( int gfamIndex, vector < vector <string> > & LeafListList );
	vector <int> getC1NodeIds(int gfamIndex);
};


#endif