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

Last modified the: 16-03-2020
by: Wandrille Duchemin

*/



#include "ReconciledTree.h"
#include "ReconciliationEvent.h"


//Variables used to get and set some node properties
/*const string spe = "S";  //name of the node property containing their species
const string ev = "Ev";  //name of the node property containing their event
const string clnum = "Cnum";//name of the node property containing their clade id
const string porder = "Po";//name of the node property containing their post order
//const string alpha = "Al";//name of the node property containing their alpha status (whether they are in a dead lineage or no)
const string ts = "TS";//name of the node property containing their time slice
*/


/*
Takes:
	- int currentNode : a node in the reconciled tree
	- &map <int, vector <string> > cladeToLeafList : associates cladeId to list of leaf names

*/
vector <string> ReconciledTree::cladeToLeafListAssociationAux(int currentNode, map <int, vector <string> > & cladeToLeafList )
{
	vector <string> leafList;

	if(isRealLeaf(currentNode))
	{
		leafList.push_back(getNodeName(currentNode));
	}
	else
	{
		vector <int> SonsIds = getSonsId(currentNode);
		for(unsigned i = 0 ; i < SonsIds.size(); i++)
		{
			vector <string> tmp = cladeToLeafListAssociationAux(SonsIds[i], cladeToLeafList );
			for(unsigned j = 0 ; j < tmp.size(); j++)
				leafList.push_back(tmp[j]);
		}
	}


	if(leafList.size()>0)
	{
		int cladeId = getNodeCladeNum(currentNode);
	
		if(cladeToLeafList.find(cladeId) == cladeToLeafList.end())
		{
			cladeToLeafList[cladeId] = leafList;
		}
	}

	return leafList;
}




/*
Auxiliary recursive function to extract a MyGeneTree instance copying the topology of the reconciled tree.

Takes:
 - Node * currentNode : a pointer to the currenNode in the reconciled tree

Returns:
	MyGeneNode * : my gene node representation of currentNode, along with the underlying subtree
*/
MyGeneNode * ReconciledTree::getMyGeneTreeTopologyAux(Node * currentNode)
{
    MyGeneNode *node;

    vector<Node *> Sons = currentNode->getSons();
    int evt = dynamic_cast<BppInteger *> ( currentNode->getNodeProperty(ev) )->getValue()  ; 
    //cout << "ReconciledTree::getMyGeneTreeTopologyAux " << Sons.size() << ' ' << evt << endl;


    if(Sons.size() == 0) // no sons -> loss or leaf
    {
    	int evt = dynamic_cast<BppInteger *> ( currentNode->getNodeProperty(ev) )->getValue()  ; 
    	if(isExtant(evt)) // leaf -> create a node that has the name of the leaf
    	{
            node = new MyGeneNode();
            node->setName( currentNode->getName() );
    	}
    	else // loss -> return NULL that will signal its parent in order not to include branches leading to losses
    	{
    		return NULL;
    	}
    }
    else // tehre are some sons
    {
    	vector <MyGeneNode *> GeneSons;

    	for(unsigned i =0; i < Sons.size(); i++) // iterate over the sons
    	{
			MyGeneNode * placeholder = getMyGeneTreeTopologyAux( Sons[i] ); // getting that subtree

			if(placeholder != NULL) //if there is something else than a branch leading to a loss -> add it
			{
				GeneSons.push_back(placeholder);
			}
    	}

    	if(GeneSons.size() == 0) // no sons that don't lead to losses -> should not occur but if it does, return NULL
    	{
    		return NULL;
    	}
    	else if(GeneSons.size() == 1) // only one son -> we return it so Nodes with only one son are contracted and only the bifurcation + leaves are kept
    	{
    		node = GeneSons[0];
    	}
    	else //2 (or more but should not happen) sons
    	{
    		node = new MyGeneNode();

	    	for(unsigned i =0; i < GeneSons.size(); i++) // iterate over the MyGeneNode sons
    		{
				node->addSon(GeneSons[i]); // adding children	
	    	}
    	}
    }

    return node;
}




int ReconciledTree::countAnyEvent( int evttype)
{
	int counter = 0;
	vector<int> NodesIds  = getNodesId();
	int NodeId;

	for (unsigned j=0; j<NodesIds.size() ; j++)
	{
		NodeId = NodesIds[j];
		if(getNodeEvent(NodeId) == evttype)
			counter++;
	}
	return counter;
}


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
int ReconciledTree::addReconciliationEvent(ReconciliationEvent RE, int NodeId, Node * currentNode, bool VERBOSE)
{
	//here we define correspondance between TERA's event notation and our event notation


	bool loss = false;//is there a loss node to create
	bool Rtoadd = false;//is there a reception node with the same cladenum to add (only in the case of TL)

	bool TStoadd = RE.hasTimeSlice();
	string evtName = RE.getevtName();


	if(VERBOSE)
		cout << "ReconciledTree::addReconciliationEvent. NodeId: " << NodeId << " event: " << evtName << endl;

	if(evtName.length()>1)
		if(RE.getevtName()[1] == 'L') // there is a loss node if the second letter of the event name is L
			loss = true;

	int evtcode;
	if(evtName == "O")
	{
		evtcode = N; //nothing happened
	}
	else if(evtName == "C")
	{
		evtcode = C; //leaf
	} 
	else if(evtName[0] == 'S') //S and SL
		evtcode = S;
	else if(evtName[0] == 'D') //D and DD
	{
		if(evtName.length()>1)//DD
			evtcode = Bout;
		else
			evtcode = D;
	}
	else //T TL TTD TFD TLTD TLFD
	{
		if(evtName.length()==1)//T -> the reception Node will be handled in the next clade
			evtcode = Sout;
		else if(evtName.length()==2)//TL
		{
			evtcode = Sout;
			Rtoadd = true;
		}
		else//TTD TFD TLTD TLFD
		{
			//check the second to last letter to know if we are From the dead or To the dead
			if(evtName[evtName.length() - 2] == 'F')// TFD or TLFD
			{
				if(evtName.length() == 4)// TLFD -> the lost is implicit 
				{
					loss = false;
					evtcode = R; 
				}
				else //TFD -> this is actually a Bifurcation out (Bout) and the reception node (R) will be handled in the next clade
					evtcode = Bout;
			}
			else//TTD or TLTD
			{
				evtcode = Sout;
			}
		}
	}
	
	if(VERBOSE)
		cout << " -> assigned event code: " << evtcode << " loss? " << loss << endl;
	//which species is assigned to the node depends on the event code
	int speciesid;

	if(evtcode == Bout)
		speciesid = -1;
	else if(evtcode == R)
	{
		//in this case we chose the child species that is not -1
		if(RE.getidXl() == -1)
			speciesid = RE.getidXr();
		else
			speciesid = RE.getidXl();
	}
	else//default case: use idX
		speciesid = RE.getidX();


	//Node * currentNode = getNode(NodeId);

	//setting the species and evt
	setNodeSpecies(currentNode,speciesid);
	setNodeEvent(currentNode, evtcode);

	if(TStoadd)
	{
		setNodeUpperBoundaryTS(currentNode,RE.gettimeSlice());//we set the time slice
		setNodeLowerBoundaryTS(currentNode,RE.gettimeSlice());//the upper and lower boundaries are the same here
	}

	//adding the loss node
	if(loss)
	{
		//the speciesid stays the same if the event is not a SL
		if(evtcode == S)//SL case
			speciesid = RE.getidXl();

		evtcode = L;

		int newId = MyNextId();//getNextId();
		Node * newLossNode  = new Node(newId);

		nbNodes++;
		//Node * currentNode  = getNode(NodeId);
		currentNode->addSon(newLossNode);

		//setting new loss node properties
		setNodeSpecies(newLossNode,speciesid);
		setNodeEvent(newLossNode, evtcode);
		setNodeCladeNum(newLossNode, -1); // the clade num of the loss node is -1 -> the empty clade

		if(TStoadd)
		{
			setNodeUpperBoundaryTS(newLossNode,RE.gettimeSlice() -1 );//we set the time slice of the loss node as its father -1.
			setNodeLowerBoundaryTS(newLossNode,RE.gettimeSlice() -1 );//the upper and lower boundaries are the same here		
		}
	}

	//finally, we check if the event is TL. If this is the case, we will then add a reception node
	if(Rtoadd)
	{
		//the speciesid is the child that is not equal to its parent
		if(RE.getidXl() == speciesid)
			speciesid = RE.getidXr();
		else
			speciesid = RE.getidXl();

		evtcode = R;

		int newId = MyNextId(); //getNextId();
		Node * newRecNode  = new Node(newId);

		nbNodes++;
		//Node * currentNode  = getNode(NodeId);
		currentNode->addSon(newRecNode);

		//setting new reception node properties
		setNodeSpecies(newRecNode,speciesid);
		setNodeEvent(newRecNode, evtcode);


		setNodeCladeNum(newRecNode, getNodeCladeNum(NodeId) ); // the clade num of the rec node is the same as its parent

		if(TStoadd)
		{
			setNodeUpperBoundaryTS(newRecNode,RE.gettimeSlice());//we set the time slice of the loss as its father.
			setNodeLowerBoundaryTS(newRecNode,RE.gettimeSlice());//the upper and lower boundaries are the same here	
		}
		


		NodeId = newId;//set NodeId as the id of the Reception node that will serve to continue the lineage
	}

	return NodeId;
}


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
int ReconciledTree::accountForPreviousSplit(ReconciliationEvent PRE, int NodeId, bool VERBOSE)
{
	if((PRE.getevtName().compare("T") != 0)&&(PRE.getevtName().compare("TFD") != 0))//neither T or TFD -> no pre-treatment to do
		return NodeId;

	if(getNodeSpecies(NodeId) == PRE.getidX())//if we are of the the same species as the parent, no need for a reception node. This is used in the case of TFD only
		return NodeId;

	int evtcode = R;

	if(VERBOSE)
		cout << "ReconciledTree::accountForPreviousSplit. NodeId: " << NodeId << " PRE: " << PRE.getevtName() << " " << PRE.getidX() << ":" << PRE.getidXl() << "," << PRE.getidXr() << endl;

	//we want the child species id that is different to the species id of its father
	int speciesid;
	if(PRE.getidX() == PRE.getidXl())
		speciesid = PRE.getidXr();
	else
		speciesid = PRE.getidXl();

	setNodeSpecies(NodeId,speciesid); // this may have already been set
	setNodeEvent(NodeId, evtcode);
	if(PRE.hasTimeSlice())
	{
		setNodeUpperBoundaryTS(NodeId,PRE.gettimeSlice());//we set the time slice of the node as its father.
		setNodeLowerBoundaryTS(NodeId,PRE.gettimeSlice());//the upper and lower boundaries are the same here
	}
		


	//creating the new node
	int newId = MyNextId();//getNextId();
	Node * newNode  = new Node(newId);

	nbNodes++;
	//branching to the existing structure
	Node * currentNode  = getNode(NodeId);
	currentNode->addSon(newNode);
	
	setNodeCladeNum(newId, getNodeCladeNum(NodeId));//setting the new node cladenum property and updating the map


	return newId;
}



/*
Creates or annotate a leaf node at NodeId

Checks if he parent is a transfer ((Sout and different species from the parent)or Bout). If it is, the function will annotate NodeId as a reception node and create a leaf node.
Otherwise just annotate the node as a leaf.

Takes:
 - NodeId (int): id of the node to annotate with species and clade already set

 */
void ReconciledTree::CreateLeafNode(int NodeId)
{
	if(getNodeEvent(NodeId)!= C)//this is not already a leaf
	{

		int FatherEventId = getNodeEvent(getFatherId(NodeId));
		int FatherSpecies = getNodeSpecies(getFatherId(NodeId));
	
		if( ( (FatherEventId == Sout)&&(FatherSpecies != getNodeSpecies(NodeId) ) ) || (FatherEventId == Bout))
		{
			//annotate NodeId as a reception
			setNodeEvent(NodeId,R);
	
			int newId = MyNextId() ;// getNextId();
			Node * newNode  = new Node(newId);
			nbNodes++;
			//branching to the existing structure
			Node * currentNode  = getNode(NodeId);
			currentNode->addSon(newNode);
	
			//annotating
			setNodeSpecies(newId,getNodeSpecies(NodeId));
			setNodeCladeNum(newId,getNodeCladeNum(NodeId));
	
			NodeId = newId;
		}
	
		//setting NodeId as a leaf
		setNodeEvent(NodeId,C);
		setNodeUpperBoundaryTS(NodeId,0);//the time slice of a leaf is 0
		setNodeLowerBoundaryTS(NodeId,0);

	}
}


/*
Adds the reconciliation and return idUl and idUr or -1,-1 if idU is a leaf

Takes:
	- CR (CladeReconciliation) : a clade reconciliation of id idU
	- VERBOSE (bool) (default = false)
Returns:
	pair<int,int> : idUl and idUr or -1,-1 if idU is a leaf
*/
pair<int,int> ReconciledTree::addCladeReconciliation(CladeReconciliation CR, bool VERBOSE)
{

    //clock_t begin;
    //clock_t end;
    //double elapsed_secs;

    //if(VERBOSE)
    //{
    //    begin = clock();
    //}



	vector<int> NodeIds = getNodeIdsFromCladeId(CR.idU);
	int NodeId;

	if(VERBOSE)
		cout << "Adding clade " << CR.idU << " -> " << NodeIds.size() << " Nodes" << endl;


	if(NodeIds.size() == 0) //node absent from the tree: create it. THIS SHOULD ONLY HAPPEN ONCE: for the root
	{
		if(!CR.isRoot())
			throw Exception("ReconciledTree::addCladeReconciliation : CladeReconciliation with clade absent from tree while not root!");
		//node creation
		NodeId = 0;//getNextId(); //should set the id to 0. as this will be the root, this should be the first node added. Furthermore, getNextId() segfaults if there is no node in the tree...

		Node * newnode  = new Node(NodeId);

		setRootNode(newnode);//setting as the new root

		setNodeCladeNum(NodeId, CR.idU);//setting the new node cladenum property and updating the map

		nbNodes=1;
	}
	else //At least one node has the current clade id. We are interested in the last one as this is the one that was created last.
	{
		NodeId = NodeIds.back();
		NodeId = accountForPreviousSplit(CR.previousEvent, NodeId, VERBOSE);
	}

	if(VERBOSE)
		cout << "done pre-treatment." << endl;

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for CRadding PRE:" << elapsed_secs << endl;
    //    begin = clock();
    //}


    clock_t Sbegin;
    
    double elapsed_secsB1=0;
	double elapsed_secsB2=0;
	double elapsed_secsB3=0;


	//setting the node properties (and supplementary loss nodes if needed): species, event

	Node * currentNode = getNode(NodeId);

	//using the first event (always at least one event)
	int current_event = 0;
	NodeId = addReconciliationEvent(CR.Events[current_event],NodeId,currentNode,VERBOSE);

	if(VERBOSE)
		cout << "  added event " << current_event << endl;

	current_event++;

	//for each event after the first
	while(current_event < CR.Events.size())
	{

		Sbegin = clock();
		//creating the new node
		//int newId = getNextId();
		//if(newId != MyNextId())
		//	cout << "idstuff "<< newId <<  "<->" << MyNextId()<<endl;
		int newId = MyNextId();


		elapsed_secsB1 += double(clock() - Sbegin) / CLOCKS_PER_SEC;
		Sbegin = clock();

		Node * newNode  = new Node(newId);
		nbNodes++;

		//branching to the existing structure
		Node * currentNode  = getNode(NodeId);
		currentNode->addSon(newNode);
		
		elapsed_secsB2 += double(clock() - Sbegin) / CLOCKS_PER_SEC;
		Sbegin = clock();



		setNodeCladeNum(newNode, CR.idU);//setting the new node cladenum property and updating the map		

		NodeId = newId; // setting as current node

		//assigning the event
		NodeId = addReconciliationEvent(CR.Events[current_event],NodeId, newNode,VERBOSE);

		elapsed_secsB3 += double(clock() - Sbegin) / CLOCKS_PER_SEC;
		


		if(VERBOSE)
			cout << "  added event " << current_event << endl;

		current_event++;

	}



//    if(VERBOSE)
//    {
//        end = clock();
//        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//        cout << "time elapsed for CRadding EVadd:" << elapsed_secs << endl;
//
//        if(elapsed_secs > 0.01)
//        {
//			cout << "time elapsed->" << CR.Events.size() << "events" << elapsed_secs;
//			for(unsigned i = 0; i < CR.Events.size();i++)
//				cout << " - " << CR.Events[i].getevtName() ;
//
//			cout <<  " - " << elapsed_secsB1 << " " << elapsed_secsB2 << " " << elapsed_secsB3 ;
//			cout << endl;
//			
//        }
//
//        begin = clock();
//    }


	pair<int,int> idsToReturn;

	//adding children, if there are some
	if(!CR.isLeaf())
	{

		Node * currentNode  = getNode(NodeId);
		int newId;
		Node * newNode;
		//adding left child
		newId = MyNextId();//getNextId();

		newNode = new Node(newId);

		nbNodes++;
		//branching to the existing structure
		
		currentNode->addSon(newNode);
		
		setNodeCladeNum(newId, CR.idUl);//setting the new node cladenum property and updating the map				
		setNodeSpecies(newId, CR.Events.back().getidXl());//setting the node species to the one specified in the last event


		//adding right child
		newId = MyNextId();//getNextId();

		newNode = new Node(newId);

		nbNodes++;
		//branching to the existing structure
		
		currentNode->addSon(newNode);
		
		setNodeCladeNum(newId, CR.idUr);//setting the new node cladenum property and updating the map				
		setNodeSpecies(newId, CR.Events.back().getidXr());//setting the node species to the one specified in the last event

		//setting the return value
		idsToReturn.first = CR.idUl;
		idsToReturn.second = CR.idUr;

	}
	else
	{

		if(getNodeEvent(NodeId) != C)//the current event is not already a leaf event
		{

			//we have a leaf node to add and annotate
			int newId = MyNextId();//getNextId();
	
			Node * newNode  = new Node(newId);
	
			nbNodes++;
			//branching to the existing structure
			Node * currentNode  = getNode(NodeId);
			currentNode->addSon(newNode);
			
			setNodeCladeNum(newId, CR.idU);//setting the new node cladenum property and updating the map		
	
			int evtcode = C;//Current (leaf)
			int speciesid;
			//choosing the correct leaf species
			//NB: this case happen when the last event is either TL, SL or TLFD
	
			if(getNodeEvent(NodeId) == R) // The last event is a reception -> the species stays the same
				speciesid = getNodeSpecies(NodeId);
			else // SL -> species is, by convention, the right child in the last reconciliation event
				speciesid = CR.Events.back().getidXr();

			//setting properties
			setNodeSpecies(newId, speciesid);
			setNodeEvent(newId, evtcode);

			setNodeUpperBoundaryTS(newId,0);//the time slice of a leaf is 0
			setNodeLowerBoundaryTS(newId,0);


			NodeId = newId;
		}

		//Trying to get a leaf Name
		if(CR.name != "")
			setNodeName(NodeId, CR.name);

		//setting the return value
		idsToReturn.first = -1;
		idsToReturn.second = -1;
	}


    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for CRadding CHadd:" << elapsed_secs << endl;
    //    begin = clock();
    //}

	return idsToReturn;
}

/*
Takes:
	- Nodeid (int): a node id
	- evtLine (string): recPhyloXML line describing a reconciliation event
*/
void ReconciledTree::setNodeDetailsFromXMLLine(int Nodeid, string evtLine , map< string , int > & speciesNameToNodeId  , bool VERBOSE)
{

	map <string, string> * properties = new map <string,string>;
	string * value = new string(""); 
	string Tname = InterpretLineXML(evtLine,properties,value);

	setNodeDetailsFromXMLLine( Nodeid,  Tname,  properties, value , speciesNameToNodeId  , VERBOSE);
}


/*
Takes:
	- Nodeid (int): a node id
	- Tname (string) : name of the tag (here, event name)
	- properties (map <string, string> *): container for informations like species, timeSlice, ...
	- value ( string * ): should be empty, technically not used here

*/
void ReconciledTree::setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value,map< string , int > & speciesNameToNodeId  , bool VERBOSE)
{
	char * pEnd; // for str -> int conversion 
	
	
	for (map<string,string>::iterator it=properties->begin(); it!=properties->end(); ++it)
	{
		if(it->first.compare("timeSlice") == 0)
		{
			int TS =  (int) strtol(it->second.c_str(), &pEnd, 10) ;
			setNodeUpperBoundaryTS(Nodeid,TS);
			setNodeLowerBoundaryTS(Nodeid,TS);
		}
		else if (it->first.compare("speciesLocation") == 0)
		{
			string species = it->second.c_str() ;
			auto speciesIt  = speciesNameToNodeId.find(species) ;
			if( speciesIt == speciesNameToNodeId.end() )
			{
				// species not found ...
				throw Exception("ReconciledTree::setNodeDetailsFromXMLLine : unknown species : " + species );
			}
			setNodeSpecies(Nodeid , speciesIt->second);

		}
		else if (it->first.compare("destinationSpecies") == 0)
			setNodeSpecies(Nodeid, (int) strtol(it->second.c_str(), &pEnd, 10) );
		else if (it->first.compare("geneName") == 0)
			setNodeName(Nodeid,it->second);
	}


/*
C 		-->	Current (leaf)
S 		-->	Speciation (in the species tree)
L 		-->	Loss (in the species tree)
D 		-->	Duplication (in the species tree)
Sout 	-->	Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
R 		-->	Transfer reception
N 		-->	no event (to account for time slices)
Bout 	-->	Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
*/

	int evtCode = -1;
	//3. finding which event is implied
	if(Tname.compare("leaf") == 0)
	{
		setNodeUpperBoundaryTS(Nodeid,0);//leaves at time 0
		setNodeLowerBoundaryTS(Nodeid,0);

		evtCode = C;
	}
	else if(Tname.compare("speciation") == 0)
	{
		evtCode = S;
	}
	else if( (Tname.compare("branchingOut") == 0) || (Tname.compare("speciationOut") == 0))
	{
		evtCode = Sout;
	}
	else if(Tname.compare("bifurcationOut") == 0)
	{
		evtCode = Bout;
		setNodeSpecies(Nodeid, -1 );
	}
	else if(Tname.compare("duplication") == 0)
	{
		evtCode = D;
	}
	else if(Tname.compare("transferLoss") == 0) // here for retrocompatibility to the original recPhyloXML format
	{
		evtCode = R;
	}
	else if(Tname.compare("transferBack") == 0)
	{
		evtCode = R;
	}
	else if( (Tname.compare("branchingOutLoss") == 0) || (Tname.compare("speciationOutLoss") == 0))
	{
		evtCode = Sout;
		//create some Loss son in the same species, with the same time slice

		int newId = MyNextId(); //getNextId();
		Node * newLossNode  = new Node(newId);
		nbNodes++;
		Node * currentNode  = getNode(Nodeid);
		currentNode->addSon(newLossNode);

		//setting new loss node properties
		setNodeSpecies(newId,getNodeSpecies(Nodeid));
		setNodeEvent(newId, L); // setting as a LOSS
		setNodeCladeNum(newId, -1); // the clade num of the loss node is -1 -> the empty clade

	}
	else if(Tname.compare("speciationLoss") == 0)
	{
		evtCode = S;
		//create some Loss son, but in a unknown species... with time slice = current -1
		/// The loss node will be added in a later post-treatment phase as it requires a species tree to determine in which species the loss is. 
		// This won't get forgotten as these are speciations nodes with only one son

	}
	else if(Tname.compare("loss") == 0)
	{
		evtCode = L;

	}

	if(evtCode == -1)
		throw Exception("ReconciledTree::setNodeDetailsFromXMLLine : unknown event : " + Tname );
	else
		setNodeEvent(Nodeid,evtCode);

	if(VERBOSE)
	{
		cout << "set details of node " << Nodeid << ";" << evtCode << ";" << getNodeSpecies(Nodeid) << ";" ;
		if(hasNodeProperty(Nodeid,uts))
			cout <<getNodeTimeSlice(Nodeid) << endl;
		else
			cout << "NOTS" << endl;
	}
		

}

/*
read and add a clade from a recPhyloXML stream

!!! This function is recursive

Takes:
	- fileIN (ifstream): streaam to a recPhyloXML file
	- Nodeid (int) : node to put the clade IN
	- VERBOSE (bool) (default = false)
*/
void ReconciledTree::addXMLClade(ifstream& fileIN,int Nodeid, map< string , int > & speciesNameToNodeId ,  bool VERBOSE) // will need a species tree
{

	int Cnum = CladeIdToNodeIds.size(); //the clade num will be some equivalent to a BFO index

	setNodeCladeNum(Nodeid , Cnum);
	//cout << "before "<< Cnum << "- after " << CladeIdToNodeIds.size() << endl;
	//cout << nbNodes << endl; 
	if(VERBOSE)
		cout << "parsing clade " << Cnum << " id " << Nodeid<< endl;

	vector <string> evtLines;

	//2. parsing of the clade
	string line;
	getline( fileIN, line );

	map <string, string> * properties = new map <string,string>;
	string * value = new string(""); 
	string Tname = InterpretLineXML(line,properties,value);

	//cout << "value:" << *value << endl;

	while(Tname.compare("/clade") != 0) //--> line signing the end of the clade
	{
		if(Tname.compare("name") == 0) // found name
		{
			if(VERBOSE)
				cout << "found name : " << *value <<  endl;
			setNodeName(Nodeid, *value);

		}
		else if(Tname.compare("eventsRec") == 0) // reading reconciliation events
		{
			if(VERBOSE)
				cout <<  "parsing events" << endl;
			
			*value = "";
			delete properties;

			getline( fileIN, line );
			properties = new map <string,string>;
			Tname = InterpretLineXML(line,properties,value);

			while(Tname.compare("/eventsRec") != 0) // here, add some nodes ...
			{
				evtLines.push_back(line); // maybe give in the already interpreted line rather than re-interpret it inside constructor...
				
				//looping
				*value = "";
				delete properties;
				getline( fileIN, line );
				properties = new map <string,string>;
				Tname = InterpretLineXML(line,properties,value);

				if(Tname.compare("") == 0)
					throw Exception("ReconciledTree::addXMLClade : error while parsing...");
			}
		}
		else if(Tname.compare("clade") == 0) // found a child --> recursion
		{
			
			//create child clade

			//creating the new node
			int newId = MyNextId(); //getNextId();

			Node * newNode  = new Node(newId);
			nbNodes++;
			//branching to the existing structure
			Node * currentNode  = getNode(Nodeid);
			currentNode->addSon(newNode);

			//recursing
			addXMLClade(fileIN, newId, speciesNameToNodeId , VERBOSE );


		}

		//reading new line
		*value = "";
		delete properties;
		getline( fileIN, line );
		properties = new map <string,string>;
		Tname = InterpretLineXML(line,properties,value);
		
	}



	//3. Event parsing
	if(evtLines.size() == 0)
	{
		throw Exception("ReconciledTree::addXMLClade : no reconciliation events found");
	}
	else if(evtLines.size()>1)
	{
		for(unsigned i = 0; i < evtLines.size() - 1; i++ ) // all events but the last
		{
			//creating parent node 
			//cout << "creating evt" << evtLines[i] << endl;

			int newId = MyNextId(); //getNextId();

			Node * newNode  = new Node(newId);

			nbNodes++;

			//branching to the existing structure
			bool ROOT = isRoot(Nodeid);

			Node * currentNode  = getNode(Nodeid);

			if(!ROOT)
			{
				Node * FNode = currentNode->getFather();
				FNode->removeSon(currentNode);
				FNode->addSon(newNode);

				newNode->addSon(currentNode);
			}
			else
			{
				currentNode->addSon(newNode);
				rootAt(newNode); // rerooting
			}			


				

			//new putting correct attributes to in the new node
			setNodeCladeNum(newId , Cnum);

			setNodeDetailsFromXMLLine(newId,evtLines[i] , speciesNameToNodeId  , VERBOSE);

		}
	}

	setNodeDetailsFromXMLLine(Nodeid,evtLines.back() , speciesNameToNodeId , VERBOSE); // setting details on the last node == the bifurcation or leaf


}

/*
Constructor that uses a stream from a recPhyloXML file
Takes
 - fileIN (ifstream): stream rom a recPhyloXML file
 - Stree (MySpeciesTree *): a pointer to the species tree
 - map< string , int > & speciesNameToNodeId
 - VERBOSE (bool) [fdefault : false]
*/
void ReconciledTree::readPhyloXMLFile(ifstream& fileIN, MySpeciesTree * Stree, map< string , int > & speciesNameToNodeId , bool VERBOSE )
{
	// presuming that the stream will yield a <clade> line forming the root of the tree
	//first find the root
	if(VERBOSE)
		cout << "Creating Reconciled_Tree." << endl;

	int rootid = 0;

	Node * newnode  = new Node(rootid);

	setRootNode(newnode);//setting as the new root

	nbNodes++;

	addXMLClade( fileIN,rootid, speciesNameToNodeId, VERBOSE ); 
	
	///post-treatment requiring species tree
	vector <int> NodesIds  = getNodesId();

	for(unsigned i = 0 ; i < NodesIds.size(); i++)
	{
		int evt = getNodeEvent(NodesIds[i]);
		vector <int> SonsIds = getSonsId(NodesIds[i]);

		if((evt == S) && (SonsIds.size() == 1)) //speciation with only one child --> speciation loss
		{
			int spF = getNodeSpecies(NodesIds[i]);
			int spC = getNodeSpecies(SonsIds[0]);

			///ugly as ****
			vector <int>  Snodes = Stree->getNodesId();
			for(unsigned j = 0 ; j < Snodes.size(); j++)
				if(Stree->getRPO(Snodes[j]) == spF)
				{
					spF = Snodes[j];
					break;
				}

			vector <int> spSons = Stree->getSonsId(spF);
			//there should be 2 sons in spSons
			if(spSons.size() != 2)
				throw Exception("ReconciledTree::readPhyloXMLFile : species in the species tree without son as a speciation.");	

			int LostIndex = 0;
			if(spC == Stree->getRPO(spSons[LostIndex]))
				LostIndex++;

			if(VERBOSE)
				cout << "Adding Loss node in species " << spSons[LostIndex] << "->" << Stree->getRPO(spSons[LostIndex]) << endl;


			//creating new loss node
			int newId = MyNextId(); //getNextId();
			Node * newLossNode  = new Node(newId);
			nbNodes++;
			Node * currentNode  = getNode(NodesIds[i]);
			currentNode->addSon(newLossNode);
	
			//setting new loss node properties
			setNodeSpecies(newId,Stree->getRPO(spSons[LostIndex]));
			setNodeEvent(newId, L); // setting as a LOSS
			setNodeCladeNum(newId, -1); // the clade num of the loss node is -1 -> the empty clade

		}
		else if((evt == S) && (SonsIds.size() == 0))
			throw Exception("ReconciledTree::readPhyloXMLFile : speciation without any son.");
	}


	/// determining timeSlice Status
	if(hasNodeProperty(getRootId(),uts))
		TimeSliceStatus = 1;
	else
		TimeSliceStatus = 0;

}

/* 
* *RECURSIVE*
* @arg int nodeid : id of the node whose subtree string is built
* @arg bool hideLosses (default: false): if true, losses and the branches leading to them wil be removed from the newick string
*
* @return string : newick representation of the subtree rooted at nodeid 
*/
string ReconciledTree::NodeString(int nodeid, bool hideLosses)
{
	string Str = "";

	vector<int> sons = getSonsId(nodeid);
	if(sons.size()> 0)
	{

		vector <string> sonStr;

		for(unsigned i = 0; i < sons.size(); i++)
		{
			string tmp =  NodeString(sons[i], hideLosses); // recursion
			if(tmp.size() >0) // don't add it the the representation is an empty str
				sonStr.push_back(tmp);
			
		}

		
		if(sonStr.size() == 0) // no son with a non-empty representation
			return Str; // return empty string

		//building Str
		Str += "(";
		for(unsigned i = 0; i < sonStr.size(); i++)
		{	
			if(i != 0)
				Str += ",";
			Str += sonStr[i];
		}
		Str += ")";
	}

	if(hasNodeName(nodeid))
		Str += getNodeName(nodeid);
	else
	{
		Str += static_cast<ostringstream*>( &(ostringstream() << nodeid) )->str() ;
	}
		

	Str += "|";

	if(!hasNodeProperty(nodeid,ev))
		Str += "NOEVENT";
	else
	{

		int evt = getNodeEvent(nodeid);
		switch(evt)
		{ 
			case C :
				Str += "Extant"; //Current (leaf)
				break;
			case S :
				Str += "Spe"; //Speciation (in the species tree)
				break;
			case L :
				if(hideLosses)// soecial option so that losses and their branch oare not in the final newick tree
				{
					return ""; // return empty string
				}
				Str += "Loss"; //Loss (in the species tree)
				break;
			case D :
				Str += "Dup"; //Duplication (in the species tree)
				break;
			case Sout :
				Str += "SpeOut"; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
				break;
			case R :
				Str += "Reception"; //Transfer reception
				break;
			case N :
				Str += "Null"; //no event (to account for time slices)
				break;
			case Bout :
				Str += "BifOut"; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
				break;
			default:
				Str += "NOEVENT";
		}
	}

	Str += "|";

	if(!hasNodeProperty(nodeid,spe))
		Str += "no species";
	else
	{
		int Number = getNodeSpecies(nodeid);
		Str += static_cast<ostringstream*>( &(ostringstream() << Number) )->str() ;
	}

	Str += "|";

	if(hasNodeProperty(nodeid,uts))
	{
		int UTS = getNodeUpperBoundaryTS(nodeid);
		int LTS = getNodeLowerBoundaryTS(nodeid);
		if(UTS != LTS)
		{
			Str += static_cast<ostringstream*>( &(ostringstream() << UTS) )->str() ;
			Str += "-";
			Str += static_cast<ostringstream*>( &(ostringstream() << LTS) )->str() ;
			
		}
		else
			Str += static_cast<ostringstream*>( &(ostringstream() << UTS) )->str() ;

	}
	else
		Str += "NOTS";

	return Str;
}

/*
Constructor from a root and a time slice status.

Intended to be used by CloneSubtree
*/
ReconciledTree::ReconciledTree(Node * root, int TSstatus): TreeTemplate<Node>(root)
{
	TimeSliceStatus = TSstatus;
}

/**
Constructor that takes a pointer to a vector of CladeReconciliation
   
Takes:
 - CladeRecs (vector<CladeReconciliation> * ) : pointer to a vector of CladeReconciliation
 - bool VERBOSE (default=false): wether the function talks or not 

 */
ReconciledTree::ReconciledTree(vector<CladeReconciliation> * CladeRecs, bool VERBOSE): TreeTemplate<Node>()
{

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    //VERBOSE=true;

    if(VERBOSE)
    {
        begin = clock();
    }


	//we do not make the assumption that the CladeReconciliation are in correct order, where correct means that the root is first and a parent always occur before its child(ren)

	//first find the root
	if(VERBOSE)
		cout << "Creating Reconciled_Tree." << endl;

	
	vector <int> CladeRecToaDD;
	int current_index = 0;

	int nbCladeRecs = CladeRecs->size();

	map<int,int> cladeNumToIndex; //mapping the cladenumber to their index in the vector *CladeRecs


	for(int i=0;i<nbCladeRecs;i++)
	{
		if(CladeRecs->at(i).isRoot())
		{
			CladeRecToaDD.push_back(CladeRecs->at(i).idU);//adding the root as the first element to do
		}
		cladeNumToIndex[CladeRecs->at(i).idU] = i;//mapping
	}

	if(CladeRecToaDD.size() == 0)//no root found -> error
		throw Exception("ReconciledTree::ReconciledTree : root reconciliation not found!");

	if(VERBOSE)
		cout << " -> map of the Clade Reconciliation done." << endl;

    if(VERBOSE)
    {
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "time elapsed for TREE CRmapping:" << elapsed_secs << endl;
        begin = clock();
    }

    int counter = 0;
	//we can now make a depth-first traversal of the CladeReconciliations
	while(! CladeRecToaDD.empty())
	{
		clock_t Sbegin = clock();
    	
    	int cladeIdToAdd = CladeRecToaDD.back();
		//CladeRecToaDD.erase(CladeRecToaDD.rbegin());//erasing the clade we just added
		CladeRecToaDD.pop_back();

		if(cladeNumToIndex.count( cladeIdToAdd ) != 0 )//the clade is present in the CladeReconciliation (it can be absent if it is a branch without any partiular event)
		{
			if(VERBOSE)
			{
				cout << "Adding clade id " << cladeIdToAdd << endl;
			}

			int current_index = cladeNumToIndex[cladeIdToAdd];

			if(VERBOSE)
			{
				cout << " -> Adding Clade reconciliation :" << endl;
				CladeRecs->at(current_index).printMe();
			}

	
			pair<int,int> new_index =  addCladeReconciliation(CladeRecs->at(current_index), VERBOSE);//add the reconciliation and return idUl and idUr or -1,-1 if idU is a leaf
	
			if(VERBOSE)
				cout << "new indexes " << new_index.first << "," << new_index.second << endl;

			if(new_index.first != -1)//adding if they exists
				CladeRecToaDD.push_back(new_index.first);
			if(new_index.second != -1)
				CladeRecToaDD.push_back(new_index.second);
	
		}
		else
		{

			if(VERBOSE)
			{
				cout << "Annotating as a leaf clade id " << cladeIdToAdd << endl;
			}


			//means that this is a leaf which have (at least) a node that has not been annotated as a leaf (however clade and species are set)
			int NodeId = getNodeIdsFromCladeId( cladeIdToAdd ).back();//retrieving the last node with that clade id
			CreateLeafNode(NodeId);//cleaner, accounts for eventual reception node tha we missed...

		}

		if(VERBOSE)
			cout << "size of CladeRecToaDD : " << CladeRecToaDD.size() << endl;


	    if(VERBOSE)
    	{
			counter++;
			clock_t Send = clock();
			elapsed_secs = double(Send - Sbegin) / CLOCKS_PER_SEC;
    	    cout << "time elapsed for CRadding " << counter << ":" << elapsed_secs << endl;	
		}
	}


    if(VERBOSE)
    {
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "time elapsed for TREE CRadding:" << elapsed_secs << " for " << counter << " clades->" << elapsed_secs / counter << endl;
        begin = clock();
    }


	//post treatment: decide ifthere are TS or no
	
	if(hasNodeProperty(getRootId(),uts))
		TimeSliceStatus = 1;
	else
		TimeSliceStatus = 0;


	if(VERBOSE)
		cout << "done" << endl;

}

/*
Constructor that uses a stream from a recPhyloXML file
Takes
 - fileIN (ifstream): stream rom a recPhyloXML file
 - Stree (MySpeciesTree *): a pointer to the species tree
 - map< string , int > & speciesNameToNodeId : associates the name of species to their node Id
 - VERBOSE (bool) [fdefault : false]
*/
ReconciledTree::ReconciledTree(ifstream& fileIN, MySpeciesTree * Stree, map< string , int > & speciesNameToNodeId, bool VERBOSE): TreeTemplate<Node>()
{

	readPhyloXMLFile( fileIN, Stree, speciesNameToNodeId , VERBOSE );
}

/** 
Constructor that uses a recPhyloXML file name
Takes
 - fileIN (ifstream): stream rom a recPhyloXML file
 - Stree (MySpeciesTree *): a pointer to the species tree
 - map< string , int > & speciesNameToNodeId : associates the name of species to their node Id
 - VERBOSE (bool) [fdefault : false]
*/
ReconciledTree::ReconciledTree(string phyloxmlFileName, MySpeciesTree * Stree , map< string , int > & speciesNameToNodeId , bool VERBOSE )
{
	ifstream fileStream(phyloxmlFileName.c_str());
	if( !fileStream.is_open() ) 
	{
		throw Exception("ReconciledTree::ReconciledTree : could not open reconciled tree file : " + phyloxmlFileName);
		exit(1);
	}


	while( !fileStream.eof() ) 
	{

		string line;
		getline( fileStream, line );
		
		map <string, string> * properties = new map <string,string>;
		string * value; 
		string Tname = InterpretLineXML(line,properties,value);

		if(Tname.compare("clade") == 0)
		{
			readPhyloXMLFile(fileStream,Stree, speciesNameToNodeId ,  VERBOSE);
			break;
		}

	}


}


/*
Returns 0 if there is no time slices (NTS), 1 if there are time slices (TS) or 2 if there are bounded time slices (BTS)
*/
int ReconciledTree::getTimeSliceStatus()
{
	return TimeSliceStatus;
}

/*
Access to NodeId with CladeId
Takes:
 - CladeId (int): id of a clade

Returns:
	vector<int>: NodeIds of the nodes with that clade OR empty vector if the clade is absent

*/
vector<int> ReconciledTree::getNodeIdsFromCladeId(int CladeId)
{
	if(CladeIdToNodeIds.count(CladeId) == 0)
		return vector<int>();
	return CladeIdToNodeIds[CladeId];
}



/*
Adds NodeId to CladeId
Takes:
 - CladeId (int): id of a clade
 - NodeId (int):  if of the node with that clade

*/
void ReconciledTree::addNodeIdToCladeIdMap(int CladeId, int NodeId)
{
	//cout << "ReconciledTree::addNodeIdToCladeIdMap " << CladeId << "," << NodeId << endl;
	if(CladeIdToNodeIds.count(CladeId) == 0)//key is absent from the map -> add it
		CladeIdToNodeIds[CladeId] = vector<int>();
	CladeIdToNodeIds[CladeId].push_back(NodeId);
}


//getter of node properties

/*
retrieve the specified node species. 

Takes:
 - NodeId (int): id of a node in the tree

Returns:
	(int): speciesid of the node
*/
int ReconciledTree::getNodeSpecies(int nodeid)
{

	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(spe) )->getValue();
	}
	//else we update the map

	NodeIdToNodeP[nodeid] = getNode(nodeid);*/


	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, spe))->getValue();
}

int ReconciledTree::getNodePostOrder(int nodeid)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(porder) )->getValue();
	}
	//else we update the map
	NodeIdToNodeP[nodeid] = getNode(nodeid);*/

	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, porder))->getValue();
}

int ReconciledTree::getNodeTimeSlice(int nodeid)
{

	if(hasNodeProperty(nodeid,ts))
	{
		/*
		map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
		if(it != NodeIdToNodeP.end())
		{
			return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(ts) )->getValue();
		}
		//else we update the map
		NodeIdToNodeP[nodeid] = getNode(nodeid);*/

		return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, ts))->getValue();
	}
	//in case the property is not set, return the upper time slice
	return getNodeUpperBoundaryTS(nodeid);
}

int ReconciledTree::getNodeUpperBoundaryTS(int nodeid)
{
	/*map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(uts) )->getValue();
	}
	//else we update the map
	NodeIdToNodeP[nodeid] = getNode(nodeid);*/

	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, uts))->getValue();
}

int ReconciledTree::getNodeLowerBoundaryTS(int nodeid)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(lts) )->getValue();
	}
	//else we update the map
	NodeIdToNodeP[nodeid] = getNode(nodeid);*/


	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, lts))->getValue();
}


int ReconciledTree::getNodeEvent(int nodeid) // retrieves the eventid of the node. eventid is an int corresponding to a reconciliation event
{

	/*map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(ev) )->getValue();
	}
	//else we update the map
	NodeIdToNodeP[nodeid] = getNode(nodeid);*/

	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, ev))->getValue();
}

int ReconciledTree::getNodeCladeNum(int nodeid) // retrieves the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{

	/*map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		return dynamic_cast<BppInteger *> ( it->second->getNodeProperty(clnum) )->getValue();
	}
	//else we update the map
	NodeIdToNodeP[nodeid] = getNode(nodeid);*/

	return dynamic_cast<BppInteger *> (getNodeProperty(nodeid, clnum))->getValue();
}

//setter of node properties
void ReconciledTree::setNodeSpecies(int nodeid, int speciesid) // set the specified node species.
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(spe, BppInteger(speciesid) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(spe, BppInteger(speciesid));
	}*/
	Node * NP = getNode(nodeid);
	NP->setNodeProperty(spe, BppInteger(speciesid));
}

void ReconciledTree::setNodeSpecies(Node * NP, int speciesid) // set the specified node species.
{
	NP->setNodeProperty(spe, BppInteger(speciesid));	
}

void ReconciledTree::setNodePostOrder(int nodeid, int postorder)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(porder, BppInteger(postorder) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(porder, BppInteger(postorder));
	}*/

	Node * NP = getNode(nodeid);
	NP->setNodeProperty(porder, BppInteger(postorder));

}

void ReconciledTree::setNodePostOrder(Node * NP, int postorder)
{
	NP->setNodeProperty(porder, BppInteger(postorder));	
}

void ReconciledTree::setNodeTimeSlice(int nodeid, int timeslice)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(ts, BppInteger(timeslice) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(ts, BppInteger(timeslice));
	}*/
	Node * NP = getNode(nodeid);

	NP->setNodeProperty(ts, BppInteger(timeslice));
}

void ReconciledTree::setNodeTimeSlice(Node * NP, int timeslice)
{
	NP->setNodeProperty(ts, BppInteger(timeslice));	
}

void ReconciledTree::setNodeUpperBoundaryTS(int nodeid, int timeslice)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(uts, BppInteger(timeslice) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(uts, BppInteger(timeslice));
	}*/
	Node * NP = getNode(nodeid);
	NP->setNodeProperty(uts, BppInteger(timeslice));
}

void ReconciledTree::setNodeUpperBoundaryTS(Node * NP, int timeslice)
{
	NP->setNodeProperty(uts, BppInteger(timeslice));
}


void ReconciledTree::setNodeLowerBoundaryTS(int nodeid, int timeslice)
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(lts, BppInteger(timeslice) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(lts, BppInteger(timeslice));
	}*/
	Node * NP = getNode(nodeid);
	NP->setNodeProperty(lts, BppInteger(timeslice));
}

void ReconciledTree::setNodeLowerBoundaryTS(Node * NP, int timeslice)
{
	NP->setNodeProperty(lts, BppInteger(timeslice));	
}


void ReconciledTree::setNodeEvent(int nodeid, int eventid) // assigns an event to the node according to the eventid
{
	/*
	map<int,Node * >::iterator it = NodeIdToNodeP.find(nodeid);
	if(it != NodeIdToNodeP.end())
	{
		it->second->setNodeProperty(ev, BppInteger(eventid) );
	}
	else
	{
		//else we update the map
		Node * NP = getNode(nodeid);
		NodeIdToNodeP[nodeid] = NP;

		NP->setNodeProperty(ev, BppInteger(eventid));
	}*/
	Node * NP = getNode(nodeid);
	NP->setNodeProperty(ev, BppInteger(eventid));
}

void ReconciledTree::setNodeEvent(Node * NP, int eventid)
{
	NP->setNodeProperty(ev, BppInteger(eventid));	
}


void ReconciledTree::setNodeCladeNum(int nodeid, int cladenum) // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{
	if(hasNodeProperty(nodeid, clnum))//property already set
	{
		if(getNodeCladeNum(nodeid) == cladenum)//already the correct property -> don't change anything
			return;
		else
			resetNodeCladeNum(nodeid); // reset the property and update the map cladeid -> nodeid
	}

	addNodeIdToCladeIdMap(cladenum, nodeid);//updating the map cladeid -> nodeid
	setNodeProperty(nodeid, clnum, BppInteger(cladenum));
}

void ReconciledTree::setNodeCladeNum(Node * NP, int cladenum) // sets the clade id of the node. Clade id is a concept used for CCPs and will typically be defined in a reference CladesAndTripartitions instance
{
	int nodeid = NP->getId();

	if(NP->hasNodeProperty(clnum))
	{
		if(getNodeCladeNum(nodeid) == cladenum)
			return;
		else
			resetNodeCladeNum(nodeid);
	}


	addNodeIdToCladeIdMap(cladenum, nodeid);//updating the map cladeid -> nodeid
	NP->setNodeProperty(clnum, BppInteger(cladenum));

}



//unset the clade id of the node. Always use this function because it also updates the CladeIdToNodeIds map
void ReconciledTree::resetNodeCladeNum(int nodeid)
{
	int CladeId = getNodeCladeNum(nodeid);

	getNode(nodeid)->deleteNodeProperty(clnum);

	///updating the map
	if(CladeIdToNodeIds.count(CladeId) != 0)//key is not absent from the map
	{
		vector<int>::iterator it;
		for (it = CladeIdToNodeIds[CladeId].begin() ; it != CladeIdToNodeIds[CladeId].end(); ++it)
    		if(*it == nodeid)
    			break;

    	if(it != CladeIdToNodeIds[CladeId].end())
	    	CladeIdToNodeIds[CladeId].erase(it);//erasing
				

		if(CladeIdToNodeIds[CladeId].size() == 0)
			CladeIdToNodeIds.erase(CladeId);
	}
	//the else case should never happen if the cladenum was set properly
}


/*
Prints node information

Takes:
 - NodeId (int): id of a node in the tree

 */
void ReconciledTree::printNode(int nodeid)
{

	cout << "Node : " << nodeid ;

	if(hasNodeName(nodeid))
		cout << " " <<getNodeName(nodeid);

	if(hasFather(nodeid))
		cout << " ; father : " << getFatherId(nodeid);
	else
		cout << " ; root ";

	if(hasNodeProperty(nodeid, spe))
		cout << " ; Species " << getNodeSpecies(nodeid);
	else
		cout << " ; no Species";

	if(hasNodeProperty(nodeid, ev))
		cout << " ; Event " << getNodeEvent(nodeid);
	else
		cout << " ; no Event";

	if(hasNodeProperty(nodeid,clnum))
		cout << " ; CladeNum " << getNodeCladeNum(nodeid);
	else
		cout << " ; no CladeNum";

	if(hasNodeProperty(nodeid,uts))
	{
		int UTS = getNodeUpperBoundaryTS(nodeid);
		int LTS = getNodeLowerBoundaryTS(nodeid);
		if(UTS != LTS)
		{
			cout << " ; Upper Time Slice " << UTS;
			cout << " ; Lower Time Slice " << LTS;
		}
		else
			cout << " ; Time Slice " << UTS;

	}
	else
		cout << " ; no TS";

//	if(!isLeaf(nodeid))
//	{
//		cout << " ; sons :" ;
//		vector <int> SonIds = getSonsId(nodeid);
//		for(unsigned i = 0; i < SonIds.size(); i++)
//			cout << SonIds[i] << " ";
//	}


	cout << endl;
}

void ReconciledTree::printMe()//prints all nodes
{

	vector <int> NodesIds = getNodesId();
	for (unsigned j=0; j<NodesIds.size() ; j++)
		printNode(NodesIds[j]);

}

string ReconciledTree::NewickString(bool hideLosses )
{
	return NodeString(getRootId(), hideLosses);
}

/*
	Removes all nodes whose event is No Event: 6
*/
void ReconciledTree::removeNoEventNodes()
{
	vector<int> NodesIds  = getNodesId();
	int NodeId;

	for (unsigned j=0; j<NodesIds.size() ; j++)
	{
		NodeId= NodesIds[j];

		if(getNodeEvent(NodeId) == N)//NoEvent node
		{
			resetNodeCladeNum(NodeId);//forgetting that the node has this cladenum

			Node * NoEventNode = getNode(NodeId);

			Node * fatherNode = NoEventNode->getFather();

			Node * sonNode = NoEventNode->getSon(0);//we want the 1st son (there should be 1 and only 1 son to a NoEvent node)

			fatherNode->removeSon(NoEventNode);
			fatherNode->addSon(sonNode);// also tells the son that its father is the fatherNode

			delete NoEventNode;
		}
	}

}
/*
Precomputes time slices at speciation and leaves nodes 

Takes:
- idXToLTS (map <int,int>) : map of species id to minimum associated time slice
*/
void ReconciledTree::PreComputeBoundedTimeSlices(map <int,int> idXToLTS)
{
	
	vector <int> IDS = getNodesId();

	for(unsigned i = 0; i < IDS.size(); i++)
	{
		int NodeId = IDS[i];
		int evt = getNodeEvent(NodeId);

		if(evt == C)// trivial case of leaves
		{
			setNodeUpperBoundaryTS(NodeId, 0);
			setNodeLowerBoundaryTS(NodeId, 0);
		}
		else if(evt == S)
		{
			int TS = idXToLTS[getNodeSpecies(NodeId)];
			setNodeUpperBoundaryTS(NodeId,TS);
			setNodeLowerBoundaryTS(NodeId,TS);
		}
	}
}

/*
!!! presumes that LTS and BTS proerties have been set for all nodes !!!
Used to fix the holes than can be left after all null nodes are deleted in BTS mode
Acts in a recursive fashion from the given node.
Check if the UTS (upper) of NodeId has a difference >1 with the LTS (lower) of its parent. If this is the case, an intermediary Null node is added

Takes:
 - int NodeId : the id of a node in the tree 
*/
void ReconciledTree::addIntermediaryNullBTS(int NodeId)
{

	vector <int> SonsIds = getSonsId(NodeId);
	for(unsigned i = 0 ; i < SonsIds.size() ; i++)
	{ // recursion
		addIntermediaryNullBTS(SonsIds[i]);
	}


	if( isRoot(NodeId) ) // the root has no parent -> nothing more to do
		return ; 


	int parentId = getFatherId(NodeId);

	//these properties shall be set and there will be errors if they aren't
	int selfUTS = getNodeUpperBoundaryTS(NodeId);
	int parentLTS = getNodeLowerBoundaryTS(parentId);

	if( (parentLTS - selfUTS) > 1 )
	{ 
		// means that there is a "time-hole" between this node and it parent.
		// we add a Null event node to guarantee that all time slices can be represented 

		int nullUTS = parentLTS - 1;
		int nullLTS = selfUTS + 1;

		int evt = N;
		int cladenum = getNodeCladeNum(NodeId);
		int species = getNodeSpecies(NodeId);


		Node * FatherNode = getNode(parentId);

		Node * currentNode = getNode(NodeId);

		int newId = MyNextId(); // getNextId();
		Node * newNode  = new Node(newId);
		nbNodes++;

		//branching to the father
		FatherNode->addSon(newNode);


		//setting the new node properties
		setNodeSpecies(newNode,species);
		setNodeEvent(newNode,evt);
		setNodeCladeNum(newNode,cladenum);

		setNodeTimeSlice(newNode, nullUTS);
		setNodeUpperBoundaryTS(newNode, nullUTS);
		setNodeLowerBoundaryTS(newNode, nullLTS);

		// rebranching to existing structure
		FatherNode->removeSon(currentNode);
		newNode->addSon(currentNode);//branching current node to the Null node before it


	}

}

/*
Recursively sets the upper and lower boundaries of a subtree.
Takes:
	- NodeId (int) : id of the root of the subtree
	- idXToUTS (map <int,int>) : map of species id to maximum associated time slice
	- idXToLTS (map <int,int>) : map of species id to minimum associated time slice
*/
void ReconciledTree::setSubtreeBoundaryTS(int NodeId, 	map <int,int> idXToUTS, map <int,int> idXToLTS)
{

	//first step : set the time slice as the current upper boundary time slice (which may be changed) so that we can keep this information

	//!!! in case there is no time slice set (example of loss nodes from XML input)
	if(! hasNodeProperty(NodeId , uts) ) // special loss cases
	{
		int currentTS = getNodeTimeSlice(getFatherId(NodeId)) - 1;
		if(currentTS < 0)
			currentTS = 0;
		//correcting the tree while we're at it
		setNodeUpperBoundaryTS(NodeId, currentTS);//set them at the current time slice
		setNodeLowerBoundaryTS(NodeId, currentTS);
	}
	
	setNodeTimeSlice(NodeId, getNodeUpperBoundaryTS(NodeId));

	//rest of the function
	int evt = getNodeEvent(NodeId);
	if(evt == C)
		return; //this is a leaf, LTS = UTS and is already set 

	bool isSpeciation = false;
	if(evt == S)
		isSpeciation = true;//in this case, LTS = UTS and is already set

	if(!isSpeciation)
	{
		//setting UTS as the UTS of its father

		int UTS;
		if(!hasFather(NodeId)) //in case the root is not a speciation
			UTS = idXToUTS[getNodeSpecies(NodeId)] ; 
		else
			UTS = getNodeUpperBoundaryTS(getFatherId(NodeId));
		


		if(evt != Bout) // Bout event are in a lost species, so we only have to check their parent in the gene tree
		{
			//check the max TS on the species branch we are in to compare
			int altUTS = idXToUTS[getNodeSpecies(NodeId)];
			//cout << NodeId << " " << evt << " : " << UTS << "<>" <<  altUTS << endl;
			if(altUTS < UTS)
				UTS = altUTS;//replace if lower
		}
		setNodeUpperBoundaryTS(NodeId,UTS);
	}

	//we can now recurse to the children

	vector <int> SonsIds = getSonsId(NodeId);

	int LTS = idXToLTS[getNodeSpecies(NodeId)];//start with a minimum value of the LTS

	for(int i = 0; i < SonsIds.size(); i++)
	{
		setSubtreeBoundaryTS(SonsIds[i],idXToUTS, idXToLTS); //RECURSION TO THE CHILDREN

		//getting the (now set) LTS of that child
		int altLTS = getNodeLowerBoundaryTS(SonsIds[i]);

		if(altLTS > LTS)//stocking it if it is the minimum
			LTS = altLTS;
	}

	if(!isSpeciation)//setting LTS is the node isn't a speciation
		setNodeLowerBoundaryTS(NodeId,LTS);

	
	setNodeTimeSlice(NodeId, getNodeUpperBoundaryTS(NodeId));

}


/*
Takes:
 - SpeciesTree (MySpeciesTree *): pointer to ultrametric species tree

Computes and set the bounded time slices for all nodes in the tree
The lower time slice are inherited from the children,
The upper time slices are inherited from the parent.

This also destroys the nodes with No Event (6) events

*/
void ReconciledTree::computeBoundedTimeSlices(MySpeciesTree * SpeciesTree)
{
//	if( getTimeSliceStatus() == 0)//first check if the time slices exists
//		cout << "ReconciledTree::computeBoundedTimeSlices : No time slices to draw from!" << endl;

	
	//1. delete NoEvent nodes
	removeNoEventNodes();

	//2. build a map between species and maximum time slice (i.e. the time slice of the previous speciation -1)

	map <int,int> idXToUTS;
	map <int,int> idXToLTS;
	vector <MySpeciesNode*> SpNodes = SpeciesTree->getNodes();

	for(unsigned j=0; j < SpNodes.size(); j++)
	{ //we want only speciation nodes: the ones whose real post order is != that of their children
		if(SpNodes[j]->isLeaf())
			continue;

		int RPO = SpNodes[j]->getInfos().realPostOrder; // getting the real species name
		int SonRPO = SpNodes[j]->getSon(0)->getInfos().realPostOrder; // species name of the (first) son

		if(RPO != SonRPO)//speciation node
		{
			int nbson = SpNodes[j]->getNumberOfSons();
			for(int i=0; i< nbson; i++)//setting time slice for all son
				idXToUTS[SpNodes[j]->getSon(i)->getInfos().realPostOrder] =  SpNodes[j]->getInfos().timeSlice - 1; 

			idXToLTS[RPO] = SpNodes[j]->getInfos().timeSlice;
		}
	}

	//2.5 precomputes the TSs at leaves and speciations
	PreComputeBoundedTimeSlices( idXToLTS);

	//3. Recursively set the UTS and LTS of nodes
	setSubtreeBoundaryTS(getRootId(), idXToUTS, idXToLTS);

	//3.5 add null event nodes when there is one of several time slices missing in a branch
	addIntermediaryNullBTS(getRootId());

	//4. setting the TimeSliceStatus to notify of BTS
	TimeSliceStatus = 2;
}

/*
Removes null event nodes and set to not time sliced
*/
void ReconciledTree::setToNonTimeSliced()
{
	if( getTimeSliceStatus() == 0)//first check if the time slices exists
		return

	//1. delete NoEvent nodes
	removeNoEventNodes();

	//2. setting the TimeSliceStatus to notify of BTS
	TimeSliceStatus = 0;

}

void ReconciledTree::SubdivideTree()
{
	//expects that node already have an exact set time slice
	if( getTimeSliceStatus() == 0)//first check if the time slices exists
		throw Exception("ReconciledTree::SubdivideTree : No time slices to draw from!");

	bool resetULTS = false;

	if( getTimeSliceStatus() == 2)// setting upper and lower time slice boundaries as equal to the TS
		resetULTS = true;



	vector<int> NodesIds  = getNodesId();
	int NodeId;

	for (unsigned j=0; j<NodesIds.size() ; j++)
	{
		NodeId= NodesIds[j];


		Node * currentNode = getNode(NodeId);
	
		int currentTS;

//		cout << NodeId << " " << currentNode->hasNodeProperty(uts) << "->" ;

		if(! currentNode->hasNodeProperty(uts) ) // special loss cases
		{
			Node * fatherNode = currentNode->getFather();

			currentTS = getNodeTimeSlice(fatherNode->getId()) - 1;
			if(currentTS < 0)
				currentTS = 0;

			//correcting the tree while we're at it
			setNodeUpperBoundaryTS(NodeId, currentTS);//set them at the current time slice
			setNodeLowerBoundaryTS(NodeId, currentTS);

		}
		else
			currentTS = getNodeTimeSlice(currentNode->getId());

//		cout << currentNode->hasNodeProperty(uts) <<endl;
	

		if(resetULTS)//reset the upper and lower time slices
		{
			setNodeUpperBoundaryTS(NodeId, currentTS);//set them at the current time slice
			setNodeLowerBoundaryTS(NodeId, currentTS);
		}

		if(!currentNode->hasFather())//this is the root node -> do nothing else
			continue;



		Node * fatherNode = currentNode->getFather();

		int fatherTS = getNodeTimeSlice(fatherNode->getId());

		int nbMissingNodes = fatherTS - currentTS - 1;


		if(nbMissingNodes > 0 )
		{
			int sp = getNodeSpecies(NodeId);
			if(getNodeEvent(NodeId) == R)
				sp = -1;

			int evt = N;
			int cladenum = getNodeCladeNum(NodeId);

			Node * newFatherNode = fatherNode;



			for(int i = 0; i< nbMissingNodes;i++)//creating the missing nodes 
			{


				int newId = MyNextId();// getNextId();

				Node * newNode  = new Node(newId);
				nbNodes++;

				//branching to the father
				newFatherNode->addSon(newNode);

				//setting the new node properties
				setNodeSpecies(newNode,sp);
				setNodeEvent(newNode,evt);
				setNodeCladeNum(newNode,cladenum);

				setNodeTimeSlice(newNode, fatherTS - i - 1);
				setNodeUpperBoundaryTS(newNode, fatherTS - i - 1);
				setNodeLowerBoundaryTS(newNode, fatherTS - i - 1);


				//incrementing
				newFatherNode = newNode;
			}

			fatherNode->removeSon(currentNode);
			newFatherNode->addSon(currentNode);//branching current node to the Null node before it
		}

	}

	//set time slice status to notify of TS
	TimeSliceStatus = 1;	
}


int ReconciledTree::getNumberOfDuplication()
{
	return countAnyEvent(D);
}

int ReconciledTree::getNumberOfLoss()
{
	return countAnyEvent(L);
}

int ReconciledTree::getNumberOfTransfer()
{
	return countAnyEvent(R);
}


/*
Checks if the node have the same species

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node * ): a node with at least the species property set

Returns:
	(bool) : true if the two nodes have the same species; false otherwise

*/
bool ReconciledTree::haveSameSpecies(Node  * n1, Node * n2)
{
	//basic checks
	if((!n1->hasNodeProperty(spe)) || (!n2->hasNodeProperty(spe)))
		throw Exception("ReconciledTree::haveSameSpecies : given nodes don't have a species (S property)");

	int n1s = dynamic_cast<BppInteger *> (n1->getNodeProperty(spe))->getValue();
	int n2s = dynamic_cast<BppInteger *> (n2->getNodeProperty(spe))->getValue();

	if(n1s != n2s)
		return false;

	return true;
}

/*
Checks if the node have compatible time slices

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node * ): a node with at least the species property set

Returns:
	(bool) : true if the two nodes have compatible time slices; false otherwise

*/
bool ReconciledTree::areTSCompatible(Node  * n1, Node * n2)
{

	if(TimeSliceStatus == 0) // no TS -> same species is sufficient
		return true;

	if(TimeSliceStatus == 1)
	{
		int n1ts = dynamic_cast<BppInteger *> (n1->getNodeProperty(uts))->getValue();
		int n2ts = dynamic_cast<BppInteger *> (n2->getNodeProperty(uts))->getValue();
		if(n1ts == n2ts)
			return true;
		return false;
	}

	int n1uts = dynamic_cast<BppInteger *> (n1->getNodeProperty(uts))->getValue();
	int n2lts = dynamic_cast<BppInteger *> (n2->getNodeProperty(lts))->getValue();

	//BTS case: TimeSliceStatus == 2
	if(n1uts < n2lts)
		return false;


	int n1lts = dynamic_cast<BppInteger *> (n1->getNodeProperty(lts))->getValue();
	int n2uts = dynamic_cast<BppInteger *> (n2->getNodeProperty(uts))->getValue();

	if(n1lts > n2uts)
		return false;

	return true;

}


/*
Checks if the node are compatible (ie. same species and comparable timeslaice).
Presumes that both nodes have the same timeslice status as the tree

Takes:
 - n1 (Node *): a node with at least the species property set
 - n2 (Node *): a node with at least the species property set

Returns:
	(bool) : true if the two nodes are compatible; false otherwise

*/
bool ReconciledTree::areCompatible(Node  * n1, Node * n2)
{
	//first species
	if(!haveSameSpecies( n1,  n2))
		return false;
	return areTSCompatible(n1, n2);
}

/*
Checks if the node are compatible (ie. same species and comparable timeslice).
Presumes that both nodes have the same timeslice status as the tree

Takes:
 - id1 (int): a node id with at least the species property set
 - id2 (int): a node id with at least the species property set

Returns:
	(bool) : true if the two nodes are compatible; false otherwise

*/
bool ReconciledTree::areCompatible(int id1, int id2)
{
	Node * n1 = getNode(id1);
	Node * n2 = getNode(id2);
	return areCompatible(n1, n2);
}


/*
Checks if AncId is the ancestor of SonId

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
	(bool): true if AncId is an ancestor of SonId
*/
bool ReconciledTree::isAncestor(int AncId, int SonId)
{
	//solving the easiest cases
	if(AncId  == SonId)
		return false;

	int RootId = getRootId();

	if(AncId == RootId) // the potential ancestor is the root so it is ancestor to SonId
		return true;

	if(SonId == RootId) // The descendant is the root -> impossible
		return false;

	Node * SonNode = getNode(SonId);

	Node * FatherNode = SonNode->getFather();

	int FatherId = FatherNode->getId();

	bool ok= false; //storing result

	while(FatherId != RootId)
	{
		if(FatherId == AncId) // found the Ancestor
		{
			ok = true;
			break;
		}
		Node * SonNode = FatherNode;
		FatherNode = SonNode->getFather();		
		FatherId = FatherNode->getId();
	}
	return ok;
}


/*
Get the path from ancestor to son. (If AncId is not an ancestor of SonId, will return the path to the root)

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
	(vector <Node *>): path including start and stop
*/
vector <Node *> ReconciledTree::pathFromAncestorToSon(int AncId, int SonId)
{


	vector <Node *> path;
	path.push_back(getNode(SonId));

	int RootId = getRootId();
	int currentId = SonId;

	while((currentId != RootId)&&(currentId != AncId))
	{
		path.insert(path.begin(),path.front()->getFather());
		currentId = path.front()->getId();
	}

	return path;
}

/*
Get the path from ancestor to son. (If AncId is not an ancestor of SonId, will return the path to the root)

Takes:
 - AncId (int): id of the potential ancestor
 - SonId (int): id of the potential descendant

Returns:
	(vector <int >): path including start and stop
*/
vector <int> ReconciledTree::idPathFromAncestorToSon(int AncId, int SonId)
{
	vector <Node *> path = pathFromAncestorToSon(AncId, SonId);

	vector <int> idpath;

	for(unsigned i = 0; i < path.size(); i++)
		idpath.push_back(path[i]->getId());

	return idpath;
}


int ReconciledTree::getIdWithName(const string &name)
{
	//go through the whole tree. Does not throw error if a node has no name
	int IdToReturn = -1;

	vector <Node *> NodeToDo;
	Node * currentNode;

	NodeToDo.push_back(getRootNode());

	while(NodeToDo.size() > 0)
	{
		currentNode = NodeToDo[0];
		if(currentNode->hasName())
		{
			//cout << "  " << currentNode->getId() << " " << currentNode->getName() << " <-> " << name << endl;
			if(currentNode->getName() == name)
			{
				IdToReturn = currentNode->getId();
				break;
			}
		}
		//cout << "ReconciledTree::getIdWithName " << currentNode->getId() << " "<< currentNode->getNumberOfSons() <<endl;

		for(unsigned i = 0 ; i < currentNode->getNumberOfSons() ; i++)
		{
			//cout << i << endl;
			NodeToDo.push_back(currentNode->getSon(i));
		}
		NodeToDo.erase(NodeToDo.begin());
	}

	if(IdToReturn == -1)
		throw Exception("ReconciledTree::getIdWithName: no nodes with name \"" + name  + "\"");

	return IdToReturn;

}


bool ReconciledTree::isRealLeaf(int id)
{
	if(getNodeEvent(id) != C)
		return false;
	return true;
}

bool ReconciledTree::isExtant(int evtcode)
{
	if(evtcode != C)
		return false;
	return true;
}

bool ReconciledTree::isSpeciation(int evtcode)
{
	if(evtcode != S)
		return false;
	return true;
}

bool ReconciledTree::isLoss(int evtcode)
{
	if(evtcode != L)
		return false;
	return true;
}

bool ReconciledTree::isDup(int evtcode)
{
	if(evtcode != D)
		return false;
	return true;
}

bool ReconciledTree::isSout(int evtcode)
{
	if(evtcode != Sout)
		return false;
	return true;
}

bool ReconciledTree::isRec(int evtcode)
{
	if(evtcode != R)
		return false;
	return true;
}

bool ReconciledTree::isNull(int evtcode)
{
	if(evtcode != N)
		return false;
	return true;
}

bool ReconciledTree::isBout(int evtcode)
{
	if(evtcode != Bout)
		return false;
	return true;
}


ReconciledTree * ReconciledTree::cloneSubtree(int newRootId)
{
	Node * newRoot = TreeTemplateTools::cloneSubtree<Node>(*this, newRootId);
	return new ReconciledTree(newRoot, getTimeSliceStatus());
}

vector <string> ReconciledTree::getRealLeavesNames()
{
	vector <string>  names;

	vector <int> LeavesId = getLeavesId();

	for(size_t i = 0; i < LeavesId.size(); i++)
	{
		if( isRealLeaf(LeavesId[i]) )
			names.push_back(getNodeName(LeavesId[i]));

	}
	return names;
}


/*
recursively extract a MyGeneTree instance copying the topology of the reconciled tree.
*/
MyGeneTree  ReconciledTree::getMyGeneTreeTopology()
{
	MyGeneNode * node =  getMyGeneTreeTopologyAux(getRootNode());
	return MyGeneTree( *node );
}


map <int, vector <string> > ReconciledTree::cladeToLeafListAssociation()
{
	map <int, vector <string> >	cladeToLeafList;
	cladeToLeafListAssociationAux(getRootId(), cladeToLeafList );

	return cladeToLeafList;
}


/*
Takes:
	- int nodeId: the id of a node in the tree

Returns:
	(vector <int> ) the list of the closest descendant of nodeId that are either leaves or speciations
*/
vector <int> ReconciledTree::getLeafOrSpeciationDescendants(int nodeId)
{
	vector <int> descendants;

	if(isLeaf(nodeId)) // leaf -> no descendant to report
		return descendants;

	vector <int> SonsId = getSonsId(nodeId);

	for(unsigned i = 0 ; i < SonsId.size(); i++)
	{
		int ev = getNodeEvent(SonsId[i]);
		if( (isSpeciation(ev)) || (isExtant(ev)) )
			descendants.push_back(SonsId[i]); // speciation or leaf found -> no recursion
		else
		{ // not found -> recursion
			vector <int> SonsDescendants = getLeafOrSpeciationDescendants(SonsId[i]);
			for(unsigned j = 0 ; j < SonsDescendants.size(); j++) // adding the detected descendant to our list
				descendants.push_back(SonsDescendants[j]);
		}
	}
	return descendants;
}

/*
	Takes:
		- string name : name of the single leaf
		- int speciesId : species id of the single leaf
*/
void ReconciledTree::makeSingleLeafTree( string name, int speciesId )
{
	assert( nbNodes == 0 );

	Node * newnode  = new Node(0);
	setRootNode(newnode);//setting as the new root
	setNodeCladeNum(0, 0);//setting the new node cladenum property and updating the map

	nbNodes++;

	//setting NodeId as a leaf
	setNodeSpecies(newnode,speciesId);
	newnode->setName(name);
	setNodeEvent(newnode,C);
	setNodeUpperBoundaryTS(newnode,0);//the time slice of a leaf is 0
	setNodeLowerBoundaryTS(newnode,0);

}