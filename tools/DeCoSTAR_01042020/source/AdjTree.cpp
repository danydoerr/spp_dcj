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

Last modified the: 13-07-2017
by: Wandrille Duchemin

*/

#include "AdjTree.h"


pair <string, string> AdjTree::leafNamesfromAname(string Aname)
{
	pair <string,string> leafNames;

	size_t pos = Aname.find_first_of("-");
	if( pos != string::npos) // the "-" was found
	{
		try
		{
			leafNames.first = Aname.substr(0,pos);
			leafNames.second = Aname.substr(pos +1 , Aname.size() - pos - 1);
		}
		catch(exception & e)
		{
			cerr << "AdjTree::leafNamesfromAname. Troubles converting the name : " << Aname << " into two ids"<<endl;
			cerr << e.what() << endl;
			exit(-1);
		}
		
	}
	else
	{
		cerr << "AdjTree::leafNamesfromAname. Troubles converting the name : " << Aname << " into two ids"<<endl;
		exit(1);
	}

	return leafNames;
}
/*
Takes:
 - int currentId : an id in the adj tree
 - int gfamIndex : wether we are interested in the ids of gfam1 or gfam2
 - &vector < vector <string> > LeafListList : a list that will contains the leaf list clade representation of all nodes that are C1 (as they are in the tree they are C1) and whose clade size is > 1
Returns;
	(vector <string>) : leaf representation of the clade of currentId
*/
vector <string> AdjTree::getC1LeafListAux( int currentId, int gfamIndex, vector < vector <string> > & LeafListList )
{
	//cout << "AdjTree::getC1LeafListAux " << currentId << " " << gfamIndex << "->" << LeafListList.size() << endl;
	//cout << "AdjTree::getC1LeafListAux " << currentId << " evt:" << getNodeEvent(currentId)<<endl;
	vector <string> LeafList;
	if(isLeaf(currentId))
	{
		if(isExtant(getNodeEvent(currentId)))
		{
			pair<string, string> leafNames =  leafNamesfromAname(getNodeName(currentId));
			if(gfamIndex == 0)
				LeafList.push_back( leafNames.first );
			else
				LeafList.push_back( leafNames.second );
		}
		//cout << "leaf ->" << LeafList.size() << " names"<<endl;
		return LeafList;
	}

	int nbNonEmptyChilden = 0;

	vector <int> sonsId = getSonsId(currentId);
	for(unsigned i = 0 ; i < sonsId.size(); i++)
	{
		vector <string> sonsLeafList = getC1LeafListAux( sonsId[i],  gfamIndex,  LeafListList );
		//cout << "son "<<i<< " yielded clade of size "<< sonsLeafList.size() << endl;
		if(sonsLeafList.size() > 0 )
		{
			nbNonEmptyChilden++;
			for(unsigned j = 0 ; j < sonsLeafList.size() ; j++)
			{
				LeafList.push_back(sonsLeafList[j]);
			}
		}
	}
	//cout << currentId <<" nbNonEmptyChilden "<< nbNonEmptyChilden << endl;
	//cout << currentId <<" LeafList "<< LeafList.size() ;
	//for(unsigned i = 0 ; i < LeafList.size(); i++)
	//	cout << " " << LeafList[i];
	//cout << endl;
	if(nbNonEmptyChilden > 1) // more than one non empty child -> new clade -> add it 
	{
		LeafListList.push_back(LeafList);
	}
	return LeafList;
}


/*
Takes:
	- n (Node *): pointer to a node in the tree

Returns:
	(int): number of AdjBreak events in the tree

*/
int AdjTree::countNbBreakAux(Node * n)
{
	int nbBreak = 0;

	
	int evt = dynamic_cast<BppInteger *> (n->getNodeProperty(ev))->getValue();


	int id = n->getId();
	if(isLivingAdjBreak(id,evt))
		nbBreak += 1;

		

	vector <Node *> sons = n->getSons();
	for(unsigned i = 0 ; i < sons.size(); i++)
		nbBreak += countNbBreakAux(sons[i]);

	return nbBreak;
}

/*
Takes:
	- n (Node *): pointer to a node in the tree

Returns:
	(int): number of Gain events in the tree

*/
int AdjTree::countNbGainAux(Node * n)
{
	int nbGain = 0;

	bool iscoev = false;

	if(n->hasNodeProperty(coev))
		if(dynamic_cast<BppInteger *> (n->getNodeProperty(coev))->getValue())
			iscoev = true;


	vector <Node *> sons = n->getSons();

	if((sons.size() > 1) && (!iscoev)) // a non-coevent that has more than 1 son -> 1 break per son above 1
		iscoev += (sons.size() - 1);


	for(unsigned i = 0 ; i < sons.size(); i++)
		nbGain += countNbGainAux(sons[i]); // recurse

	return nbGain;
}


void AdjTree::printNode(Node * n)
{
	if( n->hasName())
		cout << n->getName() << " ";

	cout << "ev : ";
	string prop = ev;
	if(!n->hasNodeProperty(prop))
		cout << "no event";
	else
		cout << dynamic_cast<BppInteger *> (n->getNodeProperty(prop))->getValue();

	cout << " (";
	prop = coev;
	if(!n->hasNodeProperty(prop))
		cout << "0";
	else
		cout << dynamic_cast<BppInteger *> (n->getNodeProperty(prop))->getValue();
	cout << ")";

	cout << " species: ";
	prop = spe;
	if(!n->hasNodeProperty(prop))
		cout << "no species";
	else
		cout << dynamic_cast<BppInteger *> (n->getNodeProperty(prop))->getValue();
	

	cout << " nodes :";
	prop = nodeid1;
	if(!n->hasNodeProperty(prop))
		cout << "no nodes";
	else
		cout << dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid1))->getValue() << "," << dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid2))->getValue();

	cout << endl;

}

void AdjTree::printNodeRec(Node * n,int level)
{
	for(unsigned i = 0; i < level ; i++)
		cout << "  ";
	cout << "|-";
	printNode(n);

	for(unsigned i = 0 ; i < n->getNumberOfSons();i++)
		printNodeRec(n->getSon(i), level + 1);
}

string AdjTree::NodeString(Node * n, bool hideLosses)
{
	string Str = "";


	int nbsons = n->getNumberOfSons();
	if(nbsons> 0)
	{

		vector <string> sonStr;

		for(unsigned i = 0; i < nbsons; i++)
		{
			string tmp =  NodeString(n->getSon(i), hideLosses); // recursion
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

	//we only want the name of real leaves
	bool isRealLeaf = false;
	if(n->hasNodeProperty(ev))
		if( isExtant( dynamic_cast<BppInteger *> (n->getNodeProperty(ev))->getValue() ))
			isRealLeaf = true;


	//printing node properties
	if( n->hasName() )
	{
		Str += n->getName();

	}
	else
	{
		if(!n->hasNodeProperty(nodeid1))
			Str += "NONAME";
		else
		{
			int Number = dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid1))->getValue();
			Str += static_cast<ostringstream*>( &(ostringstream() << Number) )->str() ;
			Str += "-" ;
			Number = dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid2))->getValue();
			Str += static_cast<ostringstream*>( &(ostringstream() << Number) )->str() ;
		}
	}





	string evtStr = "|";

	if(!n->hasNodeProperty(ev))
		evtStr += "NOEVENT";
	else
	{
		bool COEVENT = false;
		if(n->hasNodeProperty(coev))
			if(dynamic_cast<BppInteger *> (n->getNodeProperty(coev))->getValue())
				COEVENT = true;

		int evt = dynamic_cast<BppInteger *> (n->getNodeProperty(ev))->getValue();
		switch(evt)
		{ 
			case C :
				evtStr += "Extant"; //Current (leaf) <- no event text, for DeCo retro-compatibility <- no more DeCo retro comp
				break;
			case S :
				evtStr += "Spe"; //Speciation (in the species tree) <-changed for deCo retro-compatibility <- no more DeCo retro comp
				break;
			case L :
				if(hideLosses)
					return ""; // return empty str if losses are hidden

				if(COEVENT)
					evtStr += "co-";
				evtStr += "Loss"; //Loss (in the species tree)
				break;
			case D :
				if(COEVENT)
					evtStr += "co-";
				evtStr += "Dup"; //Duplication (in the species tree)
				break;
			case Sout :
				if(COEVENT)
					evtStr += "co-";
				evtStr += "SpeOut"; //Speciation to an extinct/unsampled lineage (otherwise called SpeciationOut)
				break;
			case R :
				if(COEVENT)
					evtStr += "co-";
				evtStr += "Reception"; //Transfer reception
				break;
			case N :
				if(COEVENT)
					evtStr += "co-";
				evtStr += "Null"; //no event (to account for time slices)
				break;
			case Bout :
				if(COEVENT)
					evtStr += "co-";
				evtStr += "BifOut"; //Bifurcation in an extinct/unsampled lineage (otherwise called BifurcationOut)
				break;
			case Abreak :
				evtStr += "Breakage";//adjacency break <- changed fofr DeCo retro-compatibility
				break;
			default: // unknown event
				evtStr += "NOEVENT";
		}
	}

	Str += evtStr;

	Str += "|";

	if(!n->hasNodeProperty(spe))
		Str += "NOSPECIES";
	else
	{
		int Number = dynamic_cast<BppInteger *> (n->getNodeProperty(spe))->getValue();
		Str += static_cast<ostringstream*>( &(ostringstream() << Number) )->str() ;
	}

	return Str;
}


void AdjTree::addABreakInDeadNode(int Nodeid)
{
	setNodeSpecies(Nodeid, -1 );
	setNodeCoEventStatus(Nodeid, false );
	setNodeEvent(Nodeid, Abreak);
}

/*
Takes:
	- Nodeid (int): a node id
	- evtLine (string): recPhyloXML line describing a reconciliation event
*/
void AdjTree::setNodeDetailsFromXMLLine(int Nodeid, string evtLine, bool VERBOSE)
{
	map <string, string> * properties = new map <string,string>;
	string * value = new string(""); 
	string Tname = InterpretLineXML(evtLine,properties,value);

	setNodeDetailsFromXMLLine( Nodeid,  Tname,  properties,  value, VERBOSE);
}

//same as previous; uses the interpreted xml line instead
void AdjTree::setNodeDetailsFromXMLLine(int Nodeid, string Tname, map <string, string> * properties, string * value, bool VERBOSE)
{
	char * pEnd; // for str -> int conversion 
	//cout << "in AdjTree::setNodeDetailsFromXMLLine for clade " << Nodeid << endl;
	
	for (map<string,string>::iterator it=properties->begin(); it!=properties->end(); ++it)
	{
		if (it->first.compare("speciesLocation") == 0)
			setNodeSpecies(Nodeid, (int) strtol(it->second.c_str(), &pEnd, 10) );
		else if (it->first.compare("destinationSpecies") == 0)
			setNodeSpecies(Nodeid, (int) strtol(it->second.c_str(), &pEnd, 10) );
		else if (it->first.compare("adjName") == 0)
			setNodeName(Nodeid,it->second);
		else if (it->first.compare("coevent") == 0)
			setNodeCoEventStatus(Nodeid, (int) strtol(it->second.c_str(), &pEnd, 10) );
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
Abreak  --> adjacency breakage
*/

	int evtCode = -1;
	//3. finding which event is implied
	if(Tname.compare("leaf") == 0)
	{
		evtCode = C;
	}
	else if(Tname.compare("speciation") == 0)
	{
		evtCode = S;
	}
	else if(Tname.compare("speciationOut") == 0)
	{
		evtCode = Sout;
	}
	else if(Tname.compare("branchingOut") == 0)
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
	else if(Tname.compare("transferLoss") == 0)
	{
		evtCode = R;
	}
	else if(Tname.compare("loss") == 0)
	{
		evtCode = L;
	}
	else if(Tname.compare("adjBreak") == 0)
	{
		evtCode = Abreak;
	}

	if(evtCode == -1)
		throw Exception("ReconciledTree::setNodeDetailsFromXMLLine : unknown event : " + Tname );
	else
		setNodeEvent(Nodeid,evtCode);

	if( ( hasNodeName(Nodeid) ) && ( !isExtant( evtCode ) ) && ( !isAdjBreak(evtCode) ) )
	{
		// the node name is formed of the two ids separated by "-"
		string name = getNodeName(Nodeid);

		size_t pos = name.find_first_of("-");

		if( pos != string::npos) // the "-" was found
		{
			try
			{
				string sub1 = name.substr(0,pos);
				string sub2 = name.substr(pos +1 , name.size() - pos - 1);
				int n1 = (int) strtol(sub1.c_str(), &pEnd, 10);
				int n2 = (int) strtol(sub2.c_str(), &pEnd, 10);

				setNodeNodeIds(Nodeid , pair<int,int> (n1,n2) );

			}
			catch(exception & e)
			{
				cerr << "ReconciledTree::setNodeDetailsFromXMLLine. Troubles converting the name : " << name << " into two int ids"<<endl;
				cerr << e.what() << endl;
				exit(-1);
			}

			
		}


	}


	if(VERBOSE)
	{
		cout << "set details of node " << Nodeid << ";" << evtCode << ";" ;
		if(hasNodeProperty(Nodeid,spe))
			cout  << getNodeSpecies(Nodeid) << ";";
		cout <<getNodeCoEventStatus(Nodeid) ;

		if(hasNodeProperty(Nodeid,nodeid1 ) )
		{
			pair <int,int> nids = getNodeNodeIds(Nodeid);
			cout << ";" << nids.first << "-" << nids.second ;
		}
		cout << endl ;
	}
}


/*
read and add a clade from a recPhyloXML stream

Takes:
	- fileIN (ifstream): streaam to a recPhyloXML file
	- VERBOSE (bool) (default = false)
*/
void AdjTree::addXMLClade(ifstream& fileIN, int Nodeid, bool VERBOSE) // will need a species tree
{

	if(VERBOSE)
		cout << "parsing clade " << " id " << Nodeid<< endl;

	vector <string> evtLines;

	//2. parsing of the clade
	string line;
	getline( fileIN, line );

	map <string, string> * properties = new map <string,string>;
	string * value = new string(""); 
	string Tname = InterpretLineXML(line,properties,value);

	

	while(Tname.compare("/clade") != 0) //--> line signing the end of the clade
	{

		//cout << line <<endl;

		if(Tname.compare("name") == 0) // found name
		{
			if(VERBOSE)
				cout << "found name : " << *value <<  endl;
			setNodeName(Nodeid, *value);

		}
		else if(Tname.compare("eventsAdj") == 0) // reading reconciliation events
		{
			if(VERBOSE)
				cout <<  "parsing events" << endl;
			
			*value = "";
			delete properties;

			getline( fileIN, line );
			properties = new map <string,string>;
			Tname = InterpretLineXML(line,properties,value);

			while(Tname.compare("/eventsAdj") != 0) // here, add some nodes ...
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
			int newId = getNextId();

			Node * newNode  = new Node(newId);

			//branching to the existing structure
			Node * currentNode  = getNode(Nodeid);
			currentNode->addSon(newNode);

			//recursing
			addXMLClade(fileIN, newId,  VERBOSE );


		}
		if(Tname.compare("/phylogeny") == 0) // found end of the phylogeny... should not happen but we never know
		{
			if(VERBOSE)
				cout << "found the end of the phylogeny while parsing a clade. Missing </clade> tag? Closing the clade prematurely" << endl;
			break;
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
		if(VERBOSE)
			cout << "putting Abreak in dead for node " << Nodeid << endl;
		addABreakInDeadNode(Nodeid);
		//throw Exception("AdjTree::addXMLClade : no adj events found");
	}
	else
	{
		if(evtLines.size()>1)
		{
			for(unsigned i = 0; i < evtLines.size() - 1; i++ ) // all events but the last
			{
				//creating parent node 
				
	
				int newId = getNextId();
	
				Node * newNode  = new Node(newId);
	
	
				//branching to the existing structure
				Node * currentNode  = getNode(Nodeid);
	
				bool ROOT = isRoot(Nodeid);
	
				if(!ROOT)
				{
					Node * FNode = currentNode->getFather();
					FNode->removeSon(currentNode);
					FNode->addSon(newNode);
	
					newNode->addSon(currentNode);
				}
				else //we need to re-root. the manipulation is a bit different
				{
					currentNode->addSon(newNode);
					
					rootAt(newNode); // rerooting
				}
	
				//cout << "creating evt for clade" << newId << " : " << evtLines[i] << endl;
				setNodeDetailsFromXMLLine(newId,evtLines[i], VERBOSE);
	
			}
		}
		//cout << "creating final evt for clade" << Nodeid << " : " << evtLines.back() << endl;
	
		setNodeDetailsFromXMLLine(Nodeid,evtLines.back(), VERBOSE); // setting details on the last node == the bifurcation or leaf


	}


		

}

void AdjTree::readPhyloXMLFile(ifstream& fileIN, bool VERBOSE )
{
	// presuming that the stream will yield a <clade> line froming the root of the tree
	//first find the root
	if(VERBOSE)
		cout << "Creating AdjTree." << endl;

	int rootid = 0;

	Node * newnode  = new Node(rootid);

	setRootNode(newnode);//setting as the new root


	addXMLClade( fileIN,rootid,  VERBOSE ); // recursive call to the entire tree.

}
















AdjTree::AdjTree(Node * root, int rootgain, int gfam1, int gfam2): TreeTemplate<Node>(root)
{
	RootIsGain = rootgain;
	Gfamily1 = gfam1;
	Gfamily2 = gfam2;
}



/*
Takes:
 - ifstream& fileIN : stream to a file containing an adjacency tree
 - bool rootgain : wether there is autoratimcally a gain at the root of the tree
 - int gfam1 : id of the first gene family
 - int gfam2 : id of the second gene family
 - bool VERBOSE  [default:false]

*/
AdjTree::AdjTree(ifstream& fileIN,bool rootgain, int gfam1, int gfam2, bool VERBOSE)
{
	RootIsGain = rootgain;
	Gfamily1 = gfam1;
	Gfamily2 = gfam2;

	readPhyloXMLFile( fileIN,  VERBOSE );

}


/*
Returns:
	(int): Gfamily 1

*/
int AdjTree::getGfamily1()
{
	return Gfamily1;
}

/*
Returns:
	(int): Gfamily 2

*/
int AdjTree::getGfamily2()
{
	return Gfamily2;
}

/*
Returns:
	(bool): true if the root is an adjacency gain event
*/
bool AdjTree::isRootGain()
{
	return RootIsGain;
}

/*
Takes:
	- newval (int): desired Gfamily1 value
*/
void AdjTree::setGfamily1(int newval)
{
	Gfamily1 = newval;
}

/*
Takes:
	- newval (int): desired Gfamily2 value
*/
void AdjTree::setGfamily2(int newval)
{
	Gfamily2 = newval;
}

/*
Takes:
	- newval (bool): desired RootIsGain value
*/
void AdjTree::setRootGain(bool newval)
{
	RootIsGain = newval;
}

/*
Takes:
	- id (int)

Returns:
	(int): the species the node is in
*/
int AdjTree::getNodeSpecies(int id)
{
	return dynamic_cast<BppInteger *> (getNodeProperty(id, spe))->getValue();
}

/*
Takes:
	- id (int)

Returns:
	(pair < int , int >): the nodes ids corresponding to this node in respectively Gfamily1/2 tree.
*/
pair <int,int> AdjTree::getNodeNodeIds(int id)
{
	pair <int,int> res;
	res.first =  dynamic_cast<BppInteger *> (getNodeProperty(id, nodeid1))->getValue();
	res.second = dynamic_cast<BppInteger *> (getNodeProperty(id, nodeid2))->getValue();
	return res;
}

/*
Takes:
	- id (int)

Returns:
	(bool): true if the node event concern both extremity of the adjacency (i.e: it is a co-event)
*/
bool AdjTree::getNodeCoEventStatus(int id)
{
	if(!hasNodeProperty(id,coev))
		return false;
	return (dynamic_cast<BppInteger *> (getNodeProperty(id, coev))->getValue() == 1);
}

/*
Takes:
	- id (int)

Returns:
	(int): the defining event of the node
*/
int AdjTree::getNodeEvent(int id)
{
	return dynamic_cast<BppInteger *> (getNodeProperty(id, ev))->getValue();
}



/*
Takes:
	- Node * n

Returns:
	(int): the species the node is in
*/
int AdjTree::getNodeSpecies(Node * n)
{
	return dynamic_cast<BppInteger *> (n->getNodeProperty(spe))->getValue();
}

/*
Takes:
	- Node * n

Returns:
	(pair < int , int >): the nodes ids corresponding to this node in respectively Gfamily1/2 tree.
*/
pair <int,int> AdjTree::getNodeNodeIds(Node * n)
{
	pair <int,int> res;
	res.first =  dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid1))->getValue();
	res.second = dynamic_cast<BppInteger *> (n->getNodeProperty(nodeid2))->getValue();
	return res;
}

/*
Takes:
	- Node * n

Returns:
	(bool): true if the node event concern both extremity of the adjacency (i.e: it is a co-event)
*/
bool AdjTree::getNodeCoEventStatus(Node * n)
{
	if(!n->hasNodeProperty(coev))
		return false;
	return (dynamic_cast<BppInteger *> (n->getNodeProperty(coev))->getValue() == 1);
}

/*
Takes:
	- Node * n

Returns:
	(int): the defining event of the node
*/
int AdjTree::getNodeEvent(Node * n)
{
	return dynamic_cast<BppInteger *> (n->getNodeProperty(ev))->getValue();
}




/*

Takes:
	- id (int)
	- sp (int): the species the node is in
*/
void AdjTree::setNodeSpecies(int id, int sp)
{
	setNodeProperty(id, spe, BppInteger(sp));
}

/*
Takes:
	- id (int)
	- nodeids (pair < int , int >): the nodes ids corresponding to this node in respectively Gfamily1/2 tree.
*/
void AdjTree::setNodeNodeIds(int id, pair <int,int> nodeids)
{
	setNodeProperty(id, nodeid1, BppInteger(nodeids.first));
	setNodeProperty(id, nodeid2, BppInteger(nodeids.second));
}

/*
Takes:
	- id (int)
	- co (bool): true if the node event concern both extremity of the adjacency (i.e: it is a co-event)
*/
void AdjTree::setNodeCoEventStatus(int id, bool co)
{
	if(co)
		setNodeProperty(id, coev, BppInteger(1));
	else
		setNodeProperty(id, coev, BppInteger(0));
}

/*
Takes:
	- id (int)
	- e (int): the defining event of the node
*/
void AdjTree::setNodeEvent(int id, int e)
{
	setNodeProperty(id, ev, BppInteger(e));
}

/*
Takes:
	- nodeids (pair < int , int >): the nodes ids corresponding to this node in respectively Gfamily1/2 tree.
	- sp (int): the species the node is in
	- co (bool): true if the node event concern both extremity of the adjacency (i.e: it is a co-event)
	- e (int): the defining event of the node
*/
Node * AdjTree::createNode(pair <int,int> nodeids, int sp, bool co, int e)
{
	Node * n = new Node();

	n->setNodeProperty(spe, BppInteger(sp));
	n->setNodeProperty(nodeid1, BppInteger(nodeids.first));
	n->setNodeProperty(nodeid2, BppInteger(nodeids.second));
	n->setNodeProperty(coev, BppInteger(co));
	n->setNodeProperty(ev, BppInteger(e));

	//if( isLoss(e) )
	//{
	//    cout << "created node : " << sp << " " << nodeids.first << "-"<< nodeids.second << " ev: "<< e << endl;
	//}

	return n;
}


Node * AdjTree::createBreakNode()
{
	Node * n = new Node();

	n->setNodeProperty(ev, BppInteger(Abreak));

	//cout << "created Break node."<< endl;
	return n;
}

/*
Annotates the tree leaves with leaf names and species for break nodes

Takes:
 - Rtree1 (ReconciledTree *)
 - Rtree2 (ReconciledTree *)
 - int sens1 : indicates which extremity of the first gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
 - int sens2 : indicates which extremity of the second gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)

*/
void AdjTree::refine(ReconciledTree * Rtree1, ReconciledTree * Rtree2, int sens1 , int sens2 )
{
	vector <int > leaves = getLeavesId();

	for(unsigned i = 0 ; i < leaves.size(); i++)
	{
		int it = leaves[i];


		//cout << "AdjTree::refine " <<  it << endl;

		if(getNodeEvent(it) == C) //this is a real extant leaf --> create a name
		{
			pair <int,int> nids = getNodeNodeIds(it);
			pair <string,string> names("","");

			if(Rtree1->hasNodeName(nids.first))
				names.first = Rtree1->getNodeName(nids.first);

			if(Rtree2->hasNodeName(nids.second))
				names.second = Rtree2->getNodeName(nids.second);

			string name  = names.first;
			if(sens1 == -1)
				name += "_start";
			else if(sens1 == 1)
				name += "_stop";
			name += "-";
	
			name += names.second;

			if(sens2 == -1)
				name += "_start";
			else if(sens2 == 1)
				name += "_stop";

			setNodeName(it, name);
		}
		else if(getNodeEvent(it) == L) // loss Node
		{
			continue; // keep the two nids
		}
		else if(getNodeEvent(it) == Abreak) // Adj break node
		{
			setNodeName(it, static_cast<ostringstream*>( &(ostringstream() << i) )->str());//the name is the leaf id
		}

	}
}

/*
modify the names of each node to make them account for gene extremities

Takes:
	- int sens1 : indicates which extremity of the first gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- int sens2 : indicates which extremity of the second gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- Node * current : current node in the recursion
*/
void AdjTree::setNodeNamesAux(int sens1 , int sens2, Node * current)
{

	//cout << "AdjTree::setNodeNamesAux" <<  current->getId() << " as a name? " << current->hasName() << endl;
	//
	//if( current->hasName() )
	//	cout << "  name is : " << current->getName() << endl;

	if( !current->hasName() )
	{
		pair <int,int> nids = getNodeNodeIds(current);

		stringstream ss1;
		ss1 << nids.first;
		string name = ss1.str();

		if(sens1 == -1)
			name += "_start";
		else if(sens1 == 1)
			name += "_stop";
		name += "-";

		stringstream ss2;
		ss2 << nids.second;
		name += ss2.str();

		if(sens2 == -1)
			name += "_start";
		else if(sens2 == 1)
			name += "_stop";
		
		current->setName(name);

	}

	//if( current->hasName() )
	//	cout << "  name after is : " << current->getName() << endl;


	int nbSons = current->getNumberOfSons();

	for(unsigned i = 0 ; i < nbSons ; i++)
		setNodeNamesAux( sens1 , sens2, current->getSon(i) );

}

/*
modify the names of each node to make them account for gene extremities

Takes:
	- int sens1 : indicates which extremity of the first gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- int sens2 : indicates which extremity of the second gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
*/
void AdjTree::setNodeNames(int sens1 , int sens2)
{
	//cout << "AdjTree::setNodeNames "<< sens1 << " " << sens2 << endl;
	setNodeNamesAux( sens1 , sens2, getRootNode());
}


void AdjTree::printMe()
{
	printNodeRec(getRootNode(),0);
}

string AdjTree::NewickString(bool hideLosses)
{
	return NodeString(getRootNode(), hideLosses) ;
}


bool AdjTree::isExtant(int evtcode)
{
	if(evtcode != C)
		return false;
	return true;
}

bool AdjTree::isSpeciation(int evtcode)
{
	if(evtcode != S)
		return false;
	return true;
}

bool AdjTree::isLoss(int evtcode)
{
	if(evtcode != L)
		return false;
	return true;
}

bool AdjTree::isDup(int evtcode)
{
	if(evtcode != D)
		return false;
	return true;
}

bool AdjTree::isSout(int evtcode)
{
	if(evtcode != Sout)
		return false;
	return true;
}

bool AdjTree::isRec(int evtcode)
{
	if(evtcode != R)
		return false;
	return true;
}

bool AdjTree::isNull(int evtcode)
{
	if(evtcode != N)
		return false;
	return true;
}

bool AdjTree::isBout(int evtcode)
{
	if(evtcode != Bout)
		return false;
	return true;
}

bool AdjTree::isAdjBreak(int evtcode)
{
	if(evtcode != Abreak)
		return false;
	return true;
}

/*
returns true if the event is AdjBreak AND the species is not -1
*/
bool AdjTree::isLivingAdjBreak(int Nodeid, int evtcode)
{

	if(evtcode == -1) // default value -> we search for ourselves
	{
		evtcode = getNodeEvent(Nodeid);
	}

	if(! isAdjBreak(evtcode) )
	{
		return false;
	}


	if(hasNodeProperty(Nodeid,spe))
	{
		return false; // alive break nodes don't have a given species
	}
	return true;
}

/*
Returns:
	(int): number of AdjBreak events in the tree

*/
int AdjTree::countNbBreak()
{
	return countNbBreakAux(getRootNode());
}

/*
Returns:
	(int): number of Gain events in the tree

*/
int AdjTree::countNbGain()
{
	int nbGain = 0;
	if(RootIsGain)
		nbGain +=1;

	nbGain += countNbGainAux(getRootNode());

	return nbGain;
}

AdjTree * AdjTree::cloneSubtree(int newRootId)
{
	Node * newRoot = TreeTemplateTools::cloneSubtree<Node>(*this, newRootId);
	return new AdjTree(newRoot, RootIsGain, Gfamily1, Gfamily2);
}

/*
Takes:
 - int currentId : an id in the adj tree
 - int gfamIndex : wether we are interested in the ids of gfam1 or gfam2
 - &vector < vector <string> > LeafListList : a list that will contains the leaf list clade representation of all nodes that are C1 (as they are in the tree they are C1) and whose clade size is > 1

*/
void AdjTree::getC1LeafList(int gfamIndex, vector < vector <string> > & LeafListList )
{
	getC1LeafListAux(getRootId() , gfamIndex,  LeafListList );
}

/*
Takes:
 - int gfamIndex : wether we are interested in the ids of gfam1 or gfam2

Returns:
	vector <int> : vctor of nodes id in a reconciled tree corresponding to adjacencies in the adjacency tree that aren't leaves
*/
vector <int> AdjTree::getC1NodeIds(int gfamIndex)
{
	map<int,bool> mapNodeId;

	vector <int> AdjNodes = getInnerNodesId();

	for(unsigned i = 0 ; i < AdjNodes.size(); i++)
	{ // internal nodes should always have the nid property
		pair <int,int> nids = getNodeNodeIds(AdjNodes[i]);
		if(gfamIndex == 0)
			mapNodeId[nids.first] = true;
		else
			mapNodeId[nids.second] = true;
	}

	vector <int> geneIds;
	for (map<int,bool>::iterator it=mapNodeId.begin(); it!=mapNodeId.end(); ++it)
	{
		geneIds.push_back(it->first);
	}
	return geneIds;
}