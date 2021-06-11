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

This file contains a class that manages DeCo's outputs

Created the: 06-01-2016
by: Wandrille Duchemin

Last modified the: 04-09-2016
by: Wandrille Duchemin

*/

#include "DeCoOutputManager.h"

/*
Takes:
	- OUT (ofstream) : an opened stream
	- indent_level (int): number of double-spaces at the beginning of a line
*/
void DeCoOutputManager::beginLine(ofstream& OUT, int indent_level)
{
	for(unsigned i = 0; i< indent_level ; i++)
		OUT << "  ";
}
	
/*
Writes the event of the reconciled node

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- nodeid (int): id of a node
	- indent_level (int): number of double-spaces at the beginning of a line
	- hasLoss (bool): true if one of nodeids child is a loss; false otherwise
*/
void DeCoOutputManager::WritePhyloXMLRecEvent(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level, bool hasLoss)
{
	int event = Rtree->getNodeEvent(nodeid);

	if(Rtree->isNull(event)) // if the event is a Null event it won't be written
		return;

//	cout << nodeid << " "<< indent_level << endl;

	beginLine(OUT,indent_level);
	OUT << "<";

	string eventString = "";
	string SpeciesString = "speciesLocation=";

	if(Rtree->isExtant(event))
	{
		eventString = "leaf";
	}
	else if(Rtree->isSpeciation(event))
	{
		eventString	= "speciation";
	}
	else if(Rtree->isLoss(event))
	{
		eventString	= "loss";
	}
	else if(Rtree->isDup(event))
	{
		eventString	= "duplication";
	}
	else if(Rtree->isSout(event))
	{
		eventString	= "branchingOut";
	}
	else if(Rtree->isRec(event))
	{
		eventString	= "transferBack";
		SpeciesString = "destinationSpecies=";		
	}
	else if(Rtree->isBout(event))
	{
		eventString	= "bifurcationOut";
	}
	else
	{
		eventString	= "unknownEvent";
	}

	OUT << eventString;

	if(hasLoss)
		OUT << "Loss";

	if(!Rtree->isBout(event))
	{
		OUT << " ";
		OUT << SpeciesString;
		OUT << "\"";
		OUT << Rtree->getNodeSpecies(nodeid);
		OUT << "\"";
	}

	if(Rtree->isExtant(event))
		if(Rtree->hasNodeName(nodeid))
			OUT << " geneName=\"" << Rtree->getNodeName(nodeid) << "\"" ;

	//... time slice if any.
	if((Rtree->getTimeSliceStatus() != 0) && (eventString != "leaf"))
	{
		OUT << " timeSlice=\"";
		OUT <<Rtree->getNodeTimeSlice(nodeid);
		OUT << "\"";
	}

	OUT << "></";
	OUT << eventString;
	if(hasLoss)
		OUT << "Loss";

	OUT << ">\n";

	return;
}


/*
Writes the reconciled node along with its children (recursively)

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- nodeid (int): id of a node
	- indent_level (int): number of double-spaces at the beginning of a line
*/
void DeCoOutputManager::WritePhyloXMLRecTreeAux(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level)
{
	
	int event = Rtree->getNodeEvent(nodeid);
	if(Rtree->isNull(event))// ignoring any Null event.
	{
		vector <int> children = Rtree->getSonsId(nodeid);
		WritePhyloXMLRecTreeAux(OUT, Rtree, children[0], indent_level);//recursion on the first (and only) son of the null node
	}
	else
	{

		beginLine(OUT,indent_level);
		OUT << "<clade>\n";
		indent_level++;
	
		//writing the different informations:
	
		//1. name:
		beginLine(OUT,indent_level);
		OUT << "<name>";
		if(Rtree->hasNodeName(nodeid))
		{
			OUT << Rtree->getNodeName(nodeid);
		}
		else
		{
			OUT << nodeid;
		}
		OUT << "</name>\n";
	
		//2. the reconciliation event(s)
		beginLine(OUT,indent_level);
		OUT << "<eventsRec>\n";
		indent_level++;
	
		int currentId = nodeid;
		int nextId = currentId;
		bool hasLoss = false;
		bool donext = true;
	
		while( donext )
		{
			currentId = nextId; //iteration
			hasLoss = false;
	
			vector <int> children = Rtree->getSonsId(currentId);
	
			if(children.size() == 0)
				donext = false;
	
			for(unsigned i = 0 ; i < children.size(); i++)
			{
				int event = Rtree->getNodeEvent(children[i]);
				if(Rtree->isLoss(event))
				{
					hasLoss = true;
				}
				else
				{
					if(nextId  == currentId)//this is the first non-loss child of current id -> we accept it as nextid
						nextId = children[i];
					else //this is not the first non-loss child of current id -> we get out of the loop
					{
						donext = false;
						break;
					}
				}
			}
			//we write the event of the current node
			WritePhyloXMLRecEvent(OUT,  Rtree, currentId, indent_level,  hasLoss);
		}
	
		indent_level--;
		beginLine(OUT,indent_level);
		OUT << "</eventsRec>\n";
	
	
		//3. Sons, if any
		vector <int> children = Rtree->getSonsId(currentId);
			
		for(unsigned i = 0 ; i < children.size(); i++)
		{
			WritePhyloXMLRecTreeAux( OUT,  Rtree, children[i], indent_level);
		}
	
		indent_level--;
		beginLine(OUT,indent_level);
		OUT << "</clade>\n";
	}
}

/*
Writes the event of the reconciled node

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- nodeid (int): id of a node
	- indent_level (int): number of double-spaces at the beginning of a line

*/
void DeCoOutputManager::WritePhyloXMLAdjEvent(ofstream& OUT, AdjTree * Atree, int nodeid, int indent_level)
{
	int event = Atree->getNodeEvent(nodeid);

	if(Atree->isNull(event)) // if the event is a Null event it won't be written
		return;

//	cout << nodeid << " "<< indent_level << endl;
	
	bool coevent = Atree->getNodeCoEventStatus(nodeid);
	



	beginLine(OUT,indent_level);
	OUT << "<";

	string eventString = "";
	string SpeciesString = "speciesLocation=";

	if(Atree->isExtant(event))
	{
		eventString = "leaf";
	}
	else if(Atree->isSpeciation(event))
	{
		eventString	= "speciation";
	}
	else if(Atree->isLoss(event))
	{
		eventString	= "loss";
	}
	else if(Atree->isDup(event))
	{
		eventString	= "duplication";
	}
	else if(Atree->isSout(event))
	{
		eventString	= "branchingOut";
	}
	else if(Atree->isRec(event))
	{
		eventString	= "transferBack";
		SpeciesString = "destinationSpecies=";
	}
	else if(Atree->isBout(event))
	{
		eventString	= "bifurcationOut";
	}
	else if(Atree->isAdjBreak(event))
	{
		eventString	= "adjBreak";
	}
	else
	{
		eventString	= "unknownEvent";
	}


	OUT << eventString;

	if(coevent)
		OUT << " coevent=\"1\"";
	else
		OUT << " coevent=\"0\"";

	//if(Atree->isExtant(event))
	//{
		if(Atree->hasNodeName(nodeid))
			OUT << " adjName=\"" << Atree->getNodeName(nodeid) << "\"" ;
	//}
	//else if(Atree->hasNodeProperty(nodeid,nodeid1)) // to be able to retieve them in the adjacency list
	//{
	//	pair <int,int> pid = Atree->getNodeNodeIds(nodeid);
	//	OUT << " adjName=\"" << pid.first << "-" << pid.second << "\"" ;
	//}



	if((!Atree->isBout(event))&&(!Atree->isLivingAdjBreak(nodeid,event))) //add the species information but not in the two cases where it is absent or useless
	{
		OUT << " ";
		OUT << SpeciesString;
		OUT << "\"";
		OUT << Atree->getNodeSpecies(nodeid);
		OUT << "\"";
	}

	OUT << "></";

	OUT << eventString;
	OUT << ">\n";

	return;


}

/*
Writes the adj node along with its children (recursively)

Takes:
	- OUT (ofstream) : an opened stream
	- Atree (AdjTree *) : a pointer to an adjacency tree
	- nodeid (int): id of a node
	- indent_level (int): number of double-spaces at the beginning of a line
*/
void DeCoOutputManager::WritePhyloXMLAdjTreeAux(ofstream& OUT, AdjTree * Atree,int nodeid, int indent_level)
{
	beginLine(OUT,indent_level);
	OUT << "<clade>\n";
	indent_level++;

	//writing the different informations:

	//1. name:
	beginLine(OUT,indent_level);
	OUT << "<name>";
	if(Atree->hasNodeName(nodeid))
	{
		OUT << Atree->getNodeName(nodeid);

	}
	else if(Atree->hasNodeProperty(nodeid,nodeid1))
	{
		pair <int,int> pid = Atree->getNodeNodeIds(nodeid);
		OUT << pid.first << "-" << pid.second ;

	}
	else
	{
		OUT << "NONAME";

	}

	OUT << "</name>\n";

	//2. the reconciliation event(s)
	beginLine(OUT,indent_level);
	OUT << "<eventsAdj>\n";
	indent_level++;

	int currentId = nodeid;
	int nextId = currentId;
	bool donext = true;

	while( donext )
	{
		currentId = nextId; //iteration

		vector <int> children = Atree->getSonsId(currentId);

		if(children.size() != 1)
			donext = false;
		else
			nextId = children[0];

		//we write the event of the current node
		WritePhyloXMLAdjEvent(OUT,  Atree, currentId, indent_level);
	}

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</eventsAdj>\n";


	//3. Sons, if any
	vector <int> children = Atree->getSonsId(currentId);
		
	for(unsigned i = 0 ; i < children.size(); i++)
	{
		WritePhyloXMLAdjTreeAux( OUT,  Atree, children[i], indent_level);
	}

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</clade>\n";
}

/*
Writes the sp node along with its children (recursively)

Takes:
	- OUT (ofstream) : an opened stream
	- Stree (MySpeciesTree *) : a pointer to a species tree
	- nodeid (int): id of a node
	- indent_level (int): number of double-spaces at the beginning of a line
*/
void DeCoOutputManager::WritePhyloXMLSpeTreeAux(ofstream& OUT, MySpeciesTree * Stree, int nodeid, int indent_level)
{

	MySpeciesNode * N = Stree->getNode(nodeid);
	//cout << N->getInfos().postOrder << " " << N->getInfos().timeSlice << " " << N->getInfos().realPostOrder << " " << nodeid << endl;

	beginLine(OUT,indent_level);
	OUT << "<clade>\n";
	indent_level++;

	double distance = 0;

	int currentId = nodeid;
	vector <int> children = Stree->getSonsId(currentId);

	

	while(children.size() ==1)
	{
		currentId = children[0];
		//cout  << "->" << currentId;
		children = Stree->getSonsId(currentId);
		if(Stree->hasDistanceToFather(currentId))
			distance += Stree->getDistanceToFather(currentId);
	}

	

	nodeid = currentId;
	//writing the different informations:
	/*cout << nodeid <<" RPO "<< N->getInfos().realPostOrder;
	if(Stree->isl)
	{
		cout << "name " << Stree->getNodeName(nodeid);
	}
	cout << endl;*/

	//1. name:
	beginLine(OUT,indent_level);
	OUT << "<name>";
	if(Stree->isLeaf(nodeid)) //internal nodes actually have names that do not correspond to their RPO but the TS of the speciation
	{
		OUT << Stree->getNodeName(nodeid);
	}
	else
	{
		OUT << N->getInfos().realPostOrder;
	}

	OUT << "</name>\n";



	//2. distance to father
	if(distance != 0)
	{
		beginLine(OUT,indent_level);
		OUT << "<branch_length>" << distance << "</branch_length>\n" ;
	}

	if(Stree->isSubdivided())
	{
		beginLine(OUT,indent_level);
		OUT << "<timeSlice>" << N->getInfos().timeSlice << "</timeSlice>##other\n" ;
	}

	//3. id
	beginLine(OUT,indent_level);
	OUT << "<node_id>" << N->getInfos().realPostOrder << "</node_id>\n";

	for(unsigned i = 0 ; i < children.size(); i++)
	{
		WritePhyloXMLSpeTreeAux( OUT,  Stree, children[i], indent_level);
	}	

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</clade>\n";
	

}


/*
Writes the adjacencies (i.e. co-speciation) in the adjacency tree
in the format : SP F1|G1 F2|G2

Takes:
	- OUT (ofstream) : an opened stream
	- ATree (AdjTree *) : a pointer to an adjacency tree
	- int sens1 : indicates which extremity of the first gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- int sens2 : indicates which extremity of the second gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- ReconciledTree * Rtree1 [default: NULL] : recocnilied tree that will be used to get the leaves names
	- ReconciledTree * Rtree2 [default: NULL] : recocnilied tree that will be used to get the leaves names
*/
void DeCoOutputManager::WriteAdjTreeAdjacencies(ofstream& OUT, AdjTree * ATree, int sens1, int sens2, ReconciledTree * Rtree1, ReconciledTree * Rtree2)
{


	vector <int> IDS = ATree->getNodesId();

	int f1 = ATree->getGfamily1();
	int f2 = ATree->getGfamily2();

	for(unsigned i = 0; i  <  IDS.size(); i++)
	{
		int evt = ATree->getNodeEvent(IDS[i]);
		if(ATree->isSpeciation(evt) || ATree->isExtant(evt))
		{
			if(! ATree->hasNodeProperty(IDS[i],nodeid1) ) //no nid property -> ignore
				continue;

			OUT << ATree->getNodeSpecies(IDS[i]) << " ";
			pair <int,int> gids = ATree->getNodeNodeIds(IDS[i]);

			if( ( ATree->isExtant(evt) ) && (Rtree1 != NULL) )
				OUT << Rtree1->getNodeName( gids.first );
			else
				OUT << f1 << "|" << gids.first;




			OUT << " ";

			if( ( ATree->isExtant(evt) ) && (Rtree2 != NULL) )
				OUT << Rtree2->getNodeName( gids.second );
			else
				OUT << f2 << "|" << gids.second;


			if(sens1 == -1)
				OUT << " -";
			else if(sens1 == 1)
				OUT << " +";
			else if(sens2 != 0)
				OUT << " *";

			if(sens2 == -1)
				OUT << " +";
			else if(sens2 == 1)
				OUT << " -";
			else if(sens1 != 0)
				OUT << " *";


			OUT << endl;
		}
	}

}


/*
Writes the adjacencies (i.e. co-speciation) in the adjacency forest
in the format : SP F1|G1 F2|G2

Takes:
	- OUT (ofstream) : an opened stream
	- CoEvent coevent : co event to write
*/
void DeCoOutputManager::WriteCoEvent(ofstream& OUT, CoEvent coevent)
{

	//1. 1st line: type of event, species , UTS - LTS , size of event (nb of nodes)
	//getters
	int evt = coevent.getEvent() ;
	string evtstring = "unknown";

	if( (evt == 0) || (evt == 1) || (evt == 6) ) // don't write leaf or sp adjs
		return;
	else if(evt == 2)
		evtstring = "loss";
	else if(evt == 3)
		evtstring = "duplication";
	else if(evt == 4)
		evtstring = "branchingOut";
	else if(evt == 5)
		evtstring = "transferBack";
	else if(evt == 7)
		evtstring = "bifurcationOut";


	OUT << "event: " << evtstring << " ";
	OUT << "species: " << coevent.getSpecies() << " ";
	if(coevent.getUpperTS() != -1)
	{
		OUT << "upperTS: " << coevent.getUpperTS() << " ";
		OUT << "LowerTS: " << coevent.getLowerTS() << " ";
	}
	OUT << "genes: " << coevent.getNumberOfGene() << " ";
	OUT << "connex components: " << coevent.getNbConnexComponents() << endl;


	for(unsigned i = 0 ; i <  coevent.getNumberOfAdj(); i++)
	{

		int id1 = coevent.getAdj(i).n1;
		int id2 = coevent.getAdj(i).n2;

		OUT << coevent.getGene(id1).Gfam << "|" << coevent.getGene(id1).NodeId;
		OUT << " ";
		OUT << coevent.getGene(id2).Gfam << "|" << coevent.getGene(id2).NodeId;

		OUT << endl;
	}

	OUT << endl;

}

/*
Writes the whole reconciled tree in the of stream

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- int index (default: -1): index of the Gfam to write (or -1 to omit it)
*/
void DeCoOutputManager::WritePhyloXMLRecTree(ofstream& OUT, ReconciledTree * Rtree, int index)
{
	//OUT << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	//OUT <<"<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\" xsi:schemaLocation=\"../recxml.xsd\">\n";

	int indent_level = 1;
	beginLine(OUT,indent_level);
	OUT << "<recGeneTree>\n";

	indent_level++;

	beginLine(OUT,indent_level);
	OUT << "<phylogeny rooted=\"true\">\n";

	indent_level++;

	if( index != -1)
	{
		beginLine(OUT,indent_level);
		OUT << "<id>"<< index <<"</id>"<<endl;
	}

	int rootid = Rtree->getRootId();
	WritePhyloXMLRecTreeAux( OUT,  Rtree, rootid, indent_level);

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</phylogeny>\n";

	indent_level--;
	beginLine(OUT,indent_level);

	OUT << "</recGeneTree>\n";

}

/*
Writes an adjacency tree in an ofstream

Takes:
	- OUT (ofstream) : an opened stream
	- Atree (AdjTree *) : a pointer to an adjacency tree
*/
void DeCoOutputManager::WritePhyloXMLAdjTree(ofstream& OUT, AdjTree * Atree)
{
	//OUT << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	//OUT <<"<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\" xsi:schemaLocation=\"../adjxml.xsd\">\n";

	int indent_level = 3;
	beginLine(OUT,indent_level);
	OUT << "<phylogeny rooted=\"true\" gainatroot=";
	if(Atree->isRootGain())	//telling wether there is a gain at the root or not
		OUT << "\"true\"";
	else
		OUT << "\"false\"";

	OUT << ">\n";

	indent_level++;
	int rootid = Atree->getRootId();
	WritePhyloXMLAdjTreeAux( OUT,  Atree, rootid, indent_level);

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</phylogeny>\n";
	//OUT << "</phyloxml>\n";


}


/*
Writes the whole species tree in the of stream

Takes:
	- OUT (ofstream) : an opened stream
	- Stree (MySpeciesTree *) : a pointer to a species tree

*/
void DeCoOutputManager::WritePhyloXMLSpeTree(ofstream& OUT, MySpeciesTree * Stree)
{
	OUT << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	OUT <<"<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\" xsi:schemaLocation=\"../recxml.xsd\">\n";

	int indent_level = 1;
	beginLine(OUT,indent_level);
	OUT << "<phylogeny rooted=\"true\">\n";



	indent_level++;
	int rootid = Stree->getRootId();
	WritePhyloXMLSpeTreeAux( OUT,  Stree, rootid, indent_level);

	indent_level--;
	beginLine(OUT,indent_level);
	OUT << "</phylogeny>\n";
	OUT << "</phyloxml>\n";
}

/*
Writes an adjacency tree in an ofstream

Takes:
	- OUT (ofstream) : an opened stream
	- Aforest (vector <AdjTree *> * ) : a pointer to vector of adjacency tree pointers
*/
void DeCoOutputManager::WritePhyloXMLAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, int indent)
{
	//OUT << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	//OUT <<"<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\" xsi:schemaLocation=\"../adjxml.xsd\">\n";


	for(unsigned i = 0; i < Aforest->size(); i++)
	{
		AdjTree * Atree = Aforest->at(i);

		int indent_level = indent;
		beginLine(OUT,indent_level);
		OUT << "<phylogeny rooted=\"true\" gainatroot=";
		if(Atree->isRootGain())	//telling wether there is a gain at the root or not
			OUT << "\"true\"";
		else
			OUT << "\"false\"";
	
		OUT << ">\n";
	
		indent_level++;
		int rootid = Atree->getRootId();
		WritePhyloXMLAdjTreeAux( OUT,  Atree, rootid, indent_level);
	
		indent_level--;
		beginLine(OUT,indent_level);
		OUT << "</phylogeny>\n";
	}

	//OUT << "</phyloxml>\n";
}



/*
Writes the whole reconciled tree in the of stream

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- hideLosses (bool) (default: false) : if true, losses and the branches leading to them wil be removed from the newick string
	- int index (default: -1): index of the Gfam to write (or -1 to omit it)
*/
void DeCoOutputManager::WriteNewickRecTree(ofstream& OUT, ReconciledTree * Rtree, bool hideLosses, int index )
{
	if(index != -1)
		OUT << ">family " << index << endl; 

	OUT << Rtree->NewickString(hideLosses) << ";" << endl;
}

/*
Writes an adjacency tree in an ofstream

Takes:
	- OUT (ofstream) : an opened stream
	- Aforest (vector <AdjTree *> * ) : a pointer to vector of adjacency tree pointers
	- hideLosses (bool) (default: false) : if true, losses and the branches leading to them wil be removed from the newick string
*/
void DeCoOutputManager::WriteNewickAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool hideLosses)
{
	for(unsigned i = 0; i < Aforest->size(); i++)
	{
		OUT << Aforest->at(i)->NewickString(hideLosses) << ";" << endl;
	}
}

/*
Writes the whole reconciled tree in the of stream

Takes:
	- OUT (ofstream) : an opened stream
	- Rtree (ReconciledTree *) : a pointer to a reconciled tree
	- newick (bool) : if true, the tree will be written in newick format. if false it will be written in a phyloXML like format
	- hideLosses (bool) (default: false) : if true, losses and the branches leading to them wil be removed from the newick string
	- int index (default: -1): index of the Gfam to write (or -1 to omit it)
*/
void DeCoOutputManager::WriteRecTree(ofstream& OUT, ReconciledTree * Rtree, bool newick, bool hideLosses, int index)
{
	if(newick)
		WriteNewickRecTree(OUT, Rtree, hideLosses, index);
	else
		WritePhyloXMLRecTree(OUT,Rtree, index );
}


/*
Writes an adjacency tree in an ofstream

Takes:
	- OUT (ofstream) : an opened stream
	- Aforest (vector <AdjTree *> * ) : a pointer to vector of adjacency tree pointers
	- newick (bool) : if true, the trees will be written in newick format. if false they will be written in a phyloXML like format
	- hideLosses (bool) (default: false) : if true, losses and the branches leading to them wil be removed from the newick string
*/
void DeCoOutputManager::WriteAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool newick, bool hideLosses, int indent)
{
	if(newick)
		WriteNewickAdjForest(OUT, Aforest, hideLosses);
	else
		WritePhyloXMLAdjForest(OUT,Aforest, indent);

}

/*
Writes the adjacencies (i.e. co-speciation) in the adjacency forest
in the format : SP F1|G1 F2|G2

Takes:
	- OUT (ofstream) : an opened stream
	- Aforest (vector <AdjTree *> * ) : a pointer to vector of adjacency tree pointers
	- int sens1 : indicates which extremity of the first gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- int sens2 : indicates which extremity of the second gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- ReconciledTree * Rtree1 [default: NULL] : recocnilied tree that will be used to get the leaves names
	- ReconciledTree * Rtree2 [default: NULL] : recocnilied tree that will be used to get the leaves names
*/
void DeCoOutputManager::WriteAdjForestAdjacencies(ofstream& OUT, vector<AdjTree *> * Aforest, int sens1 , int sens2, ReconciledTree * Rtree1 , ReconciledTree * Rtree2 )
{
	for(unsigned i = 0 ; i <  Aforest->size() ; i++)
		WriteAdjTreeAdjacencies(OUT, Aforest->at(i), sens1, sens2, Rtree1, Rtree2);
}

/*
Writes the adjacencies (i.e. co-speciation) in the adjacency forest
in the format : SP F1|G1 F2|G2

Takes:
	- OUT (ofstream) : an opened stream
	- vector <CoEvent> * CoEventSet : a set of coevents
*/
void DeCoOutputManager::WriteCoEventSet(ofstream& OUT, vector <CoEvent> * CoEventSet)
{
	for(unsigned i = 0 ; i <  CoEventSet->size(); i++)
	{
		//cout << "DeCoOutputManager::WriteCoEventSet" << i << endl;
		WriteCoEvent(OUT, CoEventSet->at(i));
	}
	//co-event are separated by empty lines
}