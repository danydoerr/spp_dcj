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

This file contains a class for reconciliated clades

Created the: 28-10-2015
by: Wandrille Duchemin

Last modified the: 10-02-2016
by: Wandrille Duchemin

*/

#include "CladeReconciliation.h"

//utility function
bool is_digits(const string &str)
{
    return str.find_first_not_of("0123456789") == string::npos;
}




/*
Constructor of the class that uses ecceTERA classes

Takes:
 - CandT (CladesAndTripartitions *): pointer to a CCP distribution instance which only contains the ROOTED gene tree we're reconciling
 - graph (DTLGraph *): a pointer to a graph containing reconciliations (nodes can either represent gene-species mapping or event)
 - speciesTree (MySpeciesTree *): pointer to the species tree
 - recIndex (int): index of the current clade in the reconciliation
 - reconciliation (vector< vector<DTLGraph::MyGraph::Vertex> > *): pointer to the representation of the reconciliation as a vector of vector of nodes in the graph
 - previous (ReconciliationEvent) : the event that led to this clade
 - VERBOSE (bool)
 */
CladeReconciliation::CladeReconciliation(CladesAndTripartitions * CandT, DTLGraph *graph, MySpeciesTree * speciesTree, int recIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation, ReconciliationEvent previous, bool VERBOSE)
{

	//temporary variables that are needed to be fed to the getVertexIdentfiers function
	int idX;
	double cost;

	graph->getVertexIdentfiers(reconciliation->at(recIndex)[0] , idU, idX, cost); // getting IdU

	if(VERBOSE)
		cout << "Adding the clade " << idU << endl;

	//we now parse the different events
	for( size_t z=0; z<reconciliation->at(recIndex).size(); z+=2 ) 
	{

		ReconciliationEvent rev = ReconciliationEvent(CandT, graph, speciesTree, recIndex, z, reconciliation, VERBOSE);
		Events.push_back(rev); // adding the event
	}


	if(Events.back().getevtName() == "C")//The last event is a leaf
	{
		if(VERBOSE)
			cout << "The clade" << idU << " is a leaf";
		idUl = idU;
		idUr = idU;
	} else // this is not a leaf -> set idUl and idUr correctly
	{
		//We need to know there are chidlren to this clade ...
		pair <int,int> idUchildren = CandT->getCladeSplit(idU,0); // getting the children in the CandT

		idUl = idUchildren.first;
		idUr = idUchildren.second;
	}

	//setting the previous event
	previousEvent = previous;


	//setting name if the clade is a leaf
	bool isLeaf = CandT->mClades.isLeaf( idU );
	if( isLeaf )
		name = CandT->mClades.getLeafName(idU);
	else
		name = "";
}



/*
Constructor of the class that uses a line formatted according to TERA's format

Takes:
 - recline (string): a line formatted according to TERA's format
 - speciesLeafNameToId (map<string, int>) : a map which associates each species leaf name to its realpostOrder in the species tree
 - geneLeafNameToId (map<string, int>) : a map which associates each gene leaf name to its clade number

 */
CladeReconciliation::CladeReconciliation(string recline, map<string, int> speciesLeafNameToId, map<string, int> geneLeafNameToId)
{
	/*
	Structure of the line:

	idU:<previous event line>;<event line1>;<event line2>; ... :idUl,idUr\n

	idU is always a clade number, idUl and idUr might be leaf names (identified beause they are not numbers)

	If idU is the root, then it won't have a <previous event line>. As it is always a single split, it is the only line that will have no ";" in it 

	It is supposed that leaf names (species and gene) won't contains ",", ";", ":" or "@".

	*/



	char sep1 = ':';	
	char sep2 = ';';
	char sep3 = ',';

	//finding idU and idUl,idUr
	int sep1_pos = recline.find(sep1);
	int sep1_rpos = recline.rfind(sep1);
	int sep3_rpos = recline.rfind(sep3); // finding the comma between idUl and idUr

	bool terminal = false; //will be set to true if idU is a leaf

	if(sep1_pos == sep1_rpos) // only one : character -> idU is a leaf
		terminal = true;

	else if(sep3_rpos <= sep1_rpos)
		throw Exception("CladeReconciliation::CladeReconciliation; line " + recline+ " not compliant.\n");



	//getting idU
	string idUstr = recline.substr(0,sep1_pos);

	string eventStrings;

	if(!terminal)
	{

		//getting idUl and idUr
		int cutlength = recline.length()  - sep3_rpos -1;
		if(recline[recline.length() - 1] == '\n')
			cutlength--;//crop the '\n'
	
		string idUlstr =  recline.substr(sep1_rpos + 1, sep3_rpos - sep1_rpos - 1);
		string idUrstr =  recline.substr(sep3_rpos + 1, cutlength);


		//trying too assign idU
		try
		{
			idU = atoi(idUstr.c_str());
		}
		catch(Exception e)
		{
			throw Exception("CladeReconciliation::CladeReconciliation; idU " + idUstr + " could not be resolved as an int.\n");
		}



		try//trying to assign the isolated strings to int values
		{
			//determining if that is a leaf name or a id
  			if(!is_digits(idUlstr))//leaf name
	 			idUl = (geneLeafNameToId[idUlstr] ); // extracting the leafname, retrieving the correct id and setting it
	   		else
		    	idUl = atoi(idUlstr.c_str());
			//determining if that is a leaf name or a id
	   		if(!is_digits(idUrstr))//leaf name
	  			idUr = (geneLeafNameToId[idUrstr] ); // extracting the leafname, retrieving the correct id and setting it
  			else
		    	idUr = atoi(idUrstr.c_str());
	
		}
		catch(Exception e)
		{
			throw Exception("CladeReconciliation::CladeReconciliation; idU " + idUlstr + " or " + idUrstr + " could not be resolved as an int or identified as a gene leaf.\n");	
		}

		//getting the event strings
		eventStrings = recline.substr(sep1_pos + 1 , sep1_rpos - sep1_pos - 1);
	}
	else//terminal
	{
		try
		{
//			if(geneLeafNameToId.find(idUstr) == geneLeafNameToId.end() )
//				geneLeafNameToId[idUstr] = geneLeafNameToId[idUstr].size(); // if the id doesn't exists, Create it ...

			idU = (geneLeafNameToId[idUstr] ); // extracting the leafname, retrieving the correct id and setting it
			idUr = idU;
			idUl = idU;
		}
		catch(Exception e)
		{
			throw Exception("CladeReconciliation::CladeReconciliation; idU " + idUstr + " could not be resolved as an int or identified as a gene leaf.\n");	
		}

		eventStrings = recline.substr(sep1_pos + 1);
	}
		






		


	//now the event strings

	//finding the first ";" separator
	int sep2_pos = 	eventStrings.find(sep2);

	if(sep2_pos != string::npos)//what is before that first ; is the previous event
	{
		string previouseventstr = eventStrings.substr(0,sep2_pos);
		previousEvent = ReconciliationEvent(previouseventstr, speciesLeafNameToId);

		//update eventString
		eventStrings = eventStrings.substr(sep2_pos + 1 );
		sep2_pos = 	eventStrings.find(sep2);

		while(sep2_pos != string::npos)
		{
			Events.push_back( ReconciliationEvent(eventStrings.substr(0,sep2_pos), speciesLeafNameToId) );
			//update eventString
			eventStrings = eventStrings.substr(sep2_pos + 1 );
			sep2_pos = 	eventStrings.find(sep2);

		}
	}

	//last event (and only one if this is the root)
	Events.push_back( ReconciliationEvent(eventStrings, speciesLeafNameToId) );
	

}




void CladeReconciliation::printMe()
{
	cout << idU << ":";
	if(isLeaf())
		cout << name;
	else
		cout << idUl << "," << idUr;
	cout << endl;
	cout << "previous event : ";
	previousEvent.printMe();

	int i;
	for(i=0;i < Events.size();i++)
	{
		cout << " -> event " << i << " : ";
		Events[i].printMe();
	}

}

/*
Retruns true if the clade has no previous event (i.e. it is the root of a reconciled tree)
*/
bool CladeReconciliation::isRoot()
{
	return previousEvent.is_empty(); //root if the previous event is empty
}

/*
	Returns true if the child of idU is idU itself (i.e. this is a leaf clade)
*/
bool CladeReconciliation::isLeaf()
{
	return idU == idUl;
}
