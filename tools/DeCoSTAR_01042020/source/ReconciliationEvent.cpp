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

This file contains classs for reconciliation events

Created the: 28-10-2015
by: Wandrille Duchemin

Last modified the: 23-03-2016
by: Wandrille Duchemin

*/


#include "ReconciliationEvent.h"
#include <boost/tokenizer.hpp>
#include <Bpp/Exceptions.h>



/*
Takes:
 - t (MySpeciesTree *): pointer to a species tree
 - idx (int) : a node id in the 

Returns:
	(int): idX of the nearest son that is not a Null event
*/
int ReconciliationEvent::getIdXForNullEvent(MySpeciesTree * t, int idx)
{
	vector <int> sons = t->getSonsId(idx); //getting the 1st (and only) son of the null node

	int IDX = t->getRPO(idx); // getting the real species name
	int NIDX = t->getRPO(sons[0]); // getting the real species name

	//cout << "ReconciliationEvent::getIdXForNullEvent" << IDX << "->" << NIDX << endl;
	if(( IDX != NIDX) || (sons.size() != 1))//different idx or number of children != 1
	{
		if(t->isAlpha(sons[0]))
			return -1;

		return NIDX;

	}
		

	return getIdXForNullEvent(t,sons[0]);//recurse on that node
}

int ReconciliationEvent::getidX()
{
	return idX;
}

int ReconciliationEvent::getidXl()
{
	return idXl;
}

int ReconciliationEvent::getidXr()
{
	return idXr;
}

string ReconciliationEvent::getevtName()
{
	return evtName;
}

float ReconciliationEvent::getsupport()
{
	return support;
}

int ReconciliationEvent::gettimeSlice()
{
	return timeSlice;
}

void ReconciliationEvent::setidX( int idx)
{
	idX = idx;
}

void ReconciliationEvent::setidXl( int idxl)
{
	idXl = idxl;
}

void ReconciliationEvent::setidXr( int idxr)
{
	idXr = idxr;
}

void ReconciliationEvent::setevtname( string evtname)
{
	evtName = evtname;
}

void ReconciliationEvent::setsupport( float supp)
{
	support = supp;
}


void ReconciliationEvent::settimeSlice( int ts)
{
	timeSlice = ts;
}

/*
Constructor of the class that uses ecceTERA classes

Takes:
 - CandT (CladesAndTripartitions *): pointer to a CCP distribution instance
 - graph (DTLGraph *): a pointer to a graph containing reconciliations (nodes can either represent gene-species mapping or event)
 - speciesTree (MySpeciesTree *): pointer to the species tree
 - recIndex (int): index of the current clade in the reconciliation
 - evtIndex (int): index of the current event in the reconciliation of the clade
 - reconciliation (vector< vector<DTLGraph::MyGraph::Vertex> > *): pointer to the representation of the reconciliation as a vector of vector of nodes in the graph
 - VERBOSE (bool)
 */
ReconciliationEvent::ReconciliationEvent( CladesAndTripartitions * CandT , DTLGraph * graph, MySpeciesTree * speciesTree, int recIndex, int evtIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation, bool VERBOSE)
{
	//junk variables to feed getVertexIdentfiers
	int idU;//will also be used later
	double cost;

	graph->getVertexIdentfiers((*reconciliation)[recIndex][evtIndex] , idU, idX, cost); // getting idX
	//cout << "ReconciliationEvent ->" << idU << "," << idX << endl;
	string eventName = graph->getVertexName((*reconciliation)[recIndex][evtIndex + 1]);	//getting temporary event name


	int oldIdX = idX;

	timeSlice = speciesTree->getTimeSlice(idX); // getting the time slice
	if(speciesTree->isAlpha(idX))
		idX = -1;
	else //if(speciesTree->isSubdivided())//WDF
	{
		idX = speciesTree->getRPO(idX); // getting the real species name
	}

	for( size_t k=0; k<eventName.length(); k++ ) //getting the event name
	{
		if( eventName[k] == '_' )
			break;
		evtName += eventName[k];
	}


	if(VERBOSE)
		cout << "Building ReconciliationEvent with idU " << idU << " oldIdX " << oldIdX << " idX " << idX << " TS " << timeSlice << " event " << evtName << endl;
		//if(idX != -1)
		//{
		//	vector <int> sons = speciesTree->getSonsId(idX);
		//	cout << "sons of "<< idX << " : ";
		//	for(unsigned i = 0; i < sons.size(); i++)
		//		cout << sons[i] << " , ";
		//	cout << endl;
		//}

	//Now we want to set the species children ids (and also set species id to -1 if necessary)
	if( evtName[1] == 'L' )//there is a loss
	{
		//If there is a loss, this can't be the last event of the vector, so we are safe when we fetch the evtIndex + 2 node

		int idUnext, nextIdX;
		double costnext;
		graph->getVertexIdentfiers((*reconciliation)[recIndex][evtIndex + 2] , idUnext, nextIdX, costnext); // getting the species of the next reconciliation event. The lost id will be the other one


		if( evtName == "SL" ) //SPECIATION LOSS
		{//the difficulty is in determining which son has been lost
			vector <int> sons = speciesTree->getSonsId(oldIdX);
			//cout << idX << "(" << oldIdX << ") -> " << sons[0]  << "," << sons[1] << endl;
			int Sonlid = sons[0];
			if( speciesTree->isAlpha(Sonlid)) //alpha true means that the node is in a dead lineage
				idXl = -1;
			else
				idXl = speciesTree->getRPO(Sonlid);
/*				if(speciesTree->isSubdivided())//WDF
					idXl = speciesTree->getRPO(Sonlid);
				else
					idXl = Sonlid;*/

			int Sonrid = sons[1];
			if( speciesTree->isAlpha(Sonrid) ) 
				idXr = -1;
			else
				idXr = speciesTree->getRPO(Sonrid);
/*				if(speciesTree->isSubdivided())//WDF
					idXr = speciesTree->getRPO(Sonrid);
				else
					idXr = Sonrid;*/
			
			if( nextIdX == Sonlid)//->getInfos().postOrder ) 
			{
				//we want idXl to be the lost index. so we invert them here
				int tmp = idXl;
				idXl = idXr;
				idXr = tmp;
			}
		} else if (evtName == "TL" || evtName == "TLFD")
		{
			
			if(evtName == "TLFD")
				idX = -1; // we come from the dead
			idXl = idX; // the left index is the lost one, which is also the non-transfered one -> = idX

			idXr = speciesTree->getRPO(nextIdX);
/*			if(speciesTree->isSubdivided())//WDF
				idXr = speciesTree->getRPO(nextIdX);
			else
				idXr = nextIdX;*/
		} else if( evtName == "TLTD" ) 
		{
			idXl = idX;// the left index is the lost one, which is also the non-transfered one -> = idX
			idXr = -1; // we go to the dead
		} else 
		{
			cout << "Unknown event: " << evtName << endl;
			throw Exception("ReconciliationEvent::ReconciliationEvent : unknown event" );
		}
	} else if( evtName == "S" || evtName == "T" || evtName == "TFD" || evtName == "TTD" || evtName == "D" || evtName == "DD") //bifurcation events
	{

		if( evtName == "D" ) //duplication -> species stays the same
		{
			idXl = idX;
			idXr = idX;
		} else if( evtName == "DD" ) //duplication/speciation in the dead -> all species at -1
		{
			idX = -1;
			idXl = idX;
			idXr = idX;
		} else // S T TFD or TTD
		{

			//We need to know there are chidlren to this clade ...
			pair <int,int> idUchildren = CandT->getCladeSplit(idU,0); // getting the children in the CanT
	
			int idUnext;
			double costnext;//just here as useless containers
			graph->getVertexIdentfiers((*reconciliation)[idUchildren.first][0] , idUnext, idXl, costnext); // getting idXl
			graph->getVertexIdentfiers((*reconciliation)[idUchildren.second][0] , idUnext, idXr, costnext); // getting idXr
	
			//setting -1 in case of alpha value
			if( evtName == "TFD" )
					idX = -1;
				// transfer to/from dead, replace alpha (dead) by -1
			

			//cout << idXl << "," << idXr ;

			if(  speciesTree->isAlpha(idXl)) 
				idXl = -1;
			else //if(speciesTree->isSubdivided())//WDF
				idXl = speciesTree->getRPO(idXl);
			if( speciesTree->isAlpha(idXr) ) 
				idXr = -1;
			else //if(speciesTree->isSubdivided())//WDF
				idXr = speciesTree->getRPO(idXr);


			//cout << " -> "<< idXl << "," << idXr  << endl;
		}
	} else if( evtName == "C" || evtName == "O" )//leaf or null event
	{
		idXl = idX;//the event stays the same
		idXr = idX;
	} else 
	{
		cout << "Unknown event: " << evtName << endl;
		throw Exception("ReconciliationEvent::ReconciliationEvent unknown event" );
	}


	setsupport(0);

}


/*
Constructor of the class that uses a event string formatted according to TERA's format

Takes:
 - recstr (string): a string formatted according to TERA's format
 - speciesLeafNameToId (map<string, int>) : a map which associates each species leaf name to its realpostOrder in the species tree
 */
ReconciliationEvent::ReconciliationEvent(string recstr, map<string, int> speciesLeafNameToId)
{
	char sep1 = ',';
	char sep2 = '@';

	int sep2pos = recstr.find(sep2);

	string part1_substr = recstr.substr(0, sep2pos);

	try
	{
		setsupport( atof(recstr.substr(sep2pos+1).c_str()) );
	}
	catch(Exception e)
	{
		throw Exception("ReconciliationEvent::ReconciliationEvent; string " + recstr.substr(sep2pos+1) + " could not be converted into float\n");
	}
		



	//separating what is left according to commas
	boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> > tok( part1_substr, sep );
    
    boost::tokenizer<boost::char_separator<char> >::iterator iter=tok.begin();

    try
    {
    	if((*iter)[0] == '\'')//leaf name
    		setidX(speciesLeafNameToId[iter->substr(1,iter->length() - 2)] ); // extracting the leafname, retrieving the correct id and setting it
    	else
		    setidX(atoi(iter->c_str()));
	    ++iter;

	    setevtname(*iter);
	    ++iter;

    	if((*iter)[0] == '\'')//leaf name
    		setidXl(speciesLeafNameToId[iter->substr(1,iter->length() - 2)] ); // extracting the leafname, retrieving the correct id and setting it
    	else
		    setidXl(atoi(iter->c_str()));
		
	    ++iter;
    	if((*iter)[0] == '\'')//leaf name
    		setidXr(speciesLeafNameToId[iter->substr(1,iter->length() - 2)] ); // extracting the leafname, retrieving the correct id and setting it
    	else
		    setidXr(atoi(iter->c_str()));
	    
	}
	catch(Exception e)
	{
		throw Exception("ReconciliationEvent::ReconciliationEvent; tokenized string " + part1_substr + " could not be converted into a sery of type int,str,int,int on token " + *iter +"\n");
	}

	timeSlice = -1;//default time slice value

}

/*
Takes:
	- XMLline (string): a line from a recPhyloXML file

*/
ReconciliationEvent::ReconciliationEvent(string XMLline)
{
	//0. some default values:
	setsupport(0);
	settimeSlice(-1);


	map <string, string> * properties = new map <string,string>;
	string * value; 

	//1. reading the line
	string name = InterpretLineXML(XMLline, properties, value);


	int tmpIdX;

	char * pEnd; // for str -> int conversion 
	//reporting stuff
	
	for (map<string,string>::iterator it=properties->begin(); it!=properties->end(); ++it)
	{
		if(it->first.compare("timeSlice") == 0)
			settimeSlice( (int) strtol(it->second.c_str(), &pEnd, 10) );
		else if (it->first.compare("speciesLocation") == 0)
			tmpIdX =  (int) strtol(it->second.c_str(), &pEnd, 10) ;
		else if (it->first.compare("destinationSpecies") == 0)
			tmpIdX =  (int) strtol(it->second.c_str(), &pEnd, 10) ;
	}

	//3. finding which event is implied
	if(name.compare("leaf") == 0)
	{
		evtName = "C";
		setidX(tmpIdX);
		setidXr(tmpIdX);
		setidXl(tmpIdX);
		settimeSlice(0);
	}
	else if(name.compare("speciation") == 0)
	{
		evtName = "S";
		setidX(tmpIdX);
		setidXr(-2); // -2 means that we will have to ask later...
		setidXl(-2); // -2 means that we will have to ask later...
	}
	else if( (name.compare("branchingOut") == 0) || (name.compare("speciationOut") == 0))
	{
		evtName = "TTD";
		setidX(tmpIdX);
		setidXr(-2); // -2 means that we will have to ask later...
		setidXl(-2); // -2 means that we will have to ask later...	
	}
	else if(name.compare("bifurcationOut") == 0)
	{
		evtName = "DD";
		setidX(-1);
		setidXl(-1);
		setidXr(-1);
	}
	else if(name.compare("duplication") == 0)
	{
		evtName = "D";
		setidX(tmpIdX);
		setidXr(tmpIdX);
		setidXl(tmpIdX);
	}
	else if(name.compare("transferLoss") == 0)
	{
		evtName = "TLFD";
		setidX(-1);
		setidXr(tmpIdX);
		setidXl(-1);
	}
	else if( (name.compare("branchingOutLoss") == 0) || (name.compare("speciationOutLoss") == 0))
	{
		evtName = "TLTD";
		setidX(tmpIdX);
		setidXl(tmpIdX);
		setidXr(-1);
	}
	else if(name.compare("speciationLoss") == 0)
	{
		evtName = "SL";
		setidX(tmpIdX);
		setidXr(-2); // -2 means that we will have to ask later... --> this will be the hastest as we will have to ask to the species tree...
		setidXl(-2); // -2 means that we will have to ask later...	
	}


}


bool ReconciliationEvent::hasTimeSlice()
{
	return (timeSlice > -1);
}

//empty if the event string is not set
bool ReconciliationEvent::is_empty()
{
	return evtName.empty();
}

//useful function that print the object data to cout
void ReconciliationEvent::printMe()
{
	if(is_empty())
		cout << "empty event" << endl;
	else
		cout << idX << "," << evtName << "," << idXl << "," << idXr << "@" << support << endl;
}