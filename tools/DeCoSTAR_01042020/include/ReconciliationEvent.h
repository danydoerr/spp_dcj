#ifndef RECONCILIATION_EVENT_H_
#define RECONCILIATION_EVENT_H_

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

This file contains a class for reconciliation events

Created the: 28-10-2015
by: Wandrille Duchemin

Last modified the: 20-01-2016
by: Wandrille Duchemin

*/

#include "CladesAndTripartitions.h"
#include "DTLGraph.h"
#include "MySpeciesTree.h"

#include "XMLUtils.h"

#include <string>
#include <map>
#include <Bpp/Exceptions.h>

using namespace std;
using namespace bpp;


class ReconciliationEvent
{
	protected:

		int idX; //postOrder id in the species tree
		int idXl; //postOrder id in the species tree of the left child of idX
		int idXr; //postOrder id in the species tree of the right child of idX
		//for postOrder id, -1 means that it is in a dead/unsampled lineage

		string evtName; //one of : S, SL, D, T, TL, TFD, TTD, TLTD, TLFD

		float support;

		int timeSlice;

		int getIdXForNullEvent(MySpeciesTree * t, int idx);

	public:

		int getidX();
		int getidXl();
		int getidXr();
		string getevtName();
		float getsupport();
		int gettimeSlice();


		void setidX( int idx);
		void setidXl( int idxl);
		void setidXr( int idxr);
		void setevtname( string evtname);
		void setsupport( float supp);
		void settimeSlice( int ts);

		ReconciliationEvent()
		{evtName=""; timeSlice = -1;}

		ReconciliationEvent(CladesAndTripartitions * CandT, DTLGraph *graph, MySpeciesTree * speciesTree, int recIndex, int evtIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation, bool VERBOSE);

		ReconciliationEvent(string recstr, map<string, int> speciesLeafNameToId);

		ReconciliationEvent(string XMLline);

		~ReconciliationEvent()
		{}


		bool hasTimeSlice();
		bool is_empty();
		void printMe();

};




#endif