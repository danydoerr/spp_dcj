#ifndef CLADE_RECONCILIATION_H_
#define CLADE_RECONCILIATION_H_

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
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <Bpp/Exceptions.h>
#include <fstream>


#include "CladesAndTripartitions.h"
#include "DTLGraph.h"
#include "MySpeciesTree.h"
#include "ReconciliationEvent.h"

#include "XMLUtils.h"

using namespace std;
using namespace bpp;

class CladeReconciliation
{
public:

	int idU;//clade id
	int idUl;//clade id of the left child of idU
	int idUr;//clade id of the right child of idU
	ReconciliationEvent previousEvent;//split event that lead to this clade. Left empty if this is the root clade
	vector<ReconciliationEvent> Events;
	string name;


	bool isRoot();
	bool isLeaf();

	CladeReconciliation(){Events.clear();}

	//constructor for ecceTERA classes
	CladeReconciliation(CladesAndTripartitions * CandT, DTLGraph *graph, MySpeciesTree * speciesTree, int recIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation,ReconciliationEvent previous, bool VERBOSE);


	CladeReconciliation(string recline, map<string, int> speciesLeafNameToId, map<string, int> geneLeafNameToId);


	~CladeReconciliation(){}

	void printMe();

};



#endif