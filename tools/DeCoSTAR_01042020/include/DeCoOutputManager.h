#ifndef DECO_OUTPUT_MANAGER_H_
#define DECO_OUTPUT_MANAGER_H_

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

#include <fstream>

#include "ReconciledTree.h"
#include "AdjTree.h"
#include "CoEvent.h"

using namespace std;

class DeCoOutputManager
{
	protected:
		void beginLine(ofstream& OUT, int indent_level);
		
		void WritePhyloXMLRecEvent(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level, bool hasLoss);
		void WritePhyloXMLRecTreeAux(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level);

		void WritePhyloXMLAdjEvent(ofstream& OUT, AdjTree * Atree, int nodeid, int indent_level);
		void WritePhyloXMLAdjTreeAux(ofstream& OUT, AdjTree * Atree, int nodeid, int indent_level);

		void WritePhyloXMLSpeTreeAux(ofstream& OUT, MySpeciesTree * Stree, int nodeid, int indent_level);

		void WriteAdjTreeAdjacencies(ofstream& OUT, AdjTree * ATree, int sens1 = 0 , int sens2 = 0 , ReconciledTree * Rtree1 = NULL, ReconciledTree * Rtree2 = NULL);

		void WriteCoEvent(ofstream& OUT, CoEvent coevent);

	public:
		DeCoOutputManager(){}
		~DeCoOutputManager(){}

		void WritePhyloXMLRecTree(ofstream& OUT, ReconciledTree * Rtree, int index = -1);

		void WritePhyloXMLAdjTree(ofstream& OUT, AdjTree * Atree);
		void WritePhyloXMLSpeTree(ofstream& OUT, MySpeciesTree * Stree);
		void WritePhyloXMLAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, int indent = 3);

		void WriteNewickRecTree(ofstream& OUT, ReconciledTree * Rtree, bool hideLosses = false, int index = -1);
		void WriteNewickAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool hideLosses = false);

		void WriteRecTree(ofstream& OUT, ReconciledTree * Rtree, bool newick, bool hideLosses = false, int index = -1);
		void WriteAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool newick, bool hideLosses = false, int indent=3);

		void WriteAdjForestAdjacencies(ofstream& OUT, vector<AdjTree *> * Aforest, int sens1 = 0 , int sens2 = 0 , ReconciledTree * Rtree1 = NULL, ReconciledTree * Rtree2 = NULL);

		void WriteCoEventSet(ofstream& OUT, vector <CoEvent> * CoEventSet);

};

#endif