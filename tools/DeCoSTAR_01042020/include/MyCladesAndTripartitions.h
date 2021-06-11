#ifndef MYCLADESANDTRIPARTITIONS_H_
#define MYCLADESANDTRIPARTITIONS_H_
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


/**

@file
@author Wandrille Duchemin
@created 06-Nov-2015
@modified 16-Jun-2016


**/

#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "MyGeneTree.h"

#include "CladesAndTripartitions.h"



class MyCladesAndTripartitions : public CladesAndTripartitions
{
	protected:
		//pair <int,int> getNextOrRandomClade(CladeForest * CF , int idU , double probaFixed= 1.0);

		pair<int,int> getRandomClade(int idU);
		pair<int,int> getMaxClade(int idU);
		MyGeneNode *getRandomTreeAux(int &pOrd,int idU );
		MyGeneNode *getMaxTreeAux(int &pOrd,int idU );
		pair<double,vector<string> > getTreeLikelihoodAux(MyGeneNode *node, double default_proba);

		double getMaxLikelihoodAux(int idU);


		int deleteClades(vector<int> &cladesToDelete)
		{

    		while( cladesToDelete.size() > 0 ) 
    		{
        		// check if a leaf or root is marked for deletion
        		int rootIdU = mClades.getRootClade();
        		
        		BOOST_FOREACH( int idU, cladesToDelete ) 
        		{
            		if( idU == rootIdU || mClades.isLeaf( idU ) ) 
                		return -1; // error because we want to delete the root clade or a leaf clade
        		}

        		// delete clades
				cout << "deleting " << cladesToDelete.size() << " clades"  << endl;
        		vector<int> cladeIdMapping = mClades.PUBLICdeleteClades( cladesToDelete );( cladesToDelete );

        		// returns list of clades with no splits
        		cladesToDelete = redoTripartitions( cladesToDelete, cladeIdMapping );
        	}
    

        	return 0;
			
		};


	public:
		MyCladesAndTripartitions( char charSep, vector<MyGeneTree*> &geneTrees, bool verbose, bool &overflow, string &errStr, bool polytomy = false );

	    // rooted version
//	    CladesAndTripartitions( char charSep, MyTree &geneTree );

	    // Ale reader
		MyCladesAndTripartitions( char charSep, string aleFileName, string &errStr, bool verbose );


		MyGeneTree *getRandomTree();
		MyGeneTree *getMaxTree();


		double getTreeLikelihood(MyGeneNode *node, double default_proba);

		bool isCompatible(MyGeneNode * node);
		
		double getMaxLikelihood();

		void printNonCompatible( )
		{
		    vector< vector<int> > noncompatibleLists = mClades.PUBLICgetNoncompatible();
		    for(unsigned i = 0 ; i < noncompatibleLists.size(); i++)
		    {
		    	cout << i << ":";
		    	for(unsigned j = 0; j < noncompatibleLists[i].size(); j++)
		    		cout << " " << noncompatibleLists[i][j] ;
		    	cout << endl;
		    }
		
		};

		int RestrictToClade(vector <int> cladesToForce);

		int RestrictToClade(vector <vector <string> > cladesToForce);

		MyCladesAndTripartitions(MyCladesAndTripartitions * source); // new copy contructor

		//fc to clone
		int getmGeneTreeCount(){return mGeneTreeCount;}
		string getmNewickGeneTree(){return mNewickGeneTree;}
		double getmBranchLengths(int cladeNum){return mBranchLengths[cladeNum];}
};

#endif
