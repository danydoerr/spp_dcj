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

#include <stdlib.h>
#include <time.h>
#include <boost/foreach.hpp>

#include "MyCladesAndTripartitions.h"


//Private functions


/*
Takes:
 - CladeForest * CF : pointer to a cladeForest instance
 - int idU : current idU
 - double probaFixed = 1.0 : probability not to choose the split (hypothetically) described in the CladeForest

Returns:
	(pair <int,int>):  a split of that clade chosen according to the CladeForest or randomly drawn (according to its split ratio)
*//*
pair <int,int> getNextOrRandomClade(CladeForest * CF , int idU , double probaFixed)
{

	pair <int,int> chosenSplit;

	double r = ((double) rand()/ RAND_MAX); //random results
	double ratio;//to store the various ratios


	int nbSplit = getSplitCount(idU);
	int counter;

	bool foundRandom = false;
	bool preferFixed = true;

	double rFixed = ((double) rand()/ RAND_MAX); //random results
	if( rFixed > probaFixed ) // we will chose at random
	{
		preferFixed = false;
	}

	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		if(preferFixed) // we chose the clade from the cladeForest if possible
		{
			pair <int,int> tempresult= getCladeSplit(idU,counter); // check if the clade exists in the CladeForest
	
			if( ( CF->hasClade(tempresult.first) ) || ( CF->hasClade(tempresult.second) ) ) // found
			{
				chosenSplit = tempresult;
				break;
			}
		}

		if(!foundRandom)
		{
			ratio = getSplitRatio(idU,counter);
			if( r <= ratio)
			{
				chosenSplit = getCladeSplit(idU,counter);//the chosen clade is found and returned
				foundRandom = true;
			}
			r = r - ratio;//not the chosen clade, we update r
		}

		if((foundRandom) && (!preferFixed))
			break; // we found the random choice and won't take the fixed choice -> no need to continue
	}

	return chosenSplit;
}
*/

/**
Takes:
	- int idU : id of a clade

Returns:
	pair<int,int> a split of that clade drawn randomly (according to its split ratio)
**/
pair<int,int> MyCladesAndTripartitions::getRandomClade(int idU)
{

	double r = ((double) rand()/ RAND_MAX); //random results
	double ratio;//to store the various ratios


	int nbSplit = getSplitCount(idU);
	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		ratio = getSplitRatio(idU,counter);
		if( r <= ratio)
			return getCladeSplit(idU,counter);//the chosen clade is found and returned
		r = r - ratio;//not the chosen clade, we update r
	}
	return getCladeSplit(idU,nbSplit - 1);//return the last split

}

/**
Takes:
	- int idU : id of a clade

Returns:
	pair<int,int> a split of that clade that has the maximum split ratio
**/
pair<int,int> MyCladesAndTripartitions::getMaxClade(int idU)
{

	//First, we find the list of all split of clade idU with the maximum ratio
	double maxRatio = 0;
	vector<int> maxIndex;//vector to keep the indexes of the split with the max ratio

	double ratio;

	int nbSplit = getSplitCount(idU);

	//cout << idU << ":" << nbSplit << endl;
	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		ratio = getSplitRatio(idU,counter);
		if(ratio > maxRatio)
		{
			maxRatio = ratio;
			maxIndex.clear();
			maxIndex.push_back(counter);
		}
		else if(ratio == maxRatio)
		{
			maxIndex.push_back(counter);
		}
	}


	
	if(maxIndex.size() == 1)//There is only one split with the maximum ratio: return it
	{
		return getCladeSplit(idU,maxIndex[0]);
	}

	//There is more than one split: return a random one
	int r = rand() % maxIndex.size();
	return getCladeSplit(idU,maxIndex[r]);

}


/**
Takes:
	- int &pOrd : current depth-first post-order number
	- int idU : current clade in the recursion

Returns:
	MyGeneNode * 
**/
MyGeneNode *MyCladesAndTripartitions::getRandomTreeAux(int &pOrd,int idU )
{
    MyGeneNode *node = new MyGeneNode();
    pair<int,int> cladeSplit = getRandomClade( idU );
    if( cladeSplit.first != -1 ) {
        node->addSon( getRandomTreeAux( pOrd, cladeSplit.first ) );
        node->addSon( getRandomTreeAux( pOrd, cladeSplit.second ) );
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
    } else 
        node->setName( mClades.getLeafName( idU ) );

    if( mBranchLengths[idU] != -1 )
        node->setDistanceToFather( (mBranchLengths[idU] + 1)/getCladeOccurrences(idU) );

    pOrd++;

    return node;
}


/**
Takes:
	- int &pOrd : current depth-first post-order number
	- int idU : current clade in the recursion

Returns:
	MyGeneNode * 
**/
MyGeneNode * MyCladesAndTripartitions::getMaxTreeAux(int &pOrd,int idU ) 
{
    MyGeneNode *node = new MyGeneNode();
    pair<int,int> cladeSplit = getMaxClade( idU);
    //cout << idU << "->" << cladeSplit.first << "," << cladeSplit.second << endl;
    if( cladeSplit.first != -1 ) {
        node->addSon( getMaxTreeAux( pOrd, cladeSplit.first ) );
        node->addSon( getMaxTreeAux( pOrd, cladeSplit.second ) );
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
    } else 
        node->setName( mClades.getLeafName( idU ) );

    if( mBranchLengths[idU] != -1 )
        node->setDistanceToFather( (mBranchLengths[idU] + 1)/getCladeOccurrences(idU) );

    pOrd++;

    return node;

}

/**
RECURSIVE

Takes:
	- MyGeneNode * node: a tree to compute the likelihood from
	- double default_proba: proba to give to a split that is absent from the MyCladesAndTripartitions counts
	
Returns:
	pair<double,vector<string>>:
		double: likelihood of the node (and underlying subtree) according to the CCPs distribution.
		vector<string>: list of all the leaves in the subtree rooted at node
**/
pair<double,vector<string> > MyCladesAndTripartitions::getTreeLikelihoodAux(MyGeneNode *node, double default_proba)
{

	std::pair<double, vector<string> > toreturn;

	toreturn.first = 1;//set base proba to 1.

	if (node->isLeaf())
	{
		toreturn.second.push_back(node->getName());
		return toreturn;
	}
	
	int nbson,i, cladeId;

	nbson = node->getNumberOfSons();

	if(nbson ==1)//only one son -> this is a intermediary node of some sort 
		return getTreeLikelihoodAux(node->getSon(0),default_proba);// return the value for its only son
	else if(nbson > 2)//more than 2 sons: error
		throw bpp::Exception("MyCladesAndTripartitions::getTreeLikelihoodAux: does not work when the tree is polytomic!");

	//We are now in the caser where there is more than one son

	std::pair<double, vector<string> > tempresult;//to put the result from each son

	std::pair<int,int> split;//split

	//first son
	tempresult = getTreeLikelihoodAux(node->getSon(0),default_proba);//getting the result for this son

	toreturn.first *= tempresult.first;//updating proba

	BOOST_FOREACH( string leafName, tempresult.second)
		toreturn.second.push_back(leafName);//updating leaf list

	split.first = mClades.getCladeIntFromLeafList(tempresult.second);

	tempresult.second.clear();//clearing 

	//second son
	tempresult = getTreeLikelihoodAux(node->getSon(1),default_proba);//getting the result for this son

	toreturn.first *= tempresult.first;//updating proba
	BOOST_FOREACH(string leafName, tempresult.second)
		toreturn.second.push_back(leafName);//updating leaf list

	split.second = mClades.getCladeIntFromLeafList(tempresult.second);

	tempresult.second.clear();//clearing 

	//now getting current cladeId
	cladeId = mClades.getCladeIntFromLeafList(toreturn.second);
	
	if(cladeId == -1 )
		toreturn.first *=default_proba;
	else if(split.first != -1 || split.second != -1)//if the son clades are known

		for(i=0 ; i < getSplitCount(cladeId) ; i++)
		{
			pair <int,int> testsplit = getCladeSplit(cladeId,i);
			if(split.first == testsplit.first || split.second == testsplit.first)//found the correct split
			{

				if(cladeId != 0 || mRooted)// special condition for the root clade: should not use split ratio if the distribution is unrooted
					toreturn.first *= getSplitRatio(cladeId,i);
				else if(cladeId == 0)
					toreturn.first *= ((double) getCladeOccurrences(split.first))/ ((double) mGeneTreeCount); //  if we are the "root" clade but the distribution is unrooted, we should multiply by the number of time we saw the child clade (clade probability)
				break;
			}
		}
		if(i == getSplitCount(cladeId))//the split was not found
			toreturn.first *=default_proba;

	return toreturn;
}


double MyCladesAndTripartitions::getMaxLikelihoodAux(int idU)
{

	if( mClades.isLeaf(idU) ) // stopping condition
		return 1;

	
	double maxRatio = 0;
	double ratio;

	int nbSplit = getSplitCount(idU);

	//cout << idU << ":" << nbSplit << endl;
	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		pair<int,int> cladeSplit = getCladeSplit(idU,counter);

		ratio = 1;
		if(idU != 0 || mRooted)// special condition for the root clade: should not use split ratio if the distribution is unrooted
			ratio *= getSplitRatio(idU,counter);
		else if(idU == 0)
			ratio *= ((double) getCladeOccurrences(cladeSplit.first))/ ((double) mGeneTreeCount); //  if we are the "root" clade but the distribution is unrooted, we should multiply by the number of time we saw the child clade (clade probability)


		ratio *= getMaxLikelihoodAux(cladeSplit.first) * getMaxLikelihoodAux(cladeSplit.second);

		if(ratio > maxRatio)
		{
			maxRatio = ratio;
		}
	}

	return maxRatio;

}


// Public functions



/**
 * Constructor - Compute all clades and tripartitions for the given 
 * gene trees. 
 */
MyCladesAndTripartitions::MyCladesAndTripartitions( char charSep, vector<MyGeneTree*> &geneTrees, bool verbose, bool &overflow, string &errStr, bool polytomy ) : CladesAndTripartitions( charSep, geneTrees,verbose,overflow, errStr,polytomy,NULL)
{srand (time(NULL));//init the random seed
}


/**
 * Constructor - Ale reader
 * 
 */
MyCladesAndTripartitions::MyCladesAndTripartitions( char charSep, string aleFileName, string &errStr, bool verbose ) : CladesAndTripartitions( charSep, aleFileName, errStr, verbose )
{srand (time(NULL));//init the random seed
}


/**
Returns a tree where bipartitions have been drawn at random according to their split ratios
**/
MyGeneTree * MyCladesAndTripartitions::getRandomTree()
{
 	if( mClades.getLeafCount() < 3 )
 		return getMaxTree(); 

    int pOrd = 0;
    MyGeneNode *node = getRandomTreeAux( pOrd, mClades.getRootClade() );
    return new MyGeneTree( *node );
}

/**
Returns a tree where bipartition have the maximum split ratios for the encountered clades
**/
MyGeneTree * MyCladesAndTripartitions::getMaxTree()
{
	int pOrd = 0;
	MyGeneNode *node;

	//special procedures when there is only 1 or 2 leaves in the tree.
	if(mClades.getLeafCount() == 1)
	{
		//cout << " clade count " << mClades.getCladeCount() << endl;
		//printMe();
	    node = new MyGeneNode();
	    node->setName( mClades.getLeafName( 1 ) );
	}
	else if(mClades.getLeafCount() == 2)
	{
	    node = new MyGeneNode();
    
		vector<int> Scl = mClades.getSortedClades();
		for(unsigned i = 0; i < Scl.size(); i++)
		{
			if(mClades.isLeaf(Scl[i]))
			{
				MyGeneNode *sonNode = new MyGeneNode();
				sonNode->setName( mClades.getLeafName( Scl[i] ) );
				node->addSon( sonNode);
				node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
			}		
		}
    }
    else
    {
		node = getMaxTreeAux( pOrd, mClades.getRootClade() );
    }
	

    return new MyGeneTree( *node );
}

/**
RECURSIVE

Takes:
	- MyGeneNode * node: a tree to compute the likelihood from
	- double default_proba: proba to give to a split that is absent from the MyCladesAndTripartitions counts

Returns:
	double: likelihood of the node (and underlying subtree) according to the CCPs distribution.

NB: requires the tree to be binary
**/
double MyCladesAndTripartitions::getTreeLikelihood(MyGeneNode *node, double default_proba)
{
	pair<double,vector<string> > result;
	result = getTreeLikelihoodAux(node,default_proba);//recursive function
	return result.first;
}

/*
Takes:
 - node (MyGeneNode *): node of a tree

Returns:
 (bool): true if the leafs of the subtree rooted at node all corresponds to existing clades in the MyCladesAndTripartitions instance
 */
bool MyCladesAndTripartitions::isCompatible(MyGeneNode * node)
{
	if(node->isLeaf())
	{
		vector <string> leafNames;
		leafNames.push_back(node->getName());

		try
		{
			mClades.getCladeIntFromLeafList(leafNames); // will throw an error if the leaf clade does not exists
		}
		catch(exception & e) // not very elegant.
		{
			return false;
		}
	}

	bool result= true;

	//recursion
	int nbson = node->getNumberOfSons();

	for(int i = 0; i < nbson; i++)
	{
		if(!isCompatible(node->getSon(i)))
		{
			result = false;
			break;
		}
	}
	return result;
}

double MyCladesAndTripartitions::getMaxLikelihood()
{
	return getMaxLikelihoodAux(mClades.getRootClade());
}



int MyCladesAndTripartitions::RestrictToClade(vector <int> cladesToForce)
{
	//1. from clade to force to clade to delete
	// the clades to delete are the clades incompatible with the ones we want to force
	vector< vector<int> > noncompatibleLists = mClades.PUBLICgetNoncompatible();	
	map <int,bool> cladeToDeleteMAP; // we first put the clades to delete in a map to avoid doubles
	//cout << "restricting to clades" ;
	for(unsigned i = 0 ; i < cladesToForce.size(); i++)
	{
		int idU = cladesToForce[i];
		//cout << " " << idU ;
		for(unsigned j = 0 ; j < noncompatibleLists[idU].size(); j++)
		{
			cladeToDeleteMAP[ noncompatibleLists[idU][j] ] = true;
		}
	}

	for(unsigned i = 0 ; i < cladesToForce.size(); i++)// remove the clades to restrict from the clade to delete list to if they appear in the map
	{
		int idU = cladesToForce[i];
		cladeToDeleteMAP[ idU ] = false;
	}

	//cout << endl;
	vector <int> cladeToDelete;
	for (map<int,bool>::iterator it=cladeToDeleteMAP.begin(); it!=cladeToDeleteMAP.end(); ++it) 
	{
		if(it->second) // that way we don't remove a clade to restrict
	  		cladeToDelete.push_back( it->first );
  	}
	//cout << "MyCladesAndTripartitions::RestrictToClade. Deleting "<< cladeToDelete.size() << " clades"<<endl;
	//for(unsigned i = 0 ; i < cladeToDelete.size(); i++)
	//	cout << cladeToDelete[i] << " ";
	//cout << endl;
  	//2. delete the clades to delete
	return deleteClades(cladeToDelete);
}


int MyCladesAndTripartitions::RestrictToClade(vector <vector <string> > cladesToForce)
{
	vector <int> cladesToForceIDU;

	for( unsigned i = 0 ; i < cladesToForce.size() ; i++)
	{
		int idU = mClades.getCladeIntFromLeafList( cladesToForce[i] );
		if(idU >= 0 )
			cladesToForceIDU.push_back( idU );
	
		//cout << idU <<" <- ";
		//for(unsigned j = 0 ; j < cladesToForce[i].size() ; j++)
		//	cout << cladesToForce[i][j] << " ";	
		//cout << endl;
	}
	return RestrictToClade(cladesToForceIDU);
}







/*
works if source is amalgamted trees
*/
MyCladesAndTripartitions::MyCladesAndTripartitions(MyCladesAndTripartitions * source)
{

    mRooted = false;   ///< amalgamated tree is not rooted

    mGeneTreeCount = source->getmGeneTreeCount();

    mNewickGeneTree = source->getmNewickGeneTree(); ///< save newick string for ALE output

	mClades = source->mClades;


	vector <int> sortedClades = mClades.getSortedClades();

	for(unsigned i = 0 ; i < sortedClades.size(); i++)
	{
		int cladeNum = i;// as this is how things are organized in the cladesandtripartition instance, we keep this structure

		int splitCount = source->getSplitCount( cladeNum );

		mCladeSplits.push_back( vector< pair<int,int> > () );        ///< All splits for each clade (by clade id)
		mSplitsRatio.push_back( vector <double> ()  );        ///< Occurence of each split as a vector, divided by clade occurrences.
		for(unsigned splitNum = 0 ; splitNum < splitCount; splitNum++)
		{
			mCladeSplits[cladeNum].push_back( source->getCladeSplit( cladeNum, splitNum ) ) ;
			mSplitsRatio[cladeNum].push_back( source->getSplitRatio( cladeNum, splitNum ) ) ;
		}

	    mCladeOccurrences.push_back( source->getCladeOccurrences( cladeNum ) );///< Occurence of each clade, indexed by clade number.

		mBranchLengths.push_back( source->getmBranchLengths( cladeNum ) ); ///< branch length for each clade

	}

}