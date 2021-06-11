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

This file contains a class contaning all equivalence classes between two families.

Created the: 26-02-2016
by: Wandrille Duchemin

Last modified the: 10-01-2018
by: Wandrille Duchemin

*/



#include "EquivalenceClassFamily.h"





EquivalenceClassFamily::EquivalenceClassFamily(int f1, int f2) : sens1(0), sens2(0)
{
	Gfam1 = f1;
	Gfam2 = f2;
	
	EclassList.push_back( MultiRootEquivalenceClass(f1,f2)); //EclassList.push_back(new EquivalenceClass(f1,f2)); // WMODIF
}


int EquivalenceClassFamily::getGfamily1() const
{
	return Gfam1;
}

int EquivalenceClassFamily::getGfamily2() const
{
	return Gfam2;
}

/*
Takes:
	- fam1 (int) : new gene family 1 index
*/
void EquivalenceClassFamily::setGfamily1( int fam1)
{
	Gfam1 = fam1;
}

/*
Takes:
	- fam2 (int) : new gene family 2 index
*/
void EquivalenceClassFamily::setGfamily2( int fam2)
{
	Gfam2 = fam2;
}



int EquivalenceClassFamily::getSens1() const
{
	return sens1;
}
int EquivalenceClassFamily::getSens2() const
{
	return sens2;
}

void EquivalenceClassFamily::setSens1( int s1)
{
	sens1 = s1;
	for( unsigned i = 0; i < EclassList.size(); i++)
		EclassList[i].setSens1(sens1);
}

void EquivalenceClassFamily::setSens2( int s2)
{
	sens2 = s2;
	for( unsigned i = 0; i < EclassList.size(); i++)
		EclassList[i].setSens2(sens2);
}


/*
Returns:
 (int): number of Equivalence class in the Equivalence Class family
*/
int EquivalenceClassFamily::getNbEqClasses() const
{
	return EclassList.size();
}

/*
Returns:
 (int): number of adjacencies accross all Equivalence classes in the Equivalence Class family
*/
int EquivalenceClassFamily::getNbAdj()
{
	int c =0;
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		c += EclassList[i].getNbAdj(); //WMODIF
	return c;
}

/*
Returns:
 (int): number of adjacency trees accross all Equivalence classes in the Equivalence Class family
*/
int EquivalenceClassFamily::getNbAdjTrees()
{
	int c =0;
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		c += EclassList[i].getNbAdjTrees();//WMODIF
	return c;
}


/*
Returns:
 (int): number of adjacency gain accross all Equivalence classes in the Equivalence Class family
*/
int EquivalenceClassFamily::getNbAdjGain()
{
	int c =0;
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		c += EclassList[i].getNbAdjGain();//WMODIF
	return c;
}

/*
Returns:
 (int): number of adjacency break accross all Equivalence classes in the Equivalence Class family
*/
int EquivalenceClassFamily::getNbAdjBreak()
{
	int c =0;
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		c += EclassList[i].getNbAdjBreak();//WMODIF
	return c;
}

/*
Returns true only if all Equivalence classes have their adjMatrix set
*/
bool EquivalenceClassFamily::areSetAdjMatrix()
{
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		if( !EclassList[i].isetAdjMatrix() )//WMODIF
			return false;

	return true;
}

/*
Returns true only if all Equivalence classes have their adjForest set
*/
bool EquivalenceClassFamily::areSetAdjForest()
{
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		if( !EclassList[i].isetAdjForest() )//WMODIF
			return false;

	return true;
}

//hasers and finders

/*
Takes:
	- name (string): a gene name in the gene family 1

Returns:
	(bool): true if an Equivalence class has this gene; false otherwise
*/
bool EquivalenceClassFamily::hasLname1(string name)
{
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		if( EclassList[i].hasLname1(name) )//WMODIF
			return true;

	return false;
}

/*
Takes:
	- name (string): a gene name in the gene family 2

Returns:
	(bool): true if an Equivalence class has this gene; false otherwise
*/
bool EquivalenceClassFamily::hasLname2(string name)
{
	for(unsigned i = 0 ; i < EclassList.size(); i++)
		if( EclassList[i].hasLname2(name) )//WMODIF
			return true;

	return false;
}



/*
Takes:
	- name (string): a gene name in the gene family 1

Returns:
	pair <int,int> : the first number if the Equivalence class id, the second is the gene index; both are set to -1 of the gene is not present in any Equivalence Class
*/
pair <int, int> EquivalenceClassFamily::findLname1(string name)
{
	pair <int,int> P;
	P.first = -1;
	P.second = -1;

	for(unsigned i = 0 ; i < EclassList.size(); i++)
	{
		int pos = EclassList[i].findLname1(name) ;//WMODIF
		
		if(pos != -1)
		{
			P.first = i;
			P.second = pos;
			break;
		}
	}		
	return P;	
}

/*
Takes:
	- name (string): a gene name in the gene family 2

Returns:
	pair <int,int> : the first number if the Equivalence class id, the second is the gene index; both are set to -1 of the gene is not present in any Equivalence Class
*/
pair <int, int> EquivalenceClassFamily::findLname2(string name)
{
	pair <int,int> P;
	P.first = -1;
	P.second = -1;

	for(unsigned i = 0 ; i < EclassList.size(); i++)
	{
		int pos = EclassList[i].findLname2(name) ;//WMODIF
		
		if(pos != -1)
		{
			P.first = i;
			P.second = pos;
			break;
		}
	}		
	return P;	
}




//adders with checks on gfamily names

/*
Takes:
 - f1 (int): Gfamily index of name1
 - f2 (int): Gfamily index of name2
 - s1 (int) : which end of the family 1 this is
 - s2 (int) : which end of the family 2 this is

Returns:
 (bool): true if fam1 and fam2 correspond to valid indexes sens correspond too

*/
bool EquivalenceClassFamily::isCompatible(int f1, int f2, int s1, int s2)
{
	
	bool invert = false;

	if(f1 != Gfam1)
	{
		if(f1 != Gfam2)
		{ // no fam correspond
			return false;
		}
		else if( (f2 == Gfam1) && (f1 == Gfam2) )
		{ // inversion between the two
			invert = true;
		}
		else
		{ // one correspond but not the other
			return false;
		}
	}
	else if( f2 != Gfam2)
		return false;



	if(invert)
	{ // invert is possible only if the two fam are different...
		return ( (s2 == sens1) && (s1 == sens2) );
	}

	if( f1 == f2 )
		return ( ( (s2 == sens1) && (s1 == sens2) ) || ( (s1 == sens1) && (s2 == sens2) ) );
	else
		return ( (s1 == sens1) && (s2 == sens2) ) ; 
}

/*
Adds a pair of Leaf names if fam1 and fam2 correspond to valid Gfamily names.
Arbitrarily add the adjacency to the first Equivalence Class

Takes:
 - name1 (string): name of a leaf
 - name2 (string): name of a leaf
 - fam1 (int): Gfamily index of name1
 - fam2 (int): Gfamily index of name2
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fam1 and fam2 correspond to valid indexes and the names could be added
*/
bool EquivalenceClassFamily::CheckAddAdj(string name1, string name2, int fam1, int fam2, int s1 , int s2 )
{
	if( getNbEqClasses() == 0)
		throw Exception("EquivalenceClassFamily::CheckAddAdj : tried to add an adjacency while there are no EquivalenceClass");

	return EclassList[0].CheckAddAdj( name1,  name2,  fam1,  fam2 ,s1 , s2);//WMODIF
}

/*
Adds a pair of Leaf names if fams correspond to valid Gfamily names 
Arbitrarily add the adjacency to the first Equivalence Class

Takes:
 - names (pair<string,string>): pair of leaf names
 - fams (pair<int,int>): Gfamily indexs of names
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClassFamily::CheckAddAdj(pair <string, string> names, pair <int, int> fams, int s1, int s2)
{
	return CheckAddAdj( names.first, names.second , fams.first, fams.second, s1,s2);
}

/*
Adds a list of pair of Leaf names if fams correspond to valid Gfamily names 
Arbitrarily add the adjacency to the first Equivalence Class

Takes:
 - nameList1 (vector<string>): list of leaf names
 - nameList2 (vector<string>): list of leaf names
 - fams (pair<int,int>): Gfamily indexs of names
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClassFamily::CheckAddAdjList(vector <string> nameList1, vector<string> nameList2, pair <int, int> fams, int s1, int s2)
{
	if( getNbEqClasses() == 0)
		throw Exception("EquivalenceClassFamily::CheckAddAdjList : tried to add adjacencies while there are no EquivalenceClass");

	return EclassList[0].CheckAddAdjList( nameList1,  nameList2,  fams, s1,s2);	//WMODIF
}


/*
Adds a list of pair of Leaf names if fams correspond to valid Gfamily names 
Arbitrarily adds the adjacency to the first Equivalence Class

Takes:
 - nameList (vector<pair<string,string> >): list of pairs of leaf names
 - fams (pair<int,int>): Gfamily indexs of names

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClassFamily::CheckAddAdjList(vector <pair <string,string> > nameList, pair <int, int> fams, int s1, int s2)
{
	if( getNbEqClasses() == 0)
		throw Exception("EquivalenceClassFamily::CheckAddAdjList : tried to add adjacencies while there are no EquivalenceClass");

	return EclassList[0].CheckAddAdjList( nameList,  fams, s1, s2);	//WMODIF
}



/*
Takes:
	- i (int): index of the desired Equivalence class

Returns:
	(EquivalenceClass *): pointer to the i-th Equivalence class
*/
const EquivalenceClass * EquivalenceClassFamily::getEClass(int i) const
{
	if( getNbEqClasses() <= i)
		throw Exception("EquivalenceClassFamily::getEClass : asked for non-existant EquivalenceClass.");

	EquivalenceClass * ECP;//WMODIF
	return &EclassList[i];

//	return ECP;
}

/*
Takes:
	- i (int): index of the desired Equivalence class

Returns:
	(vector <AdjTree * > *): pointer to the i-th Adjacency Forest
*/
vector <AdjTree * > * EquivalenceClassFamily::getAdjForest(int i)
{
	if( getNbEqClasses() <= i)
		throw Exception("EquivalenceClassFamily::getAdjForest : asked for non-existant EquivalenceClass.");

	return EclassList[i].getAdjForest();//WMODIF
}

/*
Replace the existsing adj forest at the given index. Presume that the given index exists

Takes
 - int i : index to replace at
 - vector <AdjTree * > * newAForest : forest to put at this position
*/
void EquivalenceClassFamily::setAdjForest(int i, vector <AdjTree * > * newAForest)
{
	if( getNbEqClasses() <= i)
		throw Exception("EquivalenceClassFamily::setAdjForest : asked for non-existant EquivalenceClass.");

	EclassList[i].setAdjForest(newAForest);//WMODIF
}

/*
Destroys the Equivalence classes and put all their adjacencies in a new one
*/
void EquivalenceClassFamily::reset()
{

	// creation of the new equivalence class
	MultiRootEquivalenceClass newEClass(Gfam1, Gfam2); // WMODIF

	for(unsigned i = 0; i <  EclassList.size(); i++)
		newEClass.merge( getEClass(i) ); // adding all ajacencies in newEClass


//	for ( vector< MultiRootEquivalenceClass  >::iterator it = EclassList.begin() ; it != EclassList.end(); ++it)
//   	{
//    	delete (it);
//   	}
	EclassList.clear();


	EclassList.push_back(newEClass);

}


/*
Takes:
 - Rtree1 (ReconciledTree *) : reconciled tree of the first gene family
 - Rtree2 (ReconciledTree *) : reconciled tree of the second gene family
 - useWholeClass (bool): Refines the Equivalence Class, considering all possible adjacencies between the two gene families and not only the observed ones.
 - forceLTrefining (bool): will refine in DeCoLT's way even if both trees have no transfer
 - verbose (bool)
*/
void EquivalenceClassFamily::refine(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool useWholeClass ,bool forceLTrefining, bool verbose)
{

	if(verbose)
	{
		cout << "refining equivalence class between gene families " << Gfam1 << " and " << Gfam2 << "." << endl; 
	}

	bool refine = true;

	if(refine)
	{
		vector < EquivalenceClass * > refinedEClasses;
	
		for(unsigned i = 0; i < EclassList.size(); i++)
		{
			
	
			
			vector<EquivalenceClass *> newEqClass;
			if(useWholeClass)
				newEqClass = EclassList[i].refineEqClassWhole(Rtree1, Rtree2, verbose);
			else
				newEqClass = EclassList[i].refineEqClass(Rtree1, Rtree2, forceLTrefining, verbose); 
	
			for(unsigned j =0; j < newEqClass.size(); j++)
				refinedEClasses.push_back(newEqClass[j]);
		}
		//cleaning
		EclassList.clear();
	
		for(unsigned i = 0; i < refinedEClasses.size(); i++)
			EclassList.push_back(MultiRootEquivalenceClass( refinedEClasses[i] ) );
		 //WMODIF
		for(unsigned i = 0; i < refinedEClasses.size(); i++)
			delete refinedEClasses[i];

	}

	if(verbose)
		cout << EclassList.size() << " classes after refinement." << endl;

}


/*
Creates the AdjMatrix of all the contained Equivalence classes.
Will replace the ancient ones if it exists

Takes:
 - map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : for artdeco
 - map<int, map<string,int> > speGeneAdjNb : for artdeco
 - map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb: used by artdeco -> degree of each gene extremities, used to modulate the costs in speciesC0C1.
 - Gcost (double): cost of a Gain
 - Bcost (double): cost of a Break
 - rtree1 (ReconciledTree *): reconciled tree for the first dimension
 - rtree2 (ReconciledTree *): reconciled tree for the second dimension
 - VERBOSE (bool)
 - boltzmann (bool) (default: false): wether to use boltzmann computation or not
 - bool LossAware  : whether to use loss aware  stuff or not
 - pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies : contains the list of free adjacencies!
 - temp (double) (default: 1) : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - absencePenalty (double) (default: -1): if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [default : 10 000] : base of the log that will be used to go from adjacency confidence score to parsimony costs

Returns:
	(vector < pair < pair<string, string> , double > > ) : vector of pairs associating an adjacency to a score
*/
vector < pair < pair<string, string> , double > > EquivalenceClassFamily::createAdjMatrix(
										map < string, map <string , double> > &adjacencyScores,
										map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb, 
										map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
										double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, 
										bool VERBOSE, bool boltzmann ,
										bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
										double temp, double absencePenalty, double adjScoreLogBase , bool interactionMode)
{
	vector < pair < pair<string, string> , double > > ScoreAssociation;
	for(unsigned i = 0; i < EclassList.size(); i++)
	{
		vector <double>  scores = EclassList[i].createAdjMatrix(adjacencyScores, speciesC0C1, speGeneAdjNb, speGeneExtremitiesAdjNb, Gcost, Bcost, rtree1, rtree2, VERBOSE,  boltzmann ,LossAware, FamiliesFreeAdjacencies,  temp,  absencePenalty , adjScoreLogBase, interactionMode);

		for(unsigned j = 0 ; j < scores.size(); j++)
			ScoreAssociation.push_back( pair < pair<string, string> , double >( EclassList[i].getAdj(j), scores[j] ) );
	}
	return ScoreAssociation;
}



/*
compute the adjMatrix of all the equivalence classes, if they are set
*/
void EquivalenceClassFamily::computeAdjMatrix()
{
	for(unsigned i = 0; i < EclassList.size(); i++)
		EclassList[i].computeAdjMatrix();//WMODIF
}


/*
Version of matrix computing that substract reconciliation event scores to co-events in the adjacency matrices

Takes:
 - WDupCost (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WLossCost (double):Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WHgtCOst (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
*/
void EquivalenceClassFamily::computeAdjMatrix(double WDupCost, double WLossCost, double WHgtCost)
{
	for(unsigned i = 0; i < EclassList.size(); i++)
		EclassList[i].computeAdjMatrix(WDupCost, WLossCost, WHgtCost);//WMODIF
}


/*
Returns:
	(bool): true if the adjacency matrices of all contained Equivalence class have been computed
*/
bool EquivalenceClassFamily::iscomputedAdjMatrix()
{

	for(unsigned i = 0; i < EclassList.size(); i++)
		if(! EclassList[i].iscomputedAdjMatrix() ) //at least one isn't set //WMODIF
			return false;
	return true; //all are set
}



/*
Backtracks the AdjMatrices and populates the AdjacencyTrees

Takes:
	- rtree1 (ReconciledTree *): reconciled tree for the first dimension
 	- rtree2 (ReconciledTree *): reconciled tree for the second dimension
	- AdjacencyTreesVector (vector < vector< AdjTree *> * > *): pointer to a vector of adjacency tree pointers vector pointer that will be updated as we build more adjacency trees
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score
 	- int &overflowed : number of EC that overflows
 	- alwaysGainAtTop (bool) [default: true]: there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 	- c1proba (double) [default = 0.5]: probability to choose c1 over c0 IF (and only if) they have the same score
 	- bool doNBBT [default = true] : whether or not most parsimonious solutions will be chosen using number of backtrack counts
*/
void EquivalenceClassFamily::backtrackAdjMatrix(ReconciledTree * rtree1, ReconciledTree * rtree2, vector < vector< AdjTree *> * >  * AdjacencyTreesVector, bool stochastic, int &overflowed, bool alwaysGainAtTop , double c1proba , bool doNBBT)
{
	//cout << "EquivalenceClassFamily::backtrackAdjMatrix init" << endl;
	for(unsigned i = 0; i < EclassList.size(); i++)
	{
		//cout << "EquivalenceClassFamily::backtrackAdjMatrix "<< i << endl;
		bool overflow = false;
		AdjacencyTreesVector->push_back( new vector< AdjTree *>);
		EclassList[i].backtrackAdjMatrix( rtree1, rtree2,  AdjacencyTreesVector->back(),  stochastic, overflow, alwaysGainAtTop , c1proba , doNBBT); //WMODIF
		if(overflow)
			overflowed++;
	}
}


/*
Backtracks the AdjMatrix of the Equivalences classes and populates their AdjForests

Takes:
	- rtree1 (ReconciledTree *): reconciled tree for the first dimension
 	- rtree2 (ReconciledTree *): reconciled tree for the second dimension
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score
 	- alwaysGainAtTop (bool) [default: true]: there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 	- c1proba (double) [default = 0.5]: probability to choose c1 over c0 IF (and only if) they have the same score
  	- bool doNBBT [default = true] : whether or not most parsimonious solutions will be chosen using number of backtrack counts	
*/
void EquivalenceClassFamily::backtrackAdjMatrixForSelf(ReconciledTree * rtree1, ReconciledTree * rtree2, bool stochastic, bool alwaysGainAtTop , double c1proba , bool doNBBT)
{

	for(unsigned i = 0; i < EclassList.size(); i++)
		EclassList[i].backtrackAdjMatrixForSelf( rtree1, rtree2, stochastic, alwaysGainAtTop , c1proba , doNBBT); //WMODIF
}
		

void EquivalenceClassFamily::printMe(bool verbose)
{
	cout << "EquivalenceClassFamily between gene family " << Gfam1 ;
	if(sens1 == 1)
		cout << " (stop) ";
	else if (sens1 == -1)
		cout << " (start) ";

	cout << " and " << Gfam2;
	if(sens2 == 1)
		cout << " (stop) ";
	else if (sens2 == -1)
		cout << " (start) ";

	cout << ". " << getNbAdj() << " adjacencies in " << getNbEqClasses() << " EquivalenceClasses." << endl;
	if(verbose)
		for(unsigned i = 0; i < EclassList.size(); i++)
			EclassList[i].printMe(verbose);//WMODIF
}


void EquivalenceClassFamily::clone( const EquivalenceClassFamily * ECF)
{
	Gfam1 = ECF->getGfamily1();
	Gfam2 = ECF->getGfamily2();

	sens1 = ECF->getSens1();
	sens2 = ECF->getSens2();

	scoreA = ECF -> getScoreA();
	
	AdjIndexMap = ECF -> getAdjIndexMap();
	AdjInScoreList = ECF -> getAdjInScoreList();
	AdjOutScoreList = ECF -> getAdjOutScoreList();
	AdjSpeList = ECF -> getAdjSpeList();
	
	TmpFileName = ECF -> getTmpFile();

	for(unsigned i = 0; i  < ECF->getNbEqClasses();i++)
	{
		EclassList.push_back(MultiRootEquivalenceClass(Gfam1,Gfam2));
		EclassList.back().clone( ECF->getEClass(i) );
	}
}

/*
Takes:
	- double threshold (default = 200)

Returns:
	(vector <int>) the number of score whose absolute log 10 are above the given threshold for each equivalence class in the family
			(in both c1 and c0) 
			(ignoring the cases with worst score (infinite in parcimony))
*/
vector <int> EquivalenceClassFamily::getNumberScoreWithAbsLog10Above(double threshold )
{
	vector < int > nbAbove;
	nbAbove.reserve( getNbEqClasses() );
	for(unsigned i = 0 ; i < getNbEqClasses() ; i++ )
	{
		nbAbove.push_back( EclassList[i].getNumberScoreWithAbsLog10Above( threshold ) );
	}
	return nbAbove;
}

///////////////// MODIF: Adelme //////////////
void EquivalenceClassFamily::setScoreA(vector < pair < pair<string, string> , double > > sAssociations)
{
	
	scoreA = sAssociations;
}

vector < pair < pair<string, string> , double > > EquivalenceClassFamily::getScoreA() const
{
	return scoreA;
}

void EquivalenceClassFamily::setAdjIndexMap(map <string, map <string, int> > AdjIM)
{
	AdjIndexMap = AdjIM;
}

void EquivalenceClassFamily::setAdjInScoreList(vector < double > AdjISL)
{
	AdjInScoreList = AdjISL;
}

void EquivalenceClassFamily::setAdjOutScoreList(vector < double > AdjOSL)
{
	AdjOutScoreList = AdjOSL;
}

void EquivalenceClassFamily::setAdjSpeList(vector < int > AdjSL)
{
	AdjSpeList = AdjSL;
}

void EquivalenceClassFamily::setSetAdjMatrixFamily(bool isSet, int i){
	EclassList[i].setSetAdjMatrix(isSet);
}

map <string, map <string, int> > EquivalenceClassFamily::getAdjIndexMap() const
{
	return AdjIndexMap;
}
vector < double > EquivalenceClassFamily::getAdjInScoreList() const 
{
	return AdjInScoreList;
}

vector < double > EquivalenceClassFamily::getAdjOutScoreList() const
{
	return AdjOutScoreList;
}

vector < int > EquivalenceClassFamily::getAdjSpeList() const 
{
	return AdjSpeList;
}

void EquivalenceClassFamily::setTmpFile(string OutDir, string OutPrefix, int index)
{
	TmpFileName = OutDir + OutPrefix + "." + boost::lexical_cast<std::string>(index) + ".tmp";
}

string EquivalenceClassFamily::getTmpFile() const
{
	return TmpFileName;
}

