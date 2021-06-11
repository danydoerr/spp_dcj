#ifndef ADJ_CLASS_FAMILY_H_
#define ADJ_CLASS_FAMILY_H_

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

#include <string>
#include <vector>

#include <Bpp/Exceptions.h>

#include "ReconciledTree.h"
#include "AdjMatrix.h"
#include "AdjTree.h"
#include "EquivalenceClass.h"

#include "MultiRootEquivalenceClass.h"

class EquivalenceClassFamily
{
protected:
	int Gfam1;
	int Gfam2;

	int sens1; // -1 if start ; 0 if indifferenciated ; 1 if stop
	int sens2; // -1 if start ; 0 if indifferenciated ; 1 if stop

	//vector <MultiRootEquivalenceClass * > EclassList; 
	vector <MultiRootEquivalenceClass > EclassList; 


	// LA modif
	map <string, map <string, int> > AdjIndexMap;
	vector < double > AdjInScoreList;
	vector < double > AdjOutScoreList;
	vector < int > AdjSpeList;
	
	vector < pair < pair<string, string> , double > > scoreA;
	
	string TmpFileName;
	//end LA modif

public:
	EquivalenceClassFamily() : sens1(0), sens2(0)
	{};

	EquivalenceClassFamily(int f1, int f2);

	~EquivalenceClassFamily()
	{
		//cout <<"destroy ECF"<<endl;
		EclassList.clear();
	}


	int getSens1() const;
	int getSens2() const;

	void setSens1( int s1);
	void setSens2( int s2);

	
	int getGfamily1() const;
	int getGfamily2() const;

	void setGfamily1( int fam1);
	void setGfamily2( int fam2);

	int getNbEqClasses() const;
	int getNbAdj(); // accross all eq classes
	int getNbAdjTrees();// accross all eq classes ; check if the adj trees have been computed

	int getNbAdjGain();// accross all eq classes ; check if the adj trees have been computed
	int getNbAdjBreak();// accross all eq classes ; check if the adj trees have been computed

	bool areSetAdjMatrix();
	bool areSetAdjForest();

	//hasers and finders
	bool hasLname1(string name);
	bool hasLname2(string name);
	pair <int, int> findLname1(string name);// returns -1 if the name is absent
	pair <int, int> findLname2(string name);// returns -1 if the name is absent

	//adders with checks on gfamily names
	bool isCompatible(int f1, int f2, int s1, int s2);

	bool CheckAddAdj(string name1, string name2, int fam1, int fam2, int s1 = 0, int s2 = 0); // adds a pair of Leaf names if fam1 and fam2 correspond to valid Gfamily names 
	bool CheckAddAdj(pair <string, string> names, pair <int, int> fams, int s1 = 0, int s2 = 0); // adds a pair of Leaf names if fams correspond to valid Gfamily names 
	bool CheckAddAdjList(vector <string> nameList1, vector<string> nameList2, pair <int, int> fams, int s1 = 0, int s2 = 0);
	bool CheckAddAdjList(vector <pair <string,string> > nameList, pair <int, int> fams, int s1 = 0, int s2 = 0); // adds several pair of leaf names


	const EquivalenceClass * getEClass(int i) const;
	vector <AdjTree * > * getAdjForest(int i);

	void setAdjForest(int i, vector <AdjTree * > * newAForest);

	void reset();

	void refine(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool useWholeClass, bool forceLTrefining, bool verbose);


	///ADJMatrixes
	vector < pair < pair<string, string> , double > > createAdjMatrix(
										map < string, map <string , double> > &adjacencyScores,
										map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb, 
										map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
										double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, 
										bool VERBOSE, bool boltzmann ,
										bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
										double temp = 1 , double absencePenalty = -1, double adjScoreLogBase =10000, bool interactionMode = false);

	void computeAdjMatrix();
	void computeAdjMatrix(double WDupCost, double WLossCost, double WHgtCost);

	bool iscomputedAdjMatrix();
	void backtrackAdjMatrix(ReconciledTree * rtree1, ReconciledTree * rtree2, vector < vector< AdjTree *> * > * AdjacencyTreesVector, bool stochastic, int &overflowed,  bool alwaysGainAtTop = true, double c1proba = 0.5, bool doNBBT = true); // creates the adj forest

	void backtrackAdjMatrixForSelf(ReconciledTree * rtree1, ReconciledTree * rtree2, bool stochastic, bool alwaysGainAtTop = true, double c1proba = 0.5, bool doNBBT = true); // creates the adj forest
		

	void printMe(bool verbose);


	void clone(const EquivalenceClassFamily * ECF);

	EquivalenceClassFamily& operator=(const EquivalenceClassFamily& rhs)
	{
		//cout << "equal operator - EquivalenceClassFamily"<<endl;

		EclassList.clear();		

		Gfam1 = rhs.getGfamily1();
		Gfam2 = rhs.getGfamily2();

		sens1 = rhs.getSens1();
		sens2 = rhs.getSens2();
	
		for(unsigned i = 0; i  < rhs.getNbEqClasses();i++)
		{
			EclassList.push_back(MultiRootEquivalenceClass(Gfam1,Gfam2));
			EclassList.back().clone( rhs.getEClass(i) );
		}

	

		return *this;
	}

	void dumpStuff()
	{
		for(unsigned i = 0; i  < getNbEqClasses();i++)
		{
			EclassList[0].dumpStuff();
		}

	};

	vector <int> getNumberScoreWithAbsLog10Above(double threshold = 200); 


	
	// LA modif
	vector < pair < pair<string, string> , double > > getScoreA() const;
	void setScoreA(vector < pair < pair<string, string> , double > > sAssociations);
	
	void setSetAdjMatrixFamily(bool isSet, int i);
	
	void setAdjIndexMap(map <string, map <string, int> > AdjIM);
	void setAdjInScoreList(vector < double > AdjISL);
	void setAdjOutScoreList(vector < double > AdjOSL);
	void setAdjSpeList(vector < int > AdjSL);

	map <string, map <string, int> > getAdjIndexMap() const;
	vector < double > getAdjInScoreList() const;
	vector < double > getAdjOutScoreList() const;
	vector < int > getAdjSpeList() const;
	
	void setTmpFile(string OutDir, string OutPrefix, int index);
	
	string getTmpFile() const;
	// end LA modif


};


#endif