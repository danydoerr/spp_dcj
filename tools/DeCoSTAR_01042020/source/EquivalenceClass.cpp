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

This file contains a class for adjacency classes

Created the: 24-11-2015
by: Wandrille Duchemin

Last modified the: 10-01-2018
by: Wandrille Duchemin

*/


#include "EquivalenceClass.h"


// protected adders

// adds a Leaf name to Lnames1
void EquivalenceClass::addLname1(string name)
{
	Lnames1.push_back(name);	
} 

// adds a Leaf name to Lnames2
void EquivalenceClass::addLname2(string name)
{
	Lnames2.push_back(name);
} 

// adds a pair of Leaf names
void EquivalenceClass::addAdj(string name1, string name2)
{
	addLname1(name1);
	addLname2(name2);
}


// adds a pair of Leaf names
void EquivalenceClass::addAdj(pair <string, string> names)
{
	addAdj(names.first,names.second);
}

// adds several pair of leaf names
void EquivalenceClass::addAdjList(vector <pair <string,string> > nameList)
{
	for(unsigned i; i < nameList.size(); i++)
	{
		addAdj(nameList[i].first,nameList[i].second);
	}
}

//protected removers

/*
Removes an element from Lnames1 at a given index. 

Takes:
 - index (int): index to remove

Returns:
	(bool) false if the index is not valid
*/
bool EquivalenceClass::removeLname1(int index)
{
	if(index < 0)
		return false;
	if(index >= Lnames1.size())
		return false;

	Lnames1.erase(Lnames1.begin() + index);
	return true;
} 


/*
Removes an element from Lnames1 at a given index. 

Takes:
 - name (string): name to remove

Returns:
	(bool) false if the name is not valid
*/
bool EquivalenceClass::removeLname1(string name)
{
	int pos = findLname1(name);
	return removeLname1(pos);
} 

/*
Removes an element from Lnames2 at a given index. 

Takes:
 - index (int): index to remove

Returns:
	(bool) false if the index is not valid
*/
bool EquivalenceClass::removeLname2(int index)
{
	if(index < 0)
		return false;
	if(index >= Lnames2.size())
		return false;

	Lnames2.erase(Lnames2.begin() + index);
	return true;
} 


/*
Removes an element from Lnames2 at a given index. 

Takes:
 - name (string): name to remove

Returns:
	(bool) false if the name is not valid
*/
bool EquivalenceClass::removeLname2(string name)
{
	int pos = findLname2(name);
	return removeLname2(pos);
} 





/*

Takes:
 - NodeId1 (int): id of a node in Rtree1
 - NodeId2 (int): id of a node in Rtree2
 - Rtree1 (ReconciledTree *): a Reconciled Tree
 - Rtree2 (ReconciledTree *): a Reconciled Tree
 - Ancestor1 (int): id of the ancestor in tree 1. (Root Id in the case of LT algorithm; LCA of Lnames1 otherwise)
 - Ancestor2 (int): id of the ancestor in tree 2. (Root Id in the case of LT algorithm; LCA of Lnames2 otherwise)

Returns:
	(pair <int, int>): pair of the highest compatible (same species (and TS)) which are ancestors of NodId1 and NodeId2 in, respectively, Rtree1 and Rtree2.
*/
pair <int, int> EquivalenceClass::getHighestCompatibleNodeIds(int NodeId1, int NodeId2, ReconciledTree * Rtree1, ReconciledTree * Rtree2, int Ancestor1, int Ancestor2)
{

	vector <Node *> Path1 = Rtree1->pathFromAncestorToSon(Ancestor1, NodeId1); // Get the path from root to son.
	vector <Node *> Path2 = Rtree2->pathFromAncestorToSon(Ancestor2, NodeId2); // Get the path from root to son.

	if(Gfamily1 == Gfamily2) // NodeId1 and NodeId2 are in the same tree
	{
		while(Path1.front()->getId() == Path2.front()->getId()) // we remove the nodes up to their LCA
		{
			Path1.erase(Path1.begin());
			Path2.erase(Path2.begin());
		}
	}

	pair <int,int> PtoReturn;
	PtoReturn.first = NodeId1;
	PtoReturn.second = NodeId2;//default pair is the one given in argument

	//double for to get the highest compatible pair of ancestor
	for(unsigned i = 0; i < Path1.size(); i++)
	{
		for(unsigned j = 0; j < Path2.size(); j++)
		{
			if(Rtree1->areCompatible(Path1[i],Path2[j])) 
			{// we found the pair we want
				PtoReturn.first = Path1[i]->getId();
				PtoReturn.second = Path2[j]->getId();	
				return PtoReturn;
			}
		}
	}
	return PtoReturn;
}


pair <int, int> EquivalenceClass::getHighestCompatibleNodeIdsStoppingReception(int NodeId1, int NodeId2, ReconciledTree * Rtree1, ReconciledTree * Rtree2, int Ancestor1, int Ancestor2)
{

	vector <Node *> Path1 = Rtree1->pathFromAncestorToSon(Ancestor1, NodeId1); // Get the path from root to son. (start and stop included)
	vector <Node *> Path2 = Rtree2->pathFromAncestorToSon(Ancestor2, NodeId2); // Get the path from root to son.

	if(Gfamily1 == Gfamily2) // NodeId1 and NodeId2 are in the same tree
	{
		while(Path1.front()->getId() == Path2.front()->getId()) // we remove the nodes up to their LCA
		{
			Path1.erase(Path1.begin());
			Path2.erase(Path2.begin());
		}
	}


	//we know that the two extremities of the adjacency are compatible.
	pair <int, int> HighestCompatible;
	HighestCompatible.first = NodeId1;
	HighestCompatible.second = NodeId2;

	int current1 = Path1.size() -1 ;
	int current2 = Path2.size() -1 ;

	bool keepGoing = true;


	//cout << "EquivalenceClass::getHighestCompatibleNodeIdsStoppingReception " <<  NodeId1 << "-" << NodeId2 << endl;


	//cout << "currents : " << current1 << "-" << current2 << endl;

	while(keepGoing)
	{
		bool firstRec = false;
		bool secondRec = false;

		bool areCompat = true; // otherwise we would have stopped

		while(areCompat) // going up alon Path1
		{
			HighestCompatible.first = Path1[current1]->getId();//updating the highest compatible
			//HighestCompatible.second = Path2[current2]->getId();

			current1--; // if compatible, we go up the tree == decrecrease current1 (remember that Path1 is from ancestor to son, so by going reverse we go son to ancestor)

			//two possible stopping conditions: we have arrived to the first node of Path1 (current1 < 0) OR the event of current1 is a reception
			if(current1 < 0)
				break;

			//cout << "moving current1 : " << current1 << ";" << Rtree1->getNodeEvent(Path1[current1]->getId()) << ";" << Rtree1->getNodeSpecies(Path1[current1]->getId()) << ";" << Rtree1->getNodeTimeSlice(Path1[current1]->getId()) << ";" << Rtree1->getNodeUpperBoundaryTS(Path1[current1]->getId())  << "-" << current2 << ";" <<  Rtree2->getNodeEvent(Path2[current2]->getId()) << ";" << Rtree2->getNodeSpecies(Path2[current2]->getId()) << ";" << Rtree2->getNodeTimeSlice(Path2[current2]->getId()) << ";" << Rtree2->getNodeUpperBoundaryTS(Path2[current2]->getId())<< endl;

			int evt = Rtree1->getNodeEvent(Path1[current1]->getId());
			if(Rtree1->isRec(evt))
			{
				firstRec = true;
				break;
			}

			if(Rtree1->getTimeSliceStatus() != 2) //not BTS
				areCompat = Rtree1->areCompatible(Path1[current1],Path2[current2]);
			else //BTS -> test on species only
				areCompat = Rtree1->haveSameSpecies(Path1[current1],Path2[current2]);

			//cout << "arecompat: "<< areCompat << endl;
		}
		//current1 is a far as it can for now. make current2 try to catch it back.

		//at this point current2 is compatible with current1 + 1
		areCompat = true;

		while(areCompat) // going up along Path2
		{
			HighestCompatible.second = Path2[current2]->getId();//updating the highest compatible
			

			current2--; // if compatible, we go up the tree == decrease current2 (remember that Path2 is from ancestor to son, so by going reverse we go son to ancestor)




			if(current2 < 0)
				break;

			//cout << "moving current2 : " << current1 << ";" << Rtree1->getNodeEvent(Path1[current1 +1 ]->getId()) << ";" << Rtree1->getNodeSpecies(Path1[current1 +1 ]->getId())<< ";" << Rtree1->getNodeTimeSlice(Path1[current1 +1 ]->getId())  << "-" << current2 << ";"<< Rtree2->getNodeEvent(Path2[current2]->getId()) << ";" << Rtree2->getNodeSpecies(Path2[current2]->getId())<< ";" << Rtree2->getNodeTimeSlice(Path2[current2]->getId()) << endl;

			int evt = Rtree2->getNodeEvent(Path2[current2]->getId());
			if(Rtree2->isRec(evt))
			{
				secondRec = true;
				break;
			}

			if(Rtree1->getTimeSliceStatus() != 2) //not BTS
				areCompat = Rtree1->areCompatible(Path1[current1 + 1 ],Path2[current2]);
			else //BTS -> test on species only
				areCompat = Rtree1->haveSameSpecies(Path1[current1 + 1 ],Path2[current2]);
		}

		//now, checking stopping conditions.
		//cout << "stopped moving : " << HighestCompatible.first << ";" << Rtree1->getNodeEvent(HighestCompatible.first) << ";" << Rtree1->getNodeSpecies(HighestCompatible.first) << "-" << HighestCompatible.second << ";" << Rtree2->getNodeEvent(HighestCompatible.second) << ";" << Rtree2->getNodeSpecies(HighestCompatible.second) << endl;

		if( (current1 < 0) || (current2 < 0) )
		{//we are the end of one of the path -> stop here
			keepGoing=false;
		}
		else if( firstRec != secondRec ) // one rec is true
		{
			keepGoing=false;
		}
		else
		{
			if(Rtree1->getTimeSliceStatus() != 2)
			{
				if( !Rtree1->areCompatible(Path1[current1],Path2[current2]) ) //there are two compatible recs and we are not at the end of either path -> check if the current are not compatible
				{ // if not compatible -> stop here
					keepGoing=false;	
				} // else keep going
			}
			else//BTS -> test on species only
			{
				if( !Rtree1->haveSameSpecies(Path1[current1],Path2[current2]) ) //there are two compatible recs and we are not at the end of either path -> check if the current are not compatible
				{ // if not compatible -> stop here
					keepGoing=false;	
				} // else keep going
			}
		}

		firstRec = false;
		secondRec = false;

	}

	//cout << "EquivalenceClass::getHighestCompatibleNodeIdsStoppingReception. found some highest. Ts:"<< Rtree1->getTimeSliceStatus() << endl;
	
	//cout << "list1 :" << current1 <<endl;
	//for( int i = Path1.size() -1 ; i > current1; i--)
	//{
	//	cout << i << ";" << Path1[i]->getId() << ";" << Rtree1->getNodeEvent(Path1[i]->getId()) << ";" << Rtree1->getNodeSpecies(Path1[i]->getId()) << ";" << Rtree1->getNodeTimeSlice((Path1[i]->getId())) << endl;
	//}
	//cout << "list2 :" << current2 <<endl;
	//for( int i = Path2.size() -1 ; i > current2; i--)
	//{
	//	cout << i << ";" << Path2[i]->getId() << ";" << Rtree2->getNodeEvent(Path2[i]->getId()) << ";" << Rtree2->getNodeSpecies(Path2[i]->getId()) << ";" << Rtree2->getNodeTimeSlice((Path2[i]->getId())) << endl;
	//}


	if( Rtree1->getTimeSliceStatus() == 2 )
	{ // for BTS, we now look at real compatibility

		//double for to get the highest compatible pair of ancestor
		for(unsigned i = current1; i < Path1.size(); i++)
		{
			for(unsigned j = current2; j < Path2.size(); j++)
			{
				if(Rtree1->areCompatible(Path1[i],Path2[j])) 
				{// we found the pair we want
					HighestCompatible.first = Path1[i]->getId();
					HighestCompatible.second = Path2[j]->getId();	
					return HighestCompatible;
				}
			}
		}
	}


	return HighestCompatible;

}



EquivalenceClass::EquivalenceClass(int fam1, int fam2)
{
	Gfamily1 = fam1; // id of the first gene family
	Gfamily2 = fam2; // id of the second gene family
	sens1 = 0;
	sens2 = 0;
	AdjForest = new vector< AdjTree *>;
	SetAdjMatrix = false;
	SetAdjForest = false;
}


/*Constructor from an equivalence class that just copies the list of adjacencies and its ancestors*/
EquivalenceClass::EquivalenceClass( EquivalenceClass *EC)
{

	AdjForest = new vector< AdjTree *>;
	
	SetAdjMatrix = false;
	SetAdjForest = false;
	//cout << " creaFEC SetAdjForest "<< SetAdjForest << endl;	

	Gfamily1 = EC->getGfamily1();
	Gfamily2 = EC->getGfamily2();

	sens1 = EC->getSens1();
	sens2 = EC->getSens2();

	ancestors = EC->getAncestors();


	vector <pair <string, string> > adjs = EC->getAdjs();
	for(unsigned i=0; i < adjs.size(); i++)
	{
		Lnames1.push_back(adjs[i].first);
		Lnames2.push_back(adjs[i].second);
	}

}



int EquivalenceClass::getSens1() const
{
	return sens1;
}
int EquivalenceClass::getSens2() const
{
	return sens2;
}

void EquivalenceClass::setSens1( int s1)
{
	sens1 = s1;
}

void EquivalenceClass::setSens2( int s2)
{
	sens2 = s2;
}


//getters
int EquivalenceClass::getGfamily1() const
{
	return Gfamily1;
}

int EquivalenceClass::getGfamily2() const
{
	return Gfamily2;
}


string EquivalenceClass::getLname1(int index) const
{
	if(index <0)
		throw Exception("EquivalenceClass::getLname1 : given a negative index");

	if(index >= Lnames1.size())
		throw Exception("EquivalenceClass::getLname1 : given an index above the maximum index");

	return Lnames1[index];
}


string EquivalenceClass::getLname2(int index) const
{
	if(index <0)
		throw Exception("EquivalenceClass::getLname1 : given a negative index");

	if(index >= Lnames2.size())
		throw Exception("EquivalenceClass::getLname1 : given an index above the maximum index");

	return Lnames2[index];	
}


pair <string,string> EquivalenceClass::getAdj(int index) const
{
	pair <string,string> adj;

	adj.first = getLname1(index);
	adj.second = getLname2(index);
	return adj;
}


int EquivalenceClass::getNbAdj() const
{
	return Lnames1.size();
}


pair <int,int> EquivalenceClass::getAncestors() const
{
	return ancestors;
}

int EquivalenceClass::getAncestor( int i) const
{
	if(i == 0)
		return ancestors.first;
	return ancestors.second;
}

/*
Returns:
	(vector <pair <string, string> >): list of the adjacencies of the clas as pairs of string corresponding to pairs of leaf names from gene family 1 and 2
*/
vector <pair <string, string> > EquivalenceClass::getAdjs() const
{
	int nbadj = getNbAdj();

	vector <pair <string, string> > adjs;

	for(unsigned i = 0; i < nbadj ; i++)
	{
		adjs.push_back(getAdj(i));
	}

	return adjs;
}



//setters
void EquivalenceClass::setGfamily1( int fam1)
{
	Gfamily1 = fam1;
}

void EquivalenceClass::setGfamily2( int fam2)
{
	Gfamily2 = fam2;	
}

void EquivalenceClass::setAncestors(pair <int,int> pa)
{
	ancestors = pa;
}

void EquivalenceClass::setAncestor( int i, int a)
{
	if(i==0)
		ancestors.first = a;
	else if(i==1)
		ancestors.second = a;
}



//hasers

/*
Takes:
 - name (string): name of a leaf

Returns:
	(bool): true if Lnames1 contains name

*/
bool EquivalenceClass::hasLname1(string name)
{
	bool has = false;
	int s = Lnames1.size();
	for(unsigned i = 0; i< s ; i++)
	{
		if(Lnames1[i] == name)
		{
			has = true;
			break;
		}
	}
	return has;
}


/*
Takes:
 - name (string): name of a leaf

Returns:
	(bool): true if Lnames1 contains name

*/
bool EquivalenceClass::hasLname2(string name)
{
	bool has = false;
	int s = Lnames2.size();
	for(unsigned i = 0; i< s ; i++)
	{
		if(Lnames2[i] == name)
		{
			has = true;
			break;
		}
	}
	return has;	
}

/*
Takes:
 - name (string): name of a leaf

Returns:
	(int): index of name in Lnames1; -1 if name is absent

*/
int EquivalenceClass::findLname1(string name)
{
	int index = -1;
	int s = Lnames1.size();
	for(unsigned i = 0; i< s ; i++)
	{
		if(Lnames1[i] == name)
		{
			index = i;
			break;
		}
	}
	return index;
}

/*
Takes:
 - name (string): name of a leaf

Returns:
	(int): index of name in Lnames2; -1 if name is absent

*/
int EquivalenceClass::findLname2(string name)
{
	int index = -1;
	int s = Lnames2.size();
	for(unsigned i = 0; i< s ; i++)
	{
		if(Lnames2[i] == name)
		{
			index = i;
			break;
		}
	}
	return index;
}


/*
Checks if the adj exists. works symetrically -> search for name1 in both Lnames1 and Lnames2

Takes:
 - name1 (string): name of a leaf
 - name2 (string): name of a leaf

Returns:
	(bool): true if Lnames1 contains name1 (/ name2) and Lnames2 contains name2 (/ name1) at the same index  

*/
bool EquivalenceClass::hasAdj(string name1, string name2)
{
	int index1 = findLname1(name1);

	bool invert = false;

	if(index1 == -1)//name1 not found in Lnames1, try in Lnames2
	{
		index1 = findLname2(name1);
		invert = true;
	}

	if(index1 == -1) // name1 not found in both Lnames1 and 2 -> adj is absent
		return false;

	string potname2;
	if(invert)//use Lnames1
		potname2 = getLname1(index1);
	else //use Lnames2
		potname2 = getLname2(index1);

	if(name2 == potname2) // name2 if what is expected
		return true; 

	return false;
}



/*
Checks if the adj exists. works symetrically -> search for name1 in both Lnames1 and Lnames2

Takes:
 - names (pair <string,string>): pair of leaf name

Returns:
	(bool): true if Lnames1 contains name1 (/ name2) and Lnames2 contains name2 (/ name1) at the same index  

*/
bool EquivalenceClass::hasAdj(pair <string , string> names)
{
	return hasAdj(names.first, names.second);
}


/*
Gets the index of the adj if it exists. works symetrically -> search for name1 in both Lnames1 and Lnames2

Takes:
 - names (pair <string,string>): pair of leaf name

Returns:
	(int): index of the Adj; -1 if the adj does not exists

*/
int EquivalenceClass::findAdj(string name1, string name2)
{
	int index1 = findLname1(name1);

	bool invert = false;

	if(index1 == -1)//name1 not found in Lnames1, try in Lnames2
	{
		index1 = findLname2(name1);
		invert = true;
	}

	if(index1 == -1) // name1 not found in both Lnames1 and 2 -> adj is absent
		return -1;

	string potname2;
	if(invert)//use Lnames1
		potname2 = getLname1(index1);
	else //use Lnames2
		potname2 = getLname2(index1);

	if(name2 == potname2) // name2 if what is expected
		return index1; 

	return -1;	
}

/*
Gets the index of the adj if it exists. works symetrically -> search for name1 in both Lnames1 and Lnames2

Takes:
 - name1 (string): name of a leaf
 - name2 (string): name of a leaf

Returns:
	(int): index of the Adj; -1 if the adj does not exists

*/
int EquivalenceClass::findAdj(pair <string , string> names)
{
	return findAdj(names.first,names.second);
}


//adders with checks on gfamily names

/*
Adds a pair of Leaf names if fam1 and fam2 correspond to valid Gfamily names 

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
bool EquivalenceClass::CheckAddAdj(string name1, string name2, int fam1, int fam2, int s1 , int s2 )
{
	//cout << "EquivalenceClass::CheckAddAdj " << name1 <<","<< name2 <<" "<< fam1 <<","<< fam2 <<" "<< s1 <<","<< s2 <<endl;

	bool invert = false;

	//checking the validity the the indexes
	if(fam1 != Gfamily1)
	{
		if(fam1 !=Gfamily2)
			return false;//fam1 is neither Gfamily1 or 2 -> abort

		invert = true;
	}

	if(invert)
	{
		if(fam2 != Gfamily1)
			return false;
	}
	else
	{
		if(fam2 != Gfamily2)
			return false;
	}

	if(invert)
	{ // invert is possible only if the two fam are different...
		if( (s2 != sens1) || (s1 != sens2) )
			return false;
	}
	else
	{
		if( Gfamily1 == Gfamily2 )
		{
			if(  ( s1  == sens1 ) && (s2 == sens2) )
				invert = false;
			else if( ( s1  == sens2 ) && (s2 == sens1) ) //same family, but invert because orientation are inverted
				invert = true;
			else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
				return false;
		}
		else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
			return false;
	}

	//cout << "really invert: " << invert << endl;

	//adding.
	if(!invert)
	{
		addLname1(name1);
		addLname2(name2);
	}
	else
	{
		addLname1(name2);
		addLname2(name1);
	}

	//printMe(true);

	return true;
}

/*
Adds a pair of Leaf names if fams correspond to valid Gfamily names 

Takes:
 - names (pair<string,string>): pair of leaf names
 - fams (pair<int,int>): Gfamily indexs of names
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClass::CheckAddAdj(pair <string, string> names, pair <int, int> fams, int s1 , int s2) 
{
	return CheckAddAdj(names.first,names.second, fams.first, fams.second, s1, s2);
}

/*
Adds a list of pair of Leaf names if fams correspond to valid Gfamily names 

Takes:
 - nameList1 (vector<string>): list of leaf names
 - nameList2 (vector<string>): list of leaf names
 - fams (pair<int,int>): Gfamily indexs of names
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClass::CheckAddAdjList(vector <string> nameList1, vector<string> nameList2, pair <int, int> fams, int s1 , int s2)
{
	
	if(nameList1.size() == 0) // empty vector...
		return false;

	//checking the fams
	int fam1 = fams.first;
	int fam2 = fams.second;
	bool invert = false;


	//checking the validity the the indexes
	if(fam1 != Gfamily1)
	{
		if(fam1 != Gfamily2)
			return false;//fam1 is neither Gfamily1 or 2 -> abort
	
		invert = true;
	}
	
	if(!invert)
	{
		if(fam2 != Gfamily2)
			return false;
	}
	else
	{
		if(fam2 != Gfamily1)
			return false;
	}

	if(invert)
	{ // invert is possible only if the two fam are different...
		if( (s2 != sens1) || (s1 != sens2) )
			return false;
	}
	else
	{
		if( Gfamily1 == Gfamily2 )
		{
			if(  ( s1  == sens1 ) && (s2 == sens2) )
				invert = false;
			else if( ( s1  == sens2 ) && (s2 == sens1) ) //same family, but invert because orientation are inverted
				invert = true;
			else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
				return false;
		}
		else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
			return false;
	}



	//adding
	for(unsigned i=0; i < nameList1.size();i++)
	{
		if(!invert)
			addAdj(nameList1[i],nameList2[i]);
		else
			addAdj(nameList2[i],nameList1[i]);
	}

	return true;
}



/*
Adds a list of pair of Leaf names if fams correspond to valid Gfamily names 

Takes:
 - nameList (vector<pair<string,string> >): list of pairs of leaf names
 - fams (pair<int,int>): Gfamily indexs of names
 - int s1 = 0 : orientation of name1
 - int s2 = 0 : orientation of name2

Returns:
 (bool): true if fams correspond to valid indexes and the names could be added
*/
bool EquivalenceClass::CheckAddAdjList(vector <pair <string,string> > nameList, pair <int, int> fams, int s1 , int s2)
{

	if(nameList.size() == 0) // empty vector...
		return false;

	//checking the fams
	int fam1 = fams.first;
	int fam2 = fams.second;
	bool invert = false;

	//checking the validity the the indexes
	if(fam1 != Gfamily1)
	{
		if(fam1 !=Gfamily2)
			return false;//fam1 is neither Gfamily1 or 2 -> abort

		invert = true;
	}

	if(!invert)
	{
		if(fam2 != Gfamily2)
			return false;
	}
	else
	{
		if(fam2 != Gfamily1)
			return false;
	}

	if(invert)
	{ // invert is possible only if the two fam are different...
		if( (s2 != sens1) || (s1 != sens2) )
			return false;
	}
	else
	{
		if( Gfamily1 == Gfamily2 )
		{
			if(  ( s1  == sens1 ) && (s2 == sens2) )
				invert = false;
			else if( ( s1  == sens2 ) && (s2 == sens1) ) //same family, but invert because orientation are inverted
				invert = true;
			else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
				return false;
		}
		else if( ( s1 != sens1 ) || (s2 != sens2) ) // non inverted but not the good orientations 
			return false;
	}


	for(unsigned i=0; i < nameList.size();i++)
	{
		
		if(!invert)
			addAdj(nameList[i].first,nameList[i].second);
		else
			addAdj(nameList[i].second,nameList[i].first);
	}

	return true;
}



//removers


/*
Removes an element from Lnames1 and Lnames2 at a given index.

Takes:
 - index (int): index to remove

Returns:
	(bool) false if the index is not valid
*/
bool EquivalenceClass::removeAdj(int index)
{
	bool res;
	res = removeLname1(index);

	if(!res)
		return false;

	res = removeLname2(index);
	return res;
}



/*
Removes an adjacency given its name.
Works symetrically: will check in both Lnames1 and Lnames for name1

Takes:
 - name1 (string) : name of a leaf
 - name2 (string) : name of a leaf

Returns:
	(bool) false if the names are not valid
*/
bool EquivalenceClass::removeAdj(string name1, string name2)
{
	int pos = findAdj( name1, name2);

	if(pos == -1) //the adjacency could not be found
		return false;

	return removeAdj(pos);
}

/*
Removes an adjacency given its name.
Works symetrically: will check in both Lnames1 and Lnames for name1

Takes:
 - names (pair<string,string>) : leaf names

Returns:
	(bool) false if the names are not valid
*/
bool EquivalenceClass::removeAdj(pair<string, string> names)
{
	return removeAdj(names.first,names.second);
}




//merge and split

/*
Adds the adjs of aClass to self. Does not change aClass

Takes:
 - aClass (EquivalenceClass *)


Returns:
	(bool) false if the GeneFamily index are differents. 
*/
bool EquivalenceClass::merge( const EquivalenceClass * aClass)
{
	pair <int,int> fams;
	fams.first = aClass->getGfamily1();
	fams.second = aClass->getGfamily2();

	int s1,s2;
	s1 = aClass->getSens1();
	s2 = aClass->getSens2();
	return CheckAddAdjList(aClass->Lnames1, aClass->Lnames2, fams, s1, s2)	;
}


/*
Takes:
 - Rtree1 (ReconciledTree *): a Reconciled Tree
 - Rtree2 (ReconciledTree *): a Reconciled Tree
 - forceLTrefining (bool): will refine in DeCoLT's way even if both trees have no transfer
 - verbose (bool)

Returns:
	(vector<EquivalenceClass *>): list of pointers to Equivalence classes that constitutes the refined version of this EquivalenceClass

*/
vector<EquivalenceClass *> EquivalenceClass::refineEqClass(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool forceLTrefining, bool verbose)
{
	vector <int> Id_Genes1;
	vector <int> Id_Genes2;

	for(unsigned i = 0 ;  i <  Lnames1.size(); i++)
	{
		Id_Genes1.push_back( Rtree1->getIdWithName(Lnames1[i]) );//getting ids
		Id_Genes2.push_back( Rtree2->getIdWithName(Lnames2[i]) );
	}

	int Ancestor1 = Rtree1->getRootId();
	int Ancestor2 = Rtree2->getRootId();

	bool DeCoLTrefining = true;

	if(!forceLTrefining)
	{
		if(Rtree1->getNumberOfTransfer() == 0 && Rtree2->getNumberOfTransfer() == 0) // no transfer in both tree -> the original DeCo way of refining equivalence class could save some computations when computing the adjacency tree
		{
			if(verbose)
				cout << "No transfer detected -> Using transfer-less refining." << endl;
			Ancestor1 = TreeTools::getLastCommonAncestor(*Rtree1,Id_Genes1);
			Ancestor2 = TreeTools::getLastCommonAncestor(*Rtree2,Id_Genes2);
			DeCoLTrefining = false;
		}
	}
	else if(verbose)
		cout << "Using tranfer-aware refining." << endl;


	vector<EquivalenceClass *> RefinedClasses; //the new refined classes

	vector < pair<int, int> > AncestorsVec; //we want to keep the ancestors id at hand

	for(unsigned i = 0 ;  i <  Lnames1.size(); i++)
	{
		int NodeId1 = Id_Genes1[i];//getting ids
		int NodeId2 = Id_Genes2[i];

		if(verbose)
			cout << "treating adj " << Lnames1[i] << "(" << NodeId1 << ")" << " - " << Lnames2[i] << "(" << NodeId2 << ")" << endl;
		
		bool ancestorfound = false;
		bool invert =false;
		unsigned j;
		pair<int, int> AdjAncestors;

		if(DeCoLTrefining)
		{
			AdjAncestors = getHighestCompatibleNodeIdsStoppingReception( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2); // getting the ancestors of the adj

			for(j = 0; j < RefinedClasses.size(); j++ ) // First check if the ancestor already exists
			{
				if(AdjAncestors.first == AncestorsVec[j].first)
				{
					if(AdjAncestors.second == AncestorsVec[j].second)
					{
						ancestorfound = true;
						break;						
					}
				}
				else if(Gfamily1 == Gfamily2) // in the case where we are in the same tree, we also have to test the case where ancestor1/2 is the ancestor of node2/1
				{
					if(AdjAncestors.first == AncestorsVec[j].second)
					{
						if(AdjAncestors.second == AncestorsVec[j].first)
						{	// here we'd want to invert the adjacency so that each gene is in the correct pool. However, we can only do that if sens1 == sens2
							if(sens1 == sens2)
							{ 
								invert = true;
								ancestorfound = true;
								break;
							}
						}
					}
				}

			}
		}
		else
		{
			for(j = 0; j < RefinedClasses.size(); j++ ) // First check if the ancestor does not already exists
			{
				//checking if the pair of ancestors are ancestors to the adjacency
				if(Rtree1->isAncestor(AncestorsVec[j].first,NodeId1))
				{
					if(Rtree2->isAncestor(AncestorsVec[j].second,NodeId2))
					{
						ancestorfound = true;
						break;
					}
				}
				else if(Gfamily1 == Gfamily2) // in the case where we are in the same tree, we also have to test the case where ancestor1/2 is the ancestor of node2/1
				{
					if(Rtree1->isAncestor(AncestorsVec[j].first,NodeId2))
					{
						if(Rtree2->isAncestor(AncestorsVec[j].second,NodeId1))
						{
							// here we'd want to invert the adjacency so that each gene is in the correct pool. However, we can only do that if sens1 == sens2
							if(sens1 == sens2)
							{
								invert = true;
								ancestorfound = true;
								break;
							}
						}
					}
				}
			}
		}



		if(ancestorfound) // the ancestor was found. Its index is j
		{
			if(verbose)
				cout << "ancestor found: " << j << ". invert:" << invert << endl;

			if(!invert)
				RefinedClasses[j]->CheckAddAdj(Lnames1[i], Lnames2[i], Gfamily1, Gfamily2, sens1, sens2);
			else
				RefinedClasses[j]->CheckAddAdj(Lnames2[i], Lnames1[i], Gfamily2, Gfamily1, sens2, sens1);
		}
		else
		{//the ancestor was not found -> find it and create a new Equivalence Class
			if(verbose)
				cout << "ancestor not found." << endl;

//			AncestorsVec.push_back( getHighestCompatibleNodeIds( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2) ); // getting the ancestors of the new EquivalenceClass

			if(DeCoLTrefining)
				AncestorsVec.push_back( AdjAncestors );//getHighestCompatibleNodeIdsStoppingReception( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2) ); // getting the ancestors of the new EquivalenceClass
			else
				AncestorsVec.push_back( getHighestCompatibleNodeIds( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2) ); // getting the ancestors of the new EquivalenceClass				

			RefinedClasses.push_back( new EquivalenceClass(Gfamily1,Gfamily2));//creating new EqClass
			RefinedClasses.back()->setAncestors(AncestorsVec.back());//setting th ancestors
			RefinedClasses.back()->setSens1(sens1);
			RefinedClasses.back()->setSens2(sens2);
			RefinedClasses.back()->CheckAddAdj(Lnames1[i], Lnames2[i], Gfamily1, Gfamily2, sens1, sens2); // Adding the adj

			if(verbose)
				cout << " new class with ancestor " << AncestorsVec.back().first << " , " << AncestorsVec.back().second << endl;


		}
	}

	//if(verbose)
	//{
	//	cout << "result of refining: "<< RefinedClasses.size() << " Eclass."<< endl;
	//	for(unsigned i = 0 ; i < RefinedClasses.size() ; i++)
	//		RefinedClasses[i]->printMe(true);
	//}



	return RefinedClasses;

}

/*
Refines the Equivalence Class, considering all possible adjacencies between the two gene families.

Takes:
 - Rtree1 (ReconciledTree *): a Reconciled Tree
 - Rtree2 (ReconciledTree *): a Reconciled Tree
 - verbose (bool)

Returns:
	(vector<EquivalenceClass *>): list of pointers to Equivalence classes that constitutes the refined version of this EquivalenceClass

*/
vector<EquivalenceClass *> EquivalenceClass::refineEqClassWhole(ReconciledTree * Rtree1, ReconciledTree * Rtree2,  bool verbose)
{
	vector <int> Id_Genes1;
	vector <int> Id_Genes2;

	for(unsigned i = 0 ;  i <  Lnames1.size(); i++)
	{
		Id_Genes1.push_back( Rtree1->getIdWithName(Lnames1[i]) );//getting ids
		Id_Genes2.push_back( Rtree2->getIdWithName(Lnames2[i]) );
	}

	int Ancestor1 = Rtree1->getRootId();
	int Ancestor2 = Rtree2->getRootId();

	vector < int > Leaves1 = Rtree1->getLeavesId();
	vector < int > Leaves2 = Rtree2->getLeavesId();

	vector<EquivalenceClass *> RefinedClasses = refineEqClass(Rtree1,  Rtree2,  true,  verbose); //the new refined classes with the adjacencies. 

	vector < pair<int, int> > AncestorsVec; //we want to keep the ancestors id at hand
	for(unsigned k = 0; k < RefinedClasses.size(); k++)
		AncestorsVec.push_back(RefinedClasses[k]->getAncestors());


	//1. creating all equivalence class for all pair of compatible leaves

	for(unsigned i = 0 ; i < Leaves1.size(); i++ )
	{
		int NodeId1 = Leaves1[i];
		int s1 = Rtree1->getNodeSpecies(NodeId1);
		for(unsigned j = 0 ; j < Leaves2.size(); j++ )
		{
			if(Gfamily1 == Gfamily2)
				if(i >= j)
					continue;//ignore if this is the same gene; or avoid to check twice the same pair
			int NodeId2 = Leaves2[j];
			int s2 = Rtree2->getNodeSpecies(NodeId2);

			if(s1 != s2)
				continue;//species are different -> ignore

			//we have two leaves in the same species--> get the highest compatible node

			bool ancestorfound = false;
			bool invert =false;
			unsigned k;
			/*
			for(k = 0; k < RefinedClasses.size(); k++ ) // First check if the ancestor does not already exists
			{
				//checking if the pair of ancestors are ancestors to the adjacency
				if(Rtree1->isAncestor(AncestorsVec[k].first,NodeId1))
					if(Rtree2->isAncestor(AncestorsVec[k].second,NodeId2))
					{
						ancestorfound = true;
						break;
					}
				else if(Gfamily1 == Gfamily2) // in the case where we are in the same tree, we also have to test the case where ancestor1/2 is the ancestor of node2/1
					if(Rtree1->isAncestor(AncestorsVec[k].first,NodeId2))
						if(Rtree2->isAncestor(AncestorsVec[k].second,NodeId1))
						{
							invert = true;
							ancestorfound = true;
							break;
						}
			}*/

			pair <int,int> AdjAncestors = getHighestCompatibleNodeIdsStoppingReception( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2); // getting the ancestors of the adj

			for(k = 0; k < RefinedClasses.size(); k++ ) // First check if the ancestor already exists
			{
				if(AdjAncestors.first == AncestorsVec[k].first)
				{
					if(AdjAncestors.second == AncestorsVec[k].second)
					{
						ancestorfound = true;
						break;						
					}
				}
				else if(Gfamily1 == Gfamily2) // in the case where we are in the same tree, we also have to test the case where ancestor1/2 is the ancestor of node2/1
				{
					if(AdjAncestors.first == AncestorsVec[k].second)
					{
						if(AdjAncestors.second == AncestorsVec[k].first)
						{
							invert = true;
							ancestorfound = true;
							break;
						}
					}
				}

			}
	
			if(ancestorfound) // the ancestor was found. Its index is j
			{
				if(verbose)
					cout << "ancestor found: " << k << ". invert:" << invert << endl;
			}
			else
			{//the ancestor was not found -> find it and create a new Equivalence Class
				if(verbose)
					cout << "ancestor not found." << endl;
	
				AncestorsVec.push_back( AdjAncestors ); // getting the ancestors of the new EquivalenceClass
				RefinedClasses.push_back( new EquivalenceClass(Gfamily1,Gfamily2));//creating new EqClass
				RefinedClasses.back()->setAncestors(AncestorsVec.back());//setting th ancestors
				RefinedClasses.back()->setSens1(sens1);
				RefinedClasses.back()->setSens2(sens2);
	
				if(verbose)
					cout << " new class with ancestor " << AncestorsVec.back().first << " , " << AncestorsVec.back().second << endl;
	
			}
		}
	}

/*
	//2.putting the adjacencies in the created equivalence classes

	for(unsigned i = 0 ;  i <  Lnames1.size(); i++)
	{
		int NodeId1 = Id_Genes1[i];//getting ids
		int NodeId2 = Id_Genes2[i];

		if(verbose)
			cout << "treating adj " << Lnames1[i] << "(" << NodeId1 << ")" << " - " << Lnames2[i] << "(" << NodeId2 << ")" << endl;
		
		bool ancestorfound = false;
		bool invert =false;
		unsigned j;
		for(j = 0; j < RefinedClasses.size(); j++ ) // First check if the ancestor does not already exists
		{
			//checking if the pair of ancestors are ancestors to the adjacency
			if(Rtree1->isAncestor(AncestorsVec[j].first,NodeId1))
				if(Rtree2->isAncestor(AncestorsVec[j].second,NodeId2))
				{
					ancestorfound = true;
					break;
				}
			else if(Gfamily1 == Gfamily2) // in the case where we are in the same tree, we also have to test the case where ancestor1/2 is the ancestor of node2/1
				if(Rtree1->isAncestor(AncestorsVec[j].first,NodeId2))
					if(Rtree2->isAncestor(AncestorsVec[j].second,NodeId1))
					{
						invert = true;
						ancestorfound = true;
						break;
					}
		}

		if(ancestorfound) // the ancestor was found. Its index is j
		{
			if(verbose)
				cout << "ancestor found: " << j << ". invert:" << invert << endl;

			if(!invert)
				RefinedClasses[j]->CheckAddAdj(Lnames1[i], Lnames2[i], Gfamily1, Gfamily2);
			else
				RefinedClasses[j]->CheckAddAdj(Lnames2[i], Lnames1[i], Gfamily2, Gfamily1);
		}
		else///////should never happen... ///////
		{//the ancestor was not found -> find it and create a new Equivalence Class
			if(verbose)
				cout << "ancestor not found." << endl;

			AncestorsVec.push_back( getHighestCompatibleNodeIds( NodeId1, NodeId2, Rtree1, Rtree2,  Ancestor1,  Ancestor2) ); // getting the ancestors of the new EquivalenceClass
			RefinedClasses.push_back( new EquivalenceClass(Gfamily1,Gfamily2));//creating new EqClass
			RefinedClasses.back()->setAncestors(AncestorsVec.back());//setting th ancestors
			RefinedClasses.back()->CheckAddAdj(Lnames1[i], Lnames2[i], Gfamily1, Gfamily2); // Adding the adj

			if(verbose)
				cout << " new class with ancestor " << AncestorsVec.back().first << " , " << AncestorsVec.back().second << endl;


		}
	}
*/
	return RefinedClasses;

}



void EquivalenceClass::printMe(bool verbose)
{
	cout << "Equivalence Class between gene family " << Gfamily1 ;
	if(sens1 == 1)
		cout << " (stop) ";
	else if (sens1 == -1)
		cout << " (start) ";

	cout << " and " << Gfamily2;
	if(sens2 == 1)
		cout << " (stop) ";
	else if (sens2 == -1)
		cout << " (start) ";

	cout << ". " << getNbAdj() << " adjacencies." << endl;

	if(verbose)
	{
		cout << "Adjacency list:" << endl;
		for(unsigned i = 0; i < getNbAdj(); i++)
		{
			cout << Lnames1[i] << " - " << Lnames2[i] << endl;
		}
	}
}

/*
Creates the AdjMatrix of the Equivalence class.
Will replace the ancient one if it exists

Takes:
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
 - map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb: : used by artdeco -> degree of each gene extremities, used to modulate the costs in speciesC0C1.
 - vector <double> adjacencyScoreVec : used by ADseq -> score between 0 and 1 of the adjacencies that denote the confidence we have in them
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
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs
*/
void EquivalenceClass::createAdjMatrixAux(map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb, 
										map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
										vector <double> &adjacencyScoreVec,
										double Gcost, double Bcost, 
										ReconciledTree * rtree1, ReconciledTree * rtree2,
										bool VERBOSE, bool boltzmann ,
										bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
										double temp , double absencePenalty, double adjScoreLogBase, bool interactionMode)
{
	int a1 = getAncestor(0);
	int a2 = getAncestor(1);

	ReconciledTree * sub1 = rtree1->cloneSubtree(a1);
	ReconciledTree * sub2 = rtree2->cloneSubtree(a2);

	//cout << " EC "<< getGfamily1() << "-" << getGfamily2();
	//cout << "ancestors : " << a1<<","<< a2<<endl;
	//cout << " subtrees : "<< sub1->getNumberOfNodes() << " " << sub2->getNumberOfNodes() << endl;
	//cout << sub1->NewickString() << endl;
	//cout << sub2->NewickString() << endl;


	//cout << "nbadj " << getAdjs().size() << endl;
	vector< pair <string,string> > adjs = getAdjs();

	Amat = new AdjMatrix(speciesC0C1, speGeneAdjNb, speGeneExtremitiesAdjNb,
							adjacencyScoreVec, 	Gcost, Bcost, sub1, sub2, adjs ,
							Gfamily1, Gfamily2, LossAware, FamiliesFreeAdjacencies, 
							VERBOSE, boltzmann, temp, absencePenalty, adjScoreLogBase, 
							getSens1(), getSens2());	

	Amat->setInteractionMode(interactionMode);


	//cout << "EquivalenceClass::createAdjMatrix from " << rtree1->getNumberOfNodes() << "*" << rtree2->getNumberOfNodes() << " to "<< sub1->getNumberOfNodes() << "*" << sub2->getNumberOfNodes() << endl;

	delete sub1;
	delete sub2;

	SetAdjMatrix = true;
}

/*
Creates the AdjMatrix of the Equivalence class.
Will replace the ancient one if it exists

Takes:
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
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

Returns:
	(vector <double>) : list of scores associated with each adj in the EC
*/
vector <double> EquivalenceClass::createAdjMatrix(map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb,
										map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
										double Gcost, double Bcost, 
										ReconciledTree * rtree1, ReconciledTree * rtree2, 
										bool VERBOSE, bool boltzmann ,
										bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
										double temp, double absencePenalty , bool interactionMode)
{

	vector <double> adjacencyScoreVec;
	// by default, we put all score at 1
	for(unsigned i = 0 ; i < getNbAdj() ; i++)
		adjacencyScoreVec.push_back(1);

	
	createAdjMatrixAux( speciesC0C1,  speGeneAdjNb, speGeneExtremitiesAdjNb,  adjacencyScoreVec,
							Gcost,  Bcost, 
							rtree1, rtree2,
							VERBOSE, boltzmann ,LossAware, FamiliesFreeAdjacencies, temp , absencePenalty, 10000, interactionMode); //last arg is default adjScoreLogBase (won't be used anyway)

	return adjacencyScoreVec;
}


/*
Creates the AdjMatrix of the Equivalence class.
Will replace the ancient one if it exists

Takes:
 - map < string, map <string , double> > &adjacencyScores : map containing scores associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : used by artdeco -> cost of having an adjacency or not for different species when the adjacency was not in the input
 - map<int, map<string,int> > speGeneAdjNb : used by artdeco -> degree of each gene, used to modulate the costs in speciesC0C1.
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
 - double adjScoreLogBase : base of the log that will be used to go from adjacency confidence score to parsimony costs

Returns:
	(vector <double>) : list of scores associated with each adj in the EC
*/
vector <double> EquivalenceClass::createAdjMatrix(map < string, map <string , double> > &adjacencyScores, 
										map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb,
										map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
										double Gcost, double Bcost, 
										ReconciledTree * rtree1, ReconciledTree * rtree2, 
										bool VERBOSE, bool boltzmann ,
										bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
										double temp, double absencePenalty, double adjScoreLogBase, bool interactionMode )
{

	vector <double> adjacencyScoreVec;
	
	//iterators to find the adjacencies scores
	map < string, map <string , double> >::iterator it1;
	map <string , double>::iterator it2;

	for(unsigned i = 0 ; i < getNbAdj() ; i++)
	{
		double score = 1; // default score if the adjacency is not in the map

		bool found = false;

		it1 = adjacencyScores.find(Lnames1[i]);
		if( it1 != adjacencyScores.end() )
		{ // found the first name
			//it2 = adjacencyScores[ Lnames1[i] ].find(Lnames2[i]); // search for the second name
			it2 = it1->second.find(Lnames2[i]); // search for the second name

			if( it2 != it1->second.end())
			{ // found the second name
				found = true;
				score = it2->second;
			}
		}

		if(!found)
		{ // reverse search: first Lnames2 then Lnames1
			it1 = adjacencyScores.find(Lnames2[i]);
			if( it1 != adjacencyScores.end() )
			{ // found the first name
				//it2 = adjacencyScores[ Lnames1[i] ].find(Lnames2[i]); // search for the second name
				it2 = it1->second.find(Lnames1[i]); // search for the second name
	
				if( it2 != it1->second.end())
				{ // found the second name
					found = true;
					score = it2->second;
				}
			}
		}

		if(VERBOSE)
		{
			if(found)
				cout << "adj "<< Lnames1[i] <<"-"<< Lnames2[i] << "-> "<< score << endl;
			//else
			//	cout << "adj "<< Lnames1[i] <<"-"<< Lnames2[i] << "-> default 1" << endl;
		}

		adjacencyScoreVec.push_back(score);
	}
	createAdjMatrixAux( speciesC0C1,  speGeneAdjNb,
							speGeneExtremitiesAdjNb,
							adjacencyScoreVec,
							Gcost,  Bcost, 
							rtree1, rtree2,
							VERBOSE, boltzmann,
							LossAware, FamiliesFreeAdjacencies,
							temp , absencePenalty, adjScoreLogBase, interactionMode);

	return adjacencyScoreVec;
}




/*
compute the adjMatrix of the equivalence class, if it is set
*/
void EquivalenceClass::computeAdjMatrix()
{
	if(!SetAdjMatrix)
		throw Exception("EquivalenceClass::computeAdjMatrix : AdjMatrix is not set.");


	Amat->computeMatrix();

}


/*
Version of matrix computing that substract reconciliation event scores to co-events in the adjacency matrix

Takes:
 - WDupCost (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WLossCost (double):Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - WHgtCOst (double): Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
*/
void EquivalenceClass::computeAdjMatrix(double WDupCost, double WLossCost, double WHgtCost)
{
	if(!SetAdjMatrix)
		throw Exception("EquivalenceClass::computeAdjMatrix : AdjMatrix is not set.");


	Amat->computeMatrix(WDupCost, WLossCost, WHgtCost);

}



/*
Returns:
 (bool): true if the Adjacency Matrix is computed

*/
bool EquivalenceClass::iscomputedAdjMatrix()
{
	if(!SetAdjMatrix)
		return false; //if the matrix isn't set, then it is not computed

	return Amat->isComputed();
}

/*
Backtracks the AdjMatrix and populates the AdjacencyTrees

Takes:
	- rtree1 (ReconciledTree *): reconciled tree for the first dimension
 	- rtree2 (ReconciledTree *): reconciled tree for the second dimension
	- AdjacencyTrees (vector< AdjTree *> *): pointer to a vector of adjacency trees vector that will be update as we build more adjacency trees
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score
 	- bool &overflow : will be changed to true if there is sign of overflow
 	- alwaysGainAtTop (bool) [default: true]: there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 	- c1proba (double) [default = 0.5]: probability to choose c1 over c0 IF (and only if) they have the same score
 	- bool doNBBT [default = true] : whether or not most parsimonious solutions will be chosen using number of backtrack counts
*/
void EquivalenceClass::backtrackAdjMatrix(ReconciledTree * rtree1, ReconciledTree * rtree2, vector< AdjTree *> * AdjacencyTrees, bool stochastic,bool &overflow, bool alwaysGainAtTop , double c1proba , bool doNBBT)
{
	//cout << "EquivalenceClass::backtrackAdjMatrix" << endl;
	if(!iscomputedAdjMatrix())
		throw Exception("EquivalenceClass::backtrackAdjMatrix : AdjMatrix is not set or computed.");

	Amat->backtrack(AdjacencyTrees,stochastic, alwaysGainAtTop, c1proba, doNBBT);


	if( ( AdjacencyTrees->size() == 0 ) && ( Lnames1.size() > 0 ) )
	{
		overflow = true;
		//cout << "!!POTENTIAL OVERFLOW!!" << endl;
		//cout << "No tree was backtracked on this equivalence class that has " << Lnames1.size() <<" adjacencies."<< endl;
		//for( unsigned i = 0 ; i < Lnames1.size(); i++)
		//	cout << Lnames1[i] <<"-"<<Lnames2[i]<<" ";
		//cout << endl;
		//cout << "If you haven't put any score on them (in the adjacency file) that would justify their absence, then you are probably suffering from overflow."<<endl;
		//cout << "To avoid this problem, choose a different boltzmann.temperature (usually one that is closer to 1)."<<endl;
	}


	for(unsigned i = 0 ; i < AdjacencyTrees->size(); i++)
	{

		AdjacencyTrees->at(i)->refine(rtree1,rtree2, sens1,sens2);
		
		AdjacencyTrees->at(i)->setNodeNames(sens1, sens2);

		AdjacencyTrees->at(i)->setGfamily1(Gfamily1);
		
		AdjacencyTrees->at(i)->setGfamily2(Gfamily2);
		
	}
}

/*
Backtracks the AdjMatrix and populates the AdjForest of the Equivalence class

Takes:
	- rtree1 (ReconciledTree *): reconciled tree for the first dimension
 	- rtree2 (ReconciledTree *): reconciled tree for the second dimension
 	- stochastic (bool): true if the backtrack is to be stochastic; false if the backtrack choses the solution with the best score
 	- alwaysGainAtTop (bool) [default: true]: there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 	- c1proba (double) [default = 0.5]: probability to choose c1 over c0 IF (and only if) they have the same score
 	- bool doNBBT [default = true] : whether or not most parsimonious solutions will be chosen using number of backtrack counts
 	- verbose
*/
void EquivalenceClass::backtrackAdjMatrixForSelf(ReconciledTree * rtree1, ReconciledTree * rtree2, bool stochastic, bool alwaysGainAtTop , double c1proba , bool doNBBT, double verbose)
{
	if(!iscomputedAdjMatrix())
		throw Exception("EquivalenceClass::backtrackAdjMatrixForSelf : AdjMatrix is not set or computed.");




	Amat->backtrack(AdjForest, stochastic, alwaysGainAtTop, c1proba, doNBBT);
	SetAdjForest = true;

	for(unsigned i = 0 ; i < AdjForest->size(); i++)
	{
		AdjForest->at(i)->refine(rtree1,rtree2 , sens1, sens2);
		AdjForest->at(i)->setNodeNames(sens1, sens2);
		AdjForest->at(i)->setGfamily1(Gfamily1);
		AdjForest->at(i)->setGfamily2(Gfamily2);
	}

	if(verbose)
	{
		cout <<  "******" << endl;
		cout << "got :" << AdjForest->size() << " trees" << endl;
	}


	///code for tree validity verification
	/*
	vector <string> seenLeaf; 

	int adjCount = 0;
	for(unsigned i = 0; i < AdjForest->size() ; i++ )
	{
		vector <int> Lid = AdjForest->at(i)->getLeavesId();
		for(unsigned j = 0; j < Lid.size() ; j++ )
		{
			int evt = AdjForest->at(i)->getNodeEvent(Lid[j]);
			if(AdjForest->at(i)->isExtant(evt))
			{
				adjCount++;
				string name = AdjForest->at(i)->getNodeName(Lid[j]);

				bool found = false;
				for(unsigned leafiter = 0; leafiter < seenLeaf.size(); leafiter++)
				{
					if(seenLeaf[leafiter] == name)
					{
						found = true;
						break;
					}
				}

				if(found)
					cout << "DUPLICATED LEAF: " << name << endl;
				else
					seenLeaf.push_back(name);
			}
		}
	}

	for(unsigned i = 0; i < Lnames1.size() ; i++)
	{
		string name = Lnames1[i] +"-" + Lnames2[i] ;
		bool found = false;
		for(unsigned leafiter = 0; leafiter < seenLeaf.size(); leafiter++)
		{
			if(seenLeaf[leafiter] == name)
			{
				found = true;
				break;
			}
		}
		if(!found)
			cout << "UNSEEN LEAF: " << name << endl;
	}
	*/
	/*
	if(adjCount != Lnames1.size())
	{
		cout << "counted " << adjCount << "adjs in the trees. Expected " << Lnames1.size();
		cout << "PBNBADJ" << endl;
		for(unsigned i = 0; i < AdjForest->size() ; i++ )
		{
			cout << AdjForest->at(i)->NewickString(false) << endl;
		}
		cout <<  "***" << endl;
		cout << rtree1->NewickString(false) << endl;
		cout <<  "***" << endl;
		cout << rtree2->NewickString(false) << endl;
		cout <<  "***" << endl;
		cout << endl;
	}

	if(verbose)
	{
		for(unsigned i = 0; i < Lnames1.size() ; i++ )
			cout << Lnames1[i] << "-" << Lnames2[i] << " ; ";
		cout << endl;
	
		cout << "ancestors : " <<  ancestors.first << "-" << ancestors.second << endl; //Id of the ancestors of the EquivalenceClass in the trees of respectively Gfamily1 and Gfamily2
		cout <<  "******" << endl;
	}
	*/
	///end of tree validation code
}

/*
Returns:
	(int): number of trees in AdjForest (throw exception if AdjForest is not set)
*/
int EquivalenceClass::getNbAdjTrees()
{
	if(!SetAdjForest)
		throw Exception("EquivalenceClass::getNbAdjTrees : AdjForest is not set.");
	return AdjForest->size();
}

/*
Returns:
	(int): number of Adjacency Gains in the trees in AdjForest (throw exception if AdjForest is not set)
*/
int EquivalenceClass::getNbAdjGain()
{
	if(!SetAdjForest)
		throw Exception("EquivalenceClass::getNbAdjGain : AdjForest is not set.");
	
	int nbgain = 0;
	for(unsigned i = 0 ; i < AdjForest->size(); i++)
		nbgain += AdjForest->at(i)->countNbGain();
	return nbgain;
}

int EquivalenceClass::getNbAdjBreak()
{
	if(!SetAdjForest)
		throw Exception("EquivalenceClass::getNbAdjBreak : AdjForest is not set.");
	
	int nbbreak = 0;
	for(unsigned i = 0 ; i < AdjForest->size(); i++)
		nbbreak += AdjForest->at(i)->countNbBreak();
	return nbbreak;
}



void EquivalenceClass::clone( const EquivalenceClass * EC)
{

	if(SetAdjMatrix)
		delete Amat;

	if(SetAdjForest)
	{
		for ( vector< AdjTree * >::iterator it = AdjForest->begin() ; it != AdjForest->end(); ++it)
	   	{
	    	delete (*it);
	   	}
	   	AdjForest->clear();
	   	delete AdjForest;
	}

	if(getNbAdj() >0)
	{
		Lnames1.clear();
		Lnames2.clear();
	}

	setAncestors(EC->getAncestors());//cloning ancestros too


	setGfamily1(EC->getGfamily1());
	setGfamily2(EC->getGfamily2());

	setSens1(EC->getSens1());
	setSens2(EC->getSens2());

	for(unsigned i = 0; i < EC->getNbAdj();i++)
		addAdj(EC->getAdj(i));

	SetAdjMatrix = EC->isetAdjMatrix();
	SetAdjForest = EC->isetAdjForest();

	if(SetAdjForest)
	{
		vector <AdjTree * > * Aforest = EC->getAdjForest();

		if(AdjForest != NULL)
		{
			delete AdjForest;
		}

		AdjForest = new vector <AdjTree * >;

		for(unsigned i = 0; i < Aforest->size();i++ )
		{
			AdjForest->push_back(Aforest->at(i)->cloneSubtree(Aforest->at(i)->getRootId()));
		}	
	}




	if(SetAdjMatrix)
	{
		Amat = EC->getAmat()->getClone();
	}
}

/*
Takes:
	- double threshold (default = 200)

Returns:
	(int) the number of score whose absolute log 10 are above the given threshold 
			(in both c1 and c0) 
			(ignoring the cases with worst score (infinite in parcimony))
*/
int EquivalenceClass::getNumberScoreWithAbsLog10Above(double threshold )
{
	if(!SetAdjMatrix)//the adj Matrix is not set -> no need
		return 0; 
	return Amat->getNumberScoreWithAbsLog10Above( threshold );
}