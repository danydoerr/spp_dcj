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

This file contains several classes intended to represent adjacency scores and solutions

Created the: 09-12-2015 
by: Wandrille Duchemin

Last modified the: 31-12-2015
by: Wandrille Duchemin

*/


#include "MyMatrixAdjSolution.h"





/*
Takes:
 - Vsolution (vector<AdjSolution>): AdjSolution to choose index from
 - tochoosefrom (vector<int>): index to choose from
 - nullScore (double)

Returns:
	(vector<int>) : list of indexes which have a non null score
*/
vector<int> MyMatrixAdjSolution::getIndexNonNull(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, double nullScore)
{
	vector <int> res;

	for(unsigned i = 0 ; i < tochoosefrom.size() ; i++)
	{
		if(Vsolution.at(tochoosefrom[i]).score != nullScore)
			res.push_back(tochoosefrom[i]);
	}

	return res;
}

/*
Takes:
 - Vsolution (vector<AdjSolution>): AdjSolution to choose index from
 - tochoosefrom (vector<int>): index to choose from
 - minimum (bool): true if the best score is the smallest one; false if this is the biggest one

Returns:
	(vector<int>) : list of indexes which are considered best
*/
vector<int> MyMatrixAdjSolution::getIndexOnlyBest(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, bool minimum)
{
	vector <int> res;
	double m;
	if(minimum)
		m = numeric_limits<double>::max();
	else
		m = -1;

	for(unsigned i = 0 ; i < tochoosefrom.size() ; i++)
	{
		bool newval = false;

		if(Vsolution[tochoosefrom[i]].score == m)
			res.push_back(tochoosefrom[i]);
		else if(minimum)
		{
//			cout << "plomin: " << Vsolution[tochoosefrom[i]].score << ":" << m << endl;
			if(Vsolution[tochoosefrom[i]].score < m)
				newval = true;
		}
		else
		{
//			cout << "plomax" << Vsolution[tochoosefrom[i]].score << ":" << m << endl;
			if(Vsolution[tochoosefrom[i]].score > m)
				newval = true;
		}

		if(newval)
		{
			res.clear();
			m = Vsolution.at(tochoosefrom[i]).score;
			res.push_back(tochoosefrom[i]);
		}
	}

	return res;
}

/*
Takes:
 - Vsolution (vector<AdjSolution>): AdjSolution to choose index from
 - tochoosefrom (vector<int>): indexes to choose from
 - wishedCoev (bool): true if only Coevoluting solution are to be chosen; false if only non coevolving solution are to be chosen

Returns:
	(vector<int>) : list of indexes which
*/
vector<int> MyMatrixAdjSolution::getIndexCoevStatus(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, bool wishedCoev)
{
	vector <int> res;

	for(unsigned i = 0 ; i < tochoosefrom.size() ; i++)
	{
		if(Vsolution[tochoosefrom[i]].coevent == wishedCoev)
			res.push_back(tochoosefrom[i]);
	}

	return res;
}



/**
 * This function allocates memory for a matrix of size d1 x d2.

Takes:
 - d1 (int)
 - d2 (int)

 */
void MyMatrixAdjSolution::setDim(int d1, int d2)
{
	dim1 = d1;
	dim2 = d2; 

	solutions = new vector<AdjSolution>[dim1*dim2];

}


/**
 * Return the matrix values at the given coordinantes, checking
 * that the coordinantes are valid and the cell has been assigned. 

Takes:
 - j (int): first coordinate
 - z (int): second coordinate

Returns:
 (vector<AdjSolution>): list of matrix values
 */
vector<AdjSolution> MyMatrixAdjSolution::getValue( int j, int z ) const
{
	if( j>dim1 || j<0 || z>dim2 || z<0 ) {
		cout << "coordinates=" << j << "," << z << endl;
		throw Exception("MyMatrixAdjSolution::getValue: invalid coordinate");
	}

	int i = j*dim2 + z;

	return solutions[i];
}



/**
 * Return the matrix values at the given coordinantes, without
 * checking validity.

Takes:
 - j (int): first coordinate
 - z (int): second coordinate

Returns:
 (vector<AdjSolution>): list of matrix values

 */
vector<AdjSolution> MyMatrixAdjSolution::getValueSure( int j, int z ) const
{   
	//without checking that it is not -1, 
	//because in this case it may be -1		 
	int i = j*dim2 + z;
	return solutions[i];
}

/**
 * Set the matrix values.

Takes:
 - j (int): first coordinate
 - z (int): second coordinate
 - values (vector<AdjSolution> &): list of matrix values

 */
void MyMatrixAdjSolution::setValues( int j, int z, vector<AdjSolution> &values ) 
{
	int i = j*dim2 + z;
	solutions[i] = values;
}





void MyMatrixAdjSolution::printMe()
{
    for(unsigned i =0 ; i < dim1; i++)
    {
        for(unsigned j =0 ; j < dim2; j++)
        {
            cout << getValueSure( i,j ).size() << "\t";
        }
        cout << endl;
    }
}
