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

This file contains the functions relative to cost computation for the class AdjMatrix

Created the: 12-05-2016
by: Wandrille Duchemin

Last modified the: 10-01-2018
by: Wandrille Duchemin

*/

#include "AdjMatrix.h"

double AdjMatrix::AggregateScore(vector <double> scores, int nbGain , int nbBreak )
{
	double s;

	if(scores.size() > 1)
	{

		// worstScore is always contaminant -> auto return worstScore (avoid overflow ?)
		if(scores[0] == worstScore)
			return worstScore;
		else if(scores[1] == worstScore)
			return worstScore;


		s = ((*this).*scoreAggregatorfunc)(scores[0], scores[1]);

		for(unsigned i = 2; i < scores.size(); i++)
		{
			if(scores[i] == worstScore)
				return worstScore;
			s = ((*this).*scoreAggregatorfunc)(s, scores[i]);
		}
	}
	else if(scores.size() == 1)
		s = scores[0];
	else
	{//giving base value?
		s = bestScore; //to adapt later for boltzmann refinement
	}

	//applying Gains and Breaks
	if(!useBoltzmann)
	{
		for(unsigned i = 0; i < nbGain; i++)
			s = ((*this).*scoreAggregatorfunc)(s, getGainCost());
	
		for(unsigned i = 0; i < nbBreak; i++)
			s = ((*this).*scoreAggregatorfunc)(s, getBreakCost());
	}
	else
		s = ((*this).*scoreAggregatorfunc)(s, getBoltzmannGainBreakCost(nbGain, nbBreak));

	return s;

}

/*
Takes:
 - Vsolution (vector <AdjSolution> ): a list of different solution

Returns
	(double): score corresponding to the best solution (non boltzmann case) or the sum of the scores of the solutions (boltzmann case)

*/
double AdjMatrix::compareScore( vector <AdjSolution>  Vsolution)
{
	vector <double> scorevector;


	for(unsigned i = 0; i < Vsolution.size(); i++)
	{
		scorevector.push_back( Vsolution.at(i).score );
	}

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
sets VsolutionC1 and VsolutionC0 for NodeId1 and NodeId2

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - VsolutionC1 (vector <AdjSolution> &): vector of solution that will be filled
 - VsolutionC1 (vector <AdjSolution> &): vector of solution that will be filled

*/
void AdjMatrix::computeSolution(int NodeId1, int NodeId2, vector <AdjSolution> &VsolutionC1, vector <AdjSolution> &VsolutionC0 )
{


	//First, determine if nodeid1 and nodeid2 are comparable
	bool incomparable = false;


	Node * n1 = Rtree1.getNode(NodeId1);
	Node * n2 = Rtree2.getNode(NodeId2);

	// both node are TS compatible
	int evt1 = Rtree1.getNodeEvent(NodeId1);
	int evt2 = Rtree2.getNodeEvent(NodeId2);
	

	// The choice is made that the TS compatibility is set by the reconciled trees rather than through the Adj algorithm.
	// This, for instance, allows a simple DeCo algorithm to perform on a full time sliced tree and respecting TS contraints

	if(!Rtree1.areTSCompatible(n1,n2))//if the node aren't TS compatible: their adjacency is impossible
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " are not time compatible. Setting values accordingly." << endl;
		}

		incomparable =true;
	}
	else if(!Rtree1.haveSameSpecies(n1, n2))
	{
		incomparable =true;	
		if(decoLTalgo) // in the case of Bout + Rec -> the species compatibility need not necessarily hold
		{
			if( ( ( Rtree1.isBout(evt1)) && ( Rtree1.isRec(evt2) ) ) || ( ( Rtree1.isBout(evt2)) && ( Rtree1.isRec(evt1) ) ) )
				incomparable = false;
		}
	}


	if(incomparable)
		VsolutionC1 =  SolutionC1DefaultImpossibleCase();


	// special case --> one of the vent is a loss
	bool firstLoss = Rtree1.isLoss(evt1);
	bool secondLoss = Rtree2.isLoss(evt2);
	
	if( (firstLoss) || (secondLoss) )
	{
		if( (firstLoss) && (secondLoss) )
		{
			if(!incomparable)
				VsolutionC1 = SolutionC1LossWithLoss(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0LossWithLoss(NodeId1,NodeId2);
			return;
		}
		else
		{
			if(!incomparable)
				VsolutionC1 = SolutionC1LossWithOther(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0LossWithOther(NodeId1,NodeId2);
			return;
		}
	}


	//special case --> the events are extants
	int nbSons1 = n1->getNumberOfSons();
	int nbSons2 = n2->getNumberOfSons();

	if((evt1 == evt2) && (nbSons1 == 0)) // extant with extant
	{
		if(!incomparable)
			VsolutionC1 = SolutionC1ExtantWithExtant( NodeId1, NodeId2);
		VsolutionC0 = SolutionC0ExtantWithExtant( NodeId1, NodeId2);
		return;
	}


	//out of the special cases. 
	vector<AdjSolution> tmpC1;
	vector<AdjSolution> tmpC0;

	bool inTheDead = ((Rtree1.getNodeSpecies(NodeId1) == -1) || (Rtree2.getNodeSpecies(NodeId2) == -1));
	bool isRec = ((Rtree1.isRec(evt1)) && (Rtree2.isRec(evt2))); // only true if both event are receptions
	bool isSout = ((Rtree1.isSout(evt1)) && (Rtree2.isSout(evt2))); // only true if both event are SpeciationOut
	bool isSpe1 = Rtree1.isSpeciation(evt1); // only true if both event are SpeciationOut
	bool isSpe2 = Rtree2.isSpeciation(evt2);

	bool noAsynchC0 = (((isSpe1)&&(isSpe2))&&(!incomparable)); // having double speciation compute asynchronous c0 leads to an ambiguity that triplicate leaf-forests-like solutions.
	//additionnaly, we can do asynchronous computation is the two species are incomparable

	// n1 happened before n2

	//C1
	if(!incomparable)
	{
		if(nbSons1 == 1)
		{
			tmpC1 = SolutionC1asynchronousOneChild( NodeId1, NodeId2, false,  inTheDead, isRec);
		}
		else if( (nbSons1 == 2) && (!isSpe1)) //if n1 is a speciation, it can't be before something
		{
			tmpC1 = SolutionC1asynchronousTwoChildren( NodeId1, NodeId2, false, inTheDead);
		}
		//adding to the already existing solutions
		if(VsolutionC1.size() == 0)
		{
			VsolutionC1 = tmpC1;
		}
		else
		{
			for(unsigned i = 0; i < tmpC1.size();i++)
			{
				if(tmpC1[i].score <= VsolutionC1.front().score)
					VsolutionC1.insert(VsolutionC1.begin(),tmpC1[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
				else
					VsolutionC1.push_back(tmpC1[i]);
			}
		}

	}

	//C0

	if(nbSons1 == 1)
	{
		tmpC0 = SolutionC0asynchronousOneChild( NodeId1, NodeId2, false);
	}
	else if(nbSons1 == 2)
	{
		if( !((isSpe1)&&(!incomparable)) )// n1 is a speciation -> a comparable node can't be after it
			tmpC0 = SolutionC0asynchronousTwoChildren( NodeId1, NodeId2, false);

		//if(!noAsynchC0)// having double speciation compute asynchronous c0 leads to an ambiguity that triplicate leaf-forests-like solutions.
	}


	//adding to the already existing solutions
	if(VsolutionC0.size() == 0)
	{
		VsolutionC0 = tmpC0;
	}
	else
	{
		for(unsigned i = 0; i < tmpC0.size();i++)
		{
			if(tmpC0[i].score <= VsolutionC0.front().score)
				VsolutionC0.insert(VsolutionC0.begin(),tmpC0[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
			else
				VsolutionC0.push_back(tmpC0[i]);
	
		}
	}

	tmpC1.clear();
	tmpC0.clear();


	// n2 happened before n1

	//C1
	if(!incomparable)
	{
		if(nbSons2 == 1)
		{
			tmpC1 = SolutionC1asynchronousOneChild( NodeId1, NodeId2, true,  inTheDead, isRec);
		}
		else if( (nbSons2 == 2) && ( !isSpe2 ) ) //if n2 is a speciation, it can't be before something
		{
			tmpC1 = SolutionC1asynchronousTwoChildren( NodeId1, NodeId2, true, inTheDead);
		}

		//adding to the already existing solutions
		if(VsolutionC1.size() == 0)
		{
			VsolutionC1 = tmpC1;
		}
		else
		{
			for(unsigned i = 0; i < tmpC1.size();i++)
			{
				if(tmpC1[i].score <= VsolutionC1.front().score)
					VsolutionC1.insert(VsolutionC1.begin(),tmpC1[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
				else
					VsolutionC1.push_back(tmpC1[i]);
			}
		}

	}

	//C0
	if(nbSons2 == 1)
	{
		tmpC0 = SolutionC0asynchronousOneChild( NodeId1, NodeId2, true);
	}
	else if(nbSons2 == 2)
	{
		if( !((isSpe2)&&(!incomparable)) )// n2 is a speciation -> a comparable node can't be after it
			tmpC0 = SolutionC0asynchronousTwoChildren( NodeId1, NodeId2, true);
		//else
		//	cout << "    new restriction ON"<<endl;
		//if(!noAsynchC0)// having double speciation compute asynchronous c0 leads to an ambiguity that triplicate leaf-forests-like solutions.
		//if(!noAsynchC0)// having double speciation compute asynchronous c0 leads to an ambiguity that triplicate leaf-forests-like solutions.
		//	tmpC0 = SolutionC0asynchronousTwoChildren( NodeId1, NodeId2, true);
	}


	//adding to the already existing solutions
	if(VsolutionC0.size() == 0)
	{
		VsolutionC0 = tmpC0;
	}
	else
	{
		for(unsigned i = 0; i < tmpC0.size();i++)
		{
			if(tmpC0[i].score <= VsolutionC0.front().score)
				VsolutionC0.insert(VsolutionC0.begin(),tmpC0[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
			else
				VsolutionC0.push_back(tmpC0[i]);
	
		}
	}

	tmpC1.clear();
	tmpC0.clear();

	// n1 at the same time as n2

	bool iscoevent = (evt1 == evt2);

	bool doSynchC1 = false;
	bool doSynchC0 = true; // we do synch 0 as long as we do not treat a speciation to non speciation in a comparable context

	if(!incomparable)// condition to do synchronous C1
	{ 
		doSynchC1 = true;
		if( ( !iscoevent ) && ( ( isSpe1 ) || ( isSpe2 ) ) ) // not a coevent, but one is a speciation -> can't be c1 and synchronous
		{
			doSynchC1 = false;
			doSynchC0 = false;
		}
	}


	if(nbSons1 == 1)
	{
		if(nbSons2 == 1)
		{
			if(doSynchC1)
				tmpC1 = SolutionC1synchronousOneChild( NodeId1,  NodeId2,  inTheDead, iscoevent);

			if(doSynchC0)
				tmpC0 = SolutionC0synchronousOneChild(NodeId1, NodeId2);
		}
		else if( nbSons2 == 2)
		{
			if(doSynchC1)
				tmpC1 = SolutionC1synchronousTwoAndOneChildren( NodeId1, NodeId2, inTheDead);

			if(doSynchC0)
				tmpC0 = SolutionC0synchronousTwoAndOneChildren( NodeId1, NodeId2);
		}
	}
	else if(nbSons1 == 2)
	{
		if(nbSons2 == 1)
		{
			if(doSynchC1)
				tmpC1 = SolutionC1synchronousTwoAndOneChildren( NodeId1, NodeId2, inTheDead);

			if(doSynchC0)
				tmpC0 = SolutionC0synchronousTwoAndOneChildren( NodeId1, NodeId2);
		}
		else if( nbSons2 == 2)
		{
			if(doSynchC1)
				tmpC1 = SolutionC1synchronousTwoChildren( NodeId1, NodeId2, inTheDead,  iscoevent, isSout);

			if(doSynchC0)
				tmpC0 = SolutionC0synchronousTwoChildren( NodeId1, NodeId2);
		}	
	}

	if(doSynchC1)
	{
		//adding to the already existing solutions
		if(VsolutionC1.size() == 0)
		{
			VsolutionC1 = tmpC1;
		}
		else
		{
			for(unsigned i = 0; i < tmpC1.size();i++)
			{
				if(tmpC1[i].score <= VsolutionC1.front().score)
					VsolutionC1.insert(VsolutionC1.begin(),tmpC1[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
				else
					VsolutionC1.push_back(tmpC1[i]);
			}
		}

	}

	if(doSynchC0)
	{
		//adding to the already existing solutions
		if(VsolutionC0.size() == 0)
		{
			VsolutionC0 = tmpC0;
		}
		else
		{
			for(unsigned i = 0; i < tmpC0.size();i++)
			{
				if(tmpC0[i].score <= VsolutionC0.front().score)
					VsolutionC0.insert(VsolutionC0.begin(),tmpC0[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
				else
					VsolutionC0.push_back(tmpC0[i]);
		
			}
		}
	}

	if(VsolutionC0.size() == 0)
		VsolutionC0 = SolutionC0DefaultImpossibleCase();
	//else
	//{
	//	cout << "VsolutionC0.front().score " << VsolutionC0.front().score << endl;
	//}


	if(VsolutionC1.size() == 0)
		VsolutionC1 = SolutionC1DefaultImpossibleCase();
	

}

/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with worst score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC1DefaultImpossibleCase()
{
	vector <AdjSolution> Vsolution;
	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().score = worstScore;

	return Vsolution;
}


// simple DeCo cases
/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector <AdjSolution> AdjMatrix::SolutionC1ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1ExtantWithExtant"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,true) );

	Vsolution.back().score = getC1(NodeId1,NodeId2); // the score of that solution

	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> C1:" << Vsolution.back().score << endl;

	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector <AdjSolution> AdjMatrix::SolutionC0ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0ExtantWithExtant"<< endl;


	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = getC0(NodeId1,NodeId2); // the score of that solution

	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> C0:" << Vsolution.back().score << endl;

	return Vsolution;
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
vector <AdjSolution> AdjMatrix::SolutionC1LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithOther"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
vector <AdjSolution> AdjMatrix::SolutionC0LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithOther"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;

}




/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
vector<AdjSolution> AdjMatrix::SolutionC1LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithLoss"<< endl;
	
	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,true) );

	Vsolution.back().score = WeightedLossCost; // Always the WeightedLossCost as this is either bestScore (default) or a value to add/multiply to bestscore = that value

	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
vector<AdjSolution> AdjMatrix::SolutionC0LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithLoss"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;	
}


/*
Function specific to the case where the nodes have only one child
Two snodes can only be synchronous if they have the same event, which also means that they have the same number of children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool inTheDead: true if the events are in the dead
 - bool iscoevent : true if both events are of the same nature

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1synchronousOneChild(int NodeId1, int NodeId2, bool inTheDead, bool iscoevent)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC1synchronousOneChild"<< endl;


	vector <AdjSolution> Vsolution;

	vector <int> SonsIds1 = Rtree1.getSonsId(NodeId1); //should be size 1
	vector <int> SonsIds2 = Rtree2.getSonsId(NodeId2); //should be size 1


	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1a1b1;
	Sc1a1b1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0a1b1;
	Sc0a1b1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)

	Sc1a1b1.id1 = SonsIds1[0];
	Sc1a1b1.id2 = SonsIds2[0];

	Sc0a1b1.id1 = SonsIds1[0];
	Sc0a1b1.id2 = SonsIds2[0];

	double c1a1b1 = getC1(Sc1a1b1.id1,Sc1a1b1.id2);
	double c0a1b1 = getC0(Sc0a1b1.id1,Sc0a1b1.id2);


	int nbBreak = 0;

	// 2 possible solutions
	//case 1 -> the adj is not present >> break
	if(!inTheDead)
		nbBreak = 1;


	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 1 break if not in the dead
	currentSolution->components.push_back( Sc0a1b1 );
	caseScore.push_back( c0a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 0 gain
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	caseScore.push_back( c1a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;
}




/*
Function specific to the case where the nodes have two children
Two snodes can only be synchronous if they have the same event, which also means that they have the same number of children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool inTheDead: whether the event is in the dead or not (breaks are not counted in the dead)
 - bool iscoevent : true if both events are of the same nature
 - bool isSout: true is the events are Souts
Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1synchronousTwoChildren(int NodeId1, int NodeId2, bool inTheDead, bool iscoevent, bool isSout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC1synchronousTwoChildren"<< endl;

	vector<AdjSolution> Vsolution;

	AdjSolution * currentSolution;

	//we first get the ids of the children of NodeIdDup
	vector <int> SonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> SonsId2 = Rtree2.getSonsId(NodeId2);

	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> children1: "<< SonsId1[0] << ","<< SonsId1[1] << " - children2: "<< SonsId2[0] << ","<< SonsId2[1] << endl;

	// specifically used for Sout case
	int transferredCouple = -1;

	if(iscoevent)
	{
		if(isSout)
		{
			int CorrectSonInd1 = -1;
			int CorrectSonInd2 = -1;
	
	
			int sp =  Rtree1.getNodeSpecies(NodeId1);
			for(unsigned i = 0; i < SonsId1.size(); i++)
			{
				if(Rtree1.getNodeSpecies(SonsId1[i]) == sp)
				{
					CorrectSonInd1 = i;
					break;
				}
			}
			//sanity check
			if(CorrectSonInd1 == -1)
				throw Exception("AdjMatrix::SolutionC1SoutWithSoutSynchronous : found no son of the Sout node that stayed in the same species.");
			
	
			//we want to get the Son of the Sout that stayed in the species tree for the node 2
			
			sp =  Rtree2.getNodeSpecies(NodeId2);
			for(unsigned i = 0; i < SonsId2.size(); i++)
			{
				if(Rtree2.getNodeSpecies(SonsId2[i]) == sp)
				{
					CorrectSonInd2 = i;
					break;
				}
			}
			//sanity check
			if(CorrectSonInd2 == -1)
				throw Exception("AdjMatrix::SolutionC1SoutWithSoutSynchronous : found no son of the Sout node that stayed in the same species.");
	
	
			transferredCouple = 2 * CorrectSonInd1 + CorrectSonInd2; 
	
		}
	}


	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
	{
		for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
		{
			for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
			{
				for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
				{
					currentSolution = new AdjSolution();
					currentSolution->coevent = iscoevent;


					if( c1a1b2 && c1a2b1 ) // in that case, put a1b1 and a2b2 first so they bear the potential gain
					{
						//a1b1
						currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[0]));
						if(!c1a1b1)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[0],SonsId2[0]));
						}
						else
						{
							caseScore.push_back(getC1(SonsId1[0],SonsId2[0]));
						}

						//a2b2
						currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[1]));
						if(!c1a2b2)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[1],SonsId2[1]));
						}
						else
						{
							caseScore.push_back(getC1(SonsId1[1],SonsId2[1]));
						}

						//a1b2 -> we know they are at 1
						currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[1]));
						caseScore.push_back(getC1(SonsId1[0],SonsId2[1]));

						//a2b1 -> we know they are at 1
						currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[0]));
						caseScore.push_back(getC1(SonsId1[1],SonsId2[0]));
					}
					else // default: put a1b2 and a2b1 first so they bear the potential gain
					{
						//a1b2
						currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[1]));
						if(!c1a1b2)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[0],SonsId2[1]));
						}
						else
							caseScore.push_back(getC1(SonsId1[0],SonsId2[1]));

						//a2b1
						currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[0]));
						if(!c1a2b1)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[1],SonsId2[0]));
						}
						else
						{
							caseScore.push_back(getC1(SonsId1[1],SonsId2[0]));
						}

						//a1b1
						currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[0]));
						if(!c1a1b1)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[0],SonsId2[0]));
						}
						else
						{
							caseScore.push_back(getC1(SonsId1[0],SonsId2[0]));
						}

						//a2b2
						currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[1]));
						if(!c1a2b2)
						{
							currentSolution->components.back().C1 = false;
							caseScore.push_back(getC0(SonsId1[1],SonsId2[1]));
						}
						else
						{
							caseScore.push_back(getC1(SonsId1[1],SonsId2[1]));
						}
							
					}


					int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;


					if(interactionMode)
					{ // in interaction mode, we expect a total of 4 adjacencies between children of a and b.
					  // any less than that cost a Break, except in the dead.

						currentSolution->NbGain = 0;
						if(inTheDead)
							currentSolution->NbBreak = 0;
						else
							currentSolution->NbBreak = nbadj-4;

					}
					else
					{
						//determining the number of gains and breaks
						switch( nbadj )
						{
							case 0: // none linked -> 2 break if covent, 1 if not ; unless we are in the dead (then no break)
								currentSolution->NbGain = 0;
								if(inTheDead)
									currentSolution->NbBreak = 0;
								else if(iscoevent)
									currentSolution->NbBreak = 2;
								else
									currentSolution->NbBreak = 1;
								break;
							case 1:// only one linked -> 1 break unless we are in the dead or not co-events
								currentSolution->NbGain = 0;
								if(inTheDead)
									currentSolution->NbBreak = 0;
								else if(iscoevent)
									currentSolution->NbBreak = 1;
								else
									currentSolution->NbBreak = 0;
								break;
							case 3: // three linked -> 1 gain if coevent; 2 otherwise
								if(iscoevent)
									currentSolution->NbGain = 1;
								else
									currentSolution->NbGain = 2;
								currentSolution->NbBreak = 0;	
								break;
							case 4: // all linked -> 2 gains if coevent; 3 otherwise
								if(iscoevent)
									currentSolution->NbGain = 2;
								else
									currentSolution->NbGain = 3;
								currentSolution->NbBreak = 0;	
								break;
							default: //exactly 2 adjs
								currentSolution->NbGain = 0;
								currentSolution->NbBreak = 0;
								if(iscoevent)
								{
									if(( c1a1b1 == c1a1b2 ) || (c1a1b1 == c1a2b1 )) // scenarii where 1 child has two adjs -> one gain and 1 break unless we are in the dead
									{
										currentSolution->NbGain = 1;
										if(inTheDead)
											currentSolution->NbBreak = 0;
										else
											currentSolution->NbBreak = 1;
									}
								}
								else // not a coevent -> exactly 1 gain
								{
									currentSolution->NbGain = 1;
								}
						}
					}
	
					if(isSout) // we see if we need to correct the number of break (if the transferred couple wase lost)
					{
						bool correct = false;
						switch(transferredCouple)
						{
							case 0:
								correct = (!c1a1b1);
								break;
							case 1:
								correct = (!c1a1b2);
								break;
							case 2:
								correct = (!c1a2b1);
								break;
							case 3:
								correct = (!c1a2b2);
								break;
							default:
								break;
						}

						currentSolution->NbBreak -= correct;

						//sanity check
						if(currentSolution->NbBreak < 0)
							currentSolution->NbBreak = 0;
					}

					

					//getting the new score
					currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);
					
					if(Vsolution.empty())
						Vsolution.push_back(*currentSolution);
					else
					{
						if(currentSolution->score <= Vsolution.front().score)
							Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
						else
							Vsolution.push_back(*currentSolution);
					}
					//clearing
					delete currentSolution;
					caseScore.clear();
				}
			}
		}
	}


	return Vsolution;
}








/*
Function specific to the case where one node has two children and the other only has one.

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool inTheDead: whether 

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1synchronousTwoAndOneChildren(int NodeId1, int NodeId2, bool inTheDead)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC1synchronousTwoAndOneChildren"<< endl;


	bool iscoevent = false; // never a co-event
	vector <AdjSolution> Vsolution;

	vector <int> SonsIds1;
	vector <int> SonsIds2;

 	SonsIds1 = Rtree1.getSonsId(NodeId1);
	SonsIds2 = Rtree2.getSonsId(NodeId2);


	bool NodeId2TwoChildren = false;

	if(SonsIds1.size() == 1)
		NodeId2TwoChildren = true;


	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	Sc1ab1.id1 = SonsIds1[0];
	Sc1ab1.id2 = SonsIds2[0];

	Sc0ab1.id1 = SonsIds1[0];
	Sc0ab1.id2 = SonsIds2[0];

	if(NodeId2TwoChildren) // node 2 has 2 children
	{
		Sc1ab2.id1 = SonsIds1[0];
		Sc1ab2.id2 = SonsIds2[1];

		Sc0ab2.id1 = SonsIds1[0];
		Sc0ab2.id2 = SonsIds2[1];
	}
	else // node 1 has two children
	{
		Sc1ab2.id1 = SonsIds1[1];
		Sc1ab2.id2 = SonsIds2[0];

		Sc0ab2.id1 = SonsIds1[1];
		Sc0ab2.id2 = SonsIds2[0];
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);



	int nbBreak = 0;

	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency -> 1 break if not in the dead
	if(!inTheDead)
	{
		nbBreak = 1;
		if(interactionMode)
			nbBreak++;
	}


	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 1 break if not in the dead
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	nbBreak = 0;
	//case 2 -> child 1  kept an adjacency 
	if( (interactionMode) && (!inTheDead))
		nbBreak = 1;

	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	nbBreak = 0;
	//case 3 -> child 2 kept an adjacency 
	if( (interactionMode) && (!inTheDead))
		nbBreak = 1;

	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);




	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency -> 1 gain if we are in "neighborhood" mode
	//                                                -> 0 gain if we are in "interaction" mode
	int nbGain = 1;
	if(interactionMode)
		nbGain = 0;
	currentSolution = new AdjSolution(nbGain,0,iscoevent) ;  // 0/1 gain 0 break 

	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	return Vsolution;


}





/*
Function specific to the case where the more ancient node has only one child

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool NodeId2First: true if nodeid2 is before nodeid1
 - bool inTheDead: true if we are in the dead
 - bool isRec: true if both events are rec

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1asynchronousOneChild(int NodeId1, int NodeId2, bool NodeId2First, bool inTheDead, bool isRec)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC1asynchronousOneChild "<< NodeId2First << endl;

	vector <AdjSolution> Vsolution;

	vector <int> SonsIds;

	if(!NodeId2First) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);

	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)

	if(NodeId2First) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);


	int nbBreak = 0;
	int nbGain = 0;

	// 2 possible solutions
	//case 1 -> the adj is not present << nothing happens
	if(!inTheDead)
		nbBreak = 1;

	if(isRec)
		nbGain = 1; // two asynchronous rec == automatically a gain


	currentSolution = new AdjSolution(nbGain,nbBreak,iscoevent) ;  // 0 gain 1 break if in the dead
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 0 gain
	currentSolution = new AdjSolution(nbGain,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	caseScore.push_back( c1ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;

}




/*
Function specific to the case where the more ancient node has two children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool NodeId2First: true if nodeid2 is before nodeid1
 - bool inTheDead: whether 

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1asynchronousTwoChildren(int NodeId1, int NodeId2, bool NodeId2First, bool inTheDead)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC1asynchronousTwoChildren "<< NodeId2First << endl;

	vector <AdjSolution> Vsolution;

	vector <int> SonsIds;

	if(!NodeId2First) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);

	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	if(NodeId2First) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

		Sc1ab2.id1 = NodeId1;
		Sc1ab2.id2 = SonsIds[1];

		Sc0ab2.id1 = NodeId1;
		Sc0ab2.id2 = SonsIds[1];


	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


		Sc1ab2.id1 = SonsIds[1];
		Sc1ab2.id2 = NodeId2;

		Sc0ab2.id1 = SonsIds[1];
		Sc0ab2.id2 = NodeId2;

	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);



	int nbBreak = 0;

	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency -> 1 break if not in the dead
	if(!inTheDead)
	{
		nbBreak = 1;
		if(interactionMode)
			nbBreak++;
	}

	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 1 break if not in the dead
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	nbBreak = 0;
	//case 2 -> child 1  kept an adjacency 
	if( (!inTheDead) && (interactionMode) )
		nbBreak++;

	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	nbBreak = 0;
	//case 3 -> child 2 kept an adjacency 
	if( (!inTheDead) && (interactionMode) )
		nbBreak++;	

	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);




	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency -> 1 gain if we are in "neighborhood" mode
	//                                                -> 0 gain if we are in "interaction" mode
	int nbGain = 1;
	if(interactionMode)
		nbGain = 0;
	currentSolution = new AdjSolution(nbGain,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	return Vsolution;


}


/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCase()
{
	if(verbose)
		cout << "SolutionC0DefaultImpossibleCase"<< endl;

	vector <AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().score = bestScore;

	return Vsolution;
}



/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCaseNew(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0DefaultImpossibleCaseNew"<< endl;


	vector <AdjSolution> Vsolution;

	//NodeId1 before
	vector <AdjSolution> Vsolution1 = SolutionC0DefaultImpossibleCaseAux( NodeId1, NodeId2, false);

	for(unsigned i = 0; i < Vsolution1.size();i++)
	{
		Vsolution.push_back(Vsolution1[i]);
		//cout << "AdjMatrix::SolutionC0DefaultImpossibleCaseNew " << NodeId1 << "-" << NodeId2 << " n1 before " << Vsolution1[i].score << endl;
	}




	//NodeId2 before
	vector <AdjSolution> Vsolution2 = SolutionC0DefaultImpossibleCaseAux( NodeId1, NodeId2, true );

	for(unsigned i = 0; i < Vsolution2.size();i++)
	{
		if(Vsolution.size() != 0)
		{
			if(Vsolution2[i].score <= Vsolution.front().score)
			{
				Vsolution.insert(Vsolution.begin(),Vsolution2[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
			}
		}
		else
		{
			Vsolution.push_back(Vsolution2[i]);
		}
		//cout << "AdjMatrix::SolutionC0DefaultImpossibleCaseNew " << NodeId1 << "-" << NodeId2 << " n2 before " << Vsolution2[i].score << endl;
	}

	if(Vsolution.size() == 0) // leaf|loss  with  leaf|loss
		Vsolution = SolutionC0DefaultImpossibleCase();

	return Vsolution;
}



vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCaseAux(int NodeId1, int NodeId2, bool invert)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0DefaultImpossibleCaseAux"<< endl;

	vector <AdjSolution> Vsolution;
	
	vector <int> SonsIds;

	if(!invert) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);


	if(SonsIds.size() == 0) // cases of losses and leaves
		return Vsolution;



	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	if(invert) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);


	if(SonsIds.size() == 1) // ONLY one son 
	{
		// 2 possible solutions

		//case 1 -> the adj is not present << nothing happens
		currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
		currentSolution->components.push_back( Sc0ab1 );
		caseScore.push_back( c0ab1 );
		currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

		Vsolution.push_back(*currentSolution);

		delete currentSolution;
		caseScore.clear();

		//case 2 -> the adj is present present << 1 gain
		currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
		currentSolution->components.push_back( Sc1ab1 );
		caseScore.push_back( c1ab1 );
		currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

		if(currentSolution->score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(*currentSolution);


		delete currentSolution;
		caseScore.clear();

		return Vsolution;
	}


	/// two sons (more than that is not covered yet)

	if(invert) // node 2 before
	{
		Sc1ab2.id1 = NodeId1;
		Sc1ab2.id2 = SonsIds[1];

		Sc0ab2.id1 = NodeId1;
		Sc0ab2.id2 = SonsIds[1];
	}
	else // ndoe 1 before
	{
		Sc1ab2.id1 = SonsIds[1];
		Sc1ab2.id2 = NodeId2;

		Sc0ab2.id1 = SonsIds[1];
		Sc0ab2.id2 = NodeId2;
	}

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);


	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency 
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> child 1  kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 2 kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	return Vsolution;
}




/*
Function specific to the case where the more ancient node has only one child

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0synchronousOneChild(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0synchronousOneChild"<< endl;


	vector <AdjSolution> Vsolution;

	vector <int> SonsIds1 = Rtree1.getSonsId(NodeId1); //should be size 1
	vector <int> SonsIds2 = Rtree2.getSonsId(NodeId2); //should be size 1


	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1a1b1;
	Sc1a1b1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0a1b1;
	Sc0a1b1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)

	Sc1a1b1.id1 = SonsIds1[0];
	Sc1a1b1.id2 = SonsIds2[0];

	Sc0a1b1.id1 = SonsIds1[0];
	Sc0a1b1.id2 = SonsIds2[0];

	double c1a1b1 = getC1(Sc1a1b1.id1,Sc1a1b1.id2);
	double c0a1b1 = getC0(Sc0a1b1.id1,Sc0a1b1.id2);


	// 2 possible solutions
	//case 1 -> the adj is not present >> break


	currentSolution = new AdjSolution(0,0,false) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0a1b1 );
	caseScore.push_back( c0a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 0 gain
	currentSolution = new AdjSolution(1,0,false) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	caseScore.push_back( c1a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;	
}







/*
Function specific to the case where the more ancient node has only one child

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool NodeId2First: true if nodeid2 is before nodeid1


Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0asynchronousOneChild(int NodeId1, int NodeId2, bool NodeId2First)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0asynchronousOneChild"<< endl;

	vector <AdjSolution> Vsolution;

	vector <int> SonsIds;

	if(!NodeId2First) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);

	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)

	if(NodeId2First) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);



	// 2 possible solutions
	//case 1 -> the adj is not present << nothing happens
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	caseScore.push_back( c1ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;

}





/*
Function specific to the case where one node has two children and the other only has one.

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0synchronousTwoAndOneChildren(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0synchronousTwoAndOneChildren"<< endl;


	bool iscoevent = false; // never a co-event
	vector <AdjSolution> Vsolution;

	vector <int> SonsIds1;
	vector <int> SonsIds2;

 	SonsIds1 = Rtree1.getSonsId(NodeId1);
	SonsIds2 = Rtree2.getSonsId(NodeId2);


	bool NodeId2TwoChildren = false;

	if(SonsIds1.size() == 1)
		NodeId2TwoChildren = true;


	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	Sc1ab1.id1 = SonsIds1[0];
	Sc1ab1.id2 = SonsIds2[0];

	Sc0ab1.id1 = SonsIds1[0];
	Sc0ab1.id2 = SonsIds2[0];

	if(NodeId2TwoChildren) // node 2 has 2 children
	{
		Sc1ab2.id1 = SonsIds1[0];
		Sc1ab2.id2 = SonsIds2[1];

		Sc0ab2.id1 = SonsIds1[0];
		Sc0ab2.id2 = SonsIds2[1];
	}
	else // node 1 has two children
	{
		Sc1ab2.id1 = SonsIds1[1];
		Sc1ab2.id2 = SonsIds2[0];

		Sc0ab2.id1 = SonsIds1[1];
		Sc0ab2.id2 = SonsIds2[0];
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);





	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency -> 0 break
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	Vsolution.push_back(*currentSolution);

	
	delete currentSolution;
	caseScore.clear();

	//case 2 -> child 1  kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 2 kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);




	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency -> 2 gain
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	return Vsolution;

}




/*
Function specific to the case where the more ancient node has two children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool NodeId2First: true if nodeid2 is before nodeid1

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0asynchronousTwoChildren(int NodeId1, int NodeId2, bool NodeId2First)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0asynchronousTwoChildren "<< NodeId2First << endl;

	vector <AdjSolution> Vsolution;

	vector <int> SonsIds;

	if(!NodeId2First) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);

	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	if(NodeId2First) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

		Sc1ab2.id1 = NodeId1;
		Sc1ab2.id2 = SonsIds[1];

		Sc0ab2.id1 = NodeId1;
		Sc0ab2.id2 = SonsIds[1];


	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


		Sc1ab2.id1 = SonsIds[1];
		Sc1ab2.id2 = NodeId2;

		Sc0ab2.id1 = SonsIds[1];
		Sc0ab2.id2 = NodeId2;

	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);





	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency -> 0 break

	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	

	Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();

	//case 2 -> child 1  kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 2 kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency -> 2 gain
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	
	
	return Vsolution;


}


/*
Function specific to the case where the nodes have two children
Two snodes can only be synchronous if they have the same event, which also means that they have the same number of children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0synchronousTwoChildren(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "SolutionC0synchronousTwoChildren"<< endl;

	vector<AdjSolution> Vsolution;
	bool iscoevent = false;

	AdjSolution * currentSolution;

	//we first get the ids of the children of NodeIdDup
	vector <int> SonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> SonsId2 = Rtree2.getSonsId(NodeId2);

	//variables to store potential node names
	vector <string> SonsName1(2);//lenght will be 2 max
	vector <string> SonsName2(2);//same

	int corr_a1b1 = 0;
	int corr_a2b1 = 0;
	int corr_a1b2 = 0;
	int corr_a2b2 = 0;
	
	if(LossAware)
	{ //this is equivalent to the hasFreeAdj() function in AdjMatrix.cpp.  but with the exception that it detect which case should have the free gain!
		
		if(Rtree1.hasNodeName(SonsId1[0]))
		{
			SonsName1[0] = Rtree1.getNodeName(SonsId1[0]);
		}
		else
		{
			stringstream ss;
			ss << Gfam1 << '|' << SonsId1[0];
			SonsName1[0] = ss.str();
			
		}
		
		if(Rtree1.hasNodeName(SonsId1[1]))
		{
			SonsName1[1] = Rtree1.getNodeName(SonsId1[1]);
		}
		else
		{
			stringstream ss;
			ss << Gfam1 << '|' << SonsId1[1];
			SonsName1[1] = ss.str();
			

		}
		
		if(Rtree2.hasNodeName(SonsId2[0]))
		{
			SonsName2[0] = Rtree2.getNodeName(SonsId2[0]);
		}
		else
		{
			stringstream ss;
			ss << Gfam2 << '|' << SonsId2[0];
			SonsName2[0] = ss.str();
		}
		
		if(Rtree2.hasNodeName(SonsId2[1]))
		{
			SonsName2[1] = Rtree2.getNodeName(SonsId2[1]);
		}
		else
		{
			stringstream ss;
			ss << Gfam2 << '|' << SonsId2[1];
			SonsName2[1] = ss.str();
		}
	
		for(int i=0; i < currFreeAdjacencies.first.size(); i++)
		{
			//cout << currFreeAdjacencies.first[i].first << " , " << currFreeAdjacencies.first[i].second << endl;
			if (SonsName1[0] == currFreeAdjacencies.first[i].first)
			{	if (SonsName2[0] == currFreeAdjacencies.first[i].second)
					corr_a1b1 = 1;
				if (SonsName2[1] == currFreeAdjacencies.first[i].second)
					corr_a1b2 = 1;
			}
			if ( SonsName1[1] == currFreeAdjacencies.first[i].first)
			{
				if (SonsName2[0] == currFreeAdjacencies.first[i].second)
					corr_a2b1 = 1;
				if (SonsName2[1] == currFreeAdjacencies.first[i].second)
					corr_a2b2 = 1;
			}
			if (SonsName2[0] == currFreeAdjacencies.first[i].first)
			{
				if (SonsName1[0] == currFreeAdjacencies.first[i].second)
					corr_a1b1 = 1;
				if (SonsName1[1] == currFreeAdjacencies.first[i].second)
					corr_a2b1 = 1;
			}
			if (SonsName2[1] == currFreeAdjacencies.first[i].first)
			{
				if (SonsName1[0] == currFreeAdjacencies.first[i].second)
					corr_a1b2 = 1;
				if (SonsName1[1] == currFreeAdjacencies.first[i].second)
					corr_a2b2 = 1;
			}
		}
	}

	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
	{
		for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
		{
			for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
			{
				for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
				{
					
					
					int c1a1b1_corr = 0; //  corr_c1axbx indicates whether or not it should be corrected ; c1axbx_corr gives the corrected value!
					int c1a1b2_corr = 0; //  corr_c1axbx indicates whether or not it should be corrected ; c1axbx_corr gives the corrected value!
					int c1a2b1_corr = 0; //  corr_c1axbx indicates whether or not it should be corrected ; c1axbx_corr gives the corrected value!
					int c1a2b2_corr = 0; //  corr_c1axbx indicates whether or not it should be corrected ; c1axbx_corr gives the corrected value!
					
					currentSolution = new AdjSolution();
					currentSolution->coevent = iscoevent;
									
					//a1b2
					currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[1]));
					if(!c1a1b2)
					{
						currentSolution->components.back().C1 = false;
						caseScore.push_back(getC0(SonsId1[0],SonsId2[1]));
					}
					else
					{
						caseScore.push_back(getC1(SonsId1[0],SonsId2[1]));

						c1a1b2_corr = c1a1b2 - corr_a1b2;
						
					}
					//a2b1
					currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[0]));
					if(!c1a2b1)
					{
						currentSolution->components.back().C1 = false;
						caseScore.push_back(getC0(SonsId1[1],SonsId2[0]));
					}
					else
					{
						caseScore.push_back(getC1(SonsId1[1],SonsId2[0]));						
						c1a2b1_corr = c1a2b1 - corr_a2b1;						
					}

					//a1b1
					currentSolution->components.push_back(AdjScore(true,SonsId1[0],SonsId2[0]));
					if(!c1a1b1)
					{
						currentSolution->components.back().C1 = false;
						caseScore.push_back(getC0(SonsId1[0],SonsId2[0]));
					}
					else
					{
						caseScore.push_back(getC1(SonsId1[0],SonsId2[0]));
						c1a1b1_corr = c1a1b1 - corr_a1b1;
					}

					//a2b2
					currentSolution->components.push_back(AdjScore(true,SonsId1[1],SonsId2[1]));
					if(!c1a2b2)
					{
						currentSolution->components.back().C1 = false;
						caseScore.push_back(getC0(SonsId1[1],SonsId2[1]));
					}
					else
					{
						caseScore.push_back(getC1(SonsId1[1],SonsId2[1]));
						c1a2b2_corr = c1a2b2 - corr_a2b2;
						
					}


					//determining the number of gains and breaks
	
					int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
					
					
					if(verbose)
						cout << "AdjMatrix::SolutionC0synchronousTwoChildren "<< c1a1b1 << " "<< c1a1b2 <<" "<< c1a2b1 <<" "<< c1a2b2<< " -> " << nbadj << endl;
					
					if(LossAware){

						nbadj = c1a1b1_corr + c1a1b2_corr + c1a2b1_corr + c1a2b2_corr;//if it's different, correction has been applied.

						if(verbose)
							cout << "LossAware AdjMatrix::SolutionC0synchronousTwoChildren "<< c1a1b1_corr << " "<< c1a1b2_corr << " "<< c1a2b1_corr << " "<< c1a2b2_corr<< " -> "<< nbadj <<endl;
					}
					currentSolution->NbGain = nbadj;
					currentSolution->NbBreak = 0;
					
					//for(unsigned i = 0 ; i < caseScore.size(); i++)
					//	cout << caseScore[i] << " ";
					//cout << endl;
					////getting the new score
					currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

					//cout << "AdjMatrix::SolutionC0synchronousTwoChildren computed score : " << currentSolution->score << endl;
					
					if(Vsolution.empty())
						Vsolution.push_back(*currentSolution);
					else
					{
						if(currentSolution->score <= Vsolution.front().score)
							Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
						else
							Vsolution.push_back(*currentSolution);
					}
					//clearing
					delete currentSolution;
					caseScore.clear();
				}
			}
		}
	}


	return Vsolution;
}




















/*
Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2


*//*
vector<AdjSolution> AdjMatrix::SolutionC1General(int NodeId1, int NodeId2)
{
	vector<AdjSolution> Vsolution;

	//First, determine if nodeid1 and nodeid2 are comparable
	bool incomparable = false;


	Node * n1 = Rtree1.getNode(NodeId1);
	Node * n2 = Rtree2.getNode(NodeId2);

	// both node are TS compatible
	int evt1 = Rtree1.getNodeEvent(NodeId1);
	int evt2 = Rtree2.getNodeEvent(NodeId2);
	

	// The choice is made that the TS compatibility is set by the reconciled trees rather than through the Adj algorithm.
	// This, for instance, allows a simple DeCo algorithm to perform on a full time sliced tree and respecting TS contraints

	if(!Rtree1.areTSCompatible(n1,n2))//if the node aren't TS compatible: their adjacency is impossible
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " are not time compatible. Setting values accordingly." << endl;
		}

		incomparable =true;
	}
	else if(!Rtree1.haveSameSpecies(n1, n2))
	{
		incomparable =true;	
		if(decoLTalgo) // in the case of Bout + Rec -> the species compatibility need not necessarily hold
		{
			if( ( ( Rtree1.isBout(evt1)) && ( Rtree1.isRec(evt2) ) ) || ( ( Rtree1.isBout(evt2)) && ( Rtree1.isRec(evt1) ) ) )
				incomparable = false;
		}
	}

	if(incomparable)
	{
		Vsolution =  SolutionC1DefaultImpossibleCase();
	}
	else
	{
		// the two node are comparable. 


		//cases where species is the same


		// special case --> one of the vent is a loss
		bool firstLoss = Rtree1.isLoss(evt1);
		bool secondLoss = Rtree2.isLoss(evt2);
	
		if( (firstLoss) || (secondLoss) )
		{
			if( (firstLoss) && (secondLoss) )
			{
				Vsolution = SolutionC1LossWithLoss(NodeId1,NodeId2);
			}
			else
			{
				Vsolution = SolutionC1LossWithOther(NodeId1,NodeId2);
			}
		}
		else
		{
			//there are no losses
			vector<AdjSolution> tmp;

			int nbSons = n1->getNumberOfSons();

			bool inTheDead = (Rtree1.getNodeSpecies(NodeId1) == -1);


			bool isRec = false; // only true if both event are receptions

			//case 0 -> is the event are identical, they may be simultaneous
			if(evt1 == evt2) // simultaneous
			{
				if(nbSons == 0) // leaf with leaf
				{
					tmp = SolutionC1ExtantWithExtant( NodeId1, NodeId2);
				}
				else if(nbSons == 1)
				{
					isRec = Rtree1.isRec(evt1); // checking if the events are receptions

					tmp = SolutionC1synchronousOneChild( NodeId1,  NodeId2,  inTheDead, iscoevent);
				}
				else if(nbSons == 2)
				{
					bool isSout = Rtree1.isSout(evt1);

					tmp = SolutionC1synchronousTwoChildren( NodeId1,  NodeId2,  inTheDead, isSout);
				}

				//adding to the already existing solutions
				if(Vsolution.size() == 0)
				{
					Vsolution = tmp;
				}
				else
				{
					for(unsigned i = 0; i < tmp.size();i++)
					{
						if(tmp[i].score <= Vsolution.front().score)
							Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
						else
							Vsolution.push_back(tmp[i]);
				
					}
				}

			}


			//case 1 : NodeId1 happened before NodeId2

			if(nbSons == 1)
			{
				Vsolution = SolutionC1asynchronousOneChild( NodeId1, NodeId2, false,  inTheDead, isRec);
			}
			else if(nbSons == 2)
			{
				Vsolution = SolutionC1asynchronousTwoChildren( NodeId1, NodeId2, false, inTheDead);
			}



			//case 2 : NodeId2 happened before NodeId1
			nbSons = n2->getNumberOfSons();
			inTheDead = (Rtree2.getNodeSpecies(NodeId2) == -1);

			if(nbSons == 1)
			{
				tmp = SolutionC1asynchronousOneChild( NodeId1, NodeId2, true,  inTheDead, isRec);
			}
			else if(nbSons == 2)
			{
				tmp = SolutionC1asynchronousTwoChildren( NodeId1, NodeId2, true, inTheDead);
			}

			//adding to the already existing solutions
			if(Vsolution.size() == 0)
			{
				Vsolution = tmp;
			}
			else
			{
				for(unsigned i = 0; i < tmp.size();i++)
				{
					if(tmp[i].score <= Vsolution.front().score)
						Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
					else
						Vsolution.push_back(tmp[i]);
			
				}
			}



		}


	}
		

	return Vsolution;
}

*/

/*
Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2


*//*
vector<AdjSolution> AdjMatrix::SolutionC0General(int NodeId1, int NodeId2)
{
	vector<AdjSolution> Vsolution;

	//First, determine if nodeid1 and nodeid2 are comparable
	bool incomparable = false;


	Node * n1 = Rtree1.getNode(NodeId1);
	Node * n2 = Rtree2.getNode(NodeId2);

	// both node are TS compatible
	int evt1 = Rtree1.getNodeEvent(NodeId1);
	int evt2 = Rtree2.getNodeEvent(NodeId2);
	

	// The choice is made that the TS compatibility is set by the reconciled trees rather than through the Adj algorithm.
	// This, for instance, allows a simple DeCo algorithm to perform on a full time sliced tree and respecting TS contraints

	if(!Rtree1.areTSCompatible(n1,n2))//if the node aren't TS compatible: their adjacency is impossible
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " are not time compatible. Setting values accordingly." << endl;
		}

		incomparable =true;
	}
	else if(!Rtree1.haveSameSpecies(n1, n2))
	{
		incomparable =true;	
		if(decoLTalgo) // in the case of Bout + Rec -> the species compatibility need not necessarily hold
		{
			if( ( ( Rtree1.isBout(evt1)) && ( Rtree1.isRec(evt2) ) ) || ( ( Rtree1.isBout(evt2)) && ( Rtree1.isRec(evt1) ) ) )
				incomparable = false;	
		}
	}

	if(incomparable)
	{
		Vsolution = SolutionC0DefaultImpossibleCaseNew(NodeId1,  NodeId2);
	}
	else
	{
		// the two node are comparable. 

		// special case --> one of the avent is a loss
		bool firstLoss = Rtree1.isLoss(evt1);
		bool secondLoss = Rtree2.isLoss(evt2);
	
		if( (firstLoss) || (secondLoss) )
		{
			if( (firstLoss) && (secondLoss) )
			{
				Vsolution = SolutionC0LossWithLoss(NodeId1,NodeId2);
			}
			else
			{
				Vsolution = SolutionC0LossWithOther(NodeId1,NodeId2);
			}
		}
		else if( ( Rtree1.isExtant(evt1)  ) && ( Rtree1.isExtant(evt2) ) ) // case where both event are comparable leave --> check adjacency
		{
			Vsolution = SolutionC0ExtantWithExtant( NodeId1, NodeId2);
		}
		else
		{
			Vsolution = SolutionC0DefaultImpossibleCaseNew( NodeId1, NodeId2);
		}


	}
		

	return Vsolution;
}


*/