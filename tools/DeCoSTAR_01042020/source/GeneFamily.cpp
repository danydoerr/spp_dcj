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

This file contains a class for Gene Family as a Clade and Tripartition instance associated to a reconciled tree

Created the: 18-11-2015
by: Wandrille Duchemin

Last modified the: 31-08-2017
by: Wandrille Duchemin

*/

#include "GeneFamily.h"



/*
Takes:
 - geneTree (MyGeneTree *) : tree to check
 - errStr (string &): string to get the nature of the error
 - bool rooted (default = false): consider the tree as already rooted or not (changes the handling of potential trifurcation)

Returns:
	(bool): true if the tree is valid for CladesAndTripartions: more than leaves and no identical leaf names
*/
bool GeneFamily::CheckGeneTree(MyGeneTree * geneTree, string &errStr, bool rooted)
{

	string dupName;
	if( !geneTree->uniqueLeaves( dupName ) )
	{
		errStr = "Duplicated leaf name : " + dupName;
		return false;
	}



	MyGeneNode * rootNode = geneTree->getRootNode();



	if(rootNode->getNumberOfSons() > 2)
	{ // multifurcating roots are ok if we aren't in rooted mode
		if(rooted)
		{
			errStr = "A gene tree has a multifurcating root while rooted=1. Please resolve all multifurcation (possible using ecceTERA) before using DeCoSTAR.";
			return false;
		}
		else if(rootNode->getNumberOfSons() == 3)
			geneTree->rootTree(); // we random root it as this has no importance.
		else
		{
			errStr = "A gene tree has a root with more than 3 children (counted " + static_cast<ostringstream*>( &(ostringstream() << rootNode->getNumberOfSons()) )->str() +" children) . Please resolve all multifurcation (possible using ecceTERA) before using DeCoSTAR.";
			return false;
		}
	}



	vector <Node * > NodeToDo = rootNode->getSons();

	while(NodeToDo.size() > 0)
	{
		Node * current  = NodeToDo.back();
		NodeToDo.pop_back();


		size_t nbSons = current->getNumberOfSons();

		if(nbSons>2)
		{
			errStr = "A gene tree has a multifurcating node. Please resolve all multifurcation (possible using ecceTERA) before using DeCoSTAR.";
			return false;
		}

		for(size_t i = 0 ; i < nbSons ; i++)
			NodeToDo.push_back(current->getSon(i));

	}




	return true;
}


/*
Draw a reconciled tree from the DTLMatrix

Takes:
 - speciesTree (MySpeciesTree *)
 - random (bool) : if the reconciliation must be chosen randomly or not
 - dated (int) : 0 if the species tree is not subdivided, 2 if it is.
*/
void GeneFamily::drawReconciliation(MySpeciesTree * speciesTree, bool random, int dated)
{

    clock_t begin;
    clock_t end;
    double elapsed_secs;

    //if(VERBOSE)
    //{
    //    begin = clock();
    //}

	//cout << "drawing reconciliation random? "<< random <<endl;
	//basic checks

	if(!CanChangeReconciliation)
	{
		if(VERBOSE)
			cout << "Tried to reset the UnrootedTree while it is not possible to change the tree (CanChangeTree false). Does nothing." << endl;
		return;	
	}

	resetReconciliation();

	

	if( CCPDistrib->mClades.getLeafCount() == 1 )
	{

		ReconciledTree RT;  //RecTree();

		int idX = speciesTree->getNode( CCPDistrib->mClades.getSpeciesName(1) )->getId() ;
		idX =  speciesTree->getRPO(idX);
		RT.makeSingleLeafTree( CCPDistrib->mClades.getLeafName( 1 ) , idX );

		if(dated == 0)//if the sp tree is not subdivided, then the tree shouldn't
			RT.setToNonTimeSliced();

		setReconciliation(RT);

		ReconciliationSet = true;
		return;
	}

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for RECreset:" << elapsed_secs << endl;
    //    begin = clock();
    //}


	///// getting the clades and tripartition of the DTL matrix

	//backtracking to get the new ROOTED tree
	MyGeneNode * node = ReconciliationMatrix->backtrack( false );
	// create a tree to delete node structre easily
	MyGeneTree rootedTree = MyGeneTree(*node); 


	//we now have a rooted tree. we launch the procedure again to get the reeconciled tree... this is fastidious

	//transforming the lone tree into a CandT instance
	CladesAndTripartitions *cladesAndTripartitions;
	cladesAndTripartitions = new CladesAndTripartitions( CharSep, rootedTree );

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for RECtopobacktrack:" << elapsed_secs << endl;
    //    begin = clock();
    //}




	////////////////////////// RECONCILIATION GRAPH /////////////////////////////////////////
	//next: building the graph and transferring to the reconciled gene tree
	//Graph construction
	DTLGraph graph = ReconciliationMatrix->constructGraph(SUPERVERBOSE);
	graph.countReconciliationNumberAndCheck( 0, 
											false,//true, ///< Remove non-canonical vertices
											SUPERVERBOSE, ///< print debugging information
											false ); // necessary for weighting?
	//int nbSol = graph.getNumberSolutions(true);
	//cout << "got "<<nbSol << " solutons in the reconciliation graph"<< endl;
	//int seed = rand();
	//string path = "GraphAllRec" +  static_cast<ostringstream*>( &(ostringstream() << nbSol) )->str() + "_"  +  static_cast<ostringstream*>( &(ostringstream() << seed) )->str() +".txt" ;
	//graph.printAllReconciliations(path, false,
    //                              false);

	
    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for GRAPHbuilding:" << elapsed_secs << endl;
    //    begin = clock();
    //}


	vector< vector<DTLGraph::MyGraph::Vertex> > reconciliation;
	graph.getScoredReconciliation( 3, rootedTree.getNumberOfNodes(),
									reconciliation,
									random );//getting a reconciliation (problem = symmetric. random = false)


    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for TERA computations:" << elapsed_secs << endl;
    //    begin = clock();
    //}



	////////////////////////// CLADERECONCILIATION /////////////////////////////////////////
	//We now want to get from the reconciliation to a list of CladeReconciliation
	//preparing for the reconciliation traversal
	vector<int> cladeToPOrd = cladesAndTripartitions->getPostOrderMapping();
	int cladeCount = cladeToPOrd.size();
	vector<int> pOrdToClade( cladeCount );

	for( int j=0; j<cladeCount; j++ ) 
		pOrdToClade[cladeToPOrd[j]] = j;

	int rootClade = cladesAndTripartitions->mClades.getRootClade();
	map <int , int > idUToParent;
	map <int , CladeReconciliation*> idUToCladeReconciliation;
	vector <CladeReconciliation> CladeRecs;
	ReconciliationEvent previousEvent = ReconciliationEvent(); //creating empty event for the root clade
	for( int pOrd = cladeCount-1; pOrd >=0; pOrd-- ) 
	{
		int idU = pOrdToClade[pOrd];
		if(idU != rootClade)
		{
			if( idUToParent.count(idU) == 0) 
			{
				if(VERBOSE)
					cout << idU << " not set" << endl;
				throw Exception("GeneFamily::drawReconciliation : firstEvents not set" );
			}
			previousEvent = CladeRecs[idUToParent[idU]].Events.back(); //getting the last event of the parent clade
		}
		CladeReconciliation CR = CladeReconciliation(cladesAndTripartitions, &graph, speciesTree, idU, &reconciliation, previousEvent, SUPERVERBOSE);
		CladeRecs.push_back(CR);
		if(!CR.isLeaf())
		{//setting parenthood of the two children
			idUToParent[CR.idUl] = CladeRecs.size() -1 ;
			idUToParent[CR.idUr] = CladeRecs.size() -1 ;
		}
		if(SUPERVERBOSE)
			CR.printMe();					
	}


    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for CladeRec building:" << elapsed_secs << endl;
    //    begin = clock();
    //}
//

	////////////////////////// FROM CLADERECONCILIATIONS TO RECONCILED TREE /////////////////////////////////////////

	//We now want to create the reconciled tree
	RecTree = ReconciledTree(&CladeRecs,SUPERVERBOSE);

	if(dated == 0)//if the sp tree is not subdivided, then the tree shouldn't
		RecTree.setToNonTimeSliced();


	//cout << "***********************"<< endl;
	//RecTree.printMe();
	//cout << "***********************"<< endl;

	ReconciliationSet = true;


    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\t\ttime elapsed for TREEbuilding:" << elapsed_secs << endl;
    //    begin = clock();
    //}

	//VERBOSE=false;

	delete cladesAndTripartitions;


}

void GeneFamily::resetReconciliation()
{
	if(!CanChangeReconciliation)
	{
		if(VERBOSE)
			cout << "Tried to reset the Reconciliation while it is not possible to change it (CanChangeReconciliation false). Does nothing." << endl;
		return;		
	}

	if(ReconciliationSet)//only do something if the reconciliation is set. Otherwise do nothing
	{
		if(SUPERVERBOSE)
			cout << "effectively resetting the reconciliation" << endl;
		
		ReconciliationSet = false;
		//RecTree = ReconciledTree();
	}
}

void GeneFamily::resetUnrootedTree()
{
	//Also reset theDTLMatrix and ReconciledTree

	if(!CanChangeTree)
	{
		if(VERBOSE)
			cout << "Tried to reset the UnrootedTree while it is not possible to change the tree (CanChangeTree false). Does nothing." << endl;
		return;
	}


	if(TreeSet)//nothing to do if the tree is not set
	{
		//resetting the reconciliation
		resetReconciliation();
	
		//resetting the DTLMatrix
		if(DTLMatSet)
		{
			DTLMatSet = false;
			delete ReconciliationMatrix;
		}
		//ReconciliationMatrix = NULL;

		if(SUPERVERBOSE)
			cout << "effectively resetting the tree" << endl;

		TreeSet = false;
		//UnrootedTree = MyGeneTree();
	}

}







/*
Constructor of the class

Takes:
 - geneTrees (vector<MyGeneTree*> &): vector of gene tree pointers
 - charsep (string): separator between species name and gene name in the leaf names
 - verbose (bool)
 - superverbose (bool)
*/
GeneFamily::GeneFamily(vector<MyGeneTree*> &geneTrees, char charsep,bool verbose, bool superverbose)
{
	//booleans about the different states of the Gene Family
	CanChangeTree = true;//implies CanChangeReconciliation
	CanChangeReconciliation = true;
	TreeSet = false;
	ReconciliationSet = false;//implies TreeSet
	DTLMatSet = false;

	VERBOSE = verbose;
	SUPERVERBOSE = superverbose;

	CharSep = charsep;


	//Checking that the given gene trees are OK
	bool gTree_OK = true;
	string errString = "";
	for(unsigned i = 0; i < geneTrees.size(); i++)
	{		
		MyGeneTree * geneTree = geneTrees[i];

		//cout << "ploup "<< i << endl;

		gTree_OK = CheckGeneTree(geneTree, errString, false);

		if(!gTree_OK)
			break;
	}
	

	if(!gTree_OK) //this means that the gene trees in that distribution all have less than 3 leaves or duplicate leaf names
	{
		cerr << "GeneFamily::GeneFamily : unable to load gene tree list. error: " << errString << endl;
		throw Exception("GeneFamily::GeneFamily : unable to load gene tree list. error: " + errString);
	}
	
	


	bool overflow = false;
	errString = "";

	// putting them in a MyCladesAndTripartition instance which computes CCPs and bipartitions
	CCPDistrib = new MyCladesAndTripartitions( CharSep, geneTrees, VERBOSE, overflow,errString, false ); //the false at the end is because we don't allow polytomy


	

	if( overflow ) 
	{
		delete CCPDistrib;
		throw Exception("GeneFamily::GeneFamily : Too many possible gene trees due to polytomies -> aborting" );
	}
	if( errString != "" )
	{
		delete CCPDistrib;
		throw Exception("GeneFamily::GeneFamily : error while creating the CCP distribution : " + errString );
	}


	//setting the leaf name vector
	leafNames = geneTrees[0]->getLeavesNames();

	if(VERBOSE)
		cout << "Gene family created. " << geneTrees.size() << " trees." << endl;

	//MaxLkh = CCPDistrib->getMaxLikelihood();

}



/*
Constructor of the class: reads the gene tree distribution from an ALE file

Takes:
 - aleFileName (string): name of an ale file
 - charsep (string): separator between species name and gene name in the leaf names
 - verbose (bool)
 - superverbose (bool)
*/
GeneFamily::GeneFamily(string aleFileName, char charsep,bool verbose, bool superverbose)
{
	//booleans about the different states of the Gene Family
	CanChangeTree = true;//implies CanChangeReconciliation
	CanChangeReconciliation = true;
	TreeSet = false;
	ReconciliationSet = false;//implies TreeSet
	DTLMatSet = false;

	VERBOSE = verbose;
	SUPERVERBOSE = superverbose;

	CharSep = charsep;

	string errString = "";
	CCPDistrib = new MyCladesAndTripartitions( CharSep, aleFileName, errString, VERBOSE );

	if( errString != "" )
	{
		delete CCPDistrib;
		throw Exception("GeneFamily::GeneFamily : error while creating the CCP distribution : " + errString );
	}

	//add something to index the leaf names
	int i = 1;
	while(CCPDistrib->mClades.isLeaf(i))
	{
		leafNames.push_back(CCPDistrib->mClades.getLeafName(i));
		i++;
	}


	//MaxLkh = CCPDistrib->getMaxLikelihood();


	if(VERBOSE)
		cout << "Gene family created. from the ale file : " << aleFileName << endl;


}




/*
Constructor of the class
Gets directly an unrooted tree
don't allow a tree topology change according to a CCP distribution

Takes:
 - UnrootedgeneTree (MyGeneTree)
 - charsep (string): separator between species name and gene name in the leaf names
 - verbose (bool)
 - superverbose (bool)
*/
GeneFamily::GeneFamily(MyGeneTree UnrootedgeneTree, char charsep, bool verbose, bool superverbose, bool rooted)
{
	//booleans about the different states of the Gene Family
	CanChangeTree = false;//implies CanChangeReconciliation
	CanChangeReconciliation = true;
	TreeSet = true;
	ReconciliationSet = false;//implies TreeSet
	DTLMatSet = false;

	VERBOSE = verbose;
	SUPERVERBOSE = superverbose;

	CharSep = charsep;

	//checking the tree
	string errString = "";
	if(!CheckGeneTree(&UnrootedgeneTree, errString, rooted))
	{
		cerr << "GeneFamily creation : error whith the tree : " << errString << endl;
		exit(1);
	}

	UnrootedTree = UnrootedgeneTree;

	leafNames = UnrootedTree.getLeavesNames();

	if(VERBOSE)
		cout << "Gene family created with 1 gene tree." << endl;

}

/*
Constructor of the class.
Gets directly the reconciled tree
don't allow a tree topology change according to a CCP distribution
don't allow a reconciliation change using TERA's algorithm

Takes:
 - RTree (ReconciledTree)
 - verbose (bool)
 - superverbose (bool)
*/
GeneFamily::GeneFamily(ReconciledTree RTree, bool verbose, bool superverbose)
{
	//booleans about the different states of the Gene Family
	CanChangeTree = false;//implies CanChangeReconciliation
	CanChangeReconciliation = false;
	TreeSet = true;
	ReconciliationSet = true;//implies TreeSet
	DTLMatSet = false;


	VERBOSE = verbose;
	SUPERVERBOSE = superverbose;

	//here, plenty of pre-processing and checks might come
	
	RecTree = RTree;

	leafNames = RecTree.getRealLeavesNames();//ok

	MyGeneTree UnrootedgeneTree = RTree.getMyGeneTreeTopology();

	//RTree.printMe();


	setUnrootedTree( UnrootedgeneTree); // setting the new topology (also computes Lkh)

	if(VERBOSE)
		cout << "Gene Family created with a reconciled tree." << endl;

	//Technically, it would be possible to extract the gene tree structure from the reconciled tree in order to have the possibility to come back and change the reconciliation through a DTLMatrix Instance. However this is not what is wished here.
}

/*
draw an unrooted tree from the CCPDistrib; at random if random is true; otherwise uses best split

Takes:
 - random (bool): will the tree be drawn at random or not 

*/
void GeneFamily::makeUnrootedTree(bool random)
{
	if(!CanChangeTree)
	{
		if(VERBOSE)
			cout << "Tried to reset the UnrootedTree while it is not possible to change the tree (CanChangeTree false). Does nothing." << endl;
		return;
	}

	//starts by resetting the tree
	resetUnrootedTree();

	MyGeneTree * newTree;

	if(random)
		newTree = CCPDistrib->getRandomTree();
	else
		newTree = CCPDistrib->getMaxTree();

	//newTree->printNewick( "test.nwk" , true);

	UnrootedTree = *newTree;

	TreeSet = true;

	//setting the likelihood of the new tree
	TreeLikelihood = CCPDistrib->getTreeLikelihood(UnrootedTree.getRootNode(),0.0); // default proba is 0 here... maybe use something else...

	if(VERBOSE)
		cout << "Made unrooted tree." << endl;
}


/*
replace the current unrooted gene tree by the one given in argument

Takes:
 - UnrootedgeneTree (MyGeneTree): a gene tree

*/
void GeneFamily::setUnrootedTree(MyGeneTree UnrootedgeneTree)
{
	//if(!CanChangeTree)
	//{
	//	if(VERBOSE)
	//		cout << "Tried to reset the UnrootedTree while it is not possible to change the tree (CanChangeTree false). Does nothing." << endl;
	//	return;
	//}

	//starts by resetting the tree
	if(CanChangeTree)
		resetUnrootedTree();

	//checks the given tree: very basic check: >2 leaves and no duplicated leaves
	string errStr;
	if(!CheckGeneTree( &UnrootedgeneTree, errStr,false))
	{
		throw Exception("GeneFamily::setUnrootedTree : error whith the tree : " + errStr );	
	}

	if(CanChangeTree)
	{
		if(!CCPDistrib->isCompatible(UnrootedgeneTree.getRootNode()))//checking that the tree has correct leaves
			throw Exception("GeneFamily::setUnrootedTree : error whith the tree : incompatible with the CCP distribution." );	
	}

	UnrootedTree = UnrootedgeneTree;

	TreeSet = true;

	//setting the likelihood of the new tree
	if(CanChangeTree)
		TreeLikelihood = CCPDistrib->getTreeLikelihood(UnrootedTree.getRootNode(),0.000001); // default proba is 0 here... maybe use something else...
}



/*
Compute and set the reconciliation of the gene family using TERA's algorithm

Takes:
 - speciesTree (MySpeciesTree *)
 - computeT (bool) : transfers are computed if true
 - computeTL (bool): transfer loss events are computed if true
 - DupCost (double): cost of an individual duplication
 - HgtCost (double):cost of an individual transfer
 - LossCost (double): cost of an individual loss
 - maxTS (int): maximum time slice in the species tree
 - SplitWeight (double): weight of the splits (so the gene tree topology) when coparing with reconciliation score. This is used in the case where we try all amalgamation
 - SubdividedSpTree (bool): true if the species tree is subdivided
 - tryAllAmalgamation (bool): true if all amalgamation have to be tested (using a CladesAndTripartition instance)
 - bool random : wether the reconciliation must be ddrawn randomly from the reconciliation graph or not
*/
void GeneFamily::makeReconciliation(MySpeciesTree * speciesTree, bool computeT, bool computeTL, double DupCost, double HgtCost, double LossCost, int maxTS , double SplitWeight, bool SubdividedSpTree, bool tryAllAmalgamation, bool random, bool rooted)
{

    //clock_t begin;
    //clock_t end;
    //double elapsed_secs;

    //if(VERBOSE)
    //{
    //    begin = clock();
    //}


	//basic checks

	if(!CanChangeReconciliation)
	{
		if(VERBOSE)
			cout << "Tried to reset the reconciliation while it is not possible to change the tree (CanChangeReconciliation false). Does nothing." << endl;
		return;	
	}

	resetReconciliation();

	if((!tryAllAmalgamation)&&(!TreeSet)) // trying to make the reconciliation while the unrooted tree is not set -> impossible if not trying all amalgamation
	{
		throw Exception("GeneFamily::makeReconciliation: trying to make the reconciliation while the unrooted tree is not set.");
	}

	if((tryAllAmalgamation)&&(!CanChangeTree))// trying to make the reconciliation while the CandT distribution is not set.
	{
		if(!TreeSet) // impossible to work around that
			throw Exception("GeneFamily::makeReconciliation: trying to make the reconciliation while neither the CCP distribution nor the unrooted gene tree are set.");

		//if the unrooted tree is set, don't try all amalgamation
		tryAllAmalgamation = false;
	}
	if((tryAllAmalgamation)&&( CCPDistrib->mClades.getLeafCount() <= 2))
	{
		tryAllAmalgamation=false;
		makeUnrootedTree(false);

	}

	if( CCPDistrib->mClades.getLeafCount() == 1 )
	{ // special procedure when there is only 1 leaf.
		setRecScore(0);
		int dated =0;
		if(SubdividedSpTree)
			dated = 2;
		drawReconciliation( speciesTree,random , dated );
		//RecTree.printMe();
		return;
	} 

	//setting up some basic arguments
	if(!computeT)
		computeTL = false;

	bool gWGDCosts = false;//we don't use the WGD computation
	bool fixedCost = true;// reconciliation costs are fixed
	bool useBestSplits = true;//always true if we want to backtrack after

	int gMaxIterations = 1;//only iteration as we don't update costs
	bool gUpdateWeightOnly = false;//no cost update

	//creating the CandT that we will use in the reconciliation
	CladesAndTripartitions *cladesAndTripartitions;

	

	if(tryAllAmalgamation)
	{
		cladesAndTripartitions = CCPDistrib; // use the already existing MyCladesAndTripartition instance
	}
	else
	{
		//use the unrooted tree
		if( (!rooted) && ( UnrootedTree.getNumberOfLeaves() > 2 ) )
		{
			vector< MyGeneTree* >  Vtree;
			Vtree.push_back(&UnrootedTree);
	
			bool overflow = false;
			string errStr ="";
	
			cladesAndTripartitions = new CladesAndTripartitions( CharSep, Vtree, SUPERVERBOSE, overflow, errStr, false);
	
			if( overflow ) 
			{
				delete cladesAndTripartitions;
				throw Exception("GeneFamily::makeReconciliation : Too many possible gene trees due to polytomies -> aborting" );
			}
			if( errStr != "" )
			{
				delete cladesAndTripartitions;
				throw Exception("GeneFamily::makeReconciliation : error while creating the CCP distribution : " + errStr );
			}
			if(SUPERVERBOSE)
				UnrootedTree.printMe();///Wdebug
		}
		else
		{
			cladesAndTripartitions = new CladesAndTripartitions( CharSep, UnrootedTree);// Vtree, SUPERVERBOSE, overflow, errStr, false);		
		}
	}


    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\ttime elapsed for RECpreppin:" << elapsed_secs << endl;
    //    begin = clock();
    //}

	if(VERBOSE)
		cout << "1st DTLMatrix computation" << endl;

	//Creating the DTLMatrix	
	ReconciliationMatrix = new DTLMatrix( speciesTree, cladesAndTripartitions, gWGDCosts, fixedCost, computeT, computeTL, DupCost, HgtCost, LossCost, maxTS, SplitWeight, useBestSplits, 0 );
	DTLMatSet = true;
	//Computing the matrix
	int dated =0;
	if(SubdividedSpTree)
		dated = 2;

	ReconciliationMatrix->calculateMatrix( VERBOSE, gMaxIterations, gUpdateWeightOnly,dated);

	if(SUPERVERBOSE)
		cout << "1st DTLMatrix calculated" << endl;

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\ttime elapsed for REC1:" << elapsed_secs << endl;
    //    begin = clock();
    //}

			

	//if(VERBOSE)
	//{
	//	cout << "Cost of a most parsimonious reconciliation : " << getRecScore() << endl;
	//}

	//backtracking to get the new ROOTED tree
	MyGeneNode * node = ReconciliationMatrix->backtrack( false );
	// create a tree to delete node structre easily
	MyGeneTree rootedTree = MyGeneTree(*node); 

	if(tryAllAmalgamation)
	{
		UnrootedTree = rootedTree;
		UnrootedTree.unroot();
	}

	////////////////////////// SECOND RECONCILIATION /////////////////////////////////////////
	//delete ReconciliationMatrix;
	if(!tryAllAmalgamation)
		delete cladesAndTripartitions;
	else //if we did all amalgamation, we now want the likelihood of the tree
		setTreeLikelihood( CCPDistrib->getTreeLikelihood(rootedTree.getRootNode(),0.0) );

	//we now have a rooted tree. we launch the procedure again to get the reeconciled tree... this is fastidious

	//transforming the lone tree into a CandT instance
	cladesAndTripartitions = new CladesAndTripartitions( CharSep, rootedTree );

	if(VERBOSE)
		cout << "2nd DTLMatrix computation" << endl;

	//Creating the DTLMatrix
	ReconciliationMatrix = new DTLMatrix( speciesTree,
										cladesAndTripartitions,
										gWGDCosts,
										fixedCost,
										computeT,
										computeTL,
										DupCost,
										HgtCost,
										LossCost,
										maxTS,
										SplitWeight,
										useBestSplits,
                                        0 );

	//Computing the matrix
	ReconciliationMatrix->calculateMatrix( VERBOSE, gMaxIterations, gUpdateWeightOnly,dated);


	//Getting the best score
	setRecScore( ReconciliationMatrix->getBestCost());

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\ttime elapsed for REC2:" << elapsed_secs << endl;
    //    begin = clock();
    //}



	drawReconciliation(speciesTree,random , dated);

    //if(VERBOSE)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "\ttime elapsed for BACKTRACK:" << elapsed_secs << endl;
    //    begin = clock();
    //}
	delete cladesAndTripartitions;

}

/*
replace the current Reconciliation by the one given as argument

Takes:
 - RTree (ReconciledTree): a reconciled gene tree

CAUTION: no particular tests are performed to ensure that a tree of the correct family was given
*/
void GeneFamily::setReconciliation(ReconciledTree RTree)
{
	if(!CanChangeReconciliation)
	{
		if(VERBOSE)
			cout << "Tried to reset the Reconciliation while it is not possible to change it (CanChangeReconciliation false). Does nothing." << endl;
		return;
	}

	//resetReconciliation();

	MyGeneTree UnrootedgeneTree = RTree.getMyGeneTreeTopology();

	setUnrootedTree( UnrootedgeneTree); // setting the new topology (also computes Lkh)

	RecTree = RTree;

	ReconciliationSet = true;

	RecScore = -1; // default value because we don't know the RTree score...
}


/*
just re-draw a tree in the already computed DTLMatrix. Won't work if no DTLMatrix have been computed


Takes:
 - speciesTree (MySpeciesTree *)
 - SubdividedSpTree (bool): true if the species tree is subdivided
 - BTS (bool): true if the reconcilied gene tree timeslicestatus has to be set to Bounded
*/
void GeneFamily::changeReconciliation(MySpeciesTree * speciesTree, bool SubdividedSpTree, bool BTS)
{
	int dated = 0;
	if(SubdividedSpTree)
		dated = 2;

	drawReconciliation(speciesTree,true,dated);

	if(BTS)
		setReconciledTreeTimeSlicesToBTS(speciesTree); // set the tree to a bounded time slice one (BTS)

}

///Getters and setters


int GeneFamily::getReconciledTreeTimeSliceStatus()
{
	// 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Asked for Reconciled tree TimeSliceStatus while the reconciled is not done (ReconciliationSet is false). Returns dummy value -1." << endl;
		return -1;
	}

	return RecTree.getTimeSliceStatus();
}

void GeneFamily::setReconciledTreeTimeSlicesToBTS(MySpeciesTree * spTree)
{
	// 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Tried to set Reconciled tree TimeSliceStatus while the reconciled is not done (ReconciliationSet is false). Does nothing." << endl;
		return;
	}
	int currentTSS = getReconciledTreeTimeSliceStatus();

	if(currentTSS == 2)
		return; // do nothing if the tss is already the good one

	/*if(currentTSS == 0)
	{
		if(VERBOSE)
			cout << "Tried to set the Reconciled Tree Timeslice status to BTS while there is no existing timeslice to draw from. Draws TSs from the species tree." << endl;
	}*/

	RecTree.computeBoundedTimeSlices(spTree);
}


void GeneFamily::setReconciledTreeTimeSlicesToTS(MySpeciesTree * spTree)
{
	// 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Tried to set Reconciled tree TimeSliceStatus while the reconciled is not done (ReconciliationSet is false). Does nothing." << endl;
		return;
	}
	int currentTSS = getReconciledTreeTimeSliceStatus();

/*	if(currentTSS == 1)
		return; // do nothing if the tss is already the good one
*/
	if(currentTSS == 0)
	{
		if(VERBOSE)
			cout << "Tried to set the Reconciled Tree Timeslice status to TS while there is no existing timeslice to draw from. Draws TSs from the species tree. ( maybe use bounded.TS=1 )" << endl;
		setReconciledTreeTimeSlicesToBTS(spTree);
	}

	RecTree.SubdivideTree();
}


void GeneFamily::setReconciledTreeTimeSlicesToNoTS()
{
	// 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Tried to set Reconciled tree TimeSliceStatus while the reconciled is not done (ReconciliationSet is false). Does nothing." << endl;
		return;
	}
	int currentTSS = getReconciledTreeTimeSliceStatus();

	if(currentTSS == 0)
		return; // do nothing if the tss is already the good one

	RecTree.setToNonTimeSliced();
}


ReconciledTree * GeneFamily::getRecTree()
{
	ReconciledTree * p = &RecTree;
	return p;
}

ReconciledTree GeneFamily::getRecTreeCopy()
{
	return RecTree;
}


double GeneFamily::getTreeLikelihood()
{
	return TreeLikelihood;
}

double GeneFamily::getNormalizedTreeLikelihood()
{
	//double LMax = CCPDistrib->getMaxLikelihood();

	return TreeLikelihood / MaxLkh ; 
}

void GeneFamily::setTreeLikelihood(double lkh)
{
	TreeLikelihood = lkh;
}

double GeneFamily::getRecScore()
{
	return RecScore;
}

void GeneFamily::setRecScore(double rscore)
{
	RecScore = rscore;
}

// Counts events in the RecTree and uses the costs to compute the score
void GeneFamily::setRecScore(double DupCost, double HgtCost, double LossCost)
{
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Tried to set the reconciliation score while the reconciled tree is not set. Does nothing." << endl;
		return;
	}

	RecScore = 0;

	RecScore = RecTree.getNumberOfDuplication() * DupCost;
	RecScore += RecTree.getNumberOfLoss() * LossCost;
	RecScore += RecTree.getNumberOfTransfer() * HgtCost;

}

void GeneFamily::printRecTree()
{
	if(!ReconciliationSet)
	{
		if(VERBOSE)
			cout << "Tried to print the reconciled tree while it is not set. Does nothing." << endl;
		return;
	}
	cout << RecTree.NewickString() << endl;
}


bool GeneFamily::hasLeaf( string lName )
{
	bool has = false;
	for(int i =0; i< leafNames.size() ; i++)
	{
		if(leafNames[i] == lName)
		{
			has = true;
			break;
		}
	}
	return has;
}


// free CCPDistrib and ReconciliationMatrix from the memory
void GeneFamily::dumpStuff()
{
	//cout << "dumping"<< endl;
	if(CanChangeTree)
	{
		delete CCPDistrib;
		CanChangeTree = false;
	}

	if(DTLMatSet)
	{
		//cout << "dumping more "<< endl;
		delete ReconciliationMatrix;
		DTLMatSet = false;
	}

}



//functions to set up a biased clades and tripartition instance

/*
For the case where ids in the adj tree does not correspond to ids in the reconciled tree

Takes:
	- AdjTree * Atree : an adjacency tree
	- int GfamIndex : wether we are interested in the ids of gfam1 or gfam2 in the adjacency tree
Returns;
	vector < vector <string> > : vector of clades represented by list of strings
*/
vector < vector <string> > GeneFamily::getLeafListToKeepFromAtreeAlone(AdjTree * Atree, int GfamIndex)
{
	vector < vector <string> > LeafListList;

	Atree->getC1LeafList( GfamIndex, LeafListList );

	return LeafListList;
}




/*
For the case where ids in the adj tree does not correspond to ids in the reconciled tree

Takes:
	- AdjTree * Atree : an adjacency tree
	- int GfamIndex : wether we are interested in the ids of gfam1 or gfam2 in the adjacency tree
Returns;
	vector < vector <string> > : vector of clades represented by list of strings

*/
vector < vector <string> > GeneFamily::getLeafListToKeepFromAtree(AdjTree * Atree, int GfamIndex)
{
	vector <int> nodeToKeepIds = Atree->getC1NodeIds(GfamIndex);

	cout << "GeneFamily::getLeafListToKeepFromAtree with gfamindex "<< GfamIndex<<". Found "<< nodeToKeepIds.size()<< " nodes to keep"<< endl;

	map <int, vector <string> > cladeToLeaves = RecTree.cladeToLeafListAssociation();

	map <int, bool > cladeToLeavesTokeep; 
	for(unsigned i = 0 ; i < nodeToKeepIds.size() ; i++) // we filte the clades we are intrested in
	{
		int cladeId = RecTree.getNodeCladeNum(nodeToKeepIds[i]);

		if(cladeToLeavesTokeep.find(cladeId) == cladeToLeavesTokeep.end()) // new clade -> add it
		{
			if(cladeToLeaves[cladeId].size() > 1) // we aren't interested in leaf clades as they, by dfinition, always exists
			{
				cladeToLeavesTokeep[cladeId] = true;
			}
		}
	}

	vector < vector <string> > leafList;
	for (map<int,bool>::iterator it=cladeToLeavesTokeep.begin(); it!=cladeToLeavesTokeep.end(); ++it)
	{
		leafList.push_back(cladeToLeaves[it->first]);
	}

	return leafList;
}


/*
Takes:
	- vector < AdjTree * > * AdjForest : a vector of adjacency trees
	- vector<int> GfamIndex : wether the current gene family is gfam 1 (0) or gfam2 (1) in each AdjTree
	- vector <bool> isCompatible : wether the current gene family and each AdjTree are compatible (ie. ids are identical)

Returns:
	vector < vector <string> > : vector of clades represented by list of strings; WITHOUT ANY GUARANTY THAT TEH LIST IS NOT REDUNDANT
*/
vector < vector <string> > GeneFamily::getLeafListToKeepFromAForest( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible)
{
	vector < vector <string> > leafList;

	for(unsigned i = 0 ; i < AdjForest->size() ; i++)
	{
		vector < vector <string> > tmp;
		if(isCompatible[i])
			tmp = getLeafListToKeepFromAtree(AdjForest->at(i), GfamIndex[i]);
		else
			tmp = getLeafListToKeepFromAtreeAlone(AdjForest->at(i), GfamIndex[i]);

		cout << "adjTree "<< i << " yielded "<< tmp.size() << " clades"<<endl;

		for(unsigned j = 0 ; j < tmp.size(); j++)
			leafList.push_back(tmp[j]);
	}

	return leafList;
}


/*
Takes:
	- vector < AdjTree * > * AdjForest : a vector of adjacency trees
	- vector<int> GfamIndex : wether the current gene family is gfam 1 (0) or gfam2 (1) in each AdjTree
	- vector <bool> isCompatible : wether the current gene family and each AdjTree are compatible (ie. ids are identical)

Returns;
	MyCladesAndTripartitions * : biased CCP distribution
*/
MyCladesAndTripartitions * GeneFamily::setBias( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible)
{

	vector < vector <string> > leafListToRestrict = getLeafListToKeepFromAForest( AdjForest , GfamIndex, isCompatible);

	//cout << "GeneFamily::setBias :" << leafListToRestrict.size() << " restricting clades" <<endl;

	MyCladesAndTripartitions * BiasedCCPDistrib = new MyCladesAndTripartitions(CCPDistrib);
	
	//cout << "ORIGINAL"<<endl;
	//CCPDistrib->printMe();
	//cout << "COPY"<<endl;
	//BiasedCCPDistrib->printMe();


	BiasedCCPDistrib->RestrictToClade( leafListToRestrict );

	return BiasedCCPDistrib;
}
