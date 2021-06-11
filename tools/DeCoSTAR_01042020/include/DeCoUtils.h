#ifndef DECO_UTIL_H_
#define DECO_UTIL_H_
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

This file contains various functions used by DeCo

Created the: 02-03-2016
by: Wandrille Duchemin

Last modified the: 10-01-2018
by: Wandrille Duchemin

*/

#include "MyGeneTree.h"
#include "MySpeciesTree.h"
#include "MyCladesAndTripartitions.h"
#include "ReconciliationEvent.h"
#include "CladeReconciliation.h"
#include "ReconciledTree.h"
#include "GeneFamily.h"
#include "EquivalenceClass.h"
#include "EquivalenceClassFamily.h"
#include "AdjMatrix.h"
#include "DeCoOutputManager.h"
#include "CoEvent.h"

#include "CladesAndTripartitions.h"
#include "DTLMatrix.h"
#include "DTLGraph.h"


#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>




#define AFOREST vector <AdjTree *>
#define ECFsample vector < AFOREST *>
#define NECFsample vector < ECFsample * >


int countConnexComponents(vector< vector<int> > *adjacencyVector);

void ReadRecPhyLoXMLFile(  vector <GeneFamily *> * GeneFamilyList, string fileName , MySpeciesTree * Stree, map< string , int > & speciesNameToNodeId,  bool verbose, bool superverbose);

void AddToFile(string filename, string line);

map <string, int> makeLeafToGFMap(vector <GeneFamily *> * GeneFamilyList);
void fillLeafToSpMap(ReconciledTree * rtree, map <string,int> &LeafToSpMap);

void printNewick( MyGeneNode *node, string fileName );
void printNewick( MySpeciesNode *node, string fileName );

double computeTopoScore( vector <GeneFamily *> GeneFamilyList, double TopoWeight = 1);
double computeReconciliationScore(vector <GeneFamily *> GeneFamilyList, double ReconWeight = 1);
double computeAdjacenciesScore(vector <EquivalenceClass> * RefinedEquivalenceClasses, double AGainCost, double ABreakCost, double AdjWeight = 1);
double computeAdjacenciesScore(vector <EquivalenceClassFamily> * ECFams, double AGainCost, double ABreakCost, double AdjWeight = 1);
double computeCoEventScore(vector <CoEvent> Lcoevent, double DupCost, double LossCost, double HGTCost, double ReconWeight = 1 );
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClass> *RefinedEquivalenceClasses, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight = 1, double ReconWeight = 1, double AdjWeight = 1 );
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClassFamily> * ECFams, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight = 1, double ReconWeight = 1, double AdjWeight = 1 );


//not used anymore - could still be useful
vector < EquivalenceClass > CreateEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, bool DeCoLTrefine,bool Verbose, bool SuperVerbose);
vector < EquivalenceClass > CreateAllPairEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList,bool Verbose, bool SuperVerbose);

//new version
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, MySpeciesTree *speciesTree, map < string, int > speciesChrNb, map < int,vector <float> > &speciesC0C1, map<int,string> &species_id_name, map<int, map<string,int> > &speGeneAdjNb, float Break, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);

vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > &adjacencies, map <string,int> &LeafToGFMap, int nbGfam, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > &adjacencies, vector < pair <int,int> >  &adjacenciesOrientations, map <string,int> &LeafToGFMap, int nbGfam, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);


bool computeArtDeCoMaps(vector< pair <string,string > > &adjacencies,
                        vector< pair <int,int> > &adjacenciesOrientations,
                         map <string,int> &LeafToSpMap,
                         MySpeciesTree *speciesTree,
                         map < string, int > speciesChrNb,
                         map < int,vector <float> > &speciesC0C1,
                         map<int,string> &species_id_name, 
                         map<int, map<string,int> > &speGeneAdjNb, 
                         map<int, map<string, pair< int,int > > > &speGeneExtremitiesAdjNb,
                         float Break,bool Verbose, bool SuperVerbose, 
                         map < string, map <string , double> > * adjacencyScores, bool includeScoredAdjs, double InferableInfoLessCladeSize = 1);



void ComputeEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList ,
                                                      map < string, map <string , double> > &adjacencyScores, map<int,vector<float> > &speciesC0C1,
                                                      map<int, map<string,int> > &speGeneAdjNb, map< int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb, 
                                                      double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,
                                                      double WDupCost, double WLossCost, double WHgtCost, 
                                                      bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty,
                                                      bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
                                                      double adjScoreLogBase = 10000 , bool interactionMode = false);

void backtrackInPlaceEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);

vector < pair < pair<string, string> , double > > ComputeOneEquivalenceClassFamily( EquivalenceClassFamily * ECF, 
                                                                                          vector <GeneFamily *> * GeneFamilyList, 
                                                                                          map < string, map <string , double> > &adjacencyScores, 
                                                                                          map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb, map< int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb, 
                                                                                          double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,
                                                                                          double WDupCost, double WLossCost, double WHgtCost, 
                                                                                          bool Verbose, bool SuperVerbose, 
                                                                                          double boltzmannTemp , double absencePenalty ,
                                                                                          bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies ,
                                                                                          double adjScoreLogBase =10000, bool interactionMode = false);

vector < pair < pair<string, string> , double > > ComputeOneEquivalenceClassFamily( EquivalenceClassFamily * ECF, 
                                                                                          ReconciledTree * Rtree1, ReconciledTree * Rtree2, 
                                                                                          map < string, map <string , double> > &adjacencyScores, 
                                                                                          map<int,vector<float> > &speciesC0C1, map<int, map<string,int> > &speGeneAdjNb, map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb, 
                                                                                          double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,
                                                                                          double WDupCost, double WLossCost, double WHgtCost, 
                                                                                          bool Verbose, bool SuperVerbose, 
                                                                                          double boltzmannTemp , double absencePenalty ,
                                                                                          bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies , 
                                                                                          double adjScoreLogBase =10000, bool interactionMode = false);



void backtrackInPlaceOneEquivalenceClassFamily(EquivalenceClassFamily * ECF, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);
void backtrackInPlaceOneEquivalenceClassFamily(EquivalenceClassFamily * ECF, ReconciledTree * Rtree1, ReconciledTree * Rtree2 , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);


void backtrackOnetimeOneEquivalenceClassFamily( ECFsample * Sample, EquivalenceClassFamily * ECF , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage, int &overflowed, bool doNBBT=true);
void backtrackNtimesOneEquivalenceClassFamily(NECFsample * Samples, EquivalenceClassFamily * ECF, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);
void backtrackNtimesNequivalenceClassFamilies(vector <NECFsample * > * AllSamples, vector <EquivalenceClassFamily > * ECFams, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);

void ComputeAndBacktrackEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList ,
                                                      map < string, map <string , double> > &adjacencyScores, map<int,vector<float> > &speciesC0C1,
                                                      map<int, map<string,int> > &speGeneAdjNb, 
                                                      map<int, map<string,pair<int,int> > > &speGeneExtremitiesAdjNb,
                                                      double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,
                                                      double WDupCost, double WLossCost, double WHgtCost,
                                                      bool Verbose, bool SuperVerbose, bool galwaysGain, 
                                                      double gC1Advantage,
                                                      bool LossAware, pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies,
                                                      double boltzmannTemp = 1, double absencePenalty = -1, double adjScoreLogBase =10000 , bool interactionMode = false);

void PopulatesCoeventsFromAdjForest( vector <AdjTree * > * AForest, vector <CoEvent> * CoEventSet, int EclassId , int gfam1, int gfam2, ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool ignoreTime = false);


MySpeciesTree * getSpeciesTree(string speciesFile, bool dateAsBootstrap = false);
int processSpeciesTree( MySpeciesTree *speciesTree , bool dateAsBootstrap, bool dated, bool computeT, bool computeTD , double hgtCost, double lossCost ,bool verbose );

void checkGeneTree( int counter, MyGeneTree *geneTree);

vector<string> readGeneDistributionsFile(string ListGeneFile);
void readGeneDistributions(vector <GeneFamily *> * GeneFamilyList, MySpeciesTree *speciesTree ,vector<string> geneFiles, map< string , int > & speciesNameToNodeId, 
                                    bool ale, bool reconciled, char charSep ,bool verbose, bool superverbose, bool rooted=false);


vector< pair <string,string > > readAdjacencies(string AdjacencyFile);

vector< pair <string,string > > readAdjacencies(string AdjacencyFile , map < string, map <string , double> > &adjacencyScores, vector < pair <int,int> > &adjacenciesOrientations, bool verbose = true);


string WriteSpeciestree(MySpeciesTree * speciesTree, bool newick, string prefix,  string char_sep = "_");

vector <string> WriteGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, bool newick, bool hideLosses,string prefix);

string WriteOneGeneFamilyReconciledTree( GeneFamily * GFam, bool newick,bool hideLosses, string prefix, int index); // not used anymore


void InitReconciledTreeFile(string fileName, bool newick);
void AddOneGeneFamilyReconciledTreeToFile( GeneFamily * GFam, bool newick,bool hideLosses, string fileName, int index);
void FinishReconciledTreeFile(string fileName, bool newick);


void InitAdjacencyTreeFile(string fileName, bool newick);
void AddECFAForestToFile(string fileName, EquivalenceClassFamily * ECF, bool newick,bool hideLosses , double GainCost = 3,  double BreakCost= 1);
void AddECFNsampleToFile(string fileName, EquivalenceClassFamily * ECF, NECFsample * sample, bool newick,bool hideLosses, double GainCost = 3,  double BreakCost= 1);
void FinishAdjacencyTreeFile(string fileName, bool newick);

void AddECFsampleToFile(string fileName, EquivalenceClassFamily * ECF, ECFsample * sample, int sampleIndex, bool newick,bool hideLosses,  double GainCost = 3,  double BreakCost= 1 , bool Init= false, bool Finish = false);


void AddAForestToFileNewick(string fileName, AFOREST * forest,bool hideLosses , int fam1 , int fam2 , int ECindex, int sens1 = 0 , int sens2 = 0, int sample = -1);
void AddAForestToFileXML(string fileName, AFOREST * forest,bool hideLosses , int fam1 , int fam2 , int ECindex, int sens1 = 0 , int sens2 = 0, int sample = -1);


vector <string> WriteECFamTrees(EquivalenceClassFamily * ECF, bool newick,bool hideLosses, string prefix);
vector <string> WriteECFsample(EquivalenceClassFamily * ECF, ECFsample * sample,int g1, int g2, int sampleIndex, bool newick,bool hideLosses, string prefix);

void WriteAdjacencies(vector < EquivalenceClassFamily > * ECFams, string prefix);
void WriteAdjacenciesOneECF( EquivalenceClassFamily * ECF, string prefix, ReconciledTree * Rtree1, ReconciledTree * Rtree2);

map <string,int> storeSpeciesChrNumber(string chrNumberFile);

void WriteGeneOneGF( ReconciledTree * RT, int famIndex, string fileName, map<int, string> &SpIdToSpName );

void UpdateAdjMapsFromTree( AdjTree * Atree, int NbSample, 
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList,
            					ReconciledTree * Rtree1,
            					ReconciledTree * Rtree2
            					);

void UpdateAdjMapsFromForest( vector < AdjTree * > * AForest, int NbSample, 
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList,
            					ReconciledTree * Rtree1,
            					ReconciledTree * Rtree2
            					);


void AddAdjMapsToFile(string fileName , int sens1, int sens2,
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList
								);


map <int,string> buildRPOtoSpNameMap(MySpeciesTree * Stree);

string removePrefix(string s, string prefix, bool removeOneMore =false);
void writeSpeciesCorrespondanceFile(string fileName, MySpeciesTree * speciesTree);


void ReadAdjMaps(int actualIndex, map <int , map < string, vector < pair<string, int > > > > & AdjGraph,
                        map <string, map <string, int> > & AdjIndexMap,
                        vector < double > & AdjInScoreList,
                        vector < double > & AdjOutScoreList,
                        vector < int > & AdjSpeList
                        );
                        
//map < int ,vector < pair <string, string> > > ReadFreeAdjacenciesFile(string filename);
void MakeFreeAdjacencies(map <int , map < string, vector < pair<string, int > > > > & AdjGraph,
                                    vector <GeneFamily *> * GeneFamilyList,
                                    map < int ,pair < vector < pair <string, string> >, bool > > & FreeAdjacencies);
pair < int, int> FindGroupLoss(int OrNode, int OrFam, 
                                                int spe_loss, int father_spe,
                                                int Node, int Fam, 
                                                map <int , map < string, vector < pair<string, int > > > > & AdjGraph,
                                                vector <GeneFamily *> * GeneFamilyList);

void saveECFsampleToTmpFile(EquivalenceClassFamily * ECF, ECFsample * sample,
                                           int sampleIndex, bool newick,bool hideLosses,  
                                           double GainCost ,  double BreakCost , bool Init, bool Finish);
                                           
void copyTmpToFile(string filename, EquivalenceClassFamily * ECF);

bool doesItExist( pair < vector < pair <string, string> >, bool > FamFreeAdjacencies, pair < string, string > potentialFreeAdjacency);

bool getNameToNodeIdMap( MySpeciesTree *speciesTree , map< string , int > & nameToNodeId ,bool verbose );

#endif