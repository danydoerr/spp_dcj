/*
Created by: Wandrille Duchemin

Last modified the: 10-01-2018
by: Wandrille Duchemin



@file DeCoSTAR.cpp
@author Wandrille Duchemin
@author Yoann Anselmetti
@author Celine Scornavacca
@author Edwin Jacox
@version 1
@date 20/07/2016

@section LICENCE
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

@section DESCRIPTION

*/


#include "DeCoUtils.h"
#include <unistd.h>

#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>


#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>



string version = "1.2";
string date = "10/01/18";
//1.2 all children adj transmission rules
//1.1.1 NBBT option to avoid memory hog
//1.1 loss aware
//1 -> paper release
//0.9.2 -> connex component computation
//0.9.1 -> extremities specific artdeco score and new option, no connex component computation

//name, type, default, description
map<string,bool> gBoolParams;
map<string,int> gIntParams;
map<string,double> gDoubleParams;
map<string,string> gStringParams; // string, char, and path 



const int gParameterCount = 45;
string const gParameters[gParameterCount][4] = {
// files
 {"parameter.file", "path", "none", "a file with input parameters" },
 {"species.file", "path", "required", "species tree file (newick). NB: default behavior wants it to be ultrametric if transfers are used (use dated.species.tree=0 to circumvent)" },
 {"gene.distribution.file", "path", "required", "gene distribution files file (one file name per line)"},
 {"adjacencies.file","path","required","adjacencies file (one adjacency per line; leafnames separated by a space)"},

// basic options
 {"with.transfer","bool","true","allow transfers the reconciliation and adjacency histories reconstruction."},
 {"dated.species.tree","bool","true","the species tree is ultrametric and dates will be used to subdivide the trees in time slices."},

// format of input
 {"char.sep", "char", "_", "char separating gene names in gene tree files. One character only."},
 {"ale","bool","false","gene tree distribution are ALE files"},
 {"already.reconciled","bool","false","gene tree distribution are reconciled gene trees in recPhyloXML format. Will skip the reconciliation phase"},
 {"verbose","int","1","show progress and timing.\n\t0: nothing is reported short of error.\n\t1: basic report (default).\n\t2: various information about reconciliation, adjacency matrix and backtracking\n\t3: maximal amount of information"},
 //{"superverbose","bool","false","show lots of text about progress and timing"},

//reconciliation options
 {"dupli.cost", "double", "2", "cost of a duplication"},
 {"HGT.cost", "double", "3", "cost of an HGT"}, 
 {"loss.cost", "double", "1", "cost of a loss"},
 {"try.all.amalgamation","bool","true","try all possible amalgamation when reconciling gene trees. Otherwise the best possible tree is used"},
 {"rooted","bool","false","specify that the root of the given trees must be kept"},

//DeCo option
 {"AGain.cost","double","2","cost of an adjacency gain"},
 {"ABreak.cost","double","1","cost of an adjacency break"},
 {"C1.Advantage","double","0.5","between 0 and 1. Probability to choose C1 (presence of adjacency) over C0 (absence of adjacency) in case of a score tie at the root of an equivalence class"},

//Loss Aware options
 {"Loss.aware","bool","false","computer adjacency matrices giving a free gain to a list of specific adjacencies based on losses in their neighbors."},
 {"Loss.iteration","int","2","number of iterations of DeCoSTAR if lossAware mode is on. Minimum is 2."},

//Scaffolding options
 {"scaffolding.mode","bool","false","use scaffolding algorithm to improve extant genomes scaffolding/assembly"},
 {"chromosome.file","path","none","used with the scaffolding.mode option. A file containing the number of chromosome in each each species (one species per line, each line comprised name of the species followed by the number of chromosome, separated by a tabulation)"},
 {"adjacency.score.log.base","double","10000","used in the case where the adjacency file also contains a score between 0 and 1. Base of the logarithm applied to this score (must be above 1)."},
 {"scaffolding.propagation.index","double","1","advanced scaffolding option. must be > 0 ; max size of the clade where we can still infer new adjacencies even though they don't possess any (ie. the only adjacency is on the outgroup of this clade)"},
 {"scaffold.includes.scored.adjs","bool","false","used in the case where the adjacency file also contains a score between 0 and 1 AND scaffolding.mode is true. If true, include the adjacencies with a score < 1 in the computation of the number of contigs."},

//Boltzmann
 {"use.boltzmann","bool","false","use Boltzmann sampling for the adjacencies history computation"},
 {"boltzmann.temperature","double","0.1","Temperature to use in the Boltzmann sampling (if used)"},
 {"nb.sample","int","1","number of samples to get from the adjacency matrix."},

//output options
 {"write.newick","bool","false","use newick format rather than phyloXML-like format"},
 {"hide.losses.newick","bool","false","if true, losses and the branches leading to them wil be removed from the newick string"},
 {"write.adjacencies","bool","true","write the adjacencies inferred in ancestral and extant species"},
 {"write.genes","bool","false","write the genes inferred in ancestral and extant species"},
 {"write.adjacency.trees","bool","false","write the inferred adjacency trees."},

 {"output.dir", "string", "none", "directory for printed files"}, 
 {"output.prefix", "string", "none", "a prefix to prepend to all output files."},


//advanced
 {"all.pair.equivalence.class","bool","false","compute adjacency histories for all pair of gene families (even if they share no adjacencies)."},
 {"bounded.TS","bool","false","use bounded time slices in adjacency history computations"},
 {"always.AGain","bool","true","always put an Ajdacency Gain at the top of an equivalence class tree"},
 {"absence.penalty","double","-1","if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization"},
 {"substract.reco.to.adj","bool","false","if set to 1, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation."},
 {"count.number.backtracks","bool","false","if set to 1, adjacency scenarios will be sampled uniformly across most parsimonious histories. potential memory overload if the reconciled gene trees are big (hundreds of leaves each)."},


 {"Topology.weight","double","1","weight of the topology in the global score" },
 {"Reconciliation.weight","double","1","weight of the reconciliation in the global score" },
 {"Adjacency.weight","double","1","weight of the adjacency in the global score" },

 {"all.children.adj.transmission.rules","bool","false","whether or not adjacencies should be expected to be transmitted to all children in the same species or not"}
 //{"dump.and.load.limit","int","100000","Advanced parameter. Number of gene families above which DeCo will load them from a temp save file in order to save RAM. Minimum of 2."},

};


void help() {
    cout << "This progam computes and outputs the most parsimonious"
            " history of a set of given adjacencies between extant genes,"
            " according to their gene trees and their reconciliation with"
            " a species tree." << endl;
    cout << "version: " << version << endl;
    cout << "__________________________________________________________"
            "_____________________" << endl;
    cout << "___________________________________INPUT__________________"
            "_____________________" << endl;
    cout << "Ex:    ./bin/DeCoSTAR parameter.file=tests/testDeCo_boltzmann.param.txt"<<endl;
    cout << "format: parameter.name | [type,default] <description>" << endl;

    for( size_t i=0; i<gParameterCount; i++ ) {
        if( gParameters[i][3] != "" )
        {
            cout << gParameters[i][0] ;
            for (size_t j = 0; j <  50 - gParameters[i][0].size() ; j++)
                cout << " ";
            cout << "| [" << gParameters[i][1] << ","
                 << gParameters[i][2] << "] " << gParameters[i][3] << endl;
        }
    }

}

string gPathPrefix = ""; // path and prefix for all output files

// change to boost::filesystem?
static bool do_mkdir( const char *path ) {
    struct stat st;
    if( stat(path, &st) != 0 ) {
        if( mkdir( path, 0777 ) != 0 && errno != EEXIST)
            return false;
    } else if( !S_ISDIR(st.st_mode) ) { // check if it is a directory
        return false;
    }

    return true;
}
/**
 ** mkpath - ensure all directories in path exist
 ** Algorithm takes the pessimistic view and works top-down to ensure
 ** each directory in path exists, rather than optimistically creating
 ** the last element and working backwards.
 **/
int mkpath( string path ) // mode_t mode )
{
    char *copypath = strdup(path.c_str());
    char *sp;
    int status = 0;
    char *pp = copypath;
    while( status == 0 && (sp = strchr(pp, '/')) != 0 ) {
        if (sp != pp) { 
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath ); //, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    //if (status == 0)
        status = do_mkdir(path.c_str()); //, mode);
    free(copypath);
    return status;
}

/*struct MatchPathSeparator {
        bool operator()( char ch ) const { return ch == '/'; }
};*/


/**
 * read input parameters
 */
void readParameters( 
        map<string,string> &params ) ///< input paramsters
{
    map<string,string>::iterator iterStr;
    bool problem = false;

    // check for unknown parameters
    map<string,bool> allParams;
    for( size_t i=0; i<gParameterCount; i++ ) 
        allParams[gParameters[i][0]] = true;
    for( iterStr = params.begin(); iterStr != params.end(); ++iterStr ) {
        map<string,bool>::iterator iterBool = allParams.find( iterStr->first );
        if( iterBool == allParams.end() ) {
            problem = true;
            cerr << "Unknown parameter: " << iterStr->first << endl;
        }
    }




    // check for a parameter file
    iterStr = params.find( "parameter.file" );  
    if( iterStr != params.end() ) {
        map<string,string> fileParams;
        bpp::AttributesTools::getAttributesMapFromFile( iterStr->second,
                                                       fileParams, "=" );
        // add file parameters if not already there
        map<string,string>::iterator fIter;
        for( fIter=fileParams.begin(); fIter!=fileParams.end(); ++fIter) 
        {
            iterStr = params.find( fIter->first );
            if( iterStr == params.end() ) 
                params[fIter->first] = fIter->second;
        }

    }
        

    for( size_t i=0; i<gParameterCount; i++ ) {
        string name = gParameters[i][0];
        string type = gParameters[i][1];
        string value;

        // get value 
        iterStr = params.find( name );  
        if( iterStr == params.end() ) {
            if( gParameters[i][1] == "required" ) {
                problem = true;
                cerr << name << " must be specified" << endl;
            }
            value = gParameters[i][2]; // default
        } else {
            value = iterStr->second;
        }

        if( type == "bool" ) {
            if ((value == "true") 
                || (value == "TRUE")
                || (value == "t")
                || (value == "T")
                || (value == "yes")
                || (value == "YES")
                || (value == "y")
                || (value == "Y")
                || (value == "1") ) 
            {
                gBoolParams[name] = true;
            } 
            else if ((value == "false") 
                || (value == "FALSE")
                || (value == "f")
                || (value == "F")
                || (value == "no")
                || (value == "NO")
                || (value == "n")
                || (value == "N")
                || (value == "0") ) 
            {
                gBoolParams[name] = false;
            }
            else {
                cerr << "Invalid boolean value (" << value << ") for " 
                     << name << endl;
                problem = true;
            }
        } else if( type == "string" ) {
            gStringParams[name] = value;
        } else if( type == "char" ) {
            if( value.size() != 1 ) {
                cerr << "Invalid char value (" << value << ") for " 
                     << name << endl;
                problem = true;
            } else {
                gStringParams[name] = value;
            }

        } else if( type == "int" ) {
            gIntParams[name] = bpp::TextTools::toInt( value );
        } else if( type == "double" ) {
            gDoubleParams[name] = bpp::TextTools::toDouble( value );
        } else if( type == "path" ) {
            if( value == "required" ) {
                cerr << name << " is a required parameter." << endl;
                problem = true;
            } else if( value != "none" ) {
                if( !bpp::FileTools::fileExists( value ) ) {
                    cerr << value << " does not exist." << endl;
                    problem = true;
                } else {
                    gStringParams[name] = value;
                }
            } else 
                gStringParams[name] = "none";
        } else {
            cerr << "Found type <" << type << "> for " << name << endl;
            throw bpp::Exception( "readParameters: unknown type" );
        }
    }

    if( problem )
        exit(1);

    // PARAMETER CHECKS -> forcing some dependent parameters

    if(gIntParams.find("verbose")->second < 0)
    {
        cerr << "verbose can't be <0 : forcing verbose=0" << endl;
        gIntParams["verbose"]=0;
    }
    else if(gIntParams.find("verbose")->second > 3)
    {
        cerr << "verbose can't be >3 : forcing verbose=3" << endl;
        gIntParams["verbose"]=3;
    }


    if( !gBoolParams.find("with.transfer")->second) //no transfer
    {
        if(gBoolParams.find("bounded.TS")->second)
        {
            cerr << "No transfer : forcing bounded.TS=false" << endl;
            gBoolParams["bounded.TS"] = false;
        }
        if(gBoolParams.find("dated.species.tree")->second)
        {
            cerr << "No transfer : forcing dated.species.tree=false" << endl;
            gBoolParams["dated.species.tree"] = false;  
        }
    }
    else // with transfers.
    {
        if(!gBoolParams.find("dated.species.tree")->second)
            if(gBoolParams.find("bounded.TS")->second)
            {
                cerr << "Non dated species tree : forcing bounded.TS=false" << endl;
                gBoolParams["bounded.TS"] = false;
            }
    }

    if( !gBoolParams.find("use.boltzmann")->second) //no boltzmann sampling
    {
        if(gDoubleParams.find("boltzmann.temperature")->second != 0.1 )
        {
            cerr << "No Boltzmann sampling: forcing boltzmann.temperature=0.1" << endl;
            gDoubleParams["boltzmann.temperature"] = 0.1;
        }
        //if(gIntParams.find("nb.sample")->second != 1)
        //{
        //    cout << "No Boltzmann sampling: forcing nb.sample=1" << endl;
        //    gIntParams["nb.sample"] = 1;
        //}
    }       
    else // Boltzmann sampling
    {
    //  if( gBoolParams.find("with.transfer")->second) //no transfer
    //  {
    //      cout << "Boltzmann sampling AND transfer might not yield consistent results." << endl;
    //  }

        if( gBoolParams.find("substract.reco.to.adj")->second )
        {
            cerr << "Boltzmann sampling: forcing substract.reco.to.adj=false" << endl;
            gBoolParams["substract.reco.to.adj"] = false;   
        }

        //if(gBoolParams.find("write.adjacencies")->second )
        //{
        //    cerr << "Boltzmann sampling: forcing write.adjacencies=false" << endl;
        //    gBoolParams["write.adjacencies"] = false;      
        //}
    }

    if( gBoolParams.find("ale")->second && gBoolParams.find("already.reconciled")->second)
    {
        cerr << "already.reconciled : forcing ale=false" << endl;
        gBoolParams["ale"] = false; 
    }

    if( gBoolParams.find("rooted")->second)
    {
        if( gBoolParams.find("try.all.amalgamation")->second )
        {
            cerr << "rooted=1 forces try.all.amalgamation=0" << endl;
            gBoolParams["try.all.amalgamation"] = false;
        }
    }  



    if( gBoolParams.find("substract.reco.to.adj")->second  && gDoubleParams.find("Adjacency.weight")->second == 0) //substract reco score but adj weight = 0 -> division by 0!!
    {
        cerr << "substract.reco.to.adj=true and Adjacency.weight=0. Will set Adjacency.weight to 0.000001 in Adjacency matrix computation to avoid division by 0." << endl;
        // modification on the fly
    }

    //ensuring 0 <= C1.advantage <= 1
    if(gDoubleParams.find("C1.Advantage")->second < 0 )
    {
        cerr << "Negative value of C1.Advantage: forcing C1.Advantage=0" << endl;
        gDoubleParams["C1.Advantage"] = 0;
    }
    else if (gDoubleParams.find("C1.Advantage")->second > 1 )
    {
        cerr << "value of C1.Advantage > 1: forcing C1.Advantage=1" << endl;
        gDoubleParams["C1.Advantage"] = 1;
    }

    if(gBoolParams.find("hide.losses.newick")->second)
        if(!gBoolParams.find("write.newick")->second)
        {
            cerr << "hide.losses.newick = 1 : forcing write.newick=1" << endl;
            gBoolParams["write.newick"] = true;
        }

    // create path prefix for output files
    if( gStringParams.find("output.prefix")->second != "none" ) 
        gPathPrefix = gStringParams.find("output.prefix")->second;
    string outputDir = gStringParams.find("output.dir")->second;
    if( outputDir != "none" ) {
        if( !mkpath( outputDir ) ) {
            cerr << "ERROR: " << outputDir << " is not a directory." 
                 << endl;
            exit(1);
        }
        gPathPrefix = outputDir + "/" + gPathPrefix;
    }

    //if( gIntParams.find("dump.and.load.limit")->second < 2)
    //{
    //    cout << "forcing dump.and.load.limit at the minimum value of 2."<< endl;
    //    gIntParams["dump.and.load.limit"] = 2;
    //}

    if( gDoubleParams.find("adjacency.score.log.base")->second <= 1)
    {
        cerr << "adjacency.score.log.base is a log base : must be >1. Forcing at default value of 10 000." << endl;
        gDoubleParams["adjacency.score.log.base"]= 10000;
    }

    if( gDoubleParams.find("scaffolding.propagation.index")->second <= 0)
    {
        cerr << "scaffolding.propagation.index must be above 0. Forcing at default value of 1."<<endl;
        gDoubleParams["scaffolding.propagation.index"]= 1;
    }

    if( gBoolParams.find("Loss.aware")->second )
    {
        if(gIntParams.find("Loss.iteration")->second <= 1)
        {
            cerr << "Loss.iteration must be above 1. Forcing at default value of 2."<<endl;
            gIntParams["Loss.iteration"]= 2;
        }
    }

}


void printParameters()
{

    for (map<string,bool>::iterator it=gBoolParams.begin(); it!=gBoolParams.end(); ++it)
    {
        cout << it->first ;
        for( unsigned j = 0; j < 50 - it->first.size(); j++)
            cout <<".";
        cout << it->second << endl;
    }

    for (map<string,int>::iterator it=gIntParams.begin(); it!=gIntParams.end(); ++it)
    {
        cout << it->first ;
        for( unsigned j = 0; j < 50 - it->first.size(); j++)
            cout <<".";
        cout << it->second << endl;
    }

    for (map<string,double>::iterator it=gDoubleParams.begin(); it!=gDoubleParams.end(); ++it)
    {
        cout << it->first ;
        for( unsigned j = 0; j < 50 - it->first.size() ; j++)
            cout <<".";
        cout << it->second << endl;
    }

    for (map<string,string>::iterator it=gStringParams.begin(); it!=gStringParams.end(); ++it)
    {
        cout << it->first ;
        for( unsigned j = 0; j < 50 - it->first.size() ; j++)
            cout <<".";
        cout << it->second << endl;
    }


}

/////////////////////////////////////////////////
// Main 
/////////////////////////////////////////////////

int main(int args, char ** argv)
{
    if(args == 1)
    {
        help();
        exit(0);
    }
    
    try 
    {
        // fill global parameter variables
        map<string, string> params = 
            bpp::AttributesTools::parseOptions(args, argv);
        readParameters( params );

        bool randomtree = false;    
        bool dateAsBootstrap = false; // a bit harsh -> make option?

        int VerboseLevel = gIntParams.find("verbose")->second;
        bool verbose = (VerboseLevel > 1);
        bool superverbose = (VerboseLevel > 2);

        if( VerboseLevel > 0 ) {
            cout << "******************************************************************" << endl;
            cout << "*                 DeCoSTAR, version "<< version;
            for(unsigned i = 0 ; i < 29 - version.length(); i++)
                cout << " ";
            cout << "*" << endl;
            cout << "*  Authors: W. Duchemin, Y. Anselmetti, C. Scornavacca, E. Jacox *"<< endl;
            cout << "*                 Created     02/03/16                           *" << endl;
            cout << "*                 Last Modif  "<< date <<"                           *" << endl;
            cout << "******************************************************************" << endl;
            cout << endl;
        }

        if( VerboseLevel > 1 ) 
            printParameters();

        //if( gBoolParams.find("verbose")->second ) 
        //    bpp::ApplicationTools::startTimer();



        clock_t begin;
        clock_t end;
        double elapsed_secs;

        //string WrittenFileFileName = gPathPrefix;
        //if( WrittenFileFileName[ WrittenFileFileName.size() - 1 ] != '/')
        //    WrittenFileFileName += ".";
        //WrittenFileFileName += "writtenFiles.txt";

        //ofstream ofs;
        //ofs.open(WrittenFileFileName.c_str(),ofstream::out );
        //ofs.close(); // we close it directly
        //AddToFile(WrittenFileFileName, WrittenFileFileName ); //forst, we write our own name


        ///////////////////////////////////////
        /// Reading data
        ///////////////////////////////////////

        MySpeciesTree * speciesTree = getSpeciesTree(gStringParams.find("species.file")->second, dateAsBootstrap);
        int maxTS = processSpeciesTree( speciesTree ,  dateAsBootstrap, gBoolParams.find("dated.species.tree")->second , gBoolParams.find("with.transfer")->second, gBoolParams.find("with.transfer")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second, verbose );

        //map< string , int > * nameToNodeId = new map< string , int > ();
        map< string , int > speciesNameToNodeId = map< string , int >();
        if( gBoolParams.find("already.reconciled")->second )
        {
            bool nbMatching = getNameToNodeIdMap( speciesTree , speciesNameToNodeId , verbose );
            if( not nbMatching )
            {
                cerr << "Reconciled tree input :" << endl;
                cerr << "The number of unique names in the species tree does not match the number of nodes in the species tree." << endl;
                cerr << "Ensure that each node in the species tree has a unique name (as bootstrap)!" << endl;
                exit(1);
            }

        }

        vector<string> geneFiles = readGeneDistributionsFile( gStringParams.find("gene.distribution.file")->second );
        
        vector <GeneFamily *> * GeneFamilyList = new vector <GeneFamily *>;

        readGeneDistributions( GeneFamilyList, speciesTree , geneFiles, speciesNameToNodeId,
                                gBoolParams.find("ale")->second ,  
                                gBoolParams.find("already.reconciled")->second, 
                                gStringParams.find("char.sep")->second[0] ,
                                verbose, superverbose , 
                                gBoolParams.find("rooted")->second );



        if(VerboseLevel > 0)
            cout << "read " << GeneFamilyList->size() << " gene families." << endl;



        //vector< pair <string,string > > adjacencies = readAdjacencies( gStringParams.find("adjacencies.file")->second );
        map < string, map <string , double> > adjacencyScores;
        vector < pair <int,int> > adjacenciesOrientations;
        vector< pair <string,string > > adjacencies = readAdjacencies(gStringParams.find("adjacencies.file")->second , adjacencyScores, adjacenciesOrientations, verbose);


        if(VerboseLevel > 0)
            cout << "read " << adjacencies.size() << " adjacencies." << endl;


        map <string,int> LeafToGFMap = makeLeafToGFMap(GeneFamilyList);
        
        map <string,int> LeafToSpMap;

        int nbGfam = GeneFamilyList->size();


        vector <string> ReconciledTreesFileNames;


        ////////////////////////////////////////////////////
        /// Reporting the species tree and the reconciled trees
        ////////////////////////////////////////////////////

        // this writes the species tree and add the filename to the file of written filenames
        //AddToFile(WrittenFileFileName, WriteSpeciestree( speciesTree,  gBoolParams.find("write.newick")->second , gPathPrefix) ); 
        WriteSpeciestree( speciesTree,  gBoolParams.find("write.newick")->second , gPathPrefix, gStringParams.find("char.sep")->second) ;
        if( (gBoolParams.find("write.genes")->second) || (gBoolParams.find("write.adjacencies")->second) )
        {
            //we write a species correspondance file
            string fileName = gPathPrefix;
            if(( fileName[ fileName.size() - 1 ] != '/')&&(fileName.size() > 0))
            {
                fileName += ".";
            }
            fileName += "species.txt";

            writeSpeciesCorrespondanceFile(fileName, speciesTree);

        }

        //write all rec trees in one file
        string recFileName = gPathPrefix;
        string adjTreeFileName = gPathPrefix;
        
        if(( recFileName[ recFileName.size() - 1 ] != '/')&&(recFileName.size() > 0))
        {
            recFileName += ".";
            adjTreeFileName += ".";
        }

        recFileName += "reconciliations";
        adjTreeFileName += "adjacencyTrees";

        if( gBoolParams.find("write.newick")->second ) 
        {
            recFileName += ".newick" ;
            adjTreeFileName += ".newick" ;
        }
        else
        {
            recFileName += ".xml";
            adjTreeFileName += ".xml";
        }

        InitReconciledTreeFile( recFileName, gBoolParams.find("write.newick")->second );

        if(VerboseLevel>0)
        {
            if(!gBoolParams.find("already.reconciled")->second)
                cout << "Gene family reconciliation"<< endl;
            else
                cout << "Gene family already reconciled"<< endl;
        }


        map <int, string> RPOtoSpNameMap; // for genes and adj writing
        //if( (gBoolParams.find("write.genes")->second) || (gBoolParams.find("write.adjacencies")->second) )
        //{
        //    RPOtoSpNameMap = buildRPOtoSpNameMap(speciesTree); // hop
        //}That was a false good idea to replace leaves RPO by leaves names in the genes output file


        string genesFileName = "";

        if(gBoolParams.find("write.genes")->second)
        { // touching the adjacencies.txt file in order to reset it 
            genesFileName = gPathPrefix;

            if(( genesFileName[ genesFileName.size() - 1 ] != '/')&&(genesFileName.size() > 0))
                genesFileName += ".";
            genesFileName += "genes.txt";
        
            ofstream ofs;
            ofs.open(genesFileName.c_str(),ofstream::out );
            ofs.close(); // we close it directly

            //AddToFile(WrittenFileFileName, filename );
        }



        ///////////////////////////////////////
        /// Computing MPR for each gene family
        ///////////////////////////////////////
        if(VerboseLevel > 0)
            begin = clock();

        //declaring a bunch of default args

        double gWeight = 0;   ///< split weight
        if(gBoolParams.find("try.all.amalgamation")->second)
            gWeight = gDoubleParams.find("Topology.weight")->second;

        

        for(int i = 0; i < GeneFamilyList->size(); i++)
        {
            GeneFamily * Gfamily = GeneFamilyList->at(i);

            if(!gBoolParams.find("already.reconciled")->second)//only reconcile if it is not already done
            {
                if((!gBoolParams.find("try.all.amalgamation")->second) && (!gBoolParams.find("rooted")->second)) // if not all amalgamation are tried, we have to set the gene tree
                    Gfamily->makeUnrootedTree(randomtree);
    
                Gfamily->makeReconciliation(speciesTree,
                                            gBoolParams.find("with.transfer")->second, gBoolParams.find("with.transfer")->second,  // transfer and transferLoss
                                            gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second , // costs
                                            maxTS, gWeight, gBoolParams.find("dated.species.tree")->second, //here, dated.species.tree means that species tree is subdivided
                                            gBoolParams.find("try.all.amalgamation")->second,false, gBoolParams.find("rooted")->second);

            }
            else
                Gfamily->setRecScore(gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second );

            if(gBoolParams.find("bounded.TS")->second)
                Gfamily->setReconciledTreeTimeSlicesToBTS(speciesTree); // set the tree to a bounded time slice one (BTS)
            else if(gBoolParams.find("dated.species.tree")->second)
                Gfamily->setReconciledTreeTimeSlicesToTS(speciesTree); // set the tree to a complete time sliced one (TS)
            else
                Gfamily->setReconciledTreeTimeSlicesToNoTS();//no TS



            if(VerboseLevel > 1)
            {
                cout << "Tree computed and reconciled for distribution " << i << ". Likelihood: " << Gfamily->getTreeLikelihood() << ". Reconciliation score: " << Gfamily->getRecScore() << " TS status " << Gfamily->getReconciledTreeTimeSliceStatus()<< endl;
                if(VerboseLevel > 2)
                    Gfamily->printRecTree();
            }

            if(gBoolParams.find("scaffolding.mode")->second)
                fillLeafToSpMap( Gfamily->getRecTree(), LeafToSpMap);

            /////////////////////////////////
            /// Writing the Reconciled tree
            /////////////////////////////////
            //string filename = WriteOneGeneFamilyReconciledTree( Gfamily, gBoolParams.find("write.newick")->second , gBoolParams.find("hide.losses.newick")->second, recFilePrefix, i);
            //AddToFile(WrittenFileFileName, filename ); 
            AddOneGeneFamilyReconciledTreeToFile( Gfamily, gBoolParams.find("write.newick")->second , gBoolParams.find("hide.losses.newick")->second, recFileName, i);


            if( gBoolParams.find("write.genes")->second ) // writing 
            {
                WriteGeneOneGF( Gfamily->getRecTree(), i, genesFileName, RPOtoSpNameMap);
                //cout << "write.genes "<< i << endl;
            }
                

            //if(DUMPANDLOAD)
            //{
            //    if(gBoolParams.find("write.newick")->second) // we need to write a phyloxml version
            //    {
            //        cout << "Write nEWICK and DUMPANDLOAD -> PB"<< endl;
                    //string fileName = gPathPrefix + "reconciliations.xml"
                    //AddOneGeneFamilyReconciledTreeToFile( Gfamily, false , false, recFileName, i);
                    //
                    //AddToFile(WrittenFileFileName, filename );
            //    }
                //ReconciledTreesFileNames.push_back(filename);

            //    delete Gfamily;// we destroy the geneFamily to save memory
            //}
            //else
            //{
                //freeing some memory: no need to keep the CladesAndTripartitions and the DTLMatrix instance here...
                if(!gBoolParams.find("already.reconciled")->second)//only reconcile if it is not already done
                    Gfamily->dumpStuff();
            //}
        }

        FinishReconciledTreeFile(recFileName, gBoolParams.find("write.newick")->second);


        //if(DUMPANDLOAD) // we destroy the geneFamilies to save memory
        //{
        //    // the pointers content have already been deleted
        //    GeneFamilyList->clear();
        //    delete GeneFamilyList;
        //}
        //for(map <string,int>::iterator it = LeafToSpMap.begin() ; it != LeafToSpMap.end(); ++it )
        //{
        //    cout << it->first << " -> "<< it->second << endl;
        //}



        //for(unsigned i = 0 ; i < ReconciledTreesFileNames.size() ; i++)
        //    cout << ReconciledTreesFileNames[i] << " ";
        //cout << endl;


        if(VerboseLevel > 0)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for reconciliations:" << elapsed_secs << endl;
        }


        ///////////////////////////
        /// DeCo "proper" starts here
        ///////////////////////////

        // Create map of chromosome number expected by species
        map < string, int > speciesChrNb;
        map < string, int >::iterator itChr;
        if(gBoolParams.find("scaffolding.mode")->second)
        {
            if(gStringParams.find("chromosome.file")->second != "none"){
                speciesChrNb = storeSpeciesChrNumber(gStringParams.find("chromosome.file")->second);
                if(VerboseLevel > 1){
                    cout<<endl<<endl<<"Chromosome number expected by species:"<<endl;
                    for(itChr=speciesChrNb.begin();itChr!=speciesChrNb.end();itChr++)
                    {
                        cout<<"\t"<<(*itChr).first<<"\t-> "<<(*itChr).second<<" chr"<<endl;
                    }
                }
            }
            else{
                cerr<<"ERROR: If scaffolding mode is ON, then you have to give as input a chromosome file (See param: chromosome.file)!"<<endl;
                exit(EXIT_FAILURE);
            }
        }



        ////////////////////////////////////////////////////
        /// Putting adjacencies into Equivalence Class
        ////////////////////////////////////////////////////
        
        //cout<<endl<<"Before CreateEquivalenceClassFamilies"<<endl;

        // Dictionary that will contain c1 and c0 cost for EXT/EXT case for each species if using ARt-DeCo algorithm
        map<int,vector<float> > speciesC0C1;
        map<int,vector<float> >::iterator it1;
        map<int, map<string,int> > speGeneAdjNb;
        map<int, map<string, pair< int,int > > > speGeneExtremitiesAdjNb;
        vector < EquivalenceClassFamily > * ECFams;
        map<int,string> species_id_name;

        ECFams = CreateEquivalenceClassFamilies(adjacencies, adjacenciesOrientations, LeafToGFMap, nbGfam, gBoolParams.find("with.transfer")->second , gBoolParams.find("all.pair.equivalence.class")->second, verbose, superverbose);

        if(gBoolParams.find("scaffolding.mode")->second)
        {
            if(VerboseLevel > 1)
                cout<<endl<<"Scaffolding mode enabled"<<endl;
            //ECFams = CreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, speciesTree, speciesChrNb, speciesC0C1, species_id_name, speGeneAdjNb, gDoubleParams.find("ABreak.cost")->second, gBoolParams.find("with.transfer")->second , gBoolParams.find("all.pair.equivalence.class")->second, verbose, superverbose);

            //WARNING: Use species tree that can be changed !!!
            vector<MySpeciesNode*> extantSpecies = speciesTree->getLeaves();
            for(unsigned j =0; j < extantSpecies.size() ; j++)
            {
                string name=extantSpecies[j]->getName();
                int id=speciesTree->getRPO(extantSpecies[j]->getId());
                species_id_name[id]=name;
            }

            map < string, map <string , double> > * ASP = &adjacencyScores;

            bool problem = computeArtDeCoMaps(adjacencies, adjacenciesOrientations,
                                LeafToSpMap, speciesTree, speciesChrNb, speciesC0C1,
                                species_id_name, speGeneAdjNb, speGeneExtremitiesAdjNb,
                                 gDoubleParams.find("ABreak.cost")->second,verbose, superverbose,
                                ASP ,gBoolParams.find("scaffold.includes.scored.adjs")->second,
                                gDoubleParams.find("scaffolding.propagation.index")->second
                                );

            if(problem)
            {
                cerr << "scaffolding mode: less observed contigs than expected chromosome. Switching off scaffolding.mode to avoid abberant behaviour."<<endl;
                gBoolParams["scaffolding.mode"] = false;
            }
        }
        else
        {
            //cout<<endl<<"#################\n### DeCo mode ###\n#################"<<endl;
            //ECFams = CreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, gBoolParams.find("with.transfer")->second , gBoolParams.find("all.pair.equivalence.class")->second, verbose, superverbose);
        }

        speciesChrNb.clear();
        adjacencies.clear();
        
        LeafToGFMap.clear();
        LeafToSpMap.clear();


        //cout<<"After CreateEquivalenceClassFamilies"<<endl;

        if(VerboseLevel > 1)
        {
            for(it1=speciesC0C1.begin();it1!=speciesC0C1.end();it1++)
            {
                int spe=(*it1).first;
                string species=species_id_name[spe];
                cout<<"Species "<<species<<" ("<<spe<<"):"<<endl;
                cout<<"\tC0(1,1): "<<speciesC0C1[(*it1).first][0]<<endl;
                cout<<"\tC0(0/1,1/0): "<<speciesC0C1[(*it1).first][1]<<endl;
                cout<<"\tC0(0,0): "<<speciesC0C1[(*it1).first][2]<<endl;
                cout<<"\tC1(1,1): "<<speciesC0C1[(*it1).first][3]<<endl;
                cout<<"\tC1(0/1,1/0): "<<speciesC0C1[(*it1).first][4]<<endl;
                cout<<"\tC1(0,0): "<<speciesC0C1[(*it1).first][5]<<endl;
            }
            cout<<endl;
        }


        if(VerboseLevel > 1)
        {
            for(unsigned i = 0 ; i < ECFams->size(); i++)
                ECFams->at(i).printMe(superverbose);
        }


        string adjFileName;
        if(gBoolParams.find("write.adjacencies")->second)
        { // touching the adjacencies.txt file in order to reset it 
            adjFileName = gPathPrefix;
            if(( adjFileName[ adjFileName.size() - 1 ] != '/')&&(adjFileName.size()>0))
                adjFileName += ".";
            adjFileName += "adjacencies.txt";
        
            ofstream ofs;
            ofs.open(adjFileName.c_str(),ofstream::out );
            ofs.close(); // we close it directly

            //AddToFile(WrittenFileFileName, filename );
        }




        ////////////////////////////////////////////////////
        /// Computing and backtracking Adjacency matrices
        ////////////////////////////////////////////////////
        if(gBoolParams.find("write.adjacency.trees")->second)
        {
            InitAdjacencyTreeFile( adjTreeFileName, gBoolParams.find("write.newick")->second );
        }

        if(VerboseLevel > 0)
        {
            begin = clock();
            cout << "Computing Adjacency histories" << endl;
        }

        double gAdjWeight = gDoubleParams.find("Adjacency.weight")->second;
        if(gAdjWeight == 0)
            gAdjWeight = 0.000001; // in order to avoid division by 0

       

        //NECFsample * AdjTreeSample; <- old way
        ECFsample * AdjTreeSample; // new way

        map < int, ReconciledTree * > loadedRecTrees;

        ReconciledTree * Rtree1;
        ReconciledTree * Rtree2;
        
        int overflowed = 0; // will contain the number of equivalence class (ie. adj matrix) that have a potential overflow

        int nbIter =1;//if loss.aware = false, only 1 iteration is needed.
        string TmpFolder; // for Loss.aware + write.adjacency.trees options combined.
        
        if(gBoolParams.find("Loss.aware")->second)
        {

            nbIter = gIntParams.find("Loss.iteration")->second;

            if(VerboseLevel > 1)
                cout << "Loss awareness with "<< nbIter << " iterations."<< endl;
            
            if(gBoolParams.find("write.adjacency.trees") -> second)
            {// if both loss.aware and write.adj.trees options are on, RAM will die if we save everything in it so we use tmp files instead.
                TmpFolder = gStringParams.find("output.dir")->second + "/tmp/";
                int status = mkdir((TmpFolder).c_str(),0777);
                
                //cout << status <<endl;
                if(status == -1)
                {
                    cerr << "ERROR while creating temporary folder " << TmpFolder << " error no : " << errno << endl;
                    exit(1);
                }
                
            }


        }
        
        vector <int> countingIter(nbIter);// used to display the number of matrices computed at each iteration.
        vector <double> countingTime(nbIter);
        
        
        int totSample = gIntParams.find("nb.sample")->second ;
        bool stochasticBT = ((gBoolParams.find("use.boltzmann")->second)||(totSample>1));


        //data structures used for loss aware DeCoSTAR.
        map <int , map < string, vector < pair<string, int > > > > AdjGraph;
        map < int ,pair < vector < pair <string, string> >, bool > > FreeAdjacencies;
        
        //variables for timing CPU.
        clock_t begin_iter;
        clock_t end_iter;
        double elapsed_secs_iter;
        

        for(int iter = 0; iter != nbIter; iter++)
        {

            begin_iter = clock();
            
            if(gBoolParams.find("Loss.aware")->second && VerboseLevel > 1)
                cout << "this is iteration "<< iter +1 << " over " << nbIter<<endl;

            int currCount = 0; //used to display the number of matrices computed at each iteration.


            int actualIndex = ECFams->size() - 1;

            while(actualIndex>-1)
            {

                EquivalenceClassFamily * ECF = &ECFams->at(actualIndex);
            
                int gfam1 = ECF->getGfamily1();
                int gfam2 = ECF->getGfamily2();
                
                    
                Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
                Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

                pair < vector < pair <string, string> >, bool > FamiliesFreeAdjacencies;     

                if( iter == 0 )
                { //else it's useless, because it has been done.
                    ECF->refine(Rtree1, Rtree2, gBoolParams.find("all.pair.equivalence.class")->second ,gBoolParams.find("with.transfer")->second , superverbose); 
                }
                else if(gBoolParams.find("Loss.aware")->second)
                {
                    FamiliesFreeAdjacencies = FreeAdjacencies[actualIndex];
                    //Getting the free adjacency vectors of interest for the pair of the current families
                    FreeAdjacencies[actualIndex].second = true;//either it already is computed or it's going to be computed.
                }
                
                bool CompMatrix = false; // whether or not this matrix have to be recomputed too

                if( (iter == 0) || (iter > 0 && FamiliesFreeAdjacencies.first.size() > 0 && FamiliesFreeAdjacencies.second == false))
                {
                    CompMatrix = true;
                }

                vector < pair < pair<string, string> , double > > scoreAssociations;
                
                if( CompMatrix )
                {
                    if(VerboseLevel > 1)
                        cout << "Adjacency matrix computation for the Equivalence class family " << actualIndex << endl;
                    if(VerboseLevel > 2)
                        cout << "GeneFam 1: " << gfam1 << " GeneFam 2: " << gfam2 << endl;
                        
                    scoreAssociations = ComputeOneEquivalenceClassFamily(ECF, Rtree1, Rtree2, 
                                                    adjacencyScores, speciesC0C1, speGeneAdjNb, speGeneExtremitiesAdjNb,
                                                    gDoubleParams.find("AGain.cost")->second,
                                                    gDoubleParams.find("ABreak.cost")->second,
                                                    gBoolParams.find("use.boltzmann")->second,
                                                    gBoolParams.find("substract.reco.to.adj")->second,
                                                    gDoubleParams.find("dupli.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight,
                                                    gDoubleParams.find("loss.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight,
                                                    gDoubleParams.find("HGT.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight,
                                                    verbose,superverbose,
                                                    gDoubleParams.find("boltzmann.temperature")->second,
                                                    gDoubleParams.find("absence.penalty")->second,
                                                    gBoolParams.find("Loss.aware")->second, FamiliesFreeAdjacencies,
                                                    gDoubleParams.find("adjacency.score.log.base")->second,
                                                    gBoolParams.find("all.children.adj.transmission.rules")->second
                                                    ) ;


                    if(VerboseLevel > 1)
                        cout << "Matrix computed" << endl;
                    
                    currCount++;
                }
                else
                {
                    //else there is no change from a previous run, and scoreAssociations has been computed and save previously.
                    scoreAssociations = ECF -> getScoreA();
                }
                
                if(gBoolParams.find("use.boltzmann")->second)
                { // checking if any score has its absolute log10 above 200 (so the score is either < 1e-200 or > 1e200) is order to sniff out potential overflows
                    double threshold = 200;
                    vector <int> nbAbove = ECF->getNumberScoreWithAbsLog10Above( threshold );
                    int total = 0;
                    for(unsigned i = 0 ; i < nbAbove.size();i++)
                        total += nbAbove[i];
                    if(total > 0)
                    {
                        cerr << total <<" scores with value above 1e"<<threshold << " or below 1e-"<<threshold<<" in the equivalence class family between "<<gfam1 <<" and "<<gfam2 ;
                        cerr << " -> potential overflow. Be wary of your results and try again with a boltzmann.temperature closer to 1 in case of problem."<<endl;
                    }
                }

                    
                // data structure used if we want to write adjs
                map <string, map <string, int> > AdjIndexMap;
                vector < double > AdjInScoreList;
                vector < double > AdjOutScoreList;
                vector < int > AdjSpeList;
                if(gBoolParams.find("write.adjacencies")->second)
                {
                    if(CompMatrix == true) // filling the different map
                    {
                        for(unsigned adjIndex = 0 ; adjIndex < scoreAssociations.size() ; adjIndex++)
                        {
                            AdjInScoreList.push_back( scoreAssociations[adjIndex].second );
                            AdjOutScoreList.push_back( 0 );
                            AdjIndexMap[scoreAssociations[adjIndex].first.first][scoreAssociations[adjIndex].first.second] = adjIndex;

                            int geneNodeId =  Rtree1->getIdWithName(scoreAssociations[adjIndex].first.first);
                            if(geneNodeId == -1) // sanity check
                            {
                                cerr << "ERROR : failed to find name "<< scoreAssociations[adjIndex].first.first << " in gene family "<< gfam1<<endl;
                            }
                            AdjSpeList.push_back( Rtree1->getNodeSpecies(geneNodeId) );

                            //int idx =  AdjIndexMap[scoreAssociations[adjIndex].first.first][scoreAssociations[adjIndex].first.second];
                            //cout << scoreAssociations[adjIndex].first.first << "-" << scoreAssociations[adjIndex].first.second << idx << " " << AdjSpeList[idx] << " " << AdjInScoreList[idx] << " " << AdjOutScoreList[idx]<< endl;
                        }
                    }
                    else
                    { // this was not recomputed, just use the results that were stored in the ECF
                        AdjIndexMap = ECF -> getAdjIndexMap();
                        AdjInScoreList = ECF -> getAdjInScoreList();
                        AdjOutScoreList = ECF -> getAdjOutScoreList();
                        AdjSpeList = ECF -> getAdjSpeList();        
                    }
                }



                if( CompMatrix == true)
                { // if the matrix was computed here.
                    for(unsigned sampleIndex = 1 ; sampleIndex <= totSample; sampleIndex++ ) // for each sample to do
                    {
                        
                        AdjTreeSample = new ECFsample();// new
                    
                        backtrackOnetimeOneEquivalenceClassFamily(AdjTreeSample, ECF, GeneFamilyList , 
                                                            stochasticBT,
                                                            verbose, superverbose,
                                                            gBoolParams.find("always.AGain")->second , gDoubleParams.find("C1.Advantage")->second, overflowed , gBoolParams.find("count.number.backtracks")->second);
                        
                                                
                        if(VerboseLevel > 1)
                            cout << sampleIndex << "/" << totSample << "\r";
                        
                        bool init = false;
                        bool finish = false;

                        if( sampleIndex == 1 )
                            init = true;

                        if(sampleIndex == totSample)
                            finish = true;
                        
                        if(gBoolParams.find("write.adjacency.trees")->second)
                        {
                            
                            if(iter == (nbIter - 1))// if it's the last iteration and matrice has been computed.
                            {
                                AddECFsampleToFile(adjTreeFileName, ECF, AdjTreeSample,
                                               sampleIndex, gBoolParams.find("write.newick")->second, gBoolParams.find("hide.losses.newick")->second, 
                                               gDoubleParams.find("AGain.cost")->second,
                                               gDoubleParams.find("ABreak.cost")->second, 
                                               init, finish);
                                if(gBoolParams.find("Loss.aware") ->second && finish == true)// if loss.aware option is on, tmp file that was previously filled has to be deleted.
                                {
                                    int status = remove((ECF -> getTmpFile()).c_str());//deleting the tmp file.
                                    if(status == -1)
                                    {
                                        cerr << "File removal of " << ECF -> getTmpFile() << " did not work. You might have to proceed manually."<<endl;
                                    }
                                }
                            }
                            else // this is not the last iteration -> we save the sample to a temporary file
                            {
                                // not the last iteration. save AdjTree to tmp file.
                                if(init == true) // first sample -> setting the file
                                    ECF -> setTmpFile(TmpFolder, gStringParams.find("output.prefix")->second, actualIndex);
                                
                                saveECFsampleToTmpFile(ECF, AdjTreeSample,
                                            sampleIndex, gBoolParams.find("write.newick")->second, gBoolParams.find("hide.losses.newick")->second, 
                                            gDoubleParams.find("AGain.cost")->second,
                                            gDoubleParams.find("ABreak.cost")->second, 
                                            init, finish);
                                
                            }                               
                            //saving to file, and giving filename to ECF so it knows what to look for.

                        }

                        if(gBoolParams.find("write.adjacencies")->second)
                        {
                            for(unsigned i = 0 ; i < AdjTreeSample->size() ; i++)
                            {
                                UpdateAdjMapsFromForest( AdjTreeSample->at(i), gIntParams.find("nb.sample")->second, 
                                    AdjIndexMap,
                                    AdjInScoreList,
                                    AdjOutScoreList,
                                    AdjSpeList,
                                    Rtree1,
                                    Rtree2
                                    );
                            }
                        }
                        

                        for(unsigned i = 0 ; i < AdjTreeSample->size() ; i++)
                        {
                            for(unsigned j = 0 ; j < AdjTreeSample->at(i)->size() ; j++)
                                delete AdjTreeSample->at(i)->at(j);
                            delete AdjTreeSample->at(i);
                        }
                        delete AdjTreeSample;
                        
                    }
                    if(iter != nbIter - 1)
                    {// saving the different maps used for writing adjacencies.
                        ECF -> setAdjIndexMap(AdjIndexMap);
                        ECF -> setAdjInScoreList(AdjInScoreList);
                        ECF -> setAdjOutScoreList(AdjOutScoreList);
                        ECF -> setAdjSpeList(AdjSpeList);
                    }//else it's useless, we won't need it.
                }
                else if(gBoolParams.find("write.adjacency.trees")->second && iter == nbIter -1)
                {
                    copyTmpToFile(adjTreeFileName, ECF);
                    //if last iteration, copy AdjTree from tmp file to main file.
                }


                if(VerboseLevel > 1)
                {
                    if(CompMatrix == true)
                        cout << "backtracked Equivalence class family " << actualIndex << " " << gIntParams.find("nb.sample")->second << " times." << endl;
                    else
                        cout << "Kept previouly backtracked adjacency trees results for Equivalence class family " << actualIndex << endl;
                 }
                

                if( (iter != nbIter - 1) && (gBoolParams.find("Loss.aware")->second) )
                {//If it's not the last iteration, we don't need to write the adjacencies but we need them anyway.
                    //second check should be useless in current state,as there are more than one iteration only with loss.awareness in current version (June 2017). 
                    ReadAdjMaps(actualIndex, AdjGraph,
                                AdjIndexMap,
                                AdjInScoreList,
                                AdjOutScoreList,
                                AdjSpeList); // fills the adjacency graph object
                    //cout << "added adjacencies to AdjGraph"<<endl;
                }
                else
                {//do this on the last iteration
                    if(gBoolParams.find("write.adjacencies")->second)
                    {
                        AddAdjMapsToFile(adjFileName, ECF->getSens1(), ECF->getSens2(),
                                        AdjIndexMap,
                                        AdjInScoreList,
                                        AdjOutScoreList,
                                        AdjSpeList);
                        //cout << "added adjacencies to file"<<endl;
                    }
                }

                if(iter == (nbIter - 1))//only on last iteration
                {
                    ECFams->erase(ECFams->end());// erasing the last element (the one we just treated) (erasing the last element in a vector is much more time efficient than erasing the first one)
                }
                else
                {//to free some memory we delete the adjMatrix after backtracking
                    if(CompMatrix == true)
                    { // there actually is an adj matrix to delete
                        int NbEqClasses = ECF -> getNbEqClasses();

                        for(unsigned i = 0; i < NbEqClasses ; i++)
                        {
                            delete ECF -> getEClass(i) -> getAmat();//erasing matrices

                            ECF -> setSetAdjMatrixFamily(false, i);// informating the ECF that the matrices are unset
                        }

                    }//else nothing to delete.
                }
                
                actualIndex--;
            }

            
            //Those are made to fill in the container for the free adjacencies
            if(gBoolParams.find("Loss.aware")->second)
            {
                if(iter != (nbIter - 1))//need to do that on all but the last iteration
                {
                    if(verbose > 1)
                        cout << "establishing free adj list"<<endl;
                        
                    MakeFreeAdjacencies(AdjGraph, GeneFamilyList, FreeAdjacencies);
                    
                    if(verbose > 1)
                        cout << "Free Adjacencies list established. Rerunning DeCoSTAR."<< endl;
                    AdjGraph.clear(); // reinitialize the adjacency graph.
                }//else it's useless because this is the last run.
            }
            
            //some variables for CPU timing.
            end_iter = clock();
            elapsed_secs_iter = double(end_iter - begin_iter) / CLOCKS_PER_SEC;
            countingIter[iter] = currCount;
            countingTime[iter] = elapsed_secs_iter;

        }//end of for loop
        if(gBoolParams.find("Loss.aware")-> second and VerboseLevel > 1)
        {
            for( unsigned i = 0; i < countingIter.size();i++){
                cout << "Iteration number " << i << " computed " << countingIter[i] << " matrices in " << countingTime[i] << " seconds."<<endl;
            }
        }
        

            


        if( (overflowed > 0) && (gBoolParams.find("use.boltzmann")->second) )
        {
            cout << "!!POTENTIAL OVERFLOW!!" << endl;
            cout << "No tree was backtracked on "<< overflowed <<" samples of equivalence classes."<< endl;
            
            cout << "If you haven't put any score on them (in the adjacency file) that would justify their absence, then you are probably suffering from overflow."<<endl;
            cout << "To avoid this problem, choose a different boltzmann.temperature (usually one that is closer to 1)."<<endl;
        }

        if(gBoolParams.find("write.adjacency.trees")->second)
        {
            FinishAdjacencyTreeFile(adjTreeFileName, gBoolParams.find("write.newick")->second);

            if(gBoolParams.find("Loss.aware") -> second)
            {//removing tmp folder and everything it contains.
                int status = remove((TmpFolder).c_str());
                    
                if(status == -1)
                {
                    cerr << "tmp folder removal did not work. Exited with system error " << errno<<endl;
                }
                
            }
        }

        for(map<int, ReconciledTree *>::iterator it = loadedRecTrees.begin(); it != loadedRecTrees.end(); ++it)
        {
            delete it->second;
        }        
        loadedRecTrees.clear();


        for(int i = 0; i < GeneFamilyList->size(); i++)
        {
            delete GeneFamilyList->at(i);
        }
        GeneFamilyList->clear();
        delete GeneFamilyList;


        if(VerboseLevel > 0)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for Adjacency Matrix computing and backtracking:" << elapsed_secs << endl;
            cout << "Results written in directory: " << gStringParams.find("output.dir")->second << endl;
        }


    }

    catch(exception & e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    return 0;
}