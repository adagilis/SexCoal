/*  Version:  sexchromo_Anti
 
	0. This program simulated the ARG for sex chromosomes. It allows for one locus A under sex-differential selection, and it assumes the difference between the equilibrium frequencies at A from eggs and sperm are negligible.
 	1. Changde Cheng adapted this file for EasyABC
		Thu Aug 21 11:13:11 EDT 2014
 		use my.model.input and my.model.output for parameters and stats
 
 */



// Includes from STL:
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include <string>
using std::string;



#include <ctime>

// Includes from Boost:
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
#include <boost/scoped_ptr.hpp>
using boost::scoped_ptr;
#include <boost/random.hpp>
using boost::mt19937;
//#include <boost/progress.hpp>
//using boost::progress_timer; // comment off for use together with R //
#include <limits>

// Includes for our files:
//	#include "parameters.h"		
#include "world.h"
#include "ran_mk.h"
#include "chromosome.h"
#include "sitenode.h"


// Global declarations:
boost::mt19937 gen;		
// Random generator recommended by Boost; declared globally so that it can be 
	
/**********************************************************************/
/* Changde's little popgen function */
vector<double> myPopGen(vector<Base> mySnp, int myXn ){ // first myXn chrX and rest is chrY
	vector<Base> myXseq(mySnp.begin(), mySnp.begin() + myXn), myYseq(mySnp.begin() + myXn, mySnp.end());
	vector<double> chrx,chry;
	for( int i = 0; i<4;i++){
		chrx.push_back( std::count ( myXseq.begin(), myXseq.end(), i) ); 
		chry.push_back( std::count ( myYseq.begin(), myYseq.end(), i) );
	}

	for( int i =0; i<4; i++){
		chrx[i] = chrx[i]/myXseq.size();
		chry[i] = chry[i]/myYseq.size();	
	}
 
	double myPiX = 1;
	double myPiY = 1;
	double myPiT = 1;
	double myFst = 0;
	double mySFS = 0;
		
	for( int i = 0; i<4;i++){
		myPiX -= chrx.at(i)*chrx.at(i);
		myPiY -= chry.at(i)*chry.at(i);
		myPiT -= ( chrx.at(i) + chry.at(i) )*( chrx.at(i) + chry.at(i) )/4; 
	}
	
	if(myPiT ==0 ){myFst = 0;}else{ myFst = 1 - (myPiX + myPiY)/(2*myPiT);} 
	
 vector<double> myOutStats;
	myOutStats.push_back(myPiX);
	myOutStats.push_back(myPiY);
	myOutStats.push_back(myPiT);
	myOutStats.push_back(myFst);
	return( myOutStats ); // return: pi(x), pi(y), pi(total), Fst(x,y)
}
/**********************************************************************/




int main (int argc, char *argv[]) {

 gen.seed(static_cast<unsigned int>(std::time(0)));	

	int coreid =1;
	if(argc>1) sscanf (argv[1],"%i", &coreid);	
	
 UINT totalEvents = 0;
	stringstream in;				//stream class on strings, to take in the name of parameter file//
//		in <<"param_"<<steps<<".txt";			read in parameter files named like "param_1.txt"//  
	in <<"my.model.input"<<coreid;			//read in my own parameter files names as "my.model.input"//
	Parameters params(in.str().c_str());		//read in parameters in a pre defined parameter data structure/object //
	int nPops=params.paramData->popSizeVec.size();	//get the number of populations in the simulation model //
	stringstream myOutput;				//to hold the name of output file//
	myOutput <<"my.model.output"<<coreid; 			//set the name of output file//
	ofstream fout( myOutput.str().c_str() );  	// output file stream//
		
	vector<double> myNsites = params.paramData->neut_site; //get the Nsites vector, ccd, 10/10/2014//
	int nNsites=params.paramData->neut_site.size();   // get the number of neutral sites been simulated on a chromosome        
	int nRuns=params.paramData->nRuns;		//get the number of simulation runs//
		
		// a simple stepping-stone migration matrix is formed. 			//
		// according to correspondence with Dr. Rafael Guerrero			//
		// this part of code with migration is not tested and need debug	//
		
	vector < vector <double> > mig_prob;
	mig_prob.resize(nPops);
	double m= params.paramData->mig_prob.at(0)/(2*params.paramData->popSizeVec[0]);

	for(int i = 0; i < nPops; i++) mig_prob.at(i).resize(nPops);
	for(int i = 0; i < nPops; i++){
		for(int j = 0; j < nPops; j++){
			if (nPops==1) mig_prob.at(i).at(j)=1;
			else if (i==j) mig_prob.at(i).at(j)=1 - m;					// probability of staying in current pop
			else if (i==j-1) { 
				mig_prob.at(i).at(j)= m/2;								// probability of migrating to neighboring pop at left
				if (i==0)  mig_prob.at(i).at(j)= m;						// probability of migrating if pop is at beginning of line  
			}
			else if (i==j+1) { mig_prob.at(i).at(j)=m/2;				// probability of migrating to neighboring pop at left
				if (i==nPops-1)   mig_prob.at(i).at(j)= m;				// probability of migrating if pop is at end of line
			}
			else  mig_prob.at(i).at(j)=0;								// probability of migrating to any other pop
		}
	}
		
	// a set of N maps. In each map there are pairs of the following values: 
	// <int> is the population 
	// <float> is the cumulative probability of migration, to be compared against a random uniform number [0,1] 
	//
	// Example: Set of four populations:: An individual in population 2 can migrate to pops {0,1,3} with prob {0,0.05, 0.05}   
	// so the map for population 2 will look like this: { {0, 0}, {1, 0.95}, {2, 0.9}, {3, 1} }
	// For the use of these probabilities in migration, check out comments on migrate() function
		
		vector < map< double, int> > mig;
		mig.resize(nPops);
		
		for(int i = 0; i < nPops; i++){
			map< double, int> migrate;
			migrate[mig_prob.at(i).at(i)]= i ;
			double cumulative = mig_prob.at(i).at(i);
			for(int j = 0; j < nPops; j++){
				if (i!=j && mig_prob.at(i).at(j)!=0) {
					cumulative += mig_prob.at(i).at(j);
					migrate[cumulative]= j;
				}
			}
			mig.at(i)= migrate;
		}

 // start worlds
    for (int timer=0; timer<nRuns; ++timer){
        World * world = new World(params.getpData());
 
        while(!world->simulationFinished()){
            totalEvents+= world->simulateGeneration( params.paramData->recombination_rate, mig_prob );
        }
 
		vector< pair<double, shared_ptr < ARGNode > > > allsites = world->getOutputSegments();
		std::sort( allsites.begin(), allsites.end() );
		vector<double> mySegSites;
        vector<double> myOutStat;
        int mySegSitesNumber = 0;
		for(int k=1; k< allsites.size();++k)
		{   
            shared_ptr < SiteNode > tmpSiteNode;
            tmpSiteNode.reset(  new SiteNode( allsites.at(k-1).first, allsites.at(k-1).second ) );
            mySegSitesNumber = randpoisson( (allsites.at(k).first - allsites.at(k-1).first )*tmpSiteNode->getTotalTime()*params.paramData->mutation_rate );
            for(int ss = 0; ss < mySegSitesNumber; ++ss)
            {
                mySegSites.push_back( randreal(allsites.at(k-1).first, allsites.at(k).first ) );
                vector<Base>  myOutSeq(params.paramData->initChr.size(),A);
                tmpSiteNode->snpDescendants(myOutSeq);
                /*myOutStat = myPopGen(myOutSeq,params.paramData->initNx);*/
     
                fout <<  mySegSites.back()*params.paramData->recombination_rate*double(params.paramData->popSizeVec[0])*2;
                
                /* outputs popgen statistics
                for(int i=0; i<4; i++){
                    #fout << " " << myOutStat.at(i);
                #}
                */
                 
                for(int q=0; q<myOutSeq.size();q++){
                    fout<<" "<<myOutSeq[q];
                }
        
                fout << endl;
            }
		}
  delete world;
 }
 
	fout.close();
	
	return 0;
}
