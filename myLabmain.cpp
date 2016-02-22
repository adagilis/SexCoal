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
boost::mt19937 gen;		// Random generator recommended by Boost; declared globally so that it can be 
//	 seeded in main.cpp and then used elsewhere
/**********************************************************************/

vector<double> myPopGen(vector<Base> mySnp){ //10 chrX and 10 chrY
	vector<Base> myXseq(mySnp.begin(), mySnp.begin() + 90),
				 myYseq(mySnp.begin() + 90, mySnp.end());
				 
/*
		for (int i = 0; i < myXseq.size(); ++i)
		{
	        	cout << myXseq[i] << "\t";
        		
    		}
		cout << endl;
		for (int i = 0; i < myYseq.size(); ++i)
		{
	        	cout << myYseq[i] << "\t";
        	}
		cout << endl;
		for (int i = 0; i < mySnp.size(); ++i)
		{
	        	cout << mySnp[i] << "\t";
        	}
		cout << endl;
*/				   
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
	
	for( int i = 0; i<4;i++){
		myPiX -= chrx.at(i)*chrx.at(i);
		myPiY -= chry.at(i)*chry.at(i);
		myPiT -= ( 3*chrx.at(i) + chry.at(i) )*( 3*chrx.at(i) + chry.at(i) )/16; //need change as the sample size of x and y varies
	}
	if(myPiT ==0 ){myFst = 0;}else{ myFst = 1 - (myPiX*3 + myPiY)/(4*myPiT);} 
	vector<double> myOutStats;
	myOutStats.push_back(myPiX);
	myOutStats.push_back(myPiY);
	myOutStats.push_back(myFst);
	return( myOutStats );
}
/**********************************************************************/


// ======================================================================================


int main (int argc, char *argv[]) {
    
    
    
// A static seed, useful for debuggin:
//
//gen.seed(42U);	
    
// In real use, seed the random number generator from the system clock -- don't forget to do this!:
//
	gen.seed(static_cast<unsigned int>(std::time(0)));	
    
// Things for measuring performance:
//
//	progress_timer timing;
//	timing.restart();
 	
//F	timing.restart();

	int coreid =1;
	if(argc>1) sscanf (argv[1],"%i", &coreid);	
	
    
    
//	for(int steps=0; steps<loops; ++steps){
	        UINT totalEvents = 0;
//		cout<<"Step >> "<<steps<<'\n'; 
		
		stringstream in;				//stream class on strings, to take in the name of parameter file//
//		in <<"param_"<<steps<<".txt";			//read in parameter files named like "param_1.txt"//  
		in <<"my.model.input"<<coreid;				//read in my own parameter files names as "my.model.input"//
		Parameters params(in.str().c_str());		//read in parameters in a pre defined parameter data structure/object //
		int nPops=params.paramData->popSizeVec.size();	//get the number of populations in the simulation model //
		
//		for(int n = 0; n< nPops; n++){
//			cout << "Population Sizes (" << (n+1) << "): " << params.getPopulationMaxSizes().at(n) << endl;
//		}						//print population size for each population //
		
		stringstream myOutput;				//to hold the name of output file//
		myOutput <<"my.model.output"<<coreid; 			//set the name of output file//
		ofstream fout( myOutput.str().c_str() );  	// output file stream//
		
/*        
		for(int k=0; k<params.paramData->neut_site.size(); ++k){
			fout<<" SITE_"<<params.paramData->neut_site.at(k);
		}
		fout<<'\n';					//print the location of neutral sites in the header of output file//	
*/        
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
		//
		//
		//
		
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
        //							//
        //							//
        //							//
        // Create a brave new world:				//
        //							//
        // Not sure what this exactly means, but will find out 	//
		
       // vector<double> outTime(params.paramData->neut_site.size(), 0);	//define a vector with doubles, the size of outTime
									//is the number of neutral sites, and it is initialized 
									//with zeroes. This vector will hold the coalescence time
									//for each locus
        
        for (int timer=0; timer<nRuns; ++timer){
	        for (int timer2=0; timer2<nNsites; ++timer2){
//to guard against Nsite == SDR or siteA. ccd. A better solution is to run a different simulation program at SDR and siteA.		
		if( myNsites.at(timer2) == params.paramData->Sexsite_pos || myNsites.at(timer2) == params.paramData->siteA_pos ){
			fout << " " << 2*params.paramData->freqA_inx*(1-params.paramData->freqA_inx)
				 << " " << 2*params.paramData->freqA_iny*(1-params.paramData->freqA_iny)
				 << " " << 1 
				 << endl; 
		}
		else{
		vector<double> foo(1,myNsites.at(timer2));
		//cout << myNsites.at(timer2) << endl;//debug
		//cout << foo.size() << endl;//debug
		params.setNsites( foo );
//		};
		//cout << "size" << params.paramData->neut_site.size() << endl;
		//cout << nRuns << ' ' << nNsites << ' ' << params.paramData->neut_site.at(0) << endl;//debug
            
            World * world = new World(params.getpData());
            
            while(!world->simulationFinished()){
                totalEvents+= world->simulateGeneration(mig_prob);
              }
            
//            vector< shared_ptr < ARGNode > > allNodes = world->getARGVec();	//this line looks intereseing, it may hold all the address of ARGnode information//
            
            vector< pair<double, shared_ptr < ARGNode > > > allsites = world->getOutputSegments(); // get breakpoints
            
            for(int k=0; k<params.paramData->neut_site.size();++k){//anyway, just one locus
                
 //               shared_ptr < ARGNode > tempArgNode= getNextSiteNode(params.paramData->neut_site.at(k), allNodes.at(allNodes.size()-1));

//ccd test it now
//		tempArgNode->outputARG();
					
					shared_ptr < SiteNode > tmpSiteNode;//-------------------------------->>>>>>>>>>>>>>>>>>>>
					tmpSiteNode.reset(  new SiteNode( world->getMaleMutRatio(), allsites.at(k).first, allsites.at(k).second ) );

//		SiteNode tempSiteNode = SiteNode(params.paramData->neut_site.at(k), tempArgNode);
//		fout << tmpSiteNode->getTime() << endl;

		vector<Base>  myOutSeq(params.paramData->initChr.size(),A);		

		tmpSiteNode->snpDescendants(world->getMaleMutRatio(), myOutSeq);//*????????????

		vector<double> myOutStat = myPopGen(myOutSeq); // Pi_x, Pi_y, Fst
//cout << " done pop gen " << ss << endl;					

	
			fout << myOutStat.at(0) 
			 << " " << myOutStat.at(1)
			 << " " << myOutStat.at(2) 
			 << endl; 		
	
//ccd test it now

                //outTime.at(k)+= tempArgNode->getTime()/(double)params.paramData->popSizeVec[0];
                //fout<<" "<<tempArgNode->getTime()/(double)params.paramData->popSizeVec[0] << endl;
                
            }									//need to read ARGNode file
            //fout<<'\n';
            
            // Clean up:
            //
            delete world;
        };//else
        }//nRun
	}//nNsites		
//		cout<<"Avg Events = "<<(float)totalEvents/(float)nRuns<<'\n'; 		
		
//		cout<<"Sample (Sex, A): ("<<params.paramData->initChr.at(0)->getContext().Sexsite<<","<<params.paramData->initChr.at(0)->getContext().siteA<<")-("
//		<<params.paramData->initChr.at(1)->getContext().Sexsite<<","<<params.paramData->initChr.at(1)->getContext().siteA<<")_ Avg Sample ="<<params.paramData->avg<<'\n';
//		cout<<"E[T]=";
//		for(int k=0; k<params.paramData->neut_site.size();++k){
//			cout<<" "<<outTime.at(k)/nRuns;	
//		}
//		cout<<"\nElapsed simTime ="<<timing.elapsed()/60<<"m\n\n";
        
		fout.close();
        
        
//	}
	
	return 0;
}
