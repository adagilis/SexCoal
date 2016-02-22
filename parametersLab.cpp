/**
 * Parameters.cpp
 * 
 * This reads in any and all parameters for a simulation from the file "Parameters.txt"
 * or from a file with necessary parameters for a given name
 * 
 * Modified and annotated by Changde Cheng X-2014.
 *
**/


// Includes from STL:
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;

// Includes from Boost:
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

// Includes from our files:
#include "typedefs.h"
#include "parameters.h"
#include "chromosome.h"


Parameters::Parameters(const char* insstring) {
/**
 * Constructor.
 * Reads values from a file with parameters into Parameter data class.
**/
	
	DBG("Parameters:: Constructing...")
	
	paramData.reset(new ParameterData);

	
/**
 * initail file, line and numbers for temperary storage of input.
**/
	ifstream fin( insstring);
	string line;
	int temp;
	double doubtemp;
	
/**
 * Read in different params now.
**/

/**
 * number of runs	
**/
	getline( fin, line);
	istringstream iss(line);
	while(iss >> temp){
		paramData->nRuns=temp;
	}

// for debug printing
//	cout<<"nRuns= "<<paramData->nRuns<<'\n';


/**
 * Population Size (EX: 10000).
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> temp){
		paramData->popSizeVec.push_back(temp);
	}


/**
 * Mutation rate (EX: 0.00001). Added by Changde, need further work on this.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->mutation_rate=doubtemp;//should I get it divided by effective population size?
	}

	
/**
 * Migration probability (EX: 0.1). Haven't been tested or debugged by Changde yet.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->mig_prob.push_back(doubtemp);
	}

	
/**
 * Ages of phylogenetic nodes. Assumed age of siteA first, Y second
 * in number of generations, starting by the youngest (EX: 1500 3000)
 * not tested by Changde yet.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	int nodes=0;
//	cout<<"node ages: ";
	while(iss >> doubtemp){
		paramData->ages.push_back(doubtemp);
//		cout<<paramData->ages.at(nodes)<<" ";
		++nodes;
	}
//	if(paramData->ages[1]!=0 && paramData->ages[1]!=1) cout<<"Check your ages, second number should be 0 or 1\n";


/**
 * read in sex determination region (SDR) position in Rho.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
// No rescaled by effective population size.
		paramData->Sexsite_pos =doubtemp;
	}
// for debug printing
//	cout<<"\nSDR at: "<<paramData->Sexsite_pos<<'\n';


/**
 * read in siteA position.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
// No rescaled by effective population size. =
		paramData->siteA_pos =doubtemp;
	}
// debug printing:
//	cout<<"siteA at: "<<paramData->siteA_pos<<'\n';
	

/**
 * Male to Female recombination rate.
 *.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->maleRecRatio =doubtemp;
	}
	
/**
 * Male to Female population rate.
 * .
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->malePopRatio =doubtemp;
	}
	
/**
 * Male to Female mutation rate.
 * .
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->maleMutRatio =doubtemp;
	}
		
/**
 * read in freqA in X (frequency of sex atangonisitc sites on X, [0,1]).
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->freqA_inx=doubtemp;
	}
// debug printing:
//	cout<<"freqA in X = "<<paramData->freqA_inx<<'\n';
	

/**
 * read in freqA in Y (frequency of sex antagonisitic site on Y, range [0,1]).
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->freqA_iny=doubtemp;
	}
// debug printing:
//	cout<<"freqA in Y = "<<paramData->freqA_iny<<'\n';
	

/**
 * read in freqA in pre-Y-sweep. Not used by Changde yet.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->preY_afreq=doubtemp;
	}
// debug printing:
//	cout<<"FreqA pre-Y-sweep = "<<paramData->preY_afreq<<'\n';


/**
 * Make the vector of segments from that single segment:
**/


/**
 * read in the number of carriers/chromosomes to be sampled
**/
	int nCarrier=0;
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> temp){
		nCarrier=temp;
	}
	

/**
 * Make average runs
 * 	0 = reads carriers below, 
 * 	1 = random (sex and siteA) in pop,
 * 	2 = random (siteA) in pop by SDR.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> temp){
		paramData->avg=temp;
	}
	
	
/**
 * the list of neutral sites
**/
/*	
	getline( fin, line);
	iss.clear();
	iss.str(line);
// debug printing.
//	cout<<"Neutral Sites at:";
	while(iss >> doubtemp){
// Not rescaled by effective population size.
		paramData->neut_site.push_back(doubtemp );
// debug printing
//		cout<<" "<<doubtemp;
	}
*/	
//	double mySavedSnps[378] =  {0.000542352,0.000595152,0.000674016,0.000758976,0.000783792,0.000804144,0.000819072,0.00082104,0.00084864,0.000868656,0.000896592,0.000949104,0.001031616,0.001271136,0.00128112,0.001322208,0.001381872,0.001414896,0.00141624,0.001462992,0.001553424,0.001638432,0.001761312,0.001901088,0.00190536,0.00193224,0.002033568,0.002149776,0.002241024,0.0022656,0.002319504,0.00232032,0.002389872,0.002418144,0.002472336,0.002530848,0.002548176,0.002652576,0.002675856,0.0026916,0.002703792,0.002705088,0.002726736,0.002763936,0.002786448,0.00283152,0.002941104,0.002992176,0.003031392,0.00308976,0.003121296,0.003140976,0.003245376,0.00327072,0.003289296,0.003318336,0.003331632,0.003342,0.003408912,0.0034152,0.0034344,0.003447264,0.003487632,0.00350472,0.003543216,0.003648912,0.003737424,0.00375168,0.00378288,0.00378888,0.003829392,0.003836976,0.003881904,0.003897408,0.003932736,0.00394152,0.003944256,0.004007328,0.004042608,0.004183776,0.004250448,0.004351296,0.004370064,0.004437744,0.004442736,0.004446096,0.004446816,0.004483536,0.00461496,0.004728192,0.004877616,0.00489432,0.005137152,0.00514512,0.00515136,0.005195616,0.005273328,0.005337072,0.00535584,0.005399472,0.005406144,0.005508048,0.005541552,0.005589936,0.005593152,0.005621856,0.005624736,0.005633616,0.00565056,0.005740512,0.00580056,0.005802816,0.005808864,0.00582696,0.005898144,0.005931552,0.005936016,0.005951568,0.005998512,0.006094656,0.006311424,0.006319824,0.00633768,0.0063984,0.006475872,0.006671616,0.006676416,0.006684576,0.006688992,0.006691296,0.0066972,0.006714672,0.006744768,0.006758736,0.006768912,0.006770448,0.006777456,0.006777984,0.006784464,0.006797952,0.006826992,0.006851904,0.00686208,0.00687744,0.006926448,0.007010928,0.007021152,0.0070296,0.007033344,0.007074816,0.007124256,0.007127568,0.007151808,0.007167792,0.00716904,0.007175568,0.007227744,0.00722952,0.00725136,0.007283664,0.007334256,0.007348512,0.007431744,0.007487472,0.007695264,0.007696032,0.007730208,0.007750416,0.007780272,0.007783248,0.007816848,0.007892016,0.008111568,0.008217312,0.008246352,0.008252736,0.008267376,0.008292384,0.008345664,0.0086184,0.008809248,0.008894256,0.009486768,0.009549408,0.00965976,0.009739584,0.0099852,0.010062912,0.010066464,0.010101792,0.01024008,0.010242096,0.010243872,0.010288176,0.01033056,0.010373472,0.010513248,0.010565904,0.010590768,0.010642128,0.010689408,0.0107208,0.01074288,0.01083576,0.0108948,0.010934976,0.010935984,0.0109392,0.010946112,0.010960464,0.01098888,0.011001024,0.011013312,0.011013696,0.011032848,0.011067168,0.011070288,0.01107552,0.011095344,0.011108496,0.01111512,0.011148672,0.011167824,0.011171184,0.011266704,0.011275968,0.011279616,0.01128504,0.0112956,0.011304864,0.011306256,0.01134648,0.011366448,0.011370288,0.011411376,0.011411712,0.011420256,0.011425776,0.011426928,0.011428752,0.01151592,0.01156368,0.01167672,0.01167696,0.0117348,0.011773152,0.01200912,0.012076944,0.012201936,0.012208704,0.012261984,0.012262272,0.012449664,0.012452112,0.01257072,0.012593664,0.012596256,0.012624144,0.012638112,0.012647136,0.012670656,0.012680064,0.012749664,0.012753744,0.012815376,0.01296984,0.012980928,0.013256784,0.013332,0.013369056,0.013370592,0.013609776,0.013655664,0.013724016,0.013765344,0.013912656,0.013991616,0.014104176,0.014170176,0.014405232,0.014422512,0.014423472,0.014443104,0.014459664,0.014466,0.014493504,0.01449696,0.014503536,0.014566656,0.014569488,0.014577312,0.01458024,0.01458744,0.014626368,0.01462776,0.014628816,0.014678592,0.014687136,0.01481928,0.014934864,0.014973888,0.015030864,0.015054912,0.015071568,0.015086448,0.015134064,0.015208032,0.015215184,0.015234672,0.015350064,0.015374736,0.015410064,0.015411216,0.015453072,0.015507696,0.015535968,0.015561312,0.015568512,0.01557288,0.015586416,0.01561272,0.015784368,0.015840144,0.016117248,0.016123392,0.016135392,0.016159728,0.01616808,0.01622304,0.016247856,0.016253136,0.016274976,0.016279536,0.016440192,0.016465824,0.016469328,0.01652712,0.016603104,0.016624224,0.0166344,0.016688064,0.016727952,0.016755936,0.016907088,0.01705128,0.017086944,0.017090544,0.017144112,0.017250864,0.017252496,0.017258448,0.017425248,0.017461968,0.017505456,0.01758816,0.017600976,0.017607792,0.017613504,0.017641056,0.017641728,0.017651664,0.017730816,0.017820816,0.017866608,0.01799472,0.018064608,0.018072672,0.01807632,0.018138192,0.018233136,0.018236016,0.018287424,0.01841256,0.018536928,0.018557328,0.019344576,0.019683264,0.02037024}; 
		double mySavedSnps[378] =  {0.000542352,0.000595152,0.000674016,0.000758976,0.000783792,0.000804144,0.000819072,0.00082104,0.00084864,0.000868656}; 
		paramData->neut_site.assign(&mySavedSnps[0], &mySavedSnps[0]+10);
	
//	cout<<'\n';
		if(nCarrier==0){
		int n=0;
		vector<vector<Segment> > initSegVec;
/** build the initial carriers
**/
// for some reason, the following code line doesn't work on my ubuntu laptop
// I commented it off. Changde
//		while(!fin.eof()){
		while( getline(fin, line) ){
			++n;
			initSegVec.resize(n);
			
			vector<double> carrLine;
			iss.clear();
			iss.str(line);
			while(iss >> doubtemp){
				carrLine.push_back(doubtemp);
			}
// read in context now. At this order: population, sex, siteA?
			Context initCtxt(static_cast<int> (carrLine[0]),
							 static_cast<int> (carrLine[1]),
							 static_cast<int> (carrLine[2]));
			for(int i=0; i<paramData->neut_site.size(); ++i){
				Segment initSeg;
// read in neutral sites, no test for order, but it has to be ordered from small to large
                double segtmp = paramData->neut_site[i];
				initSeg.L=segtmp; initSeg.R=segtmp;
				initSegVec.at(n-1).push_back(initSeg);
			}

// Make a prototype chromosome:

			shared_ptr< ARGNode > nulle;
			shared_ptr<Chromosome> chr( new Chromosome( initCtxt, initSegVec.at(n-1), nulle ) );
			paramData->initChr.push_back(chr);
			assert(!chr->isEmpty());
		}
	}	
	else{
		vector<double> carrLine;
		vector<Segment> initSegVec;
		
		getline(fin, line);
		iss.clear();
		iss.str(line);
		while(iss >> doubtemp){
			carrLine.push_back(doubtemp);
		}
		Context initCtxt(static_cast<int> (carrLine[0]),
						 static_cast<int> (carrLine[1]),
						 static_cast<int> (carrLine[2]));
			for(int i=0; i<paramData->neut_site.size(); ++i){
				Segment initSeg;
// read in neutral sites, no test for order, but it has to be ordered from small to large
                double segtmp = paramData->neut_site[i];
				initSeg.L=segtmp; initSeg.R=segtmp;
				initSegVec.push_back(initSeg);
			}
		
		// Make a prototype chromosome:
		//
		for(int j=0; j<nCarrier; ++j){
			
			shared_ptr< ARGNode > nulle;
			shared_ptr<Chromosome> chr( new Chromosome( initCtxt, initSegVec, nulle ) );
			paramData->initChr.push_back(chr);
		}
	}
	

	
	fin.close();

	}	
	




Parameters::~Parameters() {		// Destructor
	DBG("Parameters:: Destructing...")
}


vector<int> Parameters::getPopulationMaxSizes(){	// Returns vector of population sizes (over populations)
	return paramData->popSizeVec;
}



vector<shared_ptr<Chromosome> > Parameters::getChromVec(){	// Returns the prototype chromosome for the initial carriers
	return paramData->initChr;
}

shared_ptr<Parameters::ParameterData> Parameters::getpData(){
	return paramData;
}

void Parameters::setNsites(vector<double> myNsites){
		paramData->neut_site = myNsites;
		
		
		// Make a prototype chromosome:
		//
		
			
			vector<Segment> initSegVec;
			for(int k=0; k< myNsites.size(); ++k){
				Segment initSeg;
				initSeg.L=myNsites.at(k); initSeg.R=myNsites.at(k);
				initSegVec.push_back(initSeg);
			}
		for(int j=0; j<paramData->initChr.size(); ++j){
			paramData->initChr.at(j)->setSegmentVector(initSegVec);
		}//j
			
			


}
