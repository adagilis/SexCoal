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
 * Regional Mutation rate. Changde.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->mutation_rate=doubtemp;
	}
	
/**
 * Regional Recombination rate. Changde.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->recombination_rate=doubtemp;
  
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

// the number of epochs

	getline( fin, line);
	iss.clear();
	iss.str(line);

	while(iss >> temp){
		paramData->nEpochs = temp;
	}

// the list of breaks between epochs

	getline( fin, line);
	iss.clear();
	iss.str(line);
//	cout<<"The timeline is: [0, ";
	while(iss >> doubtemp){
		paramData->epoch_breaks.push_back(doubtemp);
//		cout<<" "<<doubtemp<<",";
	}
//	cout<<" Inf)"<<'\n';

// the list of population scaling factors for each epoch

	getline( fin, line);
	iss.clear();
	iss.str(line);
//	cout<<"The relative population sizes are: ";
	while(iss >> doubtemp){
		paramData->scaling_factors.push_back(doubtemp);
//		cout<<" "<<doubtemp;
	}
//	cout<<'\n';


/**
 * Ages of phylogenetic nodes. Assumed age of siteA first, Y second
 * in number of generations, starting by the youngest (EX: 1500 3000)
 * Changde.
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
		paramData->Sexsite_pos =doubtemp / (2* (double)paramData->popSizeVec.at(0) * paramData->recombination_rate);

	}

	


/**
 * read in siteA position in Rho.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->siteA_pos =doubtemp / (2* (double)paramData->popSizeVec.at(0) * paramData->recombination_rate);
	}

	
	

/**
 * Male to Female recombination rate.
 * Haven't used it by Changde.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->maleRecRatio =doubtemp;
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

	

/**
 * read in freqA in Y (frequency of sex antagonisitic site on Y, range [0,1]).
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->freqA_iny=doubtemp;
	}
	
	

/**
 * read in freqA in pre-Y-sweep. Not used by Changde yet.
**/
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->preY_afreq=doubtemp;
	}
	


/**
 * Make the vector of segments from that single segment:
**/


/**
 * read in the radius of SA and SDR sites
**/
 
	getline( fin, line);
	iss.clear();
	iss.str(line);
	while(iss >> doubtemp){
		paramData->myEpsi =doubtemp;
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

	getline( fin, line);
	iss.clear();
	iss.str(line);

	while(iss >> doubtemp){
		paramData->neut_site.push_back(doubtemp / (2* (double)paramData->popSizeVec.at(0)));

	}

 
  int n=0;
  int nx=0;

		vector<vector<Segment> > initSegVec;
// build the initial carriers
// Changde
		while( getline(fin, line) ){
			++n;
   initSegVec.resize(n);
			
			iss.clear();
			iss.str(line);
   vector<double> carrLine;
			while(iss >> doubtemp){
				carrLine.push_back(doubtemp);
   }
// read in context now. At this order: population, sex, siteA
			Context initCtxt(static_cast<int> (carrLine[0]),
							 static_cast<int> (carrLine[1]),
							 static_cast<int> (carrLine[2]));
			if(carrLine[1] == 0){ nx++; };
   
   Segment initSeg;
			initSeg.L=paramData->myEpsi; 
			initSeg.R=paramData->siteA_pos-paramData->myEpsi;
			initSegVec.at(n-1).push_back(initSeg);
   initSeg.L=paramData->siteA_pos+paramData->myEpsi;
			initSeg.R=1;
   initSegVec.at(n-1).push_back(initSeg);
   			
// Make a prototype chromosome:
			shared_ptr< ARGNode > nulle;
			shared_ptr<Chromosome> chr( new Chromosome( initCtxt, initSegVec.at(n-1), nulle ) );
			paramData->initChr.push_back(chr);
			assert(!chr->isEmpty());
		}
  
  paramData->initNchrs = n;
  paramData->initNx = nx;

	fin.close();

}	


Parameters::~Parameters() {		// Destructor
	DBG("Parameters:: Destructing...")
}


vector<int> Parameters::getPopulationMaxSizes(){	// Returns vector of population sizes (over populations)
	return paramData->popSizeVec;
}

vector<double> Parameters::getScalingFactor(){		// Returns vector of scaling factors.
	return paramData->scaling_factors;
}

vector<double> Parameters::getEpoBreakpoints(){		// Returns vector of epoch breakpoints.
	return paramData->epoch_breaks;
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
