#ifndef PARAMETERS_H_
#define PARAMETERS_H_

using std::vector;
// Includes from Boost:
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
#include <boost/scoped_ptr.hpp>
using boost::scoped_ptr;

// Includes from our files:
#include "typedefs.h"


// Forward declaration:
class Chromosome;

class Parameters {
	
public:
	// private:
	struct ParameterData;
	shared_ptr<ParameterData> paramData;
	
	Parameters();
	Parameters(const char* insstring);
	~Parameters();
	vector<int> getPopulationMaxSizes();
	vector<shared_ptr<Chromosome> > getChromVec();
	shared_ptr<ParameterData> getpData();
	void setNsites(vector<double>);
	vector<double> getScalingFactor();
	vector<double> getEpoBreakpoints();	
	
};

struct Parameters::ParameterData{
	//
	// Private data for Parameters (accessed by an opaque pointer)
	//
	vector<int> popSizeVec;						// Vector with population sizes (numbers of chromosomes) in each population
	vector<double> mig_prob;
	int nRuns;
	double mutation_rate;
 double recombination_rate; 
	vector<double> ages;
	double maleRecRatio;
	double freqA_inx;
	double freqA_iny;
 double preY_afreq;
	double Sexsite_pos;
	double siteA_pos;
 int initNchrs;
 int initNx;
 double myEpsi;
	int nEpochs; //Number of epochs
	vector<double> epoch_breaks; // Vector with the breakpoints between epochs, with timeline looking like: [0,t_1,...,1)
	vector<double> scaling_factors; // Vector with scaling pop sizes for t = 0,...,n
	int avg;
	vector<shared_ptr<Chromosome> >initChr;		// vector of initial carriers
	vector<double> neut_site;			// the original version of Nsites take a vector of doubles; ccd 10/10/2014

};


#endif /*PARAMETERS_H_*/


