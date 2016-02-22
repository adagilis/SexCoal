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
	
};

struct Parameters::ParameterData{
	//
	// Private data for Parameters (accessed by an opaque pointer)
	//
	vector<int> popSizeVec;						// Vector with population sizes (numbers of chromosomes) in each population
	vector<double> mig_prob;
	int nRuns;
	double mutation_rate;
	vector<double> ages;
	double maleRecRatio;
	double malePopRatio;// sex ratio
	double maleMutRatio;	
	double freqA_inx;
	double freqA_iny;
	double Sexsite_pos;
    double preY_afreq;
	double siteA_pos;
	int avg;
	vector<shared_ptr<Chromosome> >initChr;		// vector of initial carriers
	vector<double> neut_site;			// the original version of Nsites take a vector of doubles; ccd 10/10/2014

};


#endif /*PARAMETERS_H_*/


