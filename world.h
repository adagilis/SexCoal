/* world.h
 *
 * Data and methods to simulate evolution backwards.
 *
 *
 *
 */


#ifndef WORLD_
#define WORLD_


// Includes from Boost:
	#include <boost/scoped_ptr.hpp>
		using boost::scoped_ptr;
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr;
	#include <boost/weak_ptr.hpp>
		using boost::weak_ptr;

// Includes from STL:
	#include <vector>
		using std::vector;
	#include <map>
		using std::map;
#include <math.h>
using std::pair;
using std::make_pair;
		
// Includes from our files:
	#include "typedefs.h"
	#include "argnode.h"
	#include "parameters.h"

// Forward declaration:
	class Chromosome;
	
	
// ==================================================================
class World {
 private:
	struct WorldData;									// This structure contains all the private data
	scoped_ptr<WorldData> worldData;					// Pointer to private data
	bool alreadyCoalesced(int parentID);				// Determins if parentID is already on the list of parent carriers that are coalescents
	cluster_t cluster;
	
public:
	World(shared_ptr<Parameters::ParameterData> p);		
	~World();
	void makecluster (int nPops);
	bool simulationFinished();				// Returns true if the simulation is finished	
	double nGenerations();
	void Generation_pp();
	void calcTotalRecPairs();
    void calcTotalNoRecPairs();
	UINT totalNCarriers();					// Returns total number of carriers in the world
	int ctxtNcarriers(const Context ctxt); //	total number of carriers in a given context
	UINT siteNcarriers(int siteA);			//  total number of carriers with a given allelic state in site A
	UINT popNcarriers(int pop);				//	total number of carriers in a given population 
	int migrateEvent(vector < vector< double> >& mig_prob, vector<double>& rate, double total);
	int coalesceEvent(vector<double>& rate, double total);
	int recombineEvent(vector<double>& rate, double total, bool insideSA);
	int simulateGeneration(double recombinationRate, vector < vector <double> > & mig_prob);
	int migratePoisson(vector < vector< double> >& mig_prob);
	int coalescePoisson();
	int recombinePoisson();
	int checkEpoch(double waiting);
    	int timePeriod(double waiting);
	void mergeContext(bool siteORsex);
	void siteAfreqUpdate();
    void yAgefreqUpdate();
	vector< shared_ptr < ARGNode > >& getARGVec();
	Context randCtx(int pop);
	Context randCtxBySex(Context c);
	vector<double> getFreqbyPop(Context c);
	//ccd
	vector<double> getRecombBreakpoints();
	vector< pair<double, shared_ptr < ARGNode > > >& getOutputSegments();
//	vector<int> getNonEmptyChrsNumber();
	
};

struct World::WorldData		// Structure with the private data for World, accessed by an opaque pointer
{
	double generation;								// Generation number, starting at 0
	int current_epoch;								// Current epoch, starting at 0
	UINT nPops;										// Number of populations in World
	UINT nClust;
	vector<int> nChromos;	
	vector<int> popSize;
	vector< shared_ptr < ARGNode > > argNodeVec;					// Vector of pointers to all ARG nodes in the simulation
	shared_ptr< vector< vector < shared_ptr<Chromosome> > > > carriers;		
	double Sexsite_pos;
	double freqA_inx;
	double freqA_iny;
    double preY_afreq;
	double siteA_pos;
	double maleRecRatio;
	int nEpochs; //Number of epochs
	vector<double> epoch_breaks; // Vector with the breakpoints between epochs, with timeline looking like: [0,t_1,...,1)
	vector<double> scaling_factors;
	vector<double> freq;
	vector<double> totalRecPairs;
    vector<double> totalNoRecPairs;
	vector< map<double, pair<Context,Context> > >parents;
    vector< map<double, pair<Context,Context> > >NoRecParents;
    vector<double> ages;
    bool ageCheck;
    
    vector<double> recomb_breakpoints;
    vector< pair<double, shared_ptr < ARGNode > > > outputSegments;
//   vector<int> seg_copies;
//    vector< shared_ptr < ARGNode > > segNodeVec;
//   vector<int> non_empty_chrs_number;
    
    double LeftTotalSegSize;
};




#endif /*WORLD_*/

