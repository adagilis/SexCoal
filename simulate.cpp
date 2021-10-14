// simulate.cpp
//
//


// Includes from STL:

//#define RCDBGR
//#define SWEEPDBGR

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <algorithm>
using std::max;
using std::min;

// Includes from Boost:
#include <boost/foreach.hpp>
#include <boost/math/special_functions/binomial.hpp>
using boost::math::binomial_coefficient;

// Includes from our files:
#include "chromosome.h"
#include "ran_mk.h"
#include "world.h"
#include "argnode.h"
#include "sitenode.h"

//	Implimentation notes:

//		Recombination is done with interference:  there's a maximum of one recombination event per chromosome.


int World::simulateGeneration(double recombinationRate, vector < vector <double> > & mig_prob){
    //
    // this function calculates the rates of events per generation,
    // then determines which event happened
    //
    //
    
	vector <vector<double> > migmap;
	
	vector<double> mRate;
	vector<double> cRate;
	vector<double> rRate;
	vector<double> uRate;
	
	double sexToARange =worldData->siteA_pos - worldData->Sexsite_pos;
	if( sexToARange < 0 ){sexToARange = -sexToARange;}
    
	double totalM=0;
	double totalC=0;
	double totalR=0;
  	double totalU=0;
    
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
		
		Context cxt =i->first;
		unsigned int k=ctxtNcarriers(cxt);
//changde////////////////////////////////////////////////////////////////////////////////////////////
//cout << k << "\t";
//changde///////////////////////////////////////////////////////////////////////////////////////////


		if(k>0){
            
            UINT clustID=i->second;
            int pop= i->first.pop;
            
			double	clustSize = worldData->freq.at(clustID)* (double) worldData->popSize.at(pop);
            
			// calculates migration rates for each context
            
            vector<double> mfreq = getFreqbyPop(cxt);
            vector<double> migs;
            double b_mig_total=0;
            for(int l=0; l<worldData->nPops;++l){
                if(l!=pop){
                    b_mig_total += mig_prob.at(pop).at(l) * mfreq.at(l)/mfreq.at(pop);
                    migs.push_back(b_mig_total);
                }
                else migs.push_back(0);
            }
            migmap.push_back(migs);
            
            mRate.push_back(k * b_mig_total);	// rate is number of carriers in context, times probability of migration in population
            totalM += mRate[clustID];		//
            
            
            if(k > 1){														// calculates coalescence rates for each context
                unsigned int a=2;
                double ka= binomial_coefficient<double>(k,a)/(clustSize*worldData->scaling_factors.at(worldData->current_epoch)); // Modified for epoch pop size
                cRate.push_back( ka );
                totalC += ka;
                
            }
            
            
            for(int carrierID = 0; carrierID < k; ++carrierID){				// calculates recombination rates for each carrier
                
                double seg = worldData->scaling_factors.at(worldData->current_epoch)*sexToARange *recombinationRate * worldData->totalRecPairs.at(clustID); // Modified for epoch pop size
                
                rRate.push_back( seg  );
                totalR += seg;
				
            }
            
            for(int carrierID = 0; carrierID < k; ++carrierID){				// calculates recombination rates for each carrier
                
                double totLength= worldData->carriers->at(clustID).at(carrierID)->getTotalLength(worldData->Sexsite_pos,worldData->siteA_pos);
                
                double outLength=totLength-sexToARange;
                double seg = worldData->scaling_factors.at(worldData->current_epoch)*outLength *recombinationRate * worldData->totalNoRecPairs.at(clustID); // Modified for epoch pop size
                
                uRate.push_back( seg  );
                totalU += seg;
				
            }
        }
		else{vector<double> dum (1,0); migmap.push_back(dum);mRate.push_back(0);}
	}
	///////////////////
	int nEvents=0;
	double Rate= totalM+ totalC+ totalR+totalU;
   
//changde////////////////////////////////////////////////////////////////////////////////
//cout << totalM << "\t" << totalC << "\t" << totalR << "\t" << totalU << "\t";
//changde////////////////////////////////////////////////////////////////////////////////
 
    double waiting_t;
	if(Rate > 0) waiting_t= randexp(Rate);
    else if (Rate==0) waiting_t = worldData->ages[0]; //	
	double kingman= randreal(0,1);
	
	// Here, we check if the waiting time would advance us to the next epoch. If so, we reset waiting_t to the distance to next epoch, and then force no events to happen in the time period.
	int epoSwitch = 0;
	if(checkEpoch(waiting_t) == 1){
		epoSwitch = 1;
		waiting_t = worldData->epoch_breaks.at(worldData->current_epoch) - (worldData->generation);
		totalM = 0;
		totalC = 0;
		totalR = 0;
		totalU = 0;
		//cout<<"Epoch Changed\n";
		worldData->current_epoch += 1;
	}

	switch (timePeriod(waiting_t)) {
		case 0:
            //  kingman coalescent
			
			worldData->generation += waiting_t;
//changde's print out////////////////////////////////////////////////////////////////////////
//cout << worldData->generation << "\t" << waiting_t << "\t";
//changde///////////////////////////////////////////////////////////////////////////////////
            
            if (kingman < totalM/Rate) nEvents+= migrateEvent(migmap, mRate, totalM);
            else if (kingman < (totalM+totalR)/Rate) nEvents+=recombineEvent(rRate, totalR,true);
            else if (kingman < (totalM+totalR+totalU)/Rate) nEvents+=recombineEvent(uRate, totalU,false);
            else nEvents+=coalesceEvent(cRate, totalC);
            
			break;
			
            
		case 1:
            // Age of the A1/A2 polymorphism. All carriers of A2 will coalesce and become YA1
            worldData->generation = worldData->ages.at(0)+1;
            mergeContext(false);
 			worldData->ageCheck=true;
            siteAfreqUpdate();
            calcTotalRecPairs();calcTotalNoRecPairs();
			
    
			break;
			
		case 2:
            // Age of the Y chromosome.
            // All carriers with Y will coalesce and the ancestor carrier will be in context YA1.
            // site A is still polymorphic, and the X chromosomes act as autosomes
            worldData->generation = worldData->ages.at(0)+1;
            SDBG("Before mergeContext:---------------------------------")
            for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
                Context cxt =i->first;
                unsigned int k=ctxtNcarriers(cxt);
                if(k>0){
                    for(int carrierID = 0; carrierID < k; ++carrierID){
                        shared_ptr<Chromosome> chrom = worldData->carriers->at(cluster[cxt]).at(carrierID);
                        SDBG("Carrier in vector "<<cxt.pop<<cxt.Sexsite<<cxt.siteA<<" has context "<<chrom->getContext().pop<<chrom->getContext().Sexsite<<chrom->getContext().siteA)
                    };
                };};

			mergeContext(true);
			worldData->ageCheck=true;
            
            SDBG("Before freqUpdate:")
            for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
                Context cxt =i->first;
                unsigned int k=ctxtNcarriers(cxt);
                if(k>0){
                    for(int carrierID = 0; carrierID < k; ++carrierID){
                        shared_ptr<Chromosome> chrom = worldData->carriers->at(cluster[cxt]).at(carrierID);
                        SDBG("Carrier in vector "<<cxt.pop<<cxt.Sexsite<<cxt.siteA<<" has context "<<chrom->getContext().pop<<chrom->getContext().Sexsite<<chrom->getContext().siteA)
                    };
                };};
            
            yAgefreqUpdate();
            calcTotalRecPairs(); calcTotalNoRecPairs();
            
            SDBG("After freqUpdate:")
            for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
                Context cxt =i->first;
                unsigned int k=ctxtNcarriers(cxt);
                if(k>0){
                    for(int carrierID = 0; carrierID < k; ++carrierID){
                        shared_ptr<Chromosome> chrom = worldData->carriers->at(cluster[cxt]).at(carrierID);
                        SDBG("Carrier in vector "<<cxt.pop<<cxt.Sexsite<<cxt.siteA<<" has context "<<chrom->getContext().pop<<chrom->getContext().Sexsite<<chrom->getContext().siteA)
                    };
                };};

			break;
			
			            
		case 3:
            // Age of the pop expansion. 
            worldData->generation = worldData->ages.at(0)+1;
            worldData->popSize.at(0) = 14600; // old pop size
            worldData->ageCheck=true;
            calcTotalRecPairs();calcTotalNoRecPairs();
			
    
			break;
			
		default:
			cout << "Error in Simulate(): switch returns default\n";
			break;
	}
    
  	return nEvents;
}

int World::checkEpoch(double waiting){
	int epoSwitch = 0;
	double t = waiting + worldData->generation;
	double nextBreak = t+1;
	if((worldData->nEpochs-1) != (worldData->current_epoch)){
		nextBreak = worldData->epoch_breaks.at(worldData->current_epoch);
	}
	if( t > nextBreak){
		epoSwitch = 1;
	}

	return epoSwitch;
}

int World::timePeriod(double waiting){
	int toSwitch;
	double t = waiting + worldData->generation;
    
	if(   t >= worldData->ages.at(0) && worldData->ageCheck==false)
	{
            if(worldData->ages.at(1)==0) toSwitch=1;
	        else if(worldData->ages.at(1)==1) toSwitch=2;
	        else if(worldData->ages.at(1)==2) toSwitch=3;
	        else {cout<<"Error in TimePeriod(). Exiting\n"; exit(1);}
	   
    }
	else
	{ toSwitch=0;}
	
	return toSwitch;
}


int World::migrateEvent(vector < vector< double> >& mig_prob, vector<double>& rate, double total ){
	//
	// Move one carrier from its current population to a new one.
	// The function figures out in which population is the migration event, who migrates and where to.
	//
	// Each population has a rate of migration equal to ri=(Ci * Mi) ("number of carriers in that context" times "probability of migration for one carrier in that context")
	// The function calculates these rates, and divides them by the sum of them through all contexts (scales to 1).
	//	A cumulative limit is calculated, to be used then with a random [0,1).
	//
	//			0				r1			r1+r2		r1+r2+r3	  1
	//			|---------------]-------------]-------------]---------]
	//					^				^				^			^
	//			[migration in 1]	[in 2]			[in 3]		[in 4]
	//
	//
	//	This assignment is done with a map< KEY_type double, VALUE_type int> . the map makes KEY/VALUE pairs. We use the lower_bound(x) function, that
	//	returns an iterator to the first element in the map whose key does not compare less than x, i.e. it is either equal or greater.
	//
	//
	//  Modified DEC-10 by RG
	//	DBG("Migration happened in gen "<< nGenerations())
	
	// Initialize things:
	//
	
	map<double, int> whichClust;
	
	double add=0;
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {			// this loop scales the rates to one, and constructs the map
		if(ctxtNcarriers(i->first)!=0){													// used to determine in which context migration happened
			whichClust[(rate[i->second] + add) / total] = i->second;
			add += rate[i->second];
		}
	}
	
	int c = whichClust.lower_bound(randreal(0,1))->second;				// get context where mig happened, using a random number on the map
	int who= randint(0, worldData->carriers->at(c).size()-1);			// random number that will be the index of the migrant carrier
	
	shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(who);	// Get the migrant carrier
	
	
	vector<shared_ptr<Chromosome> >::iterator pos = worldData->carriers->at(c).begin() + who;
	worldData->carriers->at(c).erase(pos);
	
	//THIS FIX ASSUMES TWO POPS ONLY!!! MUST CHANGE
	int whereto=0; if( chrom->getContext().pop==0) whereto=1;
	
	
	chrom->setPopulation(whereto);													// change the context of this carrier to its new population
	
	int newC= cluster[chrom->getContext()];
	worldData->carriers->at(newC).push_back( chrom );								// add the carrier to the vector of its new population
	
	shared_ptr < ARGNode > newNode ;
	newNode.reset(new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation));
	worldData->argNodeVec.push_back(newNode);
	chrom->setDescendant(newNode);
	
	
	return 1;
}	//********************************************* end of migrateEvent

int World::coalesceEvent(vector<double>& rate, double total)
{
//	This function executes one simple coalescent event.
// 	Written by Rafael Guerrero, VI-09
//	Modified by ccd, 2014 Nov.
//	cout << "Simple coalescent occurred in gen "<< nGenerations() << endl;
	DBG("Simple coalescent occurred in gen "<< nGenerations())

	map<double, int> whichClust;
	int id=0;
	double add=0;

	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i )
	{// this loop builds the map of cumulative fractions of coalescence (binomial coeff/ total)
		COALDBG("carriers in "<< i->second<<" "<<ctxtNcarriers(i->first))
		if(ctxtNcarriers(i->first) > 1)
		{
			whichClust[(rate.at(id) + add) / total] = i->second;
			COALDBG("rate "<<(rate.at(id) + add))
			COALDBG(total)
			COALDBG("whichClust[] "<<whichClust[(rate.at(id) + add) / total] )
			add += rate.at(id);
			++id;
		}
	}
	
	if(id==0) // no clusters have more than 1 carrier.
	{
		COALDBG("no coalescence possible, all clusters have < 2 carriers")
		return 0;
	}	
	else
	{// carry out the rest of the function.
		int c = whichClust.lower_bound(randreal(0,1))->second;
		// randomly decide in which context is coalescence happening
		vector<int> carr_idx;
		for (int u=0; u< worldData->carriers->at(c).size(); ++u) carr_idx.push_back(u);
		random_shuffle(carr_idx.begin(), carr_idx.end());
		
		shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(carr_idx.at(0));
		shared_ptr<Chromosome> chrom2 = worldData->carriers->at(c).at(carr_idx.at(1));
		
		shared_ptr < ARGNode > newNode;
		newNode.reset(  new ARGNode(worldData->argNodeVec.size(), chrom, chrom2, worldData->generation) );
		
		worldData->argNodeVec.push_back(newNode);
		// Add a pointer to this node to the vector of all nodes in	the simulation
//cout << "total segs before " << worldData->LeftTotalSegSize <<" "<< chrom->getTotalSegLength() <<" "<< chrom2->getTotalSegLength() << endl;	
		double myBefore = chrom->getTotalSegLength() + chrom2->getTotalSegLength();
			
		chrom->merge( chrom2 , newNode); // merging includes setting decendant to new ARGNode.
					
		worldData->LeftTotalSegSize -= (myBefore - chrom->getTotalSegLength() );
//cout << "total segs after " << worldData->LeftTotalSegSize  <<" "<< chrom->getTotalSegLength() <<" "<< chrom2->getTotalSegLength() << endl;	
		worldData->carriers->at(c).erase(remove(worldData->carriers->at(c).begin(), worldData->carriers->at(c).end(), chrom2), worldData->carriers->at(c).end());
		DBG(totalNCarriers()<<" in "<<worldData->generation<<"\n")
/***********************************************************************/	

vector<double>::size_type outit = 0;
while ( outit < worldData->recomb_breakpoints.size() )
{
	int myCopy = 0;
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i )
	{
		UINT clustID=i->second;
		int nLocalCarriers  = worldData->carriers->at(clustID).size();
		
		for(int carrierID = 0; carrierID < nLocalCarriers; ++carrierID)
		{
			if( siteQuery( worldData->recomb_breakpoints.at( outit ), worldData->carriers->at( clustID ).at( carrierID )->getSegmentVector()  ) )
				myCopy ++;
		}
	}
	if( myCopy == 1 )
	{
		std::pair <double, shared_ptr < ARGNode > > myFinishedSegment;
		myFinishedSegment = std::make_pair ( worldData->recomb_breakpoints.at( outit ), newNode );
		worldData->outputSegments.push_back( myFinishedSegment );
		worldData->recomb_breakpoints.erase( worldData->recomb_breakpoints.begin() + outit );
	}else
	{
		++outit;
	}
//	cout << myCopy << "\t";
//if myCopy == 1, remove breakpoint from list, and put it into outlist with the argnode
}
//cout << endl;

/*	
newNode->outputNode();//ccd debug
for(UINT j=0; j<chrom->getSegmentVector().size(); j++)
cout << "  (" << chrom->getSegmentVector().at(j).L << ", " << chrom->getSegmentVector().at(j).R << ")";
cout << endl;
***********************************************************************/		
		return 1;
	}
}//************************************end of coalesceEvent()

int World::recombineEvent(vector<double>& rate, double total, bool insideSA){
	
	/*
	 // this function runs one recombination event between a carrier and a non-carrier
	 // It determines in which context, and who recombines.
	 // It adds the resulting new chromosome to the vector of carriers.
	 //	It also checks for "empty" carriers after recombination.
	 //	In this case, the carriers are not included in the Carriers vector.
	 //
	 //
	 // Given that recombination happened, each carrier has a probability of being the recombinant (Ri) equal to
	 //	the total length of its segments divided by the total length of segments in the world.
	 //
	 //	A cumulative limit is calculated, to be used then with a random [0,1).
	 //
	 // For example, in a case with four carriers:
	 //
	 //			0				R1			R1+R2		R1+R2+R3	  1
	 //			|---------------]-------------]-------------]---------]
	 //					^				^				^			^
	 //			[recombine 1]		[ 2]			[ 3]		[ 4]
	 //
	 //
	 //	This assignment is done with a map< KEY_type double, VALUE_type pair< UINT, int> > . the map makes KEY/VALUE pairs.
	 //	In the example, we would have a map with the pairs < R1 , 1 > , < R1+R2, 2 >, < R1+R2+R3, 3 >, < 1== R1+R2+R3+R4, 4 >
	 //	We use the lower_bound(x) function, that returns an iterator to the first element in the map whose key does not
	 //	compare less than x, i.e. it is either equal or greater.
	 //	In this particular map, the VALUE will be a pair <UINT, int>, the complete "address" of the carrier.
	 //	This pair is:
	 //	- the cluster[context] ( number of the context in the map of all contexts )
	 //	- the postion of the carrier in the vector of carriers in that context.
	 
	 */
//	cout << "Recombination happened in gen "<<nGenerations() << endl;
	DBG( "Recombination happened in gen "<<nGenerations())
	
	map<double, pair<UINT,int> > whichCarrier;
	
	int id=0;
	double add=0;
	for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
		UINT clustID=i->second;
		int nLocalCarriers = worldData->carriers->at(clustID).size();
		
		for(int carrierID = 0; carrierID < nLocalCarriers; ++carrierID){
			if (rate[id]!=0) whichCarrier[(rate[id] + add) / total] = make_pair(clustID, carrierID) ;
			add += rate[id];
			++id;
		}
	}
	
	pair<UINT,int> who = whichCarrier.lower_bound(randreal(0,1))->second;
	// randomly decide which carrier will recombine, using the map.
	// "who" is a pair
	// who->first is the number of the context
	// who->second is the number of the carrier in that context
	
	shared_ptr<Chromosome> chrom = worldData->carriers->at(who.first).at(who.second);			// Get the recombining carrier
    //the follosing three are the differences between inside/outside the SA region
	double carrierSum;
    pair<Context, Context> ancestors;
    double breakpoint;
    double sexToARange =worldData->siteA_pos - worldData->Sexsite_pos ;
	if( sexToARange < 0) sexToARange = -sexToARange;
    double outSegment= chrom->getTotalLength(worldData->Sexsite_pos,  worldData->siteA_pos)-sexToARange;
    
    if(insideSA){
        carrierSum = worldData->totalRecPairs.at(who.first);
        ancestors = worldData->parents.at(who.first).lower_bound(randreal(0, carrierSum))->second;
		breakpoint = randreal(0,sexToARange) + min( worldData->Sexsite_pos, worldData->siteA_pos);//compute the breakpoint; depending SDR or siteA which is smaller
// modified by changde
	}
    else{
        carrierSum = worldData->totalNoRecPairs.at(who.first);//changde
        ancestors = worldData->NoRecParents.at(who.first).lower_bound(randreal(0, carrierSum))->second;//changde
        breakpoint = chrom->calcOutBreakpoint(worldData->Sexsite_pos, worldData->siteA_pos,outSegment);	// Compute the breakpoint
	}

// to consider the case of SDR > siteA and Nsite < SDR or siteA, changde.
	if(breakpoint < worldData->Sexsite_pos){
		pair<Context, Context> ancestors_tmp = ancestors; ancestors.first = ancestors_tmp.second; ancestors.second = ancestors_tmp.first;
	}

	RCDBG("Recomb in carrier of length: "<<chrom->getTotalLength(worldData->Sexsite_pos, worldData->siteA_pos))
	RCDBG("Ancestors (pop-sex-siteA) = ("<<ancestors.first.pop<<" "<<ancestors.first.Sexsite<<" "<<ancestors.first.siteA<<") and ("<<ancestors.second.pop<<" "<<ancestors.second.Sexsite<<" "<<ancestors.second.siteA<<")")
	
	// Make a new ARG node, add it to the vector of nodes.
	shared_ptr < ARGNode > newNode;
	newNode.reset(  new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation) );
	worldData->argNodeVec.push_back(newNode);
	
	// Break the chromosome into two:
	int inSegRec = 0;//indicator of inside segment recombination
	shared_ptr<Chromosome> chrom2 = chrom->recombine(breakpoint, newNode, ancestors, inSegRec);			// Make two recombinant chromosomes
	
	vector<shared_ptr<Chromosome> >::iterator pos = worldData->carriers->at(who.first).begin() + who.second;	// an iterator to the position of chrom. needed only for the function erase()
	worldData->carriers->at(who.first).erase(pos); //erase chrom from old context
	
	int myNchrs = 0;

	if(!chrom->isEmpty()){
		worldData->carriers->at(cluster[chrom->getContext()]).push_back(chrom);	 //add chrom to new context
		myNchrs++;
	}
	
	if(!chrom2->isEmpty())
	{
		worldData->carriers->at(cluster[chrom2->getContext()]).push_back(chrom2);	//	add chrom2 only if it holds segments
		myNchrs++;
	}
	
	//ccd
	if(inSegRec == 1){
		int myCopy = 0;
		for( cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i ) {
			UINT clustID=i->second;
			int nLocalCarriers  = worldData->carriers->at(clustID).size();
		
			for(int carrierID = 0; carrierID < nLocalCarriers; ++carrierID){
				if( siteQuery( breakpoint, worldData->carriers->at( clustID ).at( carrierID )->getSegmentVector()  ) )
					myCopy ++;
			}
		}
		if( myCopy > 1 ){
			worldData->recomb_breakpoints.push_back( breakpoint );
		}
	}
	
	DBG(totalNCarriers()<<" in "<<worldData->generation<<"\n")
	return 1;
	
}
void World::mergeContext(bool siteORsex){
    
    // This function will force coalescence of carriers, representing a selective sweep of siteA or Y chromosomes (depending on the siteORsex argument
    // It will put the ancestral carrier in an originContext hard-coded here 
    
    // To merge siteA: siteORsex = false
    // To merge sexsite: siteORsex=true
    
    
	vector<shared_ptr<Chromosome> > sons;
    vector<Context> targets(2);

	Context originCtx(0,1,0); //Context {0,Y,A1}
    targets[0].setNew(0,1,1); //context to merge {0,Y,A2}
    targets[1].setNew(0,0,1); //context to merge {0,X,A2}
    
    if(siteORsex){
        originCtx.setNew(0,0,1);
        targets[0].setNew(0,1,0);
        targets[1].setNew(0,1,1);
    };
    
  
    //First, we collect all the carriers in the sons vector
    for( int iter = 0; iter < targets.size(); ++iter ) {
        Context ctx= targets[iter];
        int n=ctxtNcarriers(ctx);
        if ( n >0 ){
            for(int j=0; j<n; ++j){
                shared_ptr<Chromosome> chr=worldData->carriers->at(cluster[ctx]).at(j);
                chr->setContext(originCtx);
                sons.push_back(chr);
            }
            worldData->carriers->at(cluster[ctx]).clear();
        }
    }
    
	//now we coalesce all the carriers in the sons vector, into a single carrier "cain", later added to the origin context
	if(sons.size()>0){
		shared_ptr<Chromosome> cain=sons.at(0);
		cain->setContext(originCtx);
        
		shared_ptr < ARGNode > cainNode ;
		cainNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, worldData->generation) );
		
		for(int i=1; i<sons.size(); ++i){
			shared_ptr<Chromosome> abel=sons.at(i);
			shared_ptr < ARGNode > newNode ;
			newNode.reset(  new ARGNode(worldData->argNodeVec.size(), cain, abel, worldData->generation) );
			worldData->argNodeVec.push_back(newNode);		// Add a pointer to this node to the vector of all nodes in	the simulation
			cain->merge( abel, newNode);					// merging includes setting decendant to new ARGNode.
		}
		worldData->carriers->at(cluster[originCtx]).push_back(cain);
	}
    
}

