/**
 * World.cpp
 *
 * This class represents the world with all of its populations.
 *
 *
 */

//#define RCDBGR

// Includes from STL:
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <math.h>

// Includes from Boost:
#include <boost/foreach.hpp>

// Includes from our files:
#include "chromosome.h"
#include "ran_mk.h"
#include "world.h"
#include "argnode.h"

// ================================================================================================================


World::World(shared_ptr<Parameters::ParameterData> p){
	
	// Constructor for World
	//
	
	DBG("Population:: Constructing Population(size, shared_ptr<Chromosome>)...")
    
	int nPops = p->popSizeVec.size();
	makecluster(nPops); // builds a map of all contexts
	
	UINT nClust=cluster.size();
    
	// Wipe clean WorldData:
	//
	worldData.reset(new WorldData);
	worldData->carriers.reset( new vector < vector< shared_ptr<Chromosome> > > );
	worldData->generation = 0;
	worldData->nPops = nPops;
	worldData->nClust = nClust;
	
	worldData->maleRecRatio= p->maleRecRatio;//different recombination rate
	worldData->malePopRatio = p->malePopRatio;// sex ratio
	worldData->myYnumbRatio = 0.5*worldData->malePopRatio/(worldData->malePopRatio + 1);
	worldData->maleMutRatio = p->maleMutRatio;// different mutation rate


	
	worldData->Sexsite_pos=p->Sexsite_pos;
	worldData->siteA_pos=p->siteA_pos;
	
	
	worldData->recomb_breakpoints.push_back(p->neut_site.at(0));//assuming just one site per run
	
	worldData->freqA_inx= p->freqA_inx;
	worldData->freqA_iny= p->freqA_iny;
    worldData->preY_afreq= p->preY_afreq;
    worldData->ages = p->ages;
    worldData->ageCheck=false;
    if(worldData->ages[0]==0) worldData->ageCheck = true;
    
  	// Copy the population sizes and resize vectors to the number of contexts:
	//
	worldData->carriers->resize(nClust);
	worldData->nChromos.resize(nClust);
	worldData->freq.resize(nClust);
    worldData->popSize=p->popSizeVec;
    
	
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
        
		UINT popID=(*iter).first.pop;
		UINT siteA=(*iter).first.siteA;
		UINT Sexsite= (*iter).first.Sexsite;
		double size=0;
		
			if(siteA==1 && Sexsite==1) {
				size=p->freqA_iny *worldData->myYnumbRatio* (double)worldData->popSize.at(popID);
	            //		cout<<"size of ("<<Sexsite<<siteA<<") = "<<size<<'\n';
			}
			else if(siteA==1 && Sexsite==0) {
				size=p->freqA_inx *(1-worldData->myYnumbRatio)* (double)worldData->popSize.at(popID);
	            //	cout<<"size of ("<<Sexsite<<siteA<<") = "<<size<<'\n';
			}
	        
			else if(siteA==0 && Sexsite==1) {
				size= (1-p->freqA_iny) * worldData->myYnumbRatio *(double)worldData->popSize.at(popID);
	            //		cout<<"size of ("<<Sexsite<<siteA<<") = "<<size<<'\n';
	        }
	        
	        else if(siteA==0 && Sexsite==0) {
				size=(1-p->freqA_inx)* (1-worldData->myYnumbRatio) * (double)worldData->popSize.at(popID);
	            //		cout<<"size of ("<<Sexsite<<siteA<<") = "<<size<<'\n';
	        }
	        
	    size = size * 0.25 * (worldData->malePopRatio + 1) * (worldData->malePopRatio + 1) /worldData->malePopRatio;
		worldData->nChromos[iter->second] = size ;
        //		cout<<"nChromos of ("<<Sexsite<<siteA<<") = "<<worldData->nChromos[iter->second]<<'\n';
        
		worldData->freq.at(iter->second) = size/ (double) worldData->popSize.at(popID);
        
	}
    
	calcTotalRecPairs();
    calcTotalNoRecPairs();
	
	shared_ptr< ARGNode > nulle ;
	
	switch (p->avg){
			
		case 0: // just transfers the initCarriers
			for(int i = 0; i < p->initChr.size(); ++i){
				assert(!p->initChr.at(i)->isEmpty());
				shared_ptr<Chromosome> chr= p->initChr.at(i)->copy_values();
				shared_ptr < ARGNode > newNode ;newNode.reset( new ARGNode(worldData->argNodeVec.size(), chr) );// Make a terminal ARG node
				chr->setDescendant(newNode);
				worldData->carriers->at(cluster[chr->getContext()]).push_back( chr );					// Adds this chromosome to the carriers array
				worldData->argNodeVec.push_back(newNode);							// Add the new node to the vector of ARG nodes
			}
			break;
			
		case 1: // makes random carriers in a population
			
			for(int i = 0; i < p->initChr.size(); ++i){
				Context ctxRdm = randCtx(p->initChr.at(i)->getContext().pop);
				shared_ptr<Chromosome> chr (new Chromosome(ctxRdm, p->initChr.at(i)->getSegmentVector(), nulle));
				shared_ptr < ARGNode > newNode ;newNode.reset( new ARGNode(worldData->argNodeVec.size(), chr) );// Make a terminal ARG node
				chr->setDescendant(newNode);
				worldData->carriers->at(cluster[chr->getContext()]).push_back( chr );					// Adds this chromosome to the carriers array
				worldData->argNodeVec.push_back(newNode);							// Add the new node to the vector of ARG nodes
			}
			break;
			
		case 2: // makes random carriers of X or Y in a population
			
			for(int i = 0; i < p->initChr.size(); ++i){
				Context ctxRdm = randCtxBySex(p->initChr.at(i)->getContext());
				shared_ptr<Chromosome> chr (new Chromosome(ctxRdm, p->initChr.at(i)->getSegmentVector(), nulle));
				shared_ptr < ARGNode > newNode ;newNode.reset( new ARGNode(worldData->argNodeVec.size(), chr) );// Make a terminal ARG node
				chr->setDescendant(newNode);
				worldData->carriers->at(cluster[chr->getContext()]).push_back( chr );					// Adds this chromosome to the carriers array
				worldData->argNodeVec.push_back(newNode);							// Add the new node to the vector of ARG nodes
			}
			break;
	}
	
	
} // End constructor for World



World::~World(){
	DBG("Population:: Destructing...")
}



bool World::simulationFinished(){
	
//	if( worldData->LeftTotalSegSize <= 0 || totalNCarriers()==1) return true;
	if (worldData->recomb_breakpoints.size() == 0) return true;
	return false;
}


int World::ctxtNcarriers(const Context ctxt){
	UINT total = worldData->carriers->at(cluster[ctxt]).size();
	return total;
}

UINT World::siteNcarriers(int siteA){
	
	UINT total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (siteA==(*iter).first.siteA)
            total = worldData->carriers->at((*iter).second).size() +total;
	}
	return total;
	
}
UINT World::popNcarriers(int pop){
	
	UINT total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop)
			total = worldData->carriers->at((*iter).second).size() +total;
	}
	return total;
	
}
UINT World::totalNCarriers(){
    //
    // Returns the total number of carriers in the world
    //
	UINT total = 0;
	for(UINT i = 0; i < worldData->nClust; i++)
	{	total = worldData->carriers->at(i).size() +total;
    }
	return total;
}



double World::nGenerations(){
	return worldData->generation;
}

void World::Generation_pp (){
	worldData->generation++;
}

void World::makecluster (int nPops){
	
	int iter=0;
	for (int i=0; i<nPops; ++i){	//loop through populations
		for (int k=0; k < 2; ++k){ //loop through Sex determining sites (0=X, 1=Y)
			for (int j=0; j < 2; ++j){	//loop through alleles in siteA
                
					Context c(i, k, j);
					cluster[c]=iter;
					++iter;
				
			}
		}
	}
	
}

vector< shared_ptr < ARGNode > >& World::getARGVec(){
	return worldData->argNodeVec;
}

vector<double> World::getFreqbyPop(Context c){
	vector<double> freqs (worldData->nPops,0);
	Context geno=c;
	for (int i=0; i<worldData->nPops;++i){
		geno.pop=i;
		freqs.at(i)=worldData->freq.at(cluster[geno]);
	}
	return freqs;
}


vector<double> World::getRecombBreakpoints(){
	return worldData->recomb_breakpoints;
}

double World::getMaleMutRatio(){
		return worldData->maleMutRatio;
}


vector< pair<double, shared_ptr < ARGNode > > >& World::getOutputSegments(){
	return worldData->outputSegments;
}


Context World::randCtx(int pop){
	double draw = randreal(0, 1);
	Context c(0,0,0);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (pop==(*iter).first.pop){
			total+= (double)worldData->freq.at((*iter).second);
			if(draw<total) {
				return (*iter).first;
			}
		}
	}
	cout<<"Error in randCtx\n";
	return c;
	
}

Context World::randCtxBySex(Context c){
    
	double f= 1-worldData->myYnumbRatio;
	if (c.Sexsite==1)f=worldData->myYnumbRatio;
	
	double draw = randreal(0, f);
	double total = 0;
	for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
		if (c.Sexsite ==(*iter).first.Sexsite && c.pop==(*iter).first.pop){
			total+= (double)worldData->freq.at((*iter).second);
			if(draw<total) {
				return (*iter).first;
			}
		}
	}
	cout<<"Error in randCtxBySex\n";
	return c;
}

void World::calcTotalRecPairs (){
    RCDBG("From calcTotalRecPairs")
    
    worldData->totalRecPairs.clear();
	worldData->totalRecPairs.resize(worldData->nClust);
	worldData->parents.clear();
    worldData->parents.resize(worldData->nClust);
    
    bool preY= worldData->ages[1]==1 && worldData->ages[0]<worldData->generation;
    
    for( cluster_t::iterator target = cluster.begin(); target != cluster.end(); ++target ) {
		
		int tpop = target->first.pop;
		int tA = target->first.siteA;
		int tS = target->first.Sexsite;
		double tfreq= worldData->freq.at(target->second);
		RCDBG("Target = "<<tS<<tA)
        double sum =0;
        if(tfreq > 0.0){
            
            Context x(tpop, 0, tA);
            Context y(tpop, 1, tA);
            double aInX = worldData->freq.at(cluster[x])/(1-worldData->myYnumbRatio);
            double aInY = worldData->freq.at(cluster[y])/worldData->myYnumbRatio;
            
            double freqXinM=0.5;
           
            if(preY){freqXinM=1;aInX = worldData->freq.at(cluster[x]);};
            
            double xaInM = freqXinM*aInX;
            double yaInM = (1-freqXinM)*aInY;
            
            double mult=0.5;
            double soo;
            double maleRec =worldData->maleRecRatio;
            
            Context l_carr (0,0,0);
            Context r_carr (0,0,0);
            
            if(tS==1){ //SDR is 1 (Y), cross is always X-Y
                
                l_carr= target->first;
                r_carr= target->first;
                r_carr.Sexsite=0;
                
                double cross1 = mult* yaInM * aInX  * maleRec /(tfreq);
                if(cross1>0){
                    sum +=cross1;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 1 "<<yaInM<<" * "<<aInX<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross1)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                if(tA==1) l_carr.siteA=0; else l_carr.siteA=1;
                
                double cross2 = mult* ((1-freqXinM)-yaInM) * aInX * maleRec/ (tfreq);
                
                if(cross2>0) {
                    sum +=cross2;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 2 "<<((1-freqXinM)-yaInM)<<" * "<<aInX <<" *"<<maleRec<<" /"<<tfreq<<" = "<<cross2)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
            }
            
            
            else{//SDR is 0 (X)
                
                // First, the case when soC is MALE, cross is X-Y
                soo = (double) 1/3; if(preY){soo=0;}
                r_carr= target->first;
                r_carr.Sexsite=1;
                
                l_carr= target->first;
                
                double cross1 = mult*aInX * yaInM  * maleRec / tfreq;
                
                if(cross1>0) {
                    sum +=cross1;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
               
                RCDBG("cross 1 "<<aInX<<" *"<<yaInM<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross1)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                
                if(tA==1) l_carr.siteA=0; else l_carr.siteA=1;
                
                double cross2 = mult*(1-aInX) * yaInM * maleRec / tfreq;
                
                if(cross2>0){
                    sum +=cross2;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 2 "<<(1-aInX)<<" *"<<yaInM<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross2)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                
                //Now SoC is FEMALE, cross is X-X
                soo=(double)2/3; if(preY)soo=1.0;
                r_carr= target->first;
                l_carr= target->first;
                
                double cross3 = mult* aInX * xaInM / tfreq;
                
                if(cross3>0){
                    sum +=cross3;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 3 "<<aInX<<" * "<<xaInM<<" / "<<tfreq<<" = "<<cross3)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                
                if(tA==1) l_carr.siteA=0; else l_carr.siteA=1;
                
                double cross4 = mult* (1-aInX) * xaInM  / tfreq;
                
                if(cross4>0){
                    sum +=cross4;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 4 "<<(1-aInX)<<" * "<<xaInM<<" / "<<tfreq<<" = "<<cross4)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                double cross5 = mult* (freqXinM-xaInM) * aInX / tfreq;
                
                if(cross5>0){
                    sum +=cross5;
                    worldData->parents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 5 "<<(freqXinM-xaInM)<<" * "<<aInX<<" /"<<tfreq<<" = "<<cross5)
                RCDBG("Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
            }
            
            
        }
        worldData->totalRecPairs.at(target->second) = sum;
        
    }
    
}//end

void World::calcTotalNoRecPairs (){
    RCDBG("From calcTotalNoRecPairs:")
    worldData->totalNoRecPairs.clear();
    worldData->totalNoRecPairs.resize(worldData->nChromos.size());
    worldData->NoRecParents.clear();
    worldData->NoRecParents.resize(worldData->nChromos.size());
    bool preY= worldData->ages[1]==1 && worldData->ages[0]<worldData->generation;
    
    for( cluster_t::iterator target = cluster.begin(); target != cluster.end(); ++target ) {
        
        int tpop = target->first.pop;
        int tA = target->first.siteA;
        int tS = target->first.Sexsite;
        double tfreq= worldData->freq.at(target->second);
        
        RCDBG("Target = "<<tS<<tA)
//changde////////////////////////////////////////////////////////////////////////
//cout << "Target = "<<tS<<"\t"<<tA<<"\t";
//changde////////////////////////////////////////////////////////////////////////

        double sum =0;
        if(tfreq>0){
            Context x(tpop, 0, tA);
            Context y(tpop, 1, tA);
            double aInX = worldData->freq.at(cluster[x])/0.75;
            double aInY = worldData->freq.at(cluster[y])/0.25;
            
            double freqXinM=0.5;
            
            if(preY){freqXinM=1;aInX = worldData->freq.at(cluster[x]);}
         
            double xaInM = freqXinM*aInX;
            double yaInM = (1-freqXinM)*aInY;
            
            double mult=0.5;
            double soo;
            double maleRec =worldData->maleRecRatio;
            
            Context l_carr (0,0,0);
            Context r_carr (0,0,0);
            
            if(tS==1){ //SDR is 1 (Y), cross is always X-Y
                
                l_carr= target->first;
                r_carr= target->first;
                r_carr.Sexsite=0;
                
                double cross1 = mult* yaInM * aInX  * maleRec /(tfreq);
                sum +=cross1;
                if(cross1>0)worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                RCDBG("cross 1 "<<yaInM<<" * "<<aInX<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross1)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                
                if(tA==1) r_carr.siteA=0; else r_carr.siteA=1;
                
                double cross2 = mult* yaInM * (1- aInX) * maleRec/ (tfreq);
                sum +=cross2;
                if(cross2>0)worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                RCDBG("cross 2 "<<yaInM<<" * "<<(1- aInX) <<" *"<<maleRec<<" /"<<tfreq<<" = "<<cross2)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                
            }
            
            
            else{//SDR is 0 (X)
                
                // First, the case when soC is MALE, cross is X-Y
                soo = (double) 1/3;if(preY)soo=0;
                
                l_carr= target->first;
                r_carr= target->first;
                r_carr.Sexsite=1;
                
                double cross1 = mult*aInX * yaInM  * maleRec / tfreq;
               
                if(cross1>0){
                     sum +=cross1;
                    worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                RCDBG("cross 1 "<<aInX<<" *"<<yaInM<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross1)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                
                if(tA==1) r_carr.siteA=0; else r_carr.siteA=1;
                
                double cross2 = mult* aInX  *((1-freqXinM)- yaInM )* maleRec / tfreq;
              
                if(cross2>0){
                      sum +=cross2;
                    worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
               
                RCDBG("cross 2 "<<aInX<<" *"<<((1-freqXinM)- yaInM )<<" * "<<maleRec<<" /"<<tfreq<<" = "<<cross2)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                //Now SoC is FEMALE, cross is X-X
                soo=(double)2/3;if(preY)soo=1;
                r_carr= target->first;
                l_carr= target->first;
                
                double cross3 = mult* aInX * xaInM / tfreq;
                                if(cross3>0){
                                    sum +=cross3;
worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                
                RCDBG("cross 3 "<<aInX<<" * "<<xaInM<<" / "<<tfreq<<" = "<<cross3)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                                }
                if(tA==1) r_carr.siteA=0; else r_carr.siteA=1;
                
                double cross4 = mult* aInX * (freqXinM-xaInM ) / tfreq;
                
                if(cross4>0){sum +=cross4;
                    worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                RCDBG("cross 4 "<<aInX<<" * "<<(freqXinM-xaInM )<<" / "<<tfreq<<" = "<<cross4)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
                double cross5 = mult* xaInM * (1-aInX )/ tfreq;
                
                if(cross5>0){
                    sum +=cross5;
               worldData->NoRecParents.at(target->second)[sum]=make_pair(l_carr, r_carr);
                RCDBG("cross 5 "<<soo<<" *"<<(1-aInX)<<" * "<<xaInM<<" /"<<tfreq<<" = "<<cross5)
                RCDBG("								Parents:  L= "<<l_carr.Sexsite<<l_carr.siteA<<" R= "<<r_carr.Sexsite<<r_carr.siteA)
                }
            }
            
            
            RCDBG(sum)
        }
        worldData->totalNoRecPairs.at(target->second) = sum;

    }
}//end

void World::siteAfreqUpdate(){
    // This function will set to zero the sizes (and frequencies) of Contexts with the allele A2 (site=1)
    // It leads to neutral sex chromosome frequencies (Y=1/4, X=3/4)
    
    for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
        
        UINT popID=(*iter).first.pop;
        UINT siteA=(*iter).first.siteA;
        UINT Sexsite= (*iter).first.Sexsite;
        double size;
        
        if(siteA==1 ) {
            size=0;
        }
        else if(siteA==0 && Sexsite==1) {
            size=  worldData->myYnumbRatio *(double)worldData->popSize.at(popID);
        }
        
        else if(siteA==0 && Sexsite==0) {
            size= worldData->myYnumbRatio * (double)worldData->popSize.at(popID);
        }
        
        worldData->nChromos[iter->second] = size ;
        worldData->freq.at(iter->second) = size / (double) worldData->popSize.at(popID);
        
    }
}
void World::yAgefreqUpdate(){
    // This function will set to zero the sizes (and frequencies) of Contexts with Y in the SDR (sexsite=1)
    // It leads to autosomes with balanced polymorphism A1=1-preY_afreq, X=preY_afreq)
    
    for( cluster_t::iterator iter = cluster.begin(); iter != cluster.end(); ++iter ) {
        
        UINT popID=(*iter).first.pop;
        UINT siteA=(*iter).first.siteA;
        UINT Sexsite= (*iter).first.Sexsite;
        double size;
        
        if(Sexsite==1 ) {
            size=0;
        }
        else if(Sexsite==0 && siteA==1) {
            size=  worldData->preY_afreq *(double)worldData->popSize.at(popID);
        }
        
        else if(siteA==0 && Sexsite==0) {
            size= (1-worldData->preY_afreq) * (double)worldData->popSize.at(popID);
        }
        else size=0;
        
        worldData->nChromos[iter->second] = size ;
        worldData->freq.at(iter->second) = size / (double) worldData->popSize.at(popID);
        
    }
}
