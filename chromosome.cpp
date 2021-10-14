/*
 * chromosome.cpp
 *
 * This class is used to represent chromosomes that carry segments that
 * are ancestral to those in the sample at time 0. The data in Chromosome
 * are
 *	• the context (pop, selected site, etc.)
 *	• a vector of chromosome segments (paried left and right)
 *	• the total length of these segments
 *	• the descendant node of the chromosome
 *
 * Functionality includes operations for recombination and coalescence.
 *
 * Written Shane Pope, VIII-08;  
 *
 * modified by Mark Kirkpatrick, IX-08
 *
 * Modified and annotated by Changde Cheng, X-14; 
 */

//#define DEBUGGER

// Includes from STL:
	#include <iostream>
		using std::cout;
		using std::endl;
	#include <vector>
		using std::vector;
	#include <algorithm>
		using std::max;
		using std::min;

// Includes from Boost:
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr;
	#include <boost/weak_ptr.hpp>
		using boost::weak_ptr;
	#include <boost/foreach.hpp>

// Includes from our files:
	#include "typedefs.h"
	#include "chromosome.h"
	#include "argnode.h"
	#include "ran_mk.h"


struct Chromosome::ChromosomeData{
/*
 * This structure contains all the data for a chromosome
 */	
	Context context;
// Contains the chromomsome's SOO, SOC, population, etc.
	vector<Segment> chromosomeSegments;
// Vector of chromosome segments, each of which is a pair of floats
	shared_ptr < ARGNode > descendant;
// Pointer to the descendant ARG node of this chromosome
};

Chromosome::Chromosome(Context ctxt, vector<Segment> chromSegments,shared_ptr < ARGNode > desc ){
/*
 * Constructs a chromosome
 * Written by Shane Pope VIII-08, modified by M.K. 2-IX-08
 */
	DBG("Chromosome:: Constructing Chromosome()...")	
// Makes chromData a pointer to the private data
	chromData.reset(new ChromosomeData);
// Copy the context data:
	chromData->context = ctxt;
// Copy the segments:
	chromData->chromosomeSegments.reserve( chromSegments.size() );
	copy(chromSegments.begin(), chromSegments.end(), back_inserter(chromData->chromosomeSegments) );
// Set the descendant ARG node:
	chromData->descendant = desc;
}



Chromosome::~Chromosome(){
/*
 * Destructor
 */
	DBG("Chromosome:: Destructing...")

}

shared_ptr<Chromosome> Chromosome::copy_values(){
/* 
 * don't understand it yet(ccd)
 */	
	vector<Segment> segs=getSegmentVector();
	assert(segs.size()!=0);
	shared_ptr<Chromosome> newChromosome(new Chromosome(getContext(), segs, getDescendant()));
	return newChromosome;
}


Chromosome::Chromosome( Chromosome& chr ){
/*
 * Copy constructor, needed to put in std::vector
 * Written by Shane Pope VIII-08, modified by M.K. 2-IX-08
 */
	DBG("Chromosome:: Copy constructor called")	
	chromData.reset(new ChromosomeData);	
// Copy the context data:
	chromData->context = chr.chromData->context;
// Copy the segments:
	chromData->chromosomeSegments.reserve(
		chr.chromData->chromosomeSegments.size() 
	);
	copy(
		chr.chromData->chromosomeSegments.begin(),
		chr.chromData->chromosomeSegments.end(),
		back_inserter(chromData->chromosomeSegments)
	);
// Set the descendant ARG node:
	chromData->descendant = chr.getDescendant();
}
	
	
void Chromosome::operator = ( Chromosome& chr ){
/*
 * Equals operator;  same result as copy constructor
 * Written by Shane Pope VIII-08, modified by M.K. 2-IX-08
 */
	DBG("Chromosome:: = Operator called")	
	chromData.reset(new ChromosomeData);
//Copy the context:
	chromData->context = chr.chromData->context;
//Copy the segments:
	chromData->chromosomeSegments.reserve( chr.chromData->chromosomeSegments.size() );
	copy(chr.chromData->chromosomeSegments.begin(),
		 chr.chromData->chromosomeSegments.end(),
		 back_inserter(chromData->chromosomeSegments));	
// Set the descendant ARG node:
	chromData->descendant = chr.getDescendant();
}


/*
void Chromosome::setSOO(Sex sx) {
	chromData->context.soo = sx;
}


void Chromosome::setSOC(Sex sx) {
	chromData->context.soc = sx;
}
*/


void Chromosome::setPopulation(UINT pop) {
/*
 * Set up population context
 */
	chromData->context.pop = pop;
}


void Chromosome::setsiteA(int a) {
/*
 * Set up siteA context, should be the status of either allele A1 or A2 for
 * sex antagonistic selectoin site (siteA). 
 */
	chromData->context.siteA= a;
}


void Chromosome::setContext(Context c) {
/*
 * set context as a whole, there are two other funtions directly set pop and 
 * siteA: setsiteA() and setPopulation()
 */
	chromData->context= c;
}


void Chromosome::setDescendant(shared_ptr < ARGNode > desc) {
/*
 * Set a pointer to the descendant ARG node
 */
	chromData->descendant = desc;
}


void Chromosome::setSegmentVector(vector<Segment> mySegVec){
/*
 * set vector of chromosome segments: ccd 10/10/2014
 */
	chromData->chromosomeSegments = mySegVec;
}


/*
Sex Chromosome::getSOC() {
	return chromData->context.soc;	
}*/


UINT Chromosome::getPopulation() {
/*
 * get populatin context of that chromosome
 */
	return chromData->context.pop;	
}


int Chromosome::getsiteA() {
/*
 * get siteA context of that chromosome
 */
	return chromData->context.siteA;
}


Context Chromosome::getContext() {
/**
 * return the context of a chromosome (population, X or Y, A1 or A2)
**/
	return chromData->context;
}


vector<Segment> Chromosome::getSegmentVector(){
/**
 * Returns the pointer to a vector of chromosome segments
**/
	return chromData->chromosomeSegments;
}


double Chromosome::getTotalLength() {
/**
 * Get the total length of the chromosome
**/
	
	int nSegs = chromData->chromosomeSegments.size();
// returen zero if there is no segments at all.
	if(nSegs == 0) return 0;
/**
 * I think the following line assumes that the coordinates of the segments is 
 * is ordered as such that small ones on the left. Not sure if this is the case.
 * need to check back in other cpp files to find out.
**/
	return chromData->chromosomeSegments[nSegs-1].R - chromData->chromosomeSegments[0].L;
}



double Chromosome::getTotalSegLength() {
/**
 * Get the total length of the segments on a chromosome
**/
	
	int nSegs = chromData->chromosomeSegments.size();
// returen zero if there is no segments at all.
	if(nSegs == 0) return 0;
	
	double myOutput=0;
	for(int i = 0; i < nSegs; i++){
		myOutput += chromData->chromosomeSegments[i].R - chromData->chromosomeSegments[i].L;
	}
	
	return myOutput;//I haven't consider the case when seg.L == seg.R
}


double Chromosome::getTotalLength(double Sexpos, double siteApos) {
/**
 * Get the length for recombination considering the position of the sex determining region and selected site A
**/
	
	int nSegs = chromData->chromosomeSegments.size();
// always check if there is any segments in the vector.
	if(nSegs == 0) return 0;
	
// again assumes that segment vectors is sorted or ordered from small to large
	double border_l =min( min(Sexpos, siteApos), chromData->chromosomeSegments[0].L);
	double border_r =max( max(Sexpos, siteApos), chromData->chromosomeSegments[nSegs-1].R);

	return border_r - border_l;
}


shared_ptr < ARGNode > Chromosome::getDescendant() {
/**
 * Returns a pointer to the descendant ARG node
**/
	return chromData->descendant;
}


/***
  * I made quite some changes here. Changde
***/


double Chromosome::calcOutBreakpoint(double sexPos,double siteApos, double outLength){
/**
 * gives a breakpoint for recombination outside of the sex - siteA segment
 * Changde rewrote to allow sexPos < siteApos
**/
    
	double breakpoint;
	double my_left;
	double my_right;

	if(sexPos < siteApos){ my_left = sexPos; my_right = siteApos;}
	else{my_left = siteApos; my_right = sexPos;};
/**
 * put a breakpoint between 0.L and 0.L+outLength
 * given 0.L on the left of my_left
**/
  
	if(chromData->chromosomeSegments[0].L< my_left) {
		breakpoint = randreal(0, outLength) + chromData->chromosomeSegments[0].L;	

					
//changde//^^^^^^^^^//
//		if(breakpoint > siteApos) breakpoint = siteApos + (breakpoint-sexPos) ;//I think this is wrong. ccd//
//changde//^^^^^^^^^//
		if(breakpoint > my_left) breakpoint = my_right + (breakpoint - my_left);//this is the right version
//changde//^^^^^^^^^//
    }
	else breakpoint = randreal(0, outLength) + my_right;
    
    return breakpoint;
};


/**
 * I guess this does recombine operation
**/
shared_ptr<Chromosome> Chromosome::recombine(
	double breakpoint, 
	shared_ptr < ARGNode > descNode, 
	pair<Context, Context> ancestors,
	int& inSegRec) {

	setContext(ancestors.first);
	Context newContext=ancestors.second;
	vector<Segment> newSegs = splitSegments(breakpoint,inSegRec);

/**	
	if(newSegs.size() > 0){
		cout << newSegs.at(0).L << " " << newSegs.at(0).R << "\n"; //debug ccd
	};
*/		
// Make a new chromosome that carries the segments to the right of the breakpoint:
	shared_ptr<Chromosome> newChromosome;
// if both carriers have material, update ARGNode
	if(newSegs.size()!=0 && chromData->chromosomeSegments.size()!=0){
//cout<<'\n';
	newChromosome.reset(new Chromosome(newContext, newSegs, descNode));
// Update the descendant ARG node for this chromosome:
	chromData->descendant = descNode;
	}
// else, the old ARGNode is used and the recomb event will be "ignored". Empty chromosomes will be deleted in recombineEvent.
	else{
		newChromosome.reset(new Chromosome(newContext, newSegs, chromData->descendant));
	}
// debug printing:
// if(newSegs.size()!=0) cout<<newContext.Sexsite<<newContext.siteA<<"\n";
// if(chromData->chromosomeSegments.size()!=0) cout<<getContext().Sexsite<<getContext().siteA<<"\n";
	
	return newChromosome;
}


vector<Segment> Chromosome::splitSegments(double breakpoint, int& inSegRec){
/**
 * Erases segments on this chromosome to the right of breakpoint, and returns the vector
 *	of segments that were erased.  If the breakpoint falls within a segment, the part of
 *	that segment to the right of the breakpoint is erase, and that part of the segment
 *	is included in the vector of segments that are returned.
 *
 * M.K., IX-08
**/
	RCDBG("splitting at "<<breakpoint )
// New vector of segments that lie to right of breakpoint
	vector<Segment> newSegVec;
	int nSegs = chromData->chromosomeSegments.size();
// Loop through segments carried by this chromosome
	for(int i = 0; i < nSegs; ++i){
// Breakpoint falls to the left or within this segment		
		if(breakpoint < chromData->chromosomeSegments[i].R){
			RCDBG("		"<<breakpoint <<" < "<<chromData->chromosomeSegments[i].R)
			int within=0;
// Breakpoint falls within this segment
				if(chromData->chromosomeSegments[i].L < breakpoint)
					{
					within=1;
					inSegRec = within;//inside a segment
					RCDBG("			breakpoint within segment, "<<chromData->chromosomeSegments[i].L <<" < "<<breakpoint)	
// Add partial segment to the right of breakpoint to the new segment vector, and 
// erase it from this chromosome:

					Segment newSeg = {breakpoint, chromData->chromosomeSegments[i].R};
					RCDBG("			New seg "<<newSeg.L<<" -- "<<newSeg.R)
					newSegVec.push_back(newSeg);
					chromData->chromosomeSegments[i].R = breakpoint;
					}
				else {
					RCDBG("			adding seg "<< i)
					newSegVec.push_back(chromData->chromosomeSegments[i]);
					}
			
			for(int j = i+1; j < nSegs; j++)
// Add remaining segments to the new segment vector	
			{
				RCDBG("			adding seg "<< j)
				newSegVec.push_back(chromData->chromosomeSegments[j]);	
			}
// Erase remaining segments from this chromosome
			chromData->chromosomeSegments.resize(i + within);
			return newSegVec;
// Return										
			
		} // end breakpoint falls to the left or within this segment
	} // end loop through segments
	
	newSegVec.resize(0);		// Breakpoint falls to the right of the rightmost segment carried by this chromosome
	RCDBG(">>>>>WARNING<<<< Breakpoint falls to the right of the end of the carrier. ")
	return newSegVec;
}

void Chromosome::merge(shared_ptr <Chromosome> chrom, shared_ptr < ARGNode > descNode) {
/**
 * Coalesces this chromosome with chrom.  
 * Sets the descendant to the new descNode
**/
	
	combineDoubletVectors(chromData->chromosomeSegments, chrom->chromData->chromosomeSegments);
	chromData->descendant = descNode;
}


void Chromosome::combineDoubletVectors(vector<Segment> vector1, vector<Segment> vector2){
/**
 * Add a doublet to a doublet vector and fix overlapping etc.
**/
	vector<Segment> outputVector;
	
	if( vector1.size() == 0 && vector2.size() == 0 );	// Both vectors are empty (never happen hopefully)
	else if(vector1.size() == 0)						// vector1 is empty, copy vector2
		copy(vector2.begin(), vector2.end(), back_inserter(outputVector));
	else if(vector2.size() == 0)						// vector2 is empty, copy vector1
		copy(vector1.begin(), vector1.end(), back_inserter(outputVector));
	else{												// Neither vector is empty
		
		//This is a "state machine" with 7 states.
		//The first four states correspond to different kinds of comparisons between vector1 and vector2.
		// State 0: vector1's L compared to vector2's L
		// State 1: vector1's R compared to vector2's R
		// State 2: vector1's R compared to vector2's L
		// State 3: vector1's L compared to vector2's R 
		// State 4: vector1 has no more segments, copy over vector2
		// State 5: vector2 has no more segments, copy over vector1
		// State 6: both vectors empty, we are finished copying
		
		//Vector indices:
		UINT indexVec1 = 0;
		UINT indexVec2 = 0;
		UINT indexVecOut = 0;
		
		UINT STATE = 0;
		bool finished = false;
        
		while(!finished){
			switch (STATE){
				case 0:{ //vector1 -> L, vector2 -> L
					if(indexVec1 == vector1.size()){
						STATE = 4;
					}
					else if(indexVec2 == vector2.size()){
						STATE = 5;
					}
					else{
						Segment tempDoublet;
						outputVector.push_back(tempDoublet);
						indexVecOut = outputVector.size()-1;
							
						if(vector1.at(indexVec1).L == vector2.at(indexVec2).L){
							outputVector.at(indexVecOut).L = vector1.at(indexVec1).L;
							STATE = 1; //1 = R, 2 = L
		            	}
			            else if(vector1.at(indexVec1).L < vector2.at(indexVec2).L){
							outputVector.at(indexVecOut).L = vector1.at(indexVec1).L;
							STATE = 2; //1 = R, 2 = L
			            }
			            else{
			                outputVector.at(indexVecOut).L = vector2.at(indexVec2).L;
			                STATE = 3; //1 = L, 2 = R
			            }
					}
		            break;//case 0
	        	}
	            
		        case 1:{	//1 = R, 2 = R
		            if(vector1.at(indexVec1).R == vector2.at(indexVec2).L){
		                outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			            indexVec1++;
			            indexVec2++;
			            if(indexVec1 == vector1.size()){
			            	STATE = 4;	//vec 1 is empty
			            }
			            else if(indexVec2 == vector2.size()){
			            	STATE = 5;	//vec 2 is empty
			            }
						else{
		                	STATE = 0; //1 = L, 2 = L
						}
		            }
		            else if(vector1.at(indexVec1).R < vector2.at(indexVec2).R){
			            indexVec1++;
		                STATE = 3; //1=L, 2=R
		            }
		            else{//vector1.at(indexVec1).R > vector2.at(indexVec2).R
			            indexVec2++;
		                STATE = 2; //1=R, 2=L
		            }
		            break;//case 1
		        }
		            
		        case 2:{	//1 = R 2 = L
					if(indexVec2 == vector2.size()){
			            outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			            indexVec1++;
						STATE = 5;
					}
					else{
			            if(vector1.at(indexVec1).R == vector2.at(indexVec2).L){
			                indexVec1++;
			            	STATE = 3;
			            }
			            else if(vector1.at(indexVec1).R < vector2.at(indexVec2).L){
			                outputVector.at(indexVecOut).R = vector1.at(indexVec1).R;
			                indexVec1++;
			                STATE = 0;
			            }
			            else{
			                STATE = 1;
			        	}
					}
		            break;//case 2
		        }
		            
				case 3:{ //1 = L,  2 = R
					if(indexVec1 == vector1.size()){
			            outputVector.at(indexVecOut).R = vector2.at(indexVec2).R;
			            indexVec2++;
						STATE = 4;
					}
					else{
			            if(vector1.at(indexVec1).L == vector2.at(indexVec2).R){
			                indexVec2++;
			            	STATE = 2;
			            }
			            else if(vector1.at(indexVec1).L < vector2.at(indexVec2).R){
			                STATE = 1;
			            }
			            else{
			                outputVector.at(indexVecOut).R = vector2.at(indexVec2).R;
			                indexVec2++;
			                STATE = 0;
			        	}
					}
		            break;//case 3
		        }
		        case 4:{ //vector1 is empty, copy vector2
		            while(indexVec2 < vector2.size()){
						Segment tempDoublet;
						tempDoublet.L = vector2.at(indexVec2).L;
						tempDoublet.R = vector2.at(indexVec2).R;
						outputVector.push_back(tempDoublet);
		                indexVec2++;
		            }
		            STATE = 6; // finished
		            break;//case 4
		        }
		        case 5:{ //vector2 is empty, copy vector1
		            while(indexVec1 < vector1.size()){
						Segment tempDoublet;
						tempDoublet.L = vector1.at(indexVec1).L;
						tempDoublet.R = vector1.at(indexVec1).R;
						outputVector.push_back(tempDoublet);
		                indexVec1++;
		            }
		            STATE = 6; // finished
		            break;//case 5
		        }
		        case 6 :{ //both vectors are empty.
		            finished = true;
		            break;//case 6
				}
			} // end switch
		} // end while
	}// end else
		
	chromData->chromosomeSegments = outputVector;
}


bool Chromosome::isEmpty(){
/**
 * Return logical values for the segment size of a chromosome
**/
	int nSegs = chromData->chromosomeSegments.size();
	if(nSegs == 0) return true;
	return false;
}






