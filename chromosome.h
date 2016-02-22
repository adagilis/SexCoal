/*
 * Chromosome.h
 *
 *
 */


#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

// Includes for STL:
	#include <vector>
		using std::vector;
// Includes for Boost:
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr;
	#include <boost/scoped_ptr.hpp>
		using boost::scoped_ptr;
// Includes for our files:
	#include "typedefs.h"
	#include "argnode.h"

#include <math.h>
using std::pair;
using std::make_pair;

// Forward declaration:
	class ARGNode;


class Chromosome 
{
	
  private:
	struct ChromosomeData;								// Contains all of the chromosome's private data
	scoped_ptr<ChromosomeData> chromData;				// Pointer to the private data
	void combineDoubletVectors(vector<Segment> vector1, // Combines two segment vectors into one
							vector<Segment> vector2);
	vector<Segment> splitSegments(double breakpoint, int& inSegRec);	// Returns segments to the right of breakpoint, and
														//	deletes those segments from the input segment vector
  public:
	// Constructors and destructor:
	Chromosome();										// Empty constructor				
	Chromosome(const Context ctxt,						// The typical constructor
				vector<Segment> chromSegments,			
				shared_ptr < ARGNode > desc);			
	~Chromosome();										// Destructor

	// Copy functions:
	Chromosome( Chromosome& chr );			// Copy constructor (copies the chromosome)
	void operator = ( Chromosome& chr );	// Equals operator;  has same effect as the copy constructor	
	shared_ptr<Chromosome> copy_values();
		// Set and get data for this chromosome:
	void setSOO(Sex sx);					// Sets sex-of-origin
	void setSOC(Sex sx);					// Sets sex-of-carrier
	void setPopulation(UINT pop);			// Sets population number
	void setsiteA(int a);
	void setContext(Context c); 
	Context getContext();					// Returns whole context
	Sex getSOO();							// Returns sex-of-origin
	Sex getSOC();							// Returns sex-of-carrier
	UINT getPopulation();					// Returns population number
	int getsiteA();				// returns allele state at site A
	vector<Segment> getSegmentVector();		// Returns vector of chromosome segments
	void setSegmentVector(vector<Segment> mySegVec);		// set vector of chromosome segments
	double getTotalSegLength(); //changde: get total length of segments on a chr.	
	double getTotalLength();				// Returns total length of the chromosome, distance from L to R in segments	
	double getTotalLength(double Sexpos, double siteApos);	// Returns the maximum between the getTotalLength and the distance to the sex site. 
	double calcOutBreakpoint(double sexPos,double siteApos, double outLength);
    void setDescendant(shared_ptr < ARGNode > desc);		// Set the pointer to the descendant ARG node
	shared_ptr < ARGNode > getDescendant();				// Returns pointer to the descendant ARG node
	bool isEmpty();
	// Recombination and coalescent functions:
	void merge(shared_ptr <Chromosome> chrom, shared_ptr < ARGNode > descNode);				// Adds segments of chrom to this chromosome
	shared_ptr<Chromosome> recombine(double breakpoint, 
									shared_ptr < ARGNode > descNode,
								pair<Context, Context> ancestors,
								int& inSegRec
											);	// Deletes segments on this chromosome to the right of breakpoint, 
															//	returns a pointer to a new chromosome carrying those segments.
															//	descNode is the ARG node representing this recombination event.
};	

#endif /*CHROMOSOME_H_*/
