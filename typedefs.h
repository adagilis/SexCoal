/*
 *  typedefs.h
 *
 *  Created by Mark Kirkpatrick on 7/15/08.
 *	Modified 2-IX-08.
 *
 */
 
 
#ifndef TYPEDEFS_
#define TYPEDEFS_

#include <boost/enable_shared_from_this.hpp>

//#include "argnode.h"
#include <map>
using std::map;
// Type definitions:

typedef unsigned int UINT;


// Enumerators (used e.g. by the structure Context):

enum Sex {F, M};				// Used for Sex-of-carrier (= SOC) and Sex-of-origin (= SOO)
enum Karyotype {SS, SI, II};	// SS = standard homokaryotype, SI = heterokaryotype, II = inverted homokaryotype
enum Base {A, C, G, T};			// The four DNA bases


// Structures:

struct Segment
//
//	Represents a single segment of chromosome by its left and right map boundaries
//
{
	double L;	// Left boundary
	double R;	// Right boundary
};


struct Context 

{
public:
	//Sex soo;	// Sex-of-origin, i.e. the Sex of the parent from whom the chromosome was inherited
	//int soc;	// Sex-of-carrier, i.e. the Sex of the individual carrying the chromosome
	int	pop;	// Population number
	int siteA; // state of site A = [0,1]
	int Sexsite;
	
	Context(){}
	Context(const int D, const int E, const int C) {pop=D; Sexsite=E; siteA=C; }
	Context(const Context& other) {siteA=other.siteA; pop=other.pop; Sexsite=other.Sexsite;}
	
	bool operator<(const Context& other) const    {
		if(pop<other.pop) return true;
		if(pop==other.pop && Sexsite<other.Sexsite) return true;
		if(pop==other.pop && Sexsite==other.Sexsite && siteA < other.siteA) return true;
		//if(pop==other.pop && Sexsite==other.Sexsite && siteA == other.siteA && soc < other.soc) return true;
		return false; }

	bool operator==(const Context& other) const    {
		if(pop==other.pop && siteA==other.siteA && Sexsite == other.Sexsite) return true;
//		if(pop==other.pop && siteA==other.siteA && Sexsite == other.Sexsite && soc == other.soc) return true;
	return false; }
 
    void setNew(const int D, const int E, const int C) {pop=D; Sexsite=E; siteA=C;}

	
};

typedef map<Context, int> cluster_t;



#endif // TYPEDEFS_


// ==================================================================



//Debugging Code
#ifdef DEBUGGER
#define DBG(x) cout << x << endl;
#else
#define DBG(x)
#endif

//Current Debugging so everything isn't so messy Code
#ifdef CURRENTDEBUGGER
#define CDBG(x) cout << x << endl;
#else
#define CDBG(x)
#endif

//Recombination Debug
#ifdef RCDBGR
#define RCDBG(x) cout << x << endl;
#else
#define RCDBG(x)
#endif

//coal Debug
#ifdef COALDBGR
#define COALDBG(x) cout << x << endl;
#else
#define COALDBG(x)
#endif

// Site Query Dbg
#ifdef QDBGR
#define QDBG(x) cout << x << endl;
#else
#define QDBG(x)
#endif
    // sweeps  Dbg
#ifdef SWEEPDBGR
#define SDBG(x) cout << x << endl;
#else
#define SDBG(x)
#endif

