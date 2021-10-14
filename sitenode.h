/*
 *  sitenode.h
 *
 *  Created by Mark Kirkpatrick on 7/28/08, modified 4-IX-08 so arrays are allocated dynamically using vector<>.
 *
 *	The SiteNode class contains the data needed to represent a node in a gene tree pertaining
 *	to a site on a chromosome whose evolution is represented by an ancestral recombination graph (ARG).
 *	Definitions for the members and methods are given below.
 *
 *  The call:
 *				geneTree = SiteNode(sitePos, &argNd, NULL)
 *
 *	returns the SiteNode for the gene tree for the site sitePos using as root the ARGNode argNd.  The 
 *	call:
 *
 *				geneTree = SiteNode(sitePos, &argNd, &siteNd)
 *
 *	returns the SiteNode that is the root for the subtree descending from argNd and setting its ancestral
 *	site node to the SiteNode siteNd.
 *
 *  baseState is an array of mutations that lie between this node and the SiteNodes ancestral to it.
 *	it.  mutationList[0] is the current base at this site, even if there were no mutations along
 *	the ancestral branch.  mutationList[1], mutationList[2], ... list mutations (bases) that occurred 
 *	along that branch in reverse chronological order.
 *
 */

#ifndef SITENODE_
#define SITENODE_

#include "typedefs.h"
#include "argnode.h"


class SiteNode
{
private:
	double sitePosition;
	// The position of the site
	int nodeNumber;
	// The identifier number of this node;  agrees with the nodeNumber for the corresponding ARGNode
	double time;
	// Time (age) of the node (I changed it from int to double, ccd)
	double totaltime;
	//total time down from that node, ccd.
	Context context;
	// Context of the chromosome at this node
	vector<shared_ptr < SiteNode > > descendant;
 // Vector of pointers to descendant nodes;  descendant.size() returns the number of descendants.
 set<int> tips; // the tips from this node
	vector<Base> baseState;
	// Array with the states on the branch ancestral to this SiteNode.  baseState[0] = baseState.front() is
	// the base state of this node's ancestral node, baseState.back() is the base state at this node,
	// baseState.size() is the number of base states in this vector, baseState.size() - 1 is the 
	// number of mutations that happened between this node and its ancestor.
	// calculate total time
 //	void generateTotalTime(double &ttime);
	//generate total time from this node
	//
	// Functions that work with mutation:
	vector<Base> mutate_jc
		(Base base0, double mu,	double time);
	// base0 is the base state at the ancestral node, mu the mutation rate,
													//	and time the number of generations between this node and its ancestral node
	void putBaseStates(vector<Base> baseVec);		// Sets baseState (the vector of ancestral base states) equal to baseVec
	void addBaseState(Base a);						// add an element to baseVec
	void snpPut(double myPick, double& y, vector<Base> & myOutSeq);
	
	
public:
	// Constructors and destructor:
	SiteNode();											// Null constructor
	SiteNode(double sitePos, shared_ptr < ARGNode > argRoot);			// Typical constructor
	~SiteNode();										// Destructor
 void generateTotalTime(double &ttime);
	//generate total time from this node
	//
	// Function that puts mutations along branches:
	void mutateDescendants(Base base0, double mu);		// Put mutations on all branches of the gene tree descending from this node,
														//	given this node has base base0, a mutation rate mu
	void snpDescendants(vector<Base>& myOutSeq);								// generate a snp on genetree
	// Functions that return data about this node:
	int getNodeNumber();								// Returns the node number
	double getTime();									// Returns the time (age) of the node
	double getTotalTime();								// Returns the total time along the gene tree rooted by the node
	Context getContext();								// Returns the context of the node
	int getNDescendants();								// Return the number of descendants
 set<int> getTips(); // Get the tree tips
 void generateSubTreeLength(set<int> myTips, double &x); //Generate the total length of subtree defined by tree tips
	shared_ptr < SiteNode> getDescendant(int i);					// Return a pointer to descendant i
	Base getBaseState();								// Get the base state at this node;  equivalent to getBaseState.back()
	Base getBaseState(int i);							// Get base state i along ancestral branch

	//
	// Output data for this node and/or the whole gene tree:
	void outputNode();									// Output data for this node only
	void outputTree();									// Output the gene tree descending from this node
//	void outputTreeV2();									// a different version of outputTree ccd 10/4/2014
	void outputTreeV3( vector<Base> * myXseq, vector<Base> * myYseq );									// a different version of outputTree ccd 10/6/2014	
	

};

static bool siteQuery(double sitePos, vector<Segment> segList)
//
//	Returns TRUE if the chromosome at the ARG node *argNode carries site sitePos, otherwise returns FALSE
//
{
	for(int i=0; i<segList.size(); ++i)
	{	
		// For debugging:
		QDBG( "      siteNode::siteQuery:  "<<i<<" " << segList[i].L << " < " << sitePos << " < " << segList[i].R << " ? " )
		
		if( segList.at(i).L < sitePos && sitePos < segList.at(i).R )
		{ QDBG("					YES (belongs to seg "<<i<<")")
			return(true);
		}else if( segList.at(i).L == sitePos )//the old line causing the recomb breakpoint shared by two segments
		{ QDBG("					YES (site = seg border)")
			return(true);
		}
		
	}
	QDBG("					NO")
	return(false);
} // end siteQuery

static shared_ptr < ARGNode > getNextSiteNode(double sitePos, shared_ptr < ARGNode > argRoot)
//
//	Descends the gene tree for the site sitePos starting at the shared_ptr < ARGNode > argRoot until the next split
//	or a terminal node is reached, then returns that ARGNode.
//
{
	// Count the number and make a list of the descendant ARGNodes from *argRoot that carry our site:
	//
	int uniqueDesc=0, nSiteDesc = 0, nARGDesc = argRoot->getNDescendants();
	for(int i = 0; i < nARGDesc; i++)
		if(siteQuery(sitePos, argRoot->getDescendantSegmentVector(i)) )
		{
			nSiteDesc++;
			uniqueDesc = i;	// If there is only one descendant of this ARGNode that carries this site, uniqueDesc
			//  will be the index number of that descendant in the descendant array of argRoot.
		}
	QDBG("Query says "<<nSiteDesc<<" and unique is "<<uniqueDesc<<'\n')	
//	cout << "~~~~~~~~~~~~" << endl;
//	argRoot->outputNode();
//	cout << argRoot->getTime() << "Query says "<<nSiteDesc<<" and unique is "<<uniqueDesc<<'\n';
	// Process the ARG node depending on the number of descendants:
	//
	if(nSiteDesc > 1)	{// *argRoot is a split
		QDBG("nSite is >1, returning argRoot\n")	
		return(argRoot);
	}
	
	else if(nSiteDesc == 1)	// *argRoot is not at the root for the gene tree for this site, so we proceed 
	{					//	to the descendant ARGNode that carries this site
		QDBG("nSite is ==1, recursive step...\n")
//		cout << "nSite is ==1, recursive step...\n";
		return(getNextSiteNode(sitePos, argRoot->getDescendantNode(uniqueDesc)));
	}
	// If we reach here, nSiteDesc == 0, and this is a terminal node:
	//
	QDBG("nSite is ==0, is this a terminal node?\n")
	return(argRoot);
	
} // end getNextSiteNode

#endif 
