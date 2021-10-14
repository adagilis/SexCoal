/*
 *  siteNode.cpp
 *  Created by Mark Kirkpatrick on 7/28/08.
 *	See "siteNode.h" for documentation.
 */

#include <iostream>
	using std::cout;
	using std::endl;
#include <algorithm>
	
#include "sitenode.h"
#include "argnode.h"
#include "ran_mk.h"

SiteNode::SiteNode()	// Null constructor
{
}

SiteNode::SiteNode(double sitePos, shared_ptr < ARGNode > argRoot)		// The typical constructor
//
// Constructs the gene tree for a site at chromosome position sitePos that descends from the ARG node *argRoot
// and that has ancestral site node *ancSiteNode.  To construct the gene tree from the root of its ARG,
// call:
//			SiteNode(sitePos, *argRoot)
//
//
//
{
	// For debugging:
	QDBG("\n   Making a SiteNode at node " << argRoot->getNodeNumber() )

	// Descend the gene tree until we reach the next split at our site:
	// 
	shared_ptr < ARGNode > argSiteRoot = getNextSiteNode(sitePos, argRoot);	// Find the ARGNode that corresponds to the root for the gene tree at this site.
	
	// For debugging:
	//	if(argSiteRoot != argRoot)
	//		cout << "    (Skipped to node " << argSiteRoot->getNodeNumber() << ")" << endl;
	// Set the basic private data:
	//
	sitePosition = sitePos;
	nodeNumber = argSiteRoot->getNodeNumber();
	time = argSiteRoot->getTime();
	context = argSiteRoot->getContext();
 tips = argSiteRoot->getTips();

	// Set the descendant nodes:
	//
	int argSiteRootNDesc = argSiteRoot->getNDescendants();
	int i;
	for(i = 0; i < argSiteRootNDesc; i++)
		if(siteQuery(sitePos, argSiteRoot->getDescendantSegmentVector(i)) )				// descendant[i] of *argSiteRoot is a carrier for sitePosition, 
		{
			shared_ptr<SiteNode> d;
			d.reset(new SiteNode(sitePos, argSiteRoot->getDescendantNode(i)));
			descendant.push_back(d);															//	so set that node as a descendant of this node and make a 
		}																						//	new SiteNode for it, which then recursively descends down that 
																								//	branch of the gene tree.
	totaltime = 0; // ccd. initial totaltime
	generateTotalTime(totaltime); //ccd, calculate totaltime
} // end SiteNode::SiteNode

SiteNode::~SiteNode()		// Destructor
{
	// For debugging:
	//	cout << "  SiteNode::~SiteNode called for node " << nodeNumber << endl;
}

vector<Base> SiteNode::mutate_jc(Base base0, double mu, double t)
//
//	Returns a vector with the base states along a descendant branch from a node
//	with base state base0, given a mutation rate mu and the number of generations
//	between the nodes is t.
//
//	Uses the Jukes-Cantor model:  changes between all pairs of bases are equally likely.
//
{
	vector<Base> baseVec(1, base0);		// Make a vector of bases and set the first element to base0 
	double timeLeft = t;
	
	do{
		double waitTime = randexp(mu);	// Waiting time until the next mutation, assumed independent of the current base state
		
		if(timeLeft < waitTime)			// No more mutation, return
			return baseVec;
			
		else {							// Mutation happens
			int newBase = ( baseVec.back() + randint(1, 3) ) % 4;	// All changes assumed equally likely
			baseVec.push_back(Base(newBase));						// Converts newBase, which is an int, to a Base via a type case
		} 
			
		timeLeft -= waitTime;
		
	} while(timeLeft > 0);	// Normally a value is returned before we get here
	
	return baseVec;
}



void SiteNode::snpDescendants(vector<Base> & myOutSeq)
//	
//  usage: snpDescendants( double pick ); pick = randreal()*totaltime;
//	Puts mutations on the descendant branches of this node (that is, sets the values of the vector baseState
//	  in the descendant nodes), given base0 as the base state at the ancestral node, mu as the mutation rate
//
{
	baseState.push_back( A );//set to A
	double myPick = randreal()*totaltime;
	
/** debug */
//cout << myPick << " myPick " << endl;	
	
	double y = 0;
	
	snpPut(myPick, y, myOutSeq);		// Descend gene tree recursively to that descendant and repeat
	
}

void SiteNode::snpPut(double myPick, double& y, vector<Base> & myOutSeq)
{		
		if(descendant.size() == 0){
			myOutSeq[getNodeNumber()] = getBaseState();

		};
	for(UINT iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
	{
		descendant[iDesc]->addBaseState( getBaseState() );
		if( y < myPick )
		{	
			y += time - descendant[iDesc]->getTime();
			if( y >= myPick)
			{
				descendant[iDesc]->addBaseState( C ); // add C to basevec
				/** debug */

			}
		};
		descendant[iDesc]->snpPut(myPick, y,myOutSeq);		// Descend gene tree recursively to that descendant and repeat
	}
}

set<int> SiteNode::getTips()	// Return the identifier number of this node
{
	return(tips);
}


void SiteNode::generateSubTreeLength(set<int> myTips, double &x)
//
// for given tree tips, return the total branch length of subtree with and only with those tips
// However, must manually make sure that the tips are descedents of the Node been queried
//
{
  set<int>::iterator it;

 set<int> myintersection;
 set_intersection( tips.begin(), tips.end(),myTips.begin(),myTips.end(),std::inserter(myintersection,myintersection.begin()));
	for(UINT iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
	{
  set<int> myintersection_des;
  set<int> dtips = descendant[iDesc]->getTips();
  set_intersection( dtips.begin(),dtips.end(),myTips.begin(),myTips.end(),std::inserter(myintersection_des,myintersection_des.begin()));
		if( myintersection_des.size() > 0   )
  {
		 if( myTips.size() > myintersection_des.size() ){x += time - descendant[iDesc]->getTime(); }
		 descendant[iDesc]->generateSubTreeLength(myTips,x);		// Descend gene tree recursively to that descendant and repeat
  }
	}
}

void SiteNode::mutateDescendants(Base base0, double mu)
//
//	Puts mutations on the descendant branches of this node (that is, sets the values of the vector baseState
//	  in the descendant nodes), given base0 as the base state at the ancestral node, mu as the mutation rate
//
{
	static int nCalls;	// Equals 0 i.f.f. this is the first time the rountine is called
	if(nCalls == 0) {	// ... then this node is the root of the gene tree, so set its base state to base0
		baseState.push_back(base0);
		nCalls++;
	}	
	
	for(UINT iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
	{
		vector<Base> baseVec = mutate_jc(base0, mu, time - descendant[iDesc]->getTime());	// Make vector of base states along the branch to descendant iDesc
		descendant[iDesc]->putBaseStates(baseVec);											// Set the vector baseState of descendant[iDesc] equal to that vector
		descendant[iDesc]->mutateDescendants(descendant[iDesc]->getBaseState(), mu);		// Descend gene tree recursively to that descendant and repeat
	}
}



void SiteNode::putBaseStates(vector<Base> baseVec)
//
// Sets baseState (the vector of ancestral base states) equal to baseVec
//
{
	baseState = baseVec;
}

void SiteNode::addBaseState(Base baseA)
//
// Sets baseState (the vector of ancestral base states) equal to baseVec
//
{
	baseState.push_back( baseA );
}

int SiteNode::getNodeNumber()
{
	return(nodeNumber);
}



double SiteNode::getTime()
{
	return(time);
}


double SiteNode::getTotalTime()
{
	return(totaltime);
}


Context SiteNode::getContext()
{
	return(context);
}



int SiteNode::getNDescendants()
{
	return(descendant.size());
}



shared_ptr< SiteNode> SiteNode::getDescendant(int i)
{
	return(descendant[i]);
}



Base SiteNode::getBaseState()
{
	return(baseState.back());
}



Base SiteNode::getBaseState(int i)
{
	return(baseState.at(i));
}


void SiteNode::outputNode()
{
	cout << " SiteNode " << nodeNumber << " for site " << sitePosition << endl;
	cout << "   time = " << time << endl;
	
	if(descendant.size() == 0)
		cout << "   Terminal node." << endl;
	else
	{
		cout << "   nDescendants = " << descendant.size() << endl;
		cout << "     Node(s):  ";
		for(UINT i = 0; i < descendant.size(); i++)
			cout << descendant[i]->getNodeNumber() << "  ";
		cout << endl;
	}
	cout << "   number of base states on ancestral branch = " << baseState.size() << endl << "   ";
	
	for(UINT i=0; i < baseState.size(); i++)
		cout << "  " << baseState[i];
	if(baseState.size() > 0)
		cout << endl;
	cout << endl;
}


void SiteNode::outputTree()
{
	cout << " SiteNode " << nodeNumber << " for site " << sitePosition << endl;
	cout << "   time = " << time << endl;
	
	if(descendant.size() == 0){
		cout << "   Terminal node" << endl;
		cout << " node " << nodeNumber << " base info." << baseState.size() << endl;//ccd
	}
	else
	{
		cout << "   nDescendants = " << descendant.size() << endl;
		cout << "     Node(s):  ";
		for(UINT i = 0; i < descendant.size(); i++)
			cout << descendant[i]->getNodeNumber() << "  ";
		cout << endl;
	}
	cout << "   number of base states on ancestral branch = " << baseState.size() << endl << "   ";
	
	for(UINT i=0; i < baseState.size(); i++)
		cout << "  " << baseState[i];
	if(baseState.size() > 0)
		cout << endl;
	cout << endl;
	for(UINT i=0; i < descendant.size(); i++)
		descendant[i]->outputTree();
}



void SiteNode::outputTreeV3( vector<Base> * myXseq, vector<Base> * myYseq )//temp output tree with tip base state by ccd
{ 

	if(descendant.size() == 0){
		if( getContext().Sexsite == 0){
			myXseq->push_back(getBaseState());
		}
		else{
			myYseq->push_back(getBaseState());
		}

	};
	for(UINT i=0; i < descendant.size(); i++)
		descendant[i]->outputTreeV3( myXseq,myYseq);
}


void SiteNode::generateTotalTime(double &x)
/*
 * Get the total length of the descendant branches of the tree rooted by this node
 * 
 */
{	
	for(UINT iDesc = 0; iDesc < descendant.size(); iDesc++)	// Loop through the branches descending from this node
	{
		x += time - descendant[iDesc]->getTime();
		descendant[iDesc]->generateTotalTime(x);		// Descend gene tree recursively to that descendant and repeat
	}
	
}
