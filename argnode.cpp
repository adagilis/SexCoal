/*
 *  argnode.cpp
 *
 *  Created by Mark Kirkpatrick on 7/9/08.
 *
 */
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

using boost::shared_ptr;


#include <iostream>
	using std::cout;
	using std::endl;
	
#include "chromosome.h"
#include "argnode.h"
#include "typedefs.h"




ARGNode::ARGNode()	// Null constructor.  Used to make an (undefined) ancestor for an ARGNode
{
}


ARGNode::ARGNode(int nNode, shared_ptr<Chromosome> chr)
//
// Constructor for a terminal node.  The segments carried by this node are
//	stored in descendant.segmentVec[0]
{	
	nodeNumber = nNode;
	time = 0;
	context = chr->getContext();
 tips.insert(nNode);
}



ARGNode::ARGNode(int nNode, shared_ptr<Chromosome> chr, const double t)
//
// Constructor for a recombination node.  Also sets the current node as
//	the ancestor of the chromosome's descendant node.
{
	nodeNumber = nNode;
	time = t;
	context = chr->getContext();	
	tips = chr->getDescendant()->getTips();
 	
	Branch branch (chr->getDescendant(), chr->getSegmentVector() );
	descendant.push_back( branch );	
	chr->getDescendant()->addAncestor(this);
}



ARGNode::ARGNode(int nNode, shared_ptr<Chromosome> chr1, shared_ptr<Chromosome> chr2, const double t)
//
// Constructor for a coalescent node.  Sets the current node as
//	the ancestor of each chromosome's descendant node.
//
{
	nodeNumber = nNode;
	time = t;
	context = chr1->getContext();
	if(!(context== chr2->getContext())) cout<<"\nError ARGNode::ARGNode ===> coalescence between carriers of different contexts\n";
	Branch branch1 (chr1->getDescendant(), chr1->getSegmentVector());
	descendant.push_back( branch1 );	
	Branch branch2 (chr2->getDescendant(), chr2->getSegmentVector());
	descendant.push_back( branch2 );
	chr1->getDescendant()->addAncestor(this);
	chr2->getDescendant()->addAncestor(this);
 
 set<int> tips1 = chr1->getDescendant()->getTips();
 set<int> tips2 = chr2->getDescendant()->getTips();
 set_union( tips1.begin(), tips1.end(), tips2.begin(), tips2.end(),std::inserter(tips,tips.begin()));
}



ARGNode::~ARGNode()	// Destructor
{
}



void ARGNode::addAncestor(ARGNode * ancNode)	// Add ancNode as an ancestral node
{

	ancestor.push_back( ancNode );
	if(ancNode->getTime() < time)
		cout << endl << ">>>WARNING: ARGNode::addAncestor: ancestor at time " << ancNode->getTime() << 
						" is younger than its descendant node at time " << time << endl << endl;
}



void ARGNode::addDescendant(shared_ptr < ARGNode >descNode, vector<Segment> segVec)
{
//
// Adds descNode as a descendant ARG node, and segVec as its corresponding segment array.
//	Also sets the current node as the ancestor of that descendant.
//


	if(descNode->getTime() > time)
		cout << endl << ">>>WARNING: ARGNode::addDescendant: descendant at time " << descNode->getTime() << 
						" is older than its ancestor node at time " << time << endl << endl;
						
	Branch branch (descNode, segVec);
	descendant.push_back( branch );	
	descNode->addAncestor(this);				// Set the ancestor of descNode equal to this node
}
	

set<int> ARGNode::getTips()	// Return the identifier number of this node
{
	return(tips);
}



int ARGNode::getNodeNumber()	// Return the identifier number of this node
{
	return(nodeNumber);
}



double ARGNode::getTime()	// Return the time (= age) of this node
{
	return(time);
}



Context ARGNode::getContext()	// Return the context this node
{
	return(context);
}



int ARGNode::getNAncestors()	// Return the number of ancestral nodes from this node
{
	return(ancestor.size());
}



ARGNode* ARGNode::getAncestor(UINT i)	// Get ancestral node i
{

	if(i>ancestor.size()-1)		// Could handle this error more gracefully?
	{
		cout << endl << ">>>ERROR: ARGNode::getAncestor called for ancestor[" << i << "], but this node has only ";
		cout << ancestor.size() << " ancestors" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(ancestor[i]);
}



int ARGNode::getNDescendants()	// Return the number of descendant nodes from this node
{
	return(descendant.size());
}



shared_ptr<ARGNode> ARGNode::getDescendantNode(int i)	// Get descendant ARG node i
{

	if(i > descendant.size())	// Can we handle this error more gracefully?
	{
		cout << endl << ">>>ERROR: ARGNode::getDescendant called for descendant[" << i << "], but this node has only ";
		cout << descendant.size() << " descendants" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(descendant.at(i).node);
}



vector< Segment > ARGNode::getDescendantSegmentVector(int i)	
//
// Get the vector of chromosome segments corresponding to the branch leading to descendant i
//
{

	if(i > descendant.size())		// Is there a better way to handle this error?
	{
		cout << endl << ">>>ERROR: ARGNode::getDescendantSegmentList called for descendantSegmentList[" << i << "], but this node has only ";
		cout << descendant.size() << " descendants" << endl << endl << "Exiting..." << endl;
		exit(-1);	
	};
	return(descendant.at(i).segmentVec);
}



void ARGNode::outputNode()	// Prints out the node's data
{
	cout << " ARGNode " << nodeNumber << ":" << endl;
	cout << "  time = " << time << endl;
	cout << "  nAncestors = " << ancestor.size() << endl;
	
	if(ancestor.size() > 0)
	{
		cout << "   Ancestal node number(s): ";
		cout << ancestor[0]->getNodeNumber();
		
		for(UINT i=1; i<ancestor.size(); i++)
			cout << ", " << ancestor[i]->getNodeNumber();
		cout << endl;
	}
	cout << "  nDescendants = " << descendant.size() << endl;
	if(descendant.size() > 0)
		for(UINT i=0; i<descendant.size(); i++)
		{
			
			cout << "   descendant[" << i << "] = node " << descendant.at(i).node->getNodeNumber()  <<  endl;
			cout << "     Segment(s) descending to that node:";
			for(UINT j=0; j<descendant.at(i).segmentVec.size(); j++)
				cout << "  (" << descendant.at(i).segmentVec.at(j).L << ", " << descendant.at(i).segmentVec.at(j).R << ")";
			cout << endl;
			}
		
}



void ARGNode::outputARG()	// Prints out the ARG descending from this node
{
	cout << " ARGNode " << nodeNumber << ":" << endl;
	cout << "  time = " << time << endl;
	cout << "  nAncestors = " << ancestor.size() << endl;
	
	if(ancestor.size() > 0)
	{
		cout << "   Ancestal node number(s): ";
		cout << ancestor[0]->getNodeNumber();
		
		for(UINT i=1; i<ancestor.size(); i++)
			cout << ", " << ancestor[i]->getNodeNumber();
		cout << endl;
	}
	cout << "  nDescendants = " << descendant.size() << endl;
	if(descendant.size() > 0)
		for(UINT i=0; i<descendant.size(); i++)
		{
			
			cout << "   descendant[" << i << "] = node " << descendant.at(i).node->getNodeNumber()  <<  endl;
			cout << "     Segment(s) descending to that node:";
			for(UINT j=0; j<descendant.at(i).segmentVec.size(); j++)
				cout << "  (" << descendant.at(i).segmentVec.at(j).L << ", " << descendant.at(i).segmentVec.at(j).R << ")";
			cout << endl;
			
		}
	cout << endl;
	for(UINT i=0; i<descendant.size(); i++)
	{
		descendant.at(i).node->outputARG();
	}
}

Branch::Branch(shared_ptr<ARGNode> n, vector<Segment> segVec){
	node= n;
	segmentVec= segVec;
}

Branch::~Branch(){
	
}

Branch::Branch(const Branch& br){
	node=br.node;
	segmentVec= br.segmentVec;
	
}	

void Branch::operator =(const Branch& br){
	node=br.node;
	segmentVec= br.segmentVec;
	
}


