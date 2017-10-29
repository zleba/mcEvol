#include "SplineGen.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <cassert>

SplineGen::SplineGen()
{
	N = 0;
	minX = minY = -1e40;
	maxX = maxY =  1e40;

}

void SplineGen::AddNode(Double x, Double y)
{
	nodes.push_back( make_pair(x,y) );
	++N;

}

void SplineGen::Finalize()
{
	sort( nodes.begin(), nodes.end(), [](pair<Double,Double> a, pair<Double,Double> b) {return a.first < b.first; } );

	isInc = nodes[N-1].second >= nodes[0].second;

	//Check if the Y is decreasing
	for(int i = 1; i < N; ++i)
		if( nodes[i].second < nodes[i-1].second ) {
			cout << "Wrong ordering of nodes " <<isInc << endl;
			cout << "i =  " << i <<" "<< N<< endl;
			for(int k = max(0,i-10); k<i+10 && k < N; ++k)
				cout << k <<" "<<setprecision(15)<<  nodes[k].first<<" "<< nodes[k].second << endl;
			//cout << "nodes[i], nodes[i-1] =  "<< nodes[i].second <<" "<< nodes[i-1].second  << endl;
			assert(0);
		}

	minX = nodes[0].first;
	maxX = nodes[N-1].first;

	minY = nodes[0].second;
	maxY = nodes[N-1].second;

}


Double SplineGen::EvalY(Double x) const
{
	const Double eps = 1e-10;

	auto itU = lower_bound(nodes.begin(), nodes.end(), make_pair(x,0) ,
	       [](pair<Double,Double> a, pair<Double,Double> b)  {return a.first < b.first;}  );

	if(itU == nodes.end()) {
		if(abs(nodes.back().first - x) < eps )
			--itU;
		else {
			cout << "Problem with inversion " << endl;
			exit(1);
		}
	}

	auto itL = itU-1;

	if(itL < nodes.begin())
		return minY;

	const double &x1 = itL->first;
	const double &x2 = itU->first;

	const double &y1 = itL->second;
	const double &y2 = itU->second;


	assert( x1 <= x+eps && x-eps <= x2 );

	return ( x*(y2 - y1) + y1*x2 - y2*x1 ) / (x2 - x1);
}




Double SplineGen::EvalX(Double y) const
{
	auto itL = lower_bound(nodes.begin(), nodes.end(), make_pair(0,y),
	       [&](pair<Double,Double> a, pair<Double,Double> b)  {return isInc ^ (a.second < b.second);}  );
	auto itU = itL+1;

	if(itU == nodes.end())
		return maxX;

	const double &x1 = itL->first;
	const double &x2 = itU->first;

	const double &y1 = itL->second;
	const double &y2 = itU->second;

	assert( y1 <= y && y <= y2 );

	return ( y*(x2 - x1) + x1*y2 - x2*y1 ) / (y2 - y1);
	
}
