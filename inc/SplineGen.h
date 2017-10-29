#ifndef _SplineGen_
#define _SplineGen_


typedef double Double;

#include <vector>
#include <utility>

using namespace std;

///////////////////////////
// Class to interpolate function with linear
// Spline has generally NON-equidistant nodes
///////////////////////////

class SplineGen {
public:
	
	SplineGen();

	Double EvalY(Double x) const; //evaluate function
	Double EvalX(Double y) const; //evaluate function

	void AddNode(Double x, Double y);
	void Finalize();

private:	
	int N;               //number of intervals (nodes-1)
	Double minX, maxX;   //value of lowest and highest node
	Double minY, maxY;   //lowest and highest value
	bool isInc;

	vector< pair<Double,Double> > nodes; 

};

#endif
