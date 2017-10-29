#ifndef _Spline_
#define _Spline_


typedef double Double;

#include <vector>
#include <utility>

using namespace std;

///////////////////////////
// Class to interpolate function with linear or cubic spline
// Spline has EQUIVIDTANT nodes
///////////////////////////

template<typename T>
class Spline {
public:
	
	Spline(int _N, Double _minX, Double _maxX);
	Spline() {N = 0;}

	Double getXmin() const {return minX;}
	Double getXmax() const {return maxX;}
	int getN() const {return N;}


	virtual T func(Double) {return 0.*vals[0];}   //definition of org SLOW function
	virtual T funcD(Double) {return 0.*vals[0];}
	virtual pair<T,T> funcBoth(Double) {return make_pair(0.*vals[0],0.*vals[0]);}

	T Eval(Double x, bool lin = true); //evaluate function
	void Fill();
	void FillTogether();

	Spline GetInverseSpline(int Nnew, bool isLinear=true);
	Spline GetIntegratedSpline(bool lin = true) const;

	Spline InverseCubic(int Nnew);
	void CheckInverse(Spline<Double> &sp, bool isLin=true);

	bool isIncreasing() const;

	virtual void printNodes() const;

	void ShiftX(Double shift);
	void ShiftY(Double shift);

	void setNode(int i, T val, T valD) { vals[i] = val; valsD[i] = valD; }
protected:	
	int N;               //number of intervals (nodes-1)
	Double minX, maxX;   //value of lowest and highest node
	Double step;         //step between nodes

	vector<T> vals; //function values in nodes
	vector<T> valsD;//Dfunction values in nodes

};

#endif
