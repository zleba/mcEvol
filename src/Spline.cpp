#include "Spline.h"
#include "SplineGen.h"

#include <iostream>
#include <cstdlib>
#include <tuple>
#include <cassert>
#include <iomanip>

#include <valarray>
#include <typeinfo>


template<typename T>
Spline<T>::Spline(int _N, Double _minX, Double _maxX)
{
	N = _N;
	minX = _minX;
	maxX = _maxX;
	step = (maxX-minX) / N;

	vals.resize(N+1);
	valsD.resize(N+1);

}

template<typename T>
void Spline<T>::Fill()
{
	for(int i = 0; i < N+1; ++i) {
		Double Xval =   minX + i * step;
		vals[i]  = func( Xval );
		valsD[i] = funcD(Xval );
	}
	/*
	//app
	for(int i = 1; i < N; ++i) {
		valsD[i] = (vals[i+1]-vals[i-1])/(2.*step);
	}
	valsD[0] = (vals[1]-vals[0])/step;
	valsD[N] = (vals[N]-vals[N-1])/step;
	*/
}

template<typename T>
void Spline<T>::FillTogether()
{
	for(int i = 0; i < N+1; ++i) {
		Double Xval =   minX + i * step;
		tie(vals[i],valsD[i])  = funcBoth( Xval );
	}
}

//template<typename T>
//void print(


template<typename T>
T Spline<T>::Eval(Double x, bool lin)
{
	if(x <= minX)
		return vals[0];
	if(x >= maxX) {
		return vals[N];
	}

	Double pos = (x - minX) / step;
	int i = (int) pos;
	pos -= i;

	i = max(0, i);
	if(i >= N)
		return vals[N];
	//i = min(N-1, i);


	if(lin) {
		return vals[i]*(1-pos)  +  vals[i+1]*pos;
	}
	else {
		Double pos1 = (1-pos)*(1-pos);
		Double h00 = (1+2*pos)*pos1;
		Double h10 = pos*pos1;
		Double h01 = pos*pos*(3.-2*pos);
		Double h11 = pos*pos*(pos-1);

		T      res = (h00*vals[i]  + h01*vals[i+1] ) +
		             (h10*valsD[i] + h11*valsD[i+1])*step;

		return res;

	}
}

template<typename T>
void Spline<T>::printNodes() const
{
	cout << "Spline PrintOut (N = " << N <<") - start" <<endl;
	for(int i = 0; i < N+1; ++i) 
		cout << i <<" "<<setprecision(12)<< minX+i*step <<" "<< vals[i] << endl;
	cout << "Spline PrintOut (N = " << N <<") - end" <<endl;
}

template<>
void Spline<valarray<Double>>::printNodes() const
{
	cout << "Spline PrintOut (N = " << N <<")" <<endl;
	for(int i = 0; i < N+1; ++i) {
		cout << i <<"  "<< minX+i*step <<" ";
		for(unsigned j = 0; j < vals[i].size(); ++j)
			cout << vals[i][j] <<" ";
		cout << endl;
		cout << i <<"D "<< minX+i*step <<" ";
		for(unsigned j = 0; j < vals[i].size(); ++j)
			cout << valsD[i][j] <<" ";
		cout << endl;
	}
}





template<typename T>
bool Spline<T>::isIncreasing() const
{
	for(int i = 1; i < N+1; ++i)
		if( vals[i] < vals[i-1]) {
			return false;
			//cout << "There is problem for i " <<i<< endl;
			//exit(1);
		}

	return true;
}


//Dummy
template<>
bool Spline<valarray<Double>>::isIncreasing() const
{
	return true;
}




template<typename T>
Spline<T> Spline<T>::GetInverseSpline(int Nnew, bool isLinear)
{
	if(!isIncreasing()) {
		cout << "Spline does not homogenously grow" << endl;
		exit(1);
	}

	Double minXnew = min( vals[0], vals[N] );
	Double maxXnew = max( vals[0], vals[N] );
	Double stepnew = (maxXnew-minXnew) / Nnew;

	Spline sp(Nnew, minXnew, maxXnew);

	SplineGen spGen;

	int Nnow = N;
	Double stepNow = (maxX-minX)/Nnow; //how many points???

	//printNodes();

	for(int i = 0; i < Nnow+1; ++i) {
		Double Val = Eval(minX + i*stepNow, isLinear);
		//cout << "RADEKek "<< i <<" "<< minX + i*stepNow <<" "<<  Val <<" : "<<maxX<< endl;
		spGen.AddNode( Val,  minX + i*stepNow );
	}
	spGen.Finalize();

	for(int i = 0; i < Nnew+1; ++i) {
		//cout << "RADEK start " << i << " "<< minXnew + i*stepnew << endl;
		sp.vals[i] = spGen.EvalY( minXnew + i*stepnew );
		//cout << "RADEK end " << i << " "<< minXnew + i*stepnew << endl;
	}

	//Check if homogennous
	for(int i = 1; i < Nnew+1; ++i)
		if( sp.vals[i] < sp.vals[i-1]) {
			cout << "There is problem for i " <<i<< endl;
			exit(1);
		}

	return sp;

}

//Dummy function  GetInverseSpline(int Nnew, bool isLinear )
template<>
Spline<valarray<Double>> Spline<valarray<Double>>::GetInverseSpline(int Nnew, bool)
{

	Double minXnew = min( vals[0][0], vals[N][0] );
	Double maxXnew = max( vals[0][0], vals[N][0] );
	//Double stepnew = (maxXnew-minXnew) / Nnew;

	Spline<valarray<Double>> sp(Nnew, minXnew, maxXnew);

	return sp;
}


template<>
Spline<Double> Spline<Double>::InverseCubic(int Nnew)
{
	if(!isIncreasing()) {
		cout << "Spline does not homogenously grow" << endl;
		exit(1);
	}

	const Double eps = 1e-10;

	Double minXnew = min( vals[0], vals[N] );
	Double maxXnew = max( vals[0], vals[N] );
	Double stepnew = (maxXnew-minXnew) / Nnew;

	Spline sp(Nnew, minXnew, maxXnew);

	sp.vals[0] = minX;
	sp.vals[Nnew] = maxX;
	sp.valsD[0] = 1./valsD[0];
	sp.valsD[Nnew] = 1./valsD[N];

	for(int i = 1; i < Nnew; ++i) {
		double y = minXnew + stepnew*i;
		auto itU = lower_bound(vals.begin(), vals.end(), y,
				 [](Double a, Double b)  {return a < b;}  );

		if(itU == vals.end()) {
			if(abs(vals.back() - y) < eps )
				--itU;
			else {
				cout << "Problem with inversion " << endl;
				exit(1);
			}
		}

		auto itL = itU-1;

		if(itL < vals.begin()) {
			cout << "Danger " << endl;
			exit(1);
		}

		const double x1 = minX + (itL-vals.begin())*step;
		//const double x2 = minX + (itU-vals.begin())*step;

		const double y1 = *itL;
		const double y2 = *itU;
		const double y1D = valsD[itL-vals.begin()];
		const double y2D = valsD[itU-vals.begin()];

		assert( y1 <= y+eps && y-eps <= y2 );
		//cout << i << " "<< y1 <<" "<<y <<" "<< y2 << endl;

		//auto lin= [=](Double pos) {return y1*(1-pos) + y2*pos;};

		auto cubic= [=](Double pos) {
			Double pos1 = (1-pos)*(1-pos);
			Double h00 = (1+2*pos)*pos1;
			Double h10 = pos*pos1;
			Double h01 = pos*pos*(3.-2*pos);
			Double h11 = pos*pos*(pos-1);

		      return ( (h00*y1  + h01*y2 ) +
		             (h10*y1D + h11*y2D)*step );
						 };
		auto cubicD= [=](Double pos) {
			Double h00 = 6*pos*pos -6*pos;
			Double h10 = 3*pos*pos -4*pos+1;
			Double h01 = -6*pos*pos +6*pos;
			Double h11 = 3*pos*pos -2*pos;

		      return ( (h00*y1  + h01*y2 ) +
		             (h10*y1D + h11*y2D)*step )/step;
						 };

		Double posT = (y - y1) / (y2 - y1);

		for(int k = 0; k < 10; ++k) {
			posT = (y - cubic(posT) + (y2-y1)*posT ) / (y2-y1);
			//cout << "Iter "<<i <<" "<< posT << " "<< cubic(posT)-y<< endl;
		}
		sp.vals[i]  = x1 + step*posT;
		sp.valsD[i] = 1./cubicD(posT);

	}


	//Fill the derivative for first and last:





	return sp;
}

template<>
void Spline<Double>::CheckInverse(Spline<Double> &sp, bool isLin)
{
	Double ratSum = 0;
	int Ntry = 300;
	for(int i = 0; i < Ntry; ++i) {
		Double v = sp.minX + rand()/(RAND_MAX+0.)*(sp.maxX-sp.minX);
		Double myX = sp.Eval(v,isLin);
		Double trueV = Eval(myX, false);
		Double rat = (trueV-v)/v;
		if(abs(rat) > 1e-9)
			cout << i <<" "<<setprecision(12)<<v << " : "<<myX<<" : " << trueV <<" "<< rat<< endl;
		ratSum += abs(rat);
	}
	cout << "Mean ratio " << ratSum/Ntry << endl;
}





template<typename T>
Spline<T> Spline<T>::GetIntegratedSpline(bool lin) const
{
	Spline sp(N, minX, maxX );

	T sum = 0. * vals[0];

	sp.vals[0]  = 0;
	sp.valsD[0] = vals[0];
	for(int i = 1; i < N+1; ++i) {
		sum += (vals[i-1]+vals[i])*step/2.;
		if(!lin) sum += (valsD[i-1]-valsD[i]) * step*step/12.;
		sp.vals[i]  = sum;
		sp.valsD[i] = vals[i];
	}
	return sp;
}

template<typename T>
void Spline<T>::ShiftX(Double shift)
{
	minX += shift;
	maxX += shift;
}

template<typename T>
void Spline<T>::ShiftY(Double shift)
{
	for(int i = 0; i < N+1; ++i) 
		vals[i]  += shift;
}

//which instances are needed
template class Spline<Double>;
template class Spline<valarray<Double> >;
