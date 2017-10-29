#ifndef SplittingsIntegralSpline_H_
#define SplittingsIntegralSpline_H_

#include <functional>
#include "Spline.h"
#include <valarray>
#include "Splittings.h"
#include <cassert>
#include <iostream>

using namespace std;

class splittingGeneral : public Spline<valarray<Double>> {

public:
	splittingGeneral(int N_, Double LnQ2min, Double LnQ2max, Double zmin_,
	             function<Double(Double)> fZmax_,  function<Double(Double)> fZmaxD_ ) :
	  Spline<valarray<Double>>(N_, LnQ2min, LnQ2max),
	  zmin(zmin_), fZmax(fZmax_), fZmaxD(fZmaxD_) 
	{
		Nf = alphaSpline::GetNf( (LnQ2min+LnQ2max)/2. );
		FillTogether();
	}

	splittingGeneral() {}

	pair<valarray<Double>,valarray<Double>> funcBoth(Double LnQ2) override { 
		Double aS2pi, aS2piD;
		tie(aS2pi,aS2piD) = alphaSpline::alphaS2PiWithD(LnQ2, Nf);
		int &Order = alphaSpline::iOrder;
		Double zmax = fZmax(LnQ2);
		Double zmaxD = fZmaxD(LnQ2);

		valarray<Double> PxxI0, PxxI1, PxxI2;

		PxxI0 =  getSplittingsI0(Nf, zmin, zmax);
		if(Order >= 2) PxxI1 = getSplittingsI1(Nf, zmin, zmax);
		if(Order >= 3) PxxI2 = getSplittingsI2(Nf, zmin, zmax);

		valarray<Double> res  = aS2pi*             PxxI0;
		if(Order >= 2)   res += aS2pi*aS2pi*       PxxI1;
		if(Order >= 3)   res += aS2pi*aS2pi*aS2pi* PxxI2;

		//Calculate Derivative
		valarray<Double> resD  =   aS2piD*             PxxI0;
		if(Order >= 2)   resD += 2*aS2piD*aS2pi *      PxxI1;
		if(Order >= 3)   resD += 3*aS2piD*aS2pi*aS2pi* PxxI2;

		assert(resD[0] < 1e6);

		valarray<Double> Pxx0, Pxx1, Pxx2;

		Pxx0 = getSplittings0(Nf, zmax);
		if(Order >= 2) Pxx1 = getSplittings1(Nf, zmax);
		if(Order >= 3) Pxx2 = getSplittings2(Nf, zmax);

		/*
		assert(Pxx0[0] < 1e6);
		cout << "Pxx1[0]=" << Pxx1[0] << " "<< zmax<<" "<< Nf<< endl;
		cout << "Test " << getSplittings1(Nf, zmax)[0] << endl;
		assert(Pxx1[0] < 1e6);
		assert(Pxx2[0] < 1e6);
		*/


		               resD += aS2pi*              zmaxD*Pxx0;
		if(Order >= 2) resD += aS2pi*aS2pi*        zmaxD*Pxx1;
		if(Order >= 3) resD += aS2pi*aS2pi*aS2pi*  zmaxD*Pxx2;

		assert(resD[0] < 1e6);

		return make_pair(res, resD);

	}

	Spline<Double> sumComponents(int iStart, int iEnd)
	{
		Spline<Double> spl(N, minX, maxX);

		for(int k = 0; k < N+1; ++k) {
			Double val=0, valD=0;
			for(int i = iStart; i <= iEnd; ++i) {
				val  += vals[k][i];
				valD += valsD[k][i];
			}
			spl.setNode(k, val, valD);
		}
		return spl;
	}




static valarray<Double>  getSplittingsI0(int nf, Double zMin, Double zMax )
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg   = SplittingPgg0(zMin, zMax);
	Pqg   = 2*nf*SplittingPqg0(zMin, zMax);
	Pgq   = SplittingPgq0(zMin, zMax);
	PqIqI = SplittingPqq0(zMin, zMax);
	PqIqJ  = 0;
	PqIqbI = 0;
	PqIqbJ = 0;

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});
}

static valarray<Double>  getSplittingsI1(int nf, Double zMin, Double zMax )
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg    =  SplittingPgg1(zMin, zMax, nf);
	Pqg    =  2*nf*SplittingPqg1(zMin, zMax, nf);
	Pgq    =  SplittingPgq1(zMin, zMax, nf);
	PqIqI  =  SplittingPqIqI1(zMin, zMax, nf);
	PqIqJ  =  (nf-1)*SplittingPqIqJ1(zMin, zMax, nf);
	PqIqbI =  SplittingPqIqbI1(zMin, zMax, nf);
	PqIqbJ =  (nf-1)*SplittingPqIqJ1(zMin, zMax, nf);

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});
}

static valarray<Double>  getSplittingsI2(int nf, Double zMin, Double zMax )
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg    =  SplittingPgg2(zMin, zMax, nf);
	Pqg    =  2*nf*SplittingPqg2(zMin, zMax, nf);
	Pgq    =  SplittingPgq2(zMin, zMax, nf);
	PqIqI  =  SplittingPqIqI2(zMin, zMax, nf);
	PqIqJ  =  (nf-1)*SplittingPqIqJ2(zMin, zMax, nf);
	PqIqbI =  SplittingPqIqbI2(zMin, zMax, nf);
	PqIqbJ =  (nf-1)*SplittingPqIqbJ2(zMin, zMax, nf);

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});

}


static valarray<Double>  getSplittings0(int nf, Double z)
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg   = SplittingPgg0(z);
	Pqg   = 2*nf*SplittingPqg0(z);
	Pgq   = SplittingPgq0(z);
	PqIqI = SplittingPqq0(z);
	PqIqJ  = 0;
	PqIqbI = 0;
	PqIqbJ = 0;

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});
}

static valarray<Double>  getSplittings1(int nf, Double z)
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg    =  SplittingPgg1(z, nf);
	Pqg    =  2*nf*SplittingPqg1(z, nf);
	Pgq    =  SplittingPgq1(z, nf);
	PqIqI  =  SplittingPqIqI1(z, nf);
	PqIqJ  =  (nf-1)*SplittingPqIqJ1(z, nf);
	PqIqbI =  SplittingPqIqbI1(z, nf);
	PqIqbJ =  (nf-1)*SplittingPqIqJ1(z, nf);

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});
}

static valarray<Double>  getSplittings2(int nf, Double z)
{
	Double Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;

	//aS is divided by 2pi
	Pgg    =  SplittingPgg2(z, nf);
	Pqg    =  2*nf*SplittingPqg2(z, nf);
	Pgq    =  SplittingPgq2(z, nf);
	PqIqI  =  SplittingPqIqI2(z, nf);
	PqIqJ  =  (nf-1)*SplittingPqIqJ2(z, nf);
	PqIqbI =  SplittingPqIqbI2(z, nf);
	PqIqbJ =  (nf-1)*SplittingPqIqbJ2(z, nf);

	return valarray<Double>({Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI});

}


private:
	int Nf;
	Double zmin;
	function<Double(Double)> fZmax, fZmaxD;

};


class splittingUnintegrated : public Spline<valarray<Double>> {
public:
	splittingUnintegrated(int N_, int nf, Double zMin, Double zMax ) :
	  Spline<valarray<Double>>(N_,  atanh(2*zMin-1),  atanh(2*zMax-1))
	{
		Nf = nf;
		Fill();
	}

	splittingUnintegrated() {}

	valarray<Double> func(Double r) override { 
		Double z = (1 + tanh(r) )/2.;

		valarray<Double> Pxx0, Pxx1(0.,7), Pxx2(0.,7);

		int &Order = alphaSpline::iOrder;

		              Pxx0 = splittingGeneral::getSplittings0(Nf, z);
		if(Order >=2) Pxx1 = splittingGeneral::getSplittings1(Nf, z);
		if(Order >=3) Pxx2 = splittingGeneral::getSplittings2(Nf, z);

		valarray<Double> Pxx(7*3);

		for(int i = 0; i < 7; ++i) {
			Pxx[i]    = Pxx0[i];
			Pxx[i+7]  = Pxx1[i];
			Pxx[i+14] = Pxx2[i];
		}
		return Pxx;
	}
	Double EvalOne(int type, Double aS2pi, Double x)
	{
		auto myfun =[&](int k) {
			return aS2pi*(vals[k][type]+ aS2pi*(vals[k][type+7]+aS2pi*vals[k][type+14]));
		};

		if(x <= minX)
			return myfun(0);
		if(x >= maxX) {
			return myfun(N);
		}

		Double pos = (x - minX) / step;
		int i = (int) pos;
		pos -= i;

		i = max(0, i);
		if(i >= N)
			return myfun(N);
		//i = min(N-1, i);

		return  myfun(i)*(1-pos) + myfun(i+1)*pos;

	}
	void TestPrecision(int type, int nTry) {
		for(int i = 0; i < nTry; ++i) {
			double r = minX + (maxX-minX)/nTry*i;
			double aT = 0.1/2./M_PI;
			double approx = EvalOne(type, aT, r);
			valarray<Double> prec = func(r);
			double org = aT*(prec[type]+ aT*(prec[type+7]+aT*prec[type+14]));
			double rat = (approx-org)/org;
			cout << "RADE "<< i <<" "<< org <<" "<< rat << endl;
		}
	}

	private:
	int Nf;

};



#endif 
