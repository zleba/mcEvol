#include "sudakovSplineSingle.h"
#include "tmd.h"

#include <cmath>
#include <iostream>
#include <cassert>

#include "integration.h"

sudakovSplineSingle::sudakovSplineSingle(int _N, int Flavour, int IOrder, Double LnQ2min, Double LnQ2max) :
                 Spline(_N, LnQ2min, LnQ2max), flavour(Flavour), iOrder(IOrder)
{ 
	Double &zmin = TMD::zmin;
	Double &zmax = TMD::zmax;
	Nf = alphaSpline::GetNf( (LnQ2min+LnQ2max)/2. );

	//cout << "RADEK " << exp(minX) <<" "<< exp(maxX) << endl;

	if(zmin < 0) {
		zmin = 1e-6;
		cout << "Value of zmin set to the default value " << zmin << endl;
	}
	if(zmax < 0) {
		zmax = 1 - 1e-5;
		cout << "Value of 1-zmax set to the default value " << 1-zmax << endl;
	}


	//Init alphaS
	//alphaSpline::FixParameters(iOrder, 0.118, 5, sqrt(8315.25) );


	//LO splittings
	if(flavour == 0) {
		Pgg0 = SplittingPgg0(zmin, zmax);
		Pqg0 = SplittingPqg0(zmin, zmax);
	}
	else {
		Pqq0 = SplittingPqq0(zmin, zmax);
		Pgq0 = SplittingPgq0(zmin, zmax);
	}

	//NLO splittings
	if(iOrder >= 2) {
		if(flavour == 0) {
			Pgg1   = SplittingPgg1(zmin, zmax,Nf);
			Pqg1   = SplittingPqg1(zmin, zmax,Nf);
		}
		else {
			Pgq1       = SplittingPgq1(zmin, zmax,Nf);
			PqIqJ1     = SplittingPqIqJ1(zmin, zmax,Nf);
			PqIqI1     = SplittingPqIqI1(zmin, zmax,Nf);
			PqIqbI1    = SplittingPqIqbI1(zmin, zmax,Nf);
		}
	}

	//NNLO splittings
	if(iOrder >= 3) {
		if(flavour == 0) {
			Pgg2       = SplittingPgg2(zmin, zmax,Nf);
			Pqg2       = SplittingPqg2(zmin, zmax,Nf);
		}
		else {
			Pgq2       = SplittingPgq2(zmin, zmax,Nf);
			PqIqJ2     = SplittingPqIqJ2(zmin, zmax,Nf);
			PqIqbJ2    = SplittingPqIqbJ2(zmin, zmax,Nf);
			PqIqI2     = SplittingPqIqI2(zmin, zmax,Nf);
			PqIqbI2    = SplittingPqIqbI2(zmin, zmax,Nf);
		}
	}


	//Check that the range contains only one Nf
	for(int f = 4; f <=6; ++f) {
		if( LnQ2min < 2.*alphaSpline::LnM[f] && LnQ2max > 2*alphaSpline::LnM[f])
			assert(0);
	}
	//Check that the Nf is correct
	assert( alphaSpline::GetNf( (LnQ2min+LnQ2max)/2.) == Nf );

	Fill(); 
}


Double sudakovSplineSingle::func(Double x)
{
	Double err=0;

	//The lowest limit of integration corresponds to edge of Sudakov Spline grid
	const Double &LogQ2min = minX;

	auto integrand = [&](Double xx) { return InsideSud(xx); };

	Double res = Integral( integrand, LogQ2min, x, err );

	if(err > 1e-9)
		cout << "Error " << err << endl;

	//cout << "Err " << err << endl;
	return res;
} 

Double sudakovSplineSingle::funcD(Double)
{
	return 0.;
} 

void sudakovSplineSingle::getGluonSplittings(Double aS2pi, Double *Pgg_, Double *Pqg_) const
{
	Double &Pgg = *Pgg_;
	Double &Pqg = *Pqg_;
	//aS is in argument divided by 2pi
	Pgg = Pgg0;
	Pqg = Pqg0;
	if(iOrder >= 2) {
		Pgg   += aS2pi * Pgg1;
		Pqg   += aS2pi * Pqg1;
	}
	if(iOrder >= 3) {
		Pgg   += aS2pi*aS2pi * Pgg2;
		Pqg   += aS2pi*aS2pi * Pqg2;
	}

}

void sudakovSplineSingle::getQuarkSplittings(Double aS2pi, Double *Pgq_, Double *PqIqJ_, Double *PqIqbJ_, Double *PqIqI_, Double *PqIqbI_) const
{
	Double &Pgq    =  *Pgq_;
	Double &PqIqJ  =  *PqIqJ_;
	Double &PqIqbJ =  *PqIqbJ_;
	Double &PqIqI  =  *PqIqI_;
	Double &PqIqbI =  *PqIqbI_;

	//aS is divided by 2pi
	Pgq = Pgq0;
	PqIqI = Pqq0;
	PqIqJ  = 0;
	PqIqbJ = 0;
	PqIqbI = 0;
	if(iOrder >= 2) {
		Pgq   += aS2pi * Pgq1;
		PqIqI += aS2pi * PqIqI1;
		PqIqJ += aS2pi * PqIqJ1;
		PqIqbJ+= aS2pi * PqIqJ1; //At NLO identical to PqIqJ 
		PqIqbI+= aS2pi * PqIqbI1;
	}
	if(iOrder >= 3) {
		Pgq   += aS2pi*aS2pi * Pgq2;
		PqIqI += aS2pi*aS2pi * PqIqI2;
		PqIqJ += aS2pi*aS2pi * PqIqJ2;
		PqIqbJ+= aS2pi*aS2pi * PqIqbJ2;
		PqIqbI+= aS2pi*aS2pi * PqIqbI2;
	}

}



Double sudakovSplineSingle::InsideSud(Double x) const
{
	//int Nf = alphaSpline::GetNf(x);
	Double aS2pi = alphaSpline::alphaS( x, Nf )/(2.*M_PI);

	Double res = 0;

	if(flavour == 0) {
		Double Pgg, Pqg;
		getGluonSplittings(aS2pi, &Pgg, &Pqg);
		res = aS2pi * ( Pgg + 2*Nf*Pqg );
	}
	else {
		Double Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI;
		getQuarkSplittings(aS2pi, &Pgq, &PqIqJ, &PqIqbJ, &PqIqI, &PqIqbI);
		res = aS2pi * (Pgq + PqIqI + PqIqbI + (Nf-1)*PqIqJ + (Nf-1)*PqIqbJ );
	}

	return res;

}


void sudakovSplineSingle::printNodes() const
{
	cout << "Sudakov PrintOut (N = " << N <<")" <<endl;
	for(int i = 0; i < N; ++i) 
		cout << i <<" "<< minX+i*step <<" "<< exp((minX+i*step)) <<" "<<  exp(0.5*(minX+i*step)) <<" "<< vals[i] << endl;


}
