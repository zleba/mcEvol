#include "sudakovSpline.h"
#include "tmd.h"

#include <cmath>
#include <iostream>

#include "integration.h"

sudakovSpline::sudakovSpline(int _N, int Flavour, int IOrder, Double _minX, Double _maxX) :
               Spline<Double>(_N, _minX, _maxX), flavour(Flavour), iOrder(IOrder)
{ 
	Double &zmin = TMD::zmin;

	//cout << "RADEK " << exp(minX) <<" "<< exp(maxX) << endl;

	if(zmin < 0) {
		zmin = 1e-6;
		cout << "Value of zmin set to the default value " << zmin << endl;
	}
	if(TMD::fZmax(5) < 0) {
		//zmax = 1 - 1e-5;

		TMD::fZmax  = [](Double){return 1.-1.e-5;};
		TMD::fZmaxD = [](Double){return 0.;};

		cout << "Value of 1-zmax set to the default value " << 1-TMD::fZmax(5) << endl;
	}


	//Init alphaS
	//alphaSpline::FixParameters(iOrder, 0.118, 5, sqrt(8315.25) );

	Double zmax = TMD::fZmax((_minX+_maxX)/2.);

	//LO splittings
	Pgg0 = SplittingPgg0(zmin, zmax);
	Pqq0 = SplittingPqq0(zmin, zmax);
	Pqg0 = SplittingPqg0(zmin, zmax);
	Pgq0 = SplittingPgq0(zmin, zmax);

	//NLO splittings
	for(int nf=3; nf <=6; ++nf) {
		Pgg1[nf]   = SplittingPgg1(zmin, zmax,nf);
		Pqg1[nf]   = SplittingPqg1(zmin, zmax,nf);
		Pgq1[nf]   = SplittingPgq1(zmin, zmax,nf);
		PqIqJ1[nf] = SplittingPqIqJ1(zmin, zmax,nf);
		PqIqI1[nf] = SplittingPqIqI1(zmin, zmax,nf);
		PqIqbI1[nf]= SplittingPqIqbI1(zmin, zmax,nf);
	}

	//NNLO splittings
	for(int nf=3; nf <=6; ++nf) {
		Pgg2[nf]   = SplittingPgg2(zmin, zmax,nf);
		Pqg2[nf]   = SplittingPqg2(zmin, zmax,nf);
		Pgq2[nf]   = SplittingPgq2(zmin, zmax,nf);
		PqIqJ2[nf] = SplittingPqIqJ2(zmin, zmax,nf);
		PqIqbJ2[nf]= SplittingPqIqbJ2(zmin, zmax,nf);
		PqIqI2[nf] = SplittingPqIqI2(zmin, zmax,nf);
		PqIqbI2[nf]= SplittingPqIqbI2(zmin, zmax,nf);
	}


	LogMc2 = 2.*alphaSpline::LnM[4];
	LogMb2 = 2.*alphaSpline::LnM[5];
	LogMt2 = 2.*alphaSpline::LnM[6];


	Fill(); 
}


Double sudakovSpline::func(Double x)
{
	Double err=0;

	//The lowest limit of integration corresponds to edge of Sudakov Spline grid
	const Double &LogQ2min = minX;

	auto integrand = [&](Double xx) { return InsideSud(xx); };

	Double res  = 0;
	//Over all flavours
	if(LogQ2min < LogMc2 ) {
		if(x > LogMt2) {
			res += Integral( integrand, LogQ2min, LogMc2, err );
			res += Integral( integrand, LogMc2,   LogMb2, err );
			res += Integral( integrand, LogMb2,   LogMt2, err );
			res += Integral( integrand, LogMt2, x, err );
		}
		else if(x > LogMb2) {
			res += Integral( integrand, LogQ2min, LogMc2, err );
			res += Integral( integrand, LogMc2,   LogMb2, err );
			res += Integral( integrand, LogMb2,   x, err );
		}
		else if(x > LogMc2) {
			res += Integral( integrand, LogQ2min, LogMc2, err );
			res += Integral( integrand, LogMc2,   x, err );
		}
		else {
			res += Integral61( integrand, LogQ2min, x, err );
		}
	}
	else if(LogQ2min < LogMb2) {

		if(x > LogMt2) {
			res += Integral( integrand, LogQ2min, LogMb2, err );
			res += Integral( integrand, LogMb2,   LogMt2, err );
			res += Integral( integrand, LogMt2, x, err );
		}
		else if(x > LogMb2) {
			res += Integral( integrand, LogQ2min, LogMb2, err );
			res += Integral( integrand, LogMb2,   x, err );
		}
		else {
			res += Integral( integrand, LogQ2min, x, err );
		}
	}
	else {
		cout << "Too large starting scale " << exp(LogQ2min)<<" - not implemented "<<endl;
		cout << exp(0.5*LogMc2) << " "<< exp(0.5*LogMb2) <<" "<< exp(0.5*LogMt2) << endl;
		exit(1);
	}

	if(err > 1e-9)
		cout << "Error " << err << endl;

	//cout << "Err " << err << endl;
	return res;
} 

Double sudakovSpline::funcD(Double)
{
	return 0.;
} 


void sudakovSpline::getSplittings(Double aS, int nf, Double &Pgg, Double &Pqg, Double &Pgq, Double &PqIqJ, Double &PqIqbJ, Double &PqIqI, Double &PqIqbI) const
{
	//aS is divided by 2pi
	Pgg = Pgg0;
	Pqg = Pqg0;
	Pgq = Pgq0;
	PqIqI = Pqq0;
	PqIqJ  = 0;
	PqIqbJ = 0;
	PqIqbI = 0;
	if(iOrder >= 2) {
		//Double aS = alphaSpline::alphaS(LnQ2);
		//int nf = alphaSpline::GetNf(LnQ2);
		Pgg   += aS * Pgg1[nf];
		Pqg   += aS * Pqg1[nf];
		Pgq   += aS * Pgq1[nf];
		PqIqI += aS * PqIqI1[nf];
		PqIqJ += aS * PqIqJ1[nf];
		PqIqbJ+= aS * PqIqJ1[nf]; //At NLO identical to PqIqJ 
		PqIqbI+= aS * PqIqbI1[nf];
	}
	if(iOrder >= 3) {
		Pgg   += aS*aS * Pgg2[nf];
		Pqg   += aS*aS * Pqg2[nf];
		Pgq   += aS*aS * Pgq2[nf];
		PqIqI += aS*aS * PqIqI2[nf];
		PqIqJ += aS*aS * PqIqJ2[nf];
		PqIqbJ+= aS*aS * PqIqbJ2[nf];
		PqIqbI+= aS*aS * PqIqbI2[nf];
	}

}



Double sudakovSpline::InsideSud(Double x) const
{
	int Nf = alphaSpline::GetNf(x);
	Double aS2pi = alphaSpline::alphaS( x )/(2.*M_PI);

	Double res = 0;

	//LO part
	if(flavour == 0)
		res += aS2pi * ( Pgg0 + 2*Nf*Pqg0 );
	else
		res += aS2pi * ( Pqq0 + Pgq0 );

	if(iOrder == 1)
		return res;

	//NLO part
	if(flavour == 0)
		res += aS2pi*aS2pi * (Pgg1[Nf] + 2*Nf*Pqg1[Nf] );
	else
		res += aS2pi*aS2pi * (Pgq1[Nf] + PqIqI1[Nf] + PqIqbI1[Nf] + 2*(Nf-1)*PqIqJ1[Nf] );

	if(iOrder == 2)
		return res;

	//NNLO part
	if(flavour == 0)
		res += aS2pi*aS2pi*aS2pi * (Pgg2[Nf] + 2*Nf*Pqg2[Nf] );
	else
		res += aS2pi*aS2pi*aS2pi * (Pgq2[Nf] + PqIqI2[Nf] + PqIqbI2[Nf] + (Nf-1)*PqIqJ2[Nf] + (Nf-1)*PqIqbJ2[Nf] );


	return res;

}


void sudakovSpline::printNodes() const
{
	cout << "Sudakov PrintOut (N = " << N <<")" <<endl;
	for(int i = 0; i < N+1; ++i) 
		cout << i <<" "<< minX+i*step <<" "<< exp((minX+i*step)) <<" "<<  exp(0.5*(minX+i*step)) <<" "<< vals[i] << endl;


}
