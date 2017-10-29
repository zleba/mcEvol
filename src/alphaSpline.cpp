
#include "alphaSpline.h"

#include <iostream>
#include <iomanip>
#include <cassert>


Double alphaSpline::b[] = {};
Double alphaSpline::c[] = {};
Double alphaSpline::LnM[] = {};
Double alphaSpline::LnLambdaSave[] = {};
int   alphaSpline::iOrder = 2; //2--NLO

Double alphaSpline::b0[] = {};
Double alphaSpline::b1[] = {};
Double alphaSpline::b2[] = {};

/*
Double alphaSpline::func(Double x)
{
	return alphaS(x);
} 
Double alphaSpline::funcD(Double x)
{
	int nf = GetNf(x);
	Double L = 2.*(x-LnLambdaSave[nf]);
	return -2.*alphaS(x) / L;
} 
*/

//LO
Double alphaSpline::GetLambda1(int nf, double a)
{
	Double res =  - 1./(b0[nf]*a);
	return res;
}

//NLO
Double alphaSpline::GetLambda2(int nf, double a)
{
	Double C = b1[nf]/b0[nf];
	Double res = - 1./(b0[nf]*a) + C/b0[nf] * log(1 + 1./(C*a) );

	return res;
}


//NNLO
Double alphaSpline::GetLambda3(int nf, double a)
{
	Double B0 = b0[nf];
	Double B1 = b1[nf];
	Double B2 = b2[nf];


	Double C = B1/B0;
	Double D = sqrt(abs(4*B0*B2 - B1*B1));
	Double res = - 1./(B0*a);
	res       += C/(2*B0) * log( (B0+a*(B1+B2*a))/a/a );
	if(nf <= 5)
		res       += (B1*B1-2*B0*B2)/( B0*B0*D ) * atan( (2*B2*a+B1)/D );
	else
		res       +=-(B1*B1-2*B0*B2)/( B0*B0*D ) * atanh( (2*B2*a+B1)/D );
	
	return res;
}




///////////////////////////////////////////////////
//Fixing alphaS
///////////////////////////////////////////////////

void alphaSpline::FixParameters(int Order, Double value, int iFix, Double mZ)
{
	iOrder = Order;

	cout << "alphaS order : "<< iOrder << endl;
	if( LnM[4] == 0 || LnM[5] == 0 || LnM[6] == 0) {
		FixMasses( sqrt(2.8979862045712714), sqrt(24.904679585130495),	sqrt(28611.524722221096) );
		cout << "Quark masses set to the default values:"  << endl;
		cout << "mc="<< exp(LnM[4]) <<", mb= "<< exp(LnM[5]) <<", mt= "<< exp(LnM[6]) << endl;
	}
	else {
		cout << "Quark masses set to by user:"  << endl;
		cout << "mc="<< exp(LnM[4]) <<", mb= "<< exp(LnM[5]) <<", mt= "<< exp(LnM[6]) << endl;
	}


	b[6] = 21./6.;
	b[5] = 23./6.;
	b[4] = 25./6.;
	b[3] = 27./6.;

	if(iOrder >= 2) {
		c[6] = 39./42.;
		c[5] = 58./46.;
		c[4] = 77./50.;
		c[3] = 96./54.;
	}
	else {
		c[3]=c[4]=c[5]=c[6] = 0.;
	}

	for(int nf = 3; nf <= 6; ++nf) {
		b0[nf] = 11./2. -  1./3. * nf;
		b1[nf] = 51./2. - 19./6. * nf;
		b2[nf] = 2857./16. - 5033./144. * nf + 325./432. * nf*nf;
	}


	Double LambdaSave[7];

	// Find Lambda_4 at m_b, by requiring alphaS(nF=4,m_b) = alphaS(nF=5,m_b)

	if(mZ > 0) {
		//iFix = 5;

		Double LnLambda;
		Double a = value/(2.*M_PI);
		if(iOrder == 1)
			LnLambda = log(mZ)  +1./2*GetLambda1(iFix, a);//    - 1./(b[iFix]*a);
		else if(iOrder == 2)
			LnLambda = log(mZ)  +1./2*GetLambda2(iFix, a);// - 1./(b[iFix]*a) + c[iFix]/b[iFix] * log(1 + 1./(c[iFix]*a) );
		else
			LnLambda = log(mZ)  +1./2*GetLambda3(iFix, a);

		LambdaSave[iFix] = exp(LnLambda);
		LnLambdaSave[iFix] = LnLambda;
	}
	else {
		LambdaSave[iFix] = value;
		LnLambdaSave[iFix] = log(value);
	}

	int d = -1;
	for(int i=iFix; !(i < 3 && d==1); i-=d) {
		if(i > 6 && d == -1) { d=1; i = iFix; continue; }
		if(i == iFix) continue;
		//int d = i > iFix ? -1 : 1;
		Double mq = exp(LnM[ i>iFix ? i : i+1 ]);
		Double valueM       = alphaS(2*log(mq),i+d);

		double a = valueM/(2*M_PI);
		double LnLambda;
		if(iOrder == 1)
			LnLambda = log(mq)  +1./2*GetLambda1(i, a);//       - 1./(b[i]*a);
		else if(iOrder == 2)
			LnLambda = log(mq)  +1./2*GetLambda2(i, a);//       - 1./(b[i]*a) + c[i]/b[i] * log(1 + 1./(c[i]*a) );
		else {
			Double aCorr = a;
			if(d < 0) { //going up
				aCorr = a +  1./4 * 14/3. *a*a*a;
			}
			else { //going down
				for(int k = 0; k < 10; ++k) 
					aCorr = a - 1./4 * 14/3. *aCorr*aCorr*aCorr;
			}
			//Double aCorr = d > 0 ? a - Delta : a + Delta;
			LnLambda = log(mq)  +1./2*GetLambda3(i, aCorr);
		}

		LnLambdaSave[i] = LnLambda;
		LambdaSave[i] = exp(LnLambda);
	}
	cout << LambdaSave[4] << " "<< LambdaSave[5] << endl;

	for(int i = 3; i <= 6; ++i)
		LnLambdaSave[i] = log(LambdaSave[i]);

	cout << "AlphaS(mZ) =  " <<mZ << " "<< alphaS(2*log(mZ) ) << endl;
}


Double alphaSpline::alphaS(Double LnQ2, int nfUser)
{
	const int Niter = 10;
	//static const Double myLambda = sqrt( 91*91 * exp(- 12*M_PI/23./0.118) );

	//cout << "My lambda " << myLambda << endl;
	int nf = (nfUser==-1) ? GetNf(LnQ2) : nfUser;

	const Double &B0 = b0[nf];
	const Double &B1 = b1[nf];
	const Double &B2 = b2[nf];


	//cout << i << endl;
	Double LnLambda2 = 2.*LnLambdaSave[nf];

	Double logScale    = LnQ2 - LnLambda2;
	Double a = 1.  / (B0 * logScale);

	//cout << "RADEK before " << __LINE__ <<" "<< iOrder <<  endl;
	//NLO
	if(iOrder == 2) {
		Double C = B1/B0;
		for(int j = 0; j < Niter; ++j) {
			a = 1./ ( B0*logScale + C*log(1+ 1./(C*a) ) );
			//cout << "Iter "<< j << " "<< a << endl;
		}
	}
	//NNLO
	else if(iOrder == 3) {


		const Double C = B1/B0;
		const Double D = sqrt(abs(4*B0*B2 - B1*B1));
		//cout << "Test value " << 4*B0*B2 - B1*B1 << endl;

		const Double F = (B1*B1-2*B0*B2)/( B0*D );
		
		for(int j = 0; j < Niter; ++j) {
			Double res = B0*logScale;
			res       += C/2 * log( (B0+a*(B1+B2*a)) / (a*a) );
			if(nf <= 5)
				res       += F * atan( (2*B2*a+B1)/D );
			else
				res       +=-F * atanh( (2*B2*a+B1)/D );
			a = 1./res;
			//cout << "Iter "<< j << " "<< a << endl;
		}
	}
	else if(iOrder != 1) {
		cout << "Wrong order" << endl;
		assert(0);
	}

	//cout << "RADEK after " << __LINE__ << endl;

	return 2*M_PI * a;

}

pair<Double,Double> alphaSpline::alphaS2PiWithD(Double LnQ2, int nfUser)
{
	const int Niter = 10;
	//static const Double myLambda = sqrt( 91*91 * exp(- 12*M_PI/23./0.118) );

	//cout << "My lambda " << myLambda << endl;
	int nf = (nfUser==-1) ? GetNf(LnQ2) : nfUser;

	const Double &B0 = b0[nf];
	const Double &B1 = b1[nf];
	const Double &B2 = b2[nf];


	//cout << i << endl;
	Double LnLambda2 = 2.*LnLambdaSave[nf];

	Double logScale    = LnQ2 - LnLambda2;
	Double a = 1.  / (B0 * logScale);

	//cout << "RADEK before " << __LINE__ <<" "<< iOrder <<  endl;
	//NLO
	if(iOrder == 2) {
		Double C = B1/B0;
		for(int j = 0; j < Niter; ++j) {
			a = 1./ ( B0*logScale + C*log(1+ 1./(C*a) ) );
			//cout << "Iter "<< j << " "<< a << endl;
		}
	}
	//NNLO
	else if(iOrder == 3) {


		const Double C = B1/B0;
		const Double D = sqrt(abs(4*B0*B2 - B1*B1));
		//cout << "Test value " << 4*B0*B2 - B1*B1 << endl;

		const Double F = (B1*B1-2*B0*B2)/( B0*D );
		
		for(int j = 0; j < Niter; ++j) {
			Double res = B0*logScale;
			res       += C/2 * log( (B0+a*(B1+B2*a)) / (a*a) );
			if(nf <= 5)
				res       += F * atan( (2*B2*a+B1)/D );
			else
				res       +=-F * atanh( (2*B2*a+B1)/D );
			a = 1./res;
			//cout << "Iter "<< j << " "<< a << endl;
		}
	}
	else if(iOrder != 1) {
		cout << "Wrong order" << endl;
		assert(0);
	}

	//Calculate Derivative
	Double aD=0;
	if(iOrder == 1)      aD = B0;
	else if(iOrder == 2) aD = B0 + B1*a;
	else if(iOrder == 3) aD = B0 + B1*a + B2*a*a;

	aD *= -a*a;

	return make_pair(a, aD);

}
