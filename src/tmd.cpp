#include <iostream>

#include "sudakovSpline.h"
#include "alphaSpline.h"
#include "Sudakov.h"
#include "Splittings.h"

#include "tmd.h"
#include <random>

#include <cmath>
#include <cstdlib>
#include <numeric>
#include <cassert>
#include <iomanip>


//static TRandom *ranRoot = new TRandom3();

std::random_device rDev;
mt19937_64 *rng = new mt19937_64(rDev());



Double TMD::zmin   = -1;
function<Double(Double)> TMD::fZmax  = [](Double){return -1.;};
function<Double(Double)> TMD::fZmaxD = [](Double){return  0.;};


inline Double Rand()
{
	//return rand()/(RAND_MAX+0.);

    uniform_real_distribution<double> unif;
    //cout << unif(*rng) << endl; 
    return unif(*rng);

	//return ranRoot->Uniform();
}
int RandI(int iMax) {

    uniform_int_distribution<int> unii(0, iMax-1);
    return unii(*rng);
	//return ranRoot->Integer(iMax);
}


TMD::TMD(int N, int IOrder,  Double q2min, Double q2max) : sudSplG(N, 0, IOrder, log(q2min), log(q2max) ),
                                                           sudSplQ(N, 1, IOrder, log(q2min), log(q2max) )
{
	//ranRoot->SetSeed(0);
	iOrder = IOrder;
	sudInvG = sudSplG.GetInverseSpline(2*N);
	sudInvQ = sudSplQ.GetInverseSpline(2*N);

	LogQ2min = log(q2min);
	LogQ2max = log(q2max);


	sudMaxG = sudSplG.Eval( LogQ2max );
	sudMaxQ = sudSplQ.Eval( LogQ2max );

	cout << "Evolution order : " << iOrder << endl;


	cout <<"Zmin, Zmax " << zmin<<" "<< fZmax((q2min+q2max)/2.) << endl;
	sudakovNew.Init(5000, LogQ2min, LogQ2max, zmin, fZmax, fZmaxD);


	/*
	int Ndiv = 1024;
	double step = log(14000.*14000./2.)/(Ndiv-1.);
	for(int i =0; i <Ndiv; i++) {
		double q2 = exp(log(2.)+i*step);
		//if((i+1) %4 == 0)
		cout << i<<" "<<setprecision(10)<<(q2) << " "<< sudSplQ.func(log(q2)) << " "<< sudSplG.func(log(q2)) << endl;

	}
	exit(1);
	*/
	/*
	cout <<"Gluon 113.13 "<<  sudSplG.Eval(2*log(14000) ) << endl;
	cout <<"Quark        "<<  sudSplQ.Eval(2*log(14000) ) << endl;
	exit(0);
	*/

	sPDF.InitDistribution();

}

Double TMD::getSudakov(Double lnQ2, int fl)
{
	if(fl == 0)
		return  sudakovNew.GetSudakovForGluon(lnQ2);
	else
		return  sudakovNew.GetSudakovForQuark(lnQ2);
}

void TMD::PrintSudakovs() const
{
	cout << "Gluon sudakov " << endl;
	sudSplG.printNodes();
	cout << "Gluon inverse Sudakov " << endl;
	sudInvG.printNodes();


	cout << "Quardk sudakov " << endl;
	sudSplQ.printNodes();

}


void TMD::CompareSudakovs()
{

	Sudakov sudNew;
	sudNew.Init(3000, LogQ2min, LogQ2max, zmin,
	                [](Double){return 1-1e-5;}, [](Double){return 0;} );

	Double step = (LogQ2max-LogQ2min)/sudSplG.getN();
	
	//int category;
	//cout << sudNew.GetLnQ2forGluon(1.8, category) << " "<< 2*alphaSpline::LnM[4]<< endl;
	//cout << "Category " << category << endl;
	//exit(1);

	for(int i = 0; i < sudSplG.getN(); ++i) {
		Double logQ2 = LogQ2min + i*step;
		//Double sudGOrg = sudSplG.Eval(logQ2);
		//Double sudGNew = sudNew.GetSudakovForGluon(logQ2);
		Double sudQOrg = sudSplQ.Eval(logQ2);
		Double sudQNew = sudNew.GetSudakovForQuark(logQ2);

		//Double ratG = (sudGNew-sudGOrg)/sudGOrg;
		Double ratQ = (sudQNew-sudQOrg)/sudQOrg;
		cout << i<<" "<<exp(0.5*logQ2)<<" "<< sudQOrg<<" "<< sudQNew <<" "<<ratQ << endl;
		//cout << i<<" "<< sudQOrg<<" "<< sudQNew <<" "<<ratQ << endl;
		//if(ratG == ratG && ratG > 1e-7) {
			//cout << "Huge ratio " << ratG << endl;
			//exit(0);
		//}
	}

}

Double TMD::GetZ(Double LnQ2, Double &z) const
{
	Double zmax0 = fZmax(LnQ2);
	Double zmin0 = zmin / br.x;

	if(zmin0 > zmax0) {
		return 0;
	}

	Double rat = (1.0 - zmin0)/(1.0 - zmax0);
	Double LogRat = log(rat);

	z = 1. - (1.-zmax0) * exp(Rand()*LogRat);

	//Weight
	return LogRat*(1-z);
}

Double TMD::GetZlin(Double LnQ2, Double &z) const
{
	Double zmax0 = fZmax(LnQ2);
	Double zmin0 = zmin / br.x;

	if(zmin0 > zmax0) {
		return 0;
	}
	z = zmin0 + (zmax0-zmin0) *Rand();

	return (zmax0-zmin0);

}

void TMD::Init(int fl)
{

	br.sud = 0;
	br.LogQ2 = LogQ2min;

	if(fl == -50)
		br.weight = sPDF.genFirst(true, br.x, br.id);
	else {
		br.weight = 1.;
		br.id = fl;

        br.x = pow(10,-Rand()/4.);
		//br.x = 1.;
	}

	Double r   = 0;//sqrt(-2*log(Rand()) );
	Double phi = 2*M_PI * Rand();

	Double px = r*cos(phi);
	Double py = r*sin(phi);
	for(int i = 0; i < 3; ++i) {
		br.pT[i].setXY( px, py );
		br.pTmy[i].setXY( px, py );
	}

}



void StartPDF::InitDistribution()
{
	//Double Int[7];
	Double &dI = Int[3+1];
	Double &uI = Int[3+2];
	Double &sI = Int[3+3];

	Double &dbI = Int[3-1];
	Double &ubI = Int[3-2];
	Double &sbI = Int[3-3];

	Double &gI = Int[3+0];

	AD = 3.064320;
	ADbar= 0.1939875;
	AU = 5.107200;
	
    auto beta = [](double a, double b) {
        double err;
        return Integral61([&](double z) {return pow(z*z,a)*pow(1-z*z,b) * 2*z;}, 0, 1, err);
    };
	//auto beta = [](Double a, Double b) { return TMath::Beta(1+a, 1+b); };

	dbI = ADbar * beta(-0.1, 6);
	ubI = ADbar * beta(-0.1, 7);
	sbI = 0.2 * (dbI + ubI);
	
	dI = AD* beta(0.8,4) + dbI;
	uI = AU* beta(0.8,3) + ubI;
	sI = sbI;
	
	gI = 1.7 * beta(-0.1, 5);

	//std::function< Double(Double) > aa[10];
	//aa[0] = [](Double b) { return 2*b; } ;

	sum = accumulate(Int, Int+7, 0.0);
	for(int i = 0; i < 7; ++i)
		IntSum[i] = accumulate(Int, Int+i+1, 0.0);

	auto &d = pdf[3+1];
	auto &u = pdf[3+2];
	auto &s = pdf[3+3];

	auto &db = pdf[3-1];
	auto &ub = pdf[3-2];
	auto &sb = pdf[3-3];

	auto &g = pdf[3+0];

	ub = [=](Double x) { return ADbar*pow(x,-0.1)*pow(1-x,7); };
	db = [=](Double x) { return ADbar*pow(x,-0.1)*pow(1-x,6); };
	sb = [=](Double x) { return 0.2*(ub(x)+db(x)); };

	u  = [=](Double x) { return ub(x) + AU*pow(x,0.8)*pow(1-x,3); };
	d  = [=](Double x) { return db(x) + AD*pow(x,0.8)*pow(1-x,4); };
	s  = sb;

	g  = [](Double x) { return 1.7 * pow(x,-0.1)*pow(1-x,5); };

}

Double StartPDF::genFirst(bool DistributeEqualy, Double &x, int &id)
{
	Double weight;
	static const Double xmin = 1e-5;

	//Linear Distribution
	//x = xmin + (1.-xmin)*Rand();
	//weight = 1. - xmin;

	//Power Distribution
	static const Double a = -0.8;
	x = pow( pow(xmin,a+1) + Rand()*(1 - pow(xmin,a+1) ), 1./(a+1) );
	weight =  pow(x,-a) * ( 1 - pow(xmin,a+1) ) / (a+1);


	if(DistributeEqualy) {
		id = RandI(7) - 3;

        double pdfWgt = weight * pdf[id+3](x) * 7.;
        //x = pow(10, (ceil(log10(x) * 4) - Rand())/4.); //to do step-like input (only for testing of convolution)
        return pdfWgt;


		//return weight * pdf[id+3](x) * 7.;
	}
	else {
		Double idRand = sum * Rand();

		//int id;
		for(id = -3; id < 3; ++id)
			if(idRand < IntSum[id+3])
				break;

		return weight * pdf[id+3](x) / (Int[id+3]/sum );
	}
	
}

Vec2 RandVector(double r)
{
	Vec2 res;
	Double a, b, r2;
	do {
		a = Rand() - 0.5;
		b = Rand() - 0.5;
		r2 = a*a+b*b;
	} while(r2 > 0.25);
	
	Double fact = r/sqrt(r2);
	return Vec2(a*fact, b*fact);
}

bool TMD::Evolve()
{
	brOld = br;

	Double sudMax = (br.id == 0) ? sudMaxG : sudMaxQ;
	Double r = Rand();
	//Double Delta = -log( 1 - r*(1-exp(br.sud-sudMax ) ) );
	Double Delta = -log(r);

	if( br.sud-sudMax >= 0) {
		cout << "r=" << r<<", br.sud=" << br.sud <<",sudMax=" << sudMax << endl;
		cout << "SudMaxG "<< sudMaxG << ", sudMaxQ="<< sudMaxQ <<" id= "<<br.id <<  endl;
		exit(1);
	}

	Double sudNew = br.sud + Delta;

	if(sudNew >= sudMax) {
		sudNew = sudMax;
		//cout << "ScaleNew "<< (br.id == 0  ? sudInvG.Eval(sudNew) : sudInvQ.Eval(sudNew) ) << endl;
		//cout << "ScaleMax "<< (br.id == 0  ? sudInvG.Eval(sudMax) : sudInvQ.Eval(sudMax) ) << endl;
	}

	//Double LogQ2newOld = (br.id == 0) ? sudInvG.Eval(sudNew) : sudInvQ.Eval(sudNew);

	int cat;
	Double LogQ2new = (br.id == 0) ? sudakovNew.GetLnQ2forGluon(sudNew, cat) :
	                                 sudakovNew.GetLnQ2forQuark(sudNew, cat);

	//cout << "RADEK "<< LogQ2newOld <<" "<< LogQ2new << endl;
	/*
	if(br.id == 0) {
		int Cat;
		Double scaleOrg = sudInvG.Eval(sudNew);
		Double scaleNew = sudakovNew.GetLnQ2forGluon(sudNew, Cat);
		Double rat = (scaleNew-scaleOrg)/scaleOrg;
		if(rat > 1e-4)
			cout <<"DANGER " <<  rat << endl;
		//cout << "Gluons " <<sudNew<<" : "<< scaleOrg <<" "<< scaleNew <<" "<< Cat<< endl;
		
		
	}
	else {
		int Cat;
		Double scaleOrg = sudInvQ.Eval(sudNew);
		Double scaleNew = sudakovNew.GetLnQ2forQuark(sudNew, Cat);
		Double rat = (scaleNew-scaleOrg)/scaleOrg;
		if(rat > 1e-4)
			cout <<"DANGER " <<  rat << endl;
		//cout << "Quarks " << sudNew<<" : "<<scaleOrg <<" "<< scaleNew <<" "<<Cat<< endl;
	}
	*/


	if( LogQ2new <= br.LogQ2 && cat%2 != 1) {
		cout << "Strange " << br.LogQ2 << " "<< LogQ2new  <<" "<< cat<< endl;
		//exit(1);
	}
	//assert( LogQ2new > br.LogQ2 );

	br.LogQ2  = LogQ2new;
	br.sud    = sudNew;

	Double z;
	if(iOrder <= 2 || cat % 2 == 0)
		z = CalculateZsplitting(br);
	else {
		//cout << "I am here " << endl;
		z = CalculateZsplittingDisc(br, cat);
	}

	//Generate new pT
	Double scale = exp(0.5*br.LogQ2);

	Vec2 pTs = RandVector(1.0);
	
	//Naive qt ordering, better qt ordering, angular ordering
	Double factors[] = { scale, sqrt(1-z)*scale, (1-z)*scale };
	for(int i = 0; i < 3; ++i) {
		br.pT[i]   =   br.pTmy[i] + factors[i]*pTs;
		br.pTmy[i] = z*br.pTmy[i] + factors[i]*pTs;
	}

	return true;
}



//returns z
Double TMD::CalculateZsplitting(Branching &Br) 
{

	const int Nf = alphaSpline::GetNf(Br.LogQ2);
	//const Double aS= alphaSpline::alphaS(Br.LogQ2)/(2.*M_PI); //precise
	const Double aS= sudakovNew.getAlphaS2pi(Nf, Br.LogQ2); //approx


	//Double PqIqJ, PqIqbJ, PqIqI, PqIqbI;
	//Double Pgg, Pqg, Pgq;
	//sudSplG.getSplittings(aS, Nf, Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI);

	//cout << "Old print " << Pgg<<" "<< Pqg<<" "<< Pgq<<" "<< PqIqJ<<" "<< PqIqbJ<<" "<< PqIqI<<" "<< PqIqbI << endl;

	valarray<Double> resMy = sudakovNew.getSplittingsI(br.LogQ2);
	const Double &Pgg =   resMy[0];
	const Double &Pqg  =  resMy[1]; //to all quarks
	const Double &Pgq  =  resMy[2];
	const Double &PqIqJ  =resMy[3]; // to all other quarks
	const Double &PqIqbJ =resMy[4]; //to all other anti-quarks
	const Double &PqIqI  =resMy[5];
	const Double &PqIqbI =resMy[6];
	//valarray<Double> resOld= {Pgg, Pqg, Pgq, PqIqJ, PqIqbJ, PqIqI, PqIqbI};

	/*
	for(int i = 0; i < 7; ++i) {
		Double diff = resMy[i] - resOld[i];
		if(abs(diff) > 1e-8)
			cout << diff << endl;
	}
	*/

	//cout << "New print " << Pgg<<" "<< Pqg<<" "<< Pgq<<" "<< PqIqJ<<" "<< PqIqbJ<<" "<< PqIqI<<" "<< PqIqbI << endl;

	/*
	valarray<Double> resOld= {Pgg, 2*Nf*Pqg, Pgq, (Nf-1)*PqIqJ, (Nf-1)*PqIqbJ, PqIqI, PqIqbI};
	for(int i = 0; i < 7; ++i) {
		if(resOld[i] < 1e-5) continue;
		Double diff = (resMy[i] - resOld[i])/resOld[i];
		if(diff > 1e-4)
			cout <<"RADEK "<< i << " "<< resMy[i] << " "<< resOld[i] << endl;
	}
	*/
	//function<Double(Double)> SplittingPgg, SplittingPqg, SplittingPgq;
	//function<Double(Double)> SplittingPqIqI,SplittingPqIqbI,SplittingPqIqJ,SplittingPqIqbJ;

	/*
	if(iOrder == 1) {
		SplittingPgg=[&](Double z) { return SplittingPgg0(z); };
		SplittingPqg=[&](Double z) { return 2*Nf*SplittingPqg0(z); };
		SplittingPgq=[&](Double z) { return SplittingPgq0(z); };
		SplittingPqIqI=[&](Double z) { return SplittingPqq0(z); };
		SplittingPqIqbI=[&](Double z) {return 0.; };
		SplittingPqIqJ=[&](Double z) {return 0.; };
		SplittingPqIqbJ=[&](Double z) {return 0.; };
	}
	else if(iOrder == 2) {
		SplittingPgg=[&](Double z) { return SplittingPgg0(z)+aS*SplittingPgg1(z,Nf); };
		SplittingPqg=[&](Double z) { return 2*Nf*(SplittingPqg0(z)+aS*SplittingPqg1(z,Nf)); };
		SplittingPgq=[&](Double z) { return SplittingPgq0(z)+aS*SplittingPgq1(z,Nf); };
		SplittingPqIqI=[&](Double z) { return SplittingPqq0(z)+aS*SplittingPqIqI1(z,Nf); };
		SplittingPqIqbI=[&](Double z) {return aS*SplittingPqIqbI1(z,Nf); };
		SplittingPqIqJ=[&](Double z)  {return (Nf-1)*aS*SplittingPqIqJ1(z,Nf); };
		SplittingPqIqbJ=[&](Double z) {return (Nf-1)*aS*SplittingPqIqJ1(z,Nf); };
	}
	else if(iOrder == 3) {
		SplittingPgg=[&](Double z) { return SplittingPgg0(z)+aS*SplittingPgg1(z,Nf)+
		                                                  aS*aS*SplittingPgg2(z,Nf); };
		SplittingPqg=[&](Double z) { return 2*Nf*(SplittingPqg0(z)+aS*SplittingPqg1(z,Nf)+
		                                                  aS*aS*SplittingPqg2(z,Nf) ); };
		SplittingPgq=[&](Double z) { return SplittingPgq0(z)+aS*SplittingPgq1(z,Nf)+
		                                                    aS*aS*SplittingPgq2(z,Nf); };
		SplittingPqIqI=[&](Double z) { return SplittingPqq0(z)+aS*SplittingPqIqI1(z,Nf)+
		                                                   aS*aS*SplittingPqIqI2(z,Nf); };
		SplittingPqIqbI=[&](Double z) {return aS*SplittingPqIqbI1(z,Nf)+
		                                      aS*aS*SplittingPqIqbI2(z,Nf); };
		SplittingPqIqJ=[&](Double z)  {return (Nf-1)*(aS*SplittingPqIqJ1(z,Nf)+
		                                       aS*aS*SplittingPqIqJ2(z,Nf)); };
		SplittingPqIqbJ=[&](Double z) {return (Nf-1)*(aS*SplittingPqIqJ1(z,Nf)+
		                                       aS*aS*SplittingPqIqbJ2(z,Nf)); };
	}
	*/


	//auto SplittingPgg=[&](Double z) {return sudakovNew.getSplitting(0,Nf,aS, atanh(2*z-1));};
	//auto SplittingPqg=[&](Double z) {return sudakovNew.getSplitting(1,Nf,aS, atanh(2*z-1));};
	//auto SplittingPgq=[&](Double z) {return sudakovNew.getSplitting(2,Nf,aS, atanh(2*z-1));};
	//auto SplittingPqIqI=[&](Double z) {return sudakovNew.getSplitting(5,Nf,aS, atanh(2*z-1));};
	//auto SplittingPqIqbI=[&](Double z) {return sudakovNew.getSplitting(6,Nf,aS, atanh(2*z-1));};
	//auto SplittingPqIqJ=[&](Double z) {return sudakovNew.getSplitting(3,Nf,aS, atanh(2*z-1));};
	//auto SplittingPqIqbJ=[&](Double z) {return sudakovNew.getSplitting(4,Nf,aS, atanh(2*z-1));};

	auto woWeight = [&](function<Double(Double)> P) {
		Double zmax =  fZmax(Br.LogQ2);
		Double Max = max(P(0), P(zmax) );
		double rx, ry, fun;
		do {
			rx = Rand()*zmax;
			ry = Rand()*Max;
			fun = P(rx);
			assert(fun <= Max);
		} while( ry > fun);
		return rx;
	};

	bool unitW = false;


	Double z;

	//cout << "isNLO " << isNLO << endl;
	//From gluon
	if(Br.id == 0) {
		Double sTot = Pgg + Pqg;
		Double ggProb =  Pgg / sTot;
		Double qgProb =  Pqg / sTot;
		if( ggProb > Rand() ) {
			Br.id = 0;
			Double wz = GetZ(Br.LogQ2, z);
			Double PggWgt=sudakovNew.getSplitting(0,Nf,aS, atanh(2*z-1)); //Pgg
			if(!unitW) Br.weight *= wz * PggWgt / (ggProb*sTot); //   Pgg;
			else z = woWeight(static_cast<Double(*)(Double)>(SplittingPgg0));
		}
		else {
			Br.id = RandI(2*Nf) + 1;
			if(Br.id > Nf) Br.id = Nf - Br.id;
			Double wz = GetZlin(Br.LogQ2, z);
			Double PqgWgt = sudakovNew.getSplitting(1,Nf,aS, atanh(2*z-1));
			if(!unitW) Br.weight *= wz * PqgWgt / (qgProb*sTot);        // Pqg
			else z = woWeight(static_cast<Double(*)(Double)>(SplittingPqg0));

			//Br.sud = sudSplQ.Eval( Br.LogQ2 ); //Change in sudakov
			Br.sud = sudakovNew.GetSudakovForQuark(Br.LogQ2);
			//cout << "Quarks " << Br.sud <<" "<< mySud << endl;
		}
	}
	//From quark
	else {

		//NLO part
		Double sTot = Pgq + PqIqI + PqIqbI + PqIqJ + PqIqbJ;
		Double gqProb    = Pgq / sTot;
		Double qIqIProb  = PqIqI  / sTot;
		Double qIqbIProb = PqIqbI / sTot;
		Double qIqJProb  = PqIqJ  / sTot; // To any quark with diffrent flavour
		Double qIqbJProb = PqIqbJ / sTot; // To any quark with diffrent flavour

		Double rr = Rand();
		if( rr  < gqProb ) { //to gluon
			Br.id = 0;
			Double wz = GetZlin(Br.LogQ2, z);
			Double PgqWgt = sudakovNew.getSplitting(2,Nf,aS, atanh(2*z-1)); //Pgq
			if(!unitW) Br.weight *= wz * PgqWgt / (gqProb * sTot); 
			else z = woWeight(static_cast<Double(*)(Double)>(SplittingPgq0));
			//Br.sud = sudSplG.Eval( Br.LogQ2 ); //Change in sudakov
			Br.sud = sudakovNew.GetSudakovForGluon(Br.LogQ2);
			//cout << "Quarks " << Br.sud <<" "<< mySud << endl;
		}
		else if( rr < qIqIProb+gqProb ) { //to same quark
			Br.id = Br.id;
			Double wz = GetZ(Br.LogQ2, z);
			Double PqIqIWgt=sudakovNew.getSplitting(5,Nf,aS, atanh(2*z-1)); 
			if(!unitW) Br.weight *= wz * PqIqIWgt / ( qIqIProb * sTot);  
			else z = woWeight(static_cast<Double(*)(Double)>(SplittingPqq0));
		}
		else if( rr < qIqbIProb+qIqIProb+gqProb ) {//to same anti-quark
			Br.id = -Br.id;
			Double wz = GetZ(Br.LogQ2, z);
			Double PqIqbIWgt=sudakovNew.getSplitting(6,Nf,aS, atanh(2*z-1)); 
			Br.weight *= wz * PqIqbIWgt / (qIqbIProb * sTot);
		}
		else if( rr < qIqJProb+qIqbIProb+qIqIProb+gqProb ) { //to different flavour qI->qJ
			int idTemp = RandI(Nf-1) + 1; //from 1 to Nf-1
			if(idTemp >= abs(Br.id))
				++idTemp;
			if(Br.id > 0)
				Br.id =+idTemp;
			else
				Br.id =-idTemp;

			Double wz = GetZ(Br.LogQ2, z);
			Double PqIqJWgt=sudakovNew.getSplitting(3,Nf,aS, atanh(2*z-1)); 
			Br.weight *= wz * PqIqJWgt / (qIqJProb * sTot);

		}
		else { //to different anti-flavour qI->qbJ
			int idTemp = RandI(Nf-1) + 1; //from 1 to Nf-1
			if(idTemp >= abs(Br.id))
				++idTemp;
			if(Br.id > 0)
				Br.id =-idTemp;
			else
				Br.id =+idTemp;

			Double wz = GetZ(Br.LogQ2, z);
			Double PqIqbJWgt=sudakovNew.getSplitting(4,Nf,aS, atanh(2*z-1)); 
			Br.weight *= wz * PqIqbJWgt  / (qIqbJProb * sTot);
		}
	}
	Br.x *= z;




	/*
	Double fromSpl[] = { SplittingPggT(th), SplittingPqgT(th), SplittingPgqT(th), SplittingPqIqIT(th),
	                  SplittingPqIqbIT(th), SplittingPqIqJT(th), SplittingPqIqbJT(th) };
	Double Org[] = { SplittingPgg(z), SplittingPqg(z), SplittingPgq(z), SplittingPqIqI(z),
	                  SplittingPqIqbI(z), SplittingPqIqJ(z), SplittingPqIqbJ(z) };
	Double rat[7];
	for(int i = 0; i < 7; ++i) {
		rat[i] = (fromSpl[i]-Org[i])/Org[i];
		cout <<"Testing "<< i<<" " << rat[i] << endl;
	}
	*/

	//cout <<"TEstik "<< SplittingPgg(zTest) <<" "<<  SplittingPggT(th)  << endl;





	return z;
}



Double TMD::CalculateZsplittingDisc(Branching &Br, int category) 
{
	Double z;
	Double nf = (category+1)/2;
	Double aS2pi = sudakovNew.GetMassTresh(nf).getAlphaS2pi();

	valarray<Double> splits    =  sudakovNew.GetMassTresh(nf).GetSplittings();
	valarray<Double> splitsAbs =  sudakovNew.GetMassTresh(nf).GetSplittingsAbs();

	Double splitHG=splits[0];
	Double splitGG=splits[1];
	Double splitQQ=splits[2];
	Double splitHQ=splits[3];
	Double splitGQ=splits[4];

	Double splitHGabs=splitsAbs[0];
	Double splitGGabs=splitsAbs[1];
	Double splitQQabs=splitsAbs[2];
	Double splitHQabs=splitsAbs[3];
	Double splitGQabs=splitsAbs[4];


	Double insideSudQ = splitQQ + splitHQ + splitGQ;
	Double insideSudG = splitGG + splitHG;



	//From quark
	if(Br.id != 0) {
		double r = Rand();
		//cout << "SplittingsQ    " << splitQQ <<" "<< splitGQ << " "<< splitHQ << endl;
		//cout << "SplittingsQabs " << splitQQabs <<" "<< splitGQabs << " "<< splitHQabs << endl;
		//
		//z = zmin + (zmax-zmin)*Rand();
		Br.weight *= GetZ(Br.LogQ2, z);

		double pTot = splitQQabs + splitGQabs + splitHQabs;
		double qqProb = splitQQabs / pTot;
		double gqProb = splitGQabs / pTot;
		double hqProb = splitHQabs / pTot;

		if( r < qqProb ) {
			Br.weight *= aS2pi*aS2pi*DiscQQ(z, 0) / (qqProb*insideSudQ);
			Br.id = Br.id;
		}
		else if( r < gqProb+qqProb ) {
			Br.weight *= aS2pi*aS2pi*DiscGQ(z, 0) / (gqProb*insideSudQ);
			Br.id = 0;
			//convert to gluon sudakov
			//sudakov = 1. - (1.-sudakov)/(1.-insideSudQ)*(1.-insideSudG);
			Br.sud = sudakovNew.ConvertToGluonSud(nf, Br.sud);
		}
		else {//Quark to heavy
			Br.id = nf*(2*RandI(2)-1);
			Br.weight *= aS2pi*aS2pi*DiscHQ(z, 0) / (hqProb*insideSudQ);
			Br.sud = sudakovNew.GetSudAboveTrQ(nf);
			//cout << "Q->H " << Br.id << endl;
		}
	}
	//From gluon
	else {
		double r = Rand();
		//cout << "SplittingsG " << splitGG <<" "<< splitHG  << endl;
		//
		//z = zmin + (zmax-zmin)*Rand();
		Br.weight *= GetZ(Br.LogQ2, z);

		double pTot = splitGGabs + splitHGabs;
		double ggProb = splitGGabs/pTot;
		double hgProb = splitHGabs/pTot;

		if( r < ggProb ) {
			Br.weight *= aS2pi*aS2pi*DiscGG(z, 0) / (ggProb*insideSudG);
			Br.id = 0;
		}
		else { //To Heavy
			Br.id = nf*(2*RandI(2)-1);
			Br.weight *= aS2pi*aS2pi*DiscHG(z, 0) / (hgProb*insideSudG);
			//convert to quark sudakov
			Br.sud = sudakovNew.GetSudAboveTrQ(nf);
			//cout << "G->H " << Br.id << endl;
		}
	}

	Br.x *= z;

	return z;
}


