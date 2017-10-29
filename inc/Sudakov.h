#ifndef Sudakov_H_
#define Sudakov_H_

#include "Spline.h"
#include "alphaSpline.h"
//#include "sudakovSplineSingle.h"
#include "SplittingsIntegralSpline.h"
#include "integration.h"

#include <map>
#include <iostream>

class SudakovItem {


public:

	SudakovItem(int N, Double LnQ2min_, Double LnQ2max_, Double zmin_,
	           function<Double(Double)> fZmax_,  function<Double(Double)> fZmaxD_) :
	      splitSpline( N, LnQ2min_, LnQ2max_, zmin_, fZmax_, fZmaxD_)
	{

			//cout << "RADEK splittings all " << endl;
			//splitSpline.printNodes();
			
			sudSplG = splitSpline.sumComponents(0, 1).GetIntegratedSpline(false);

			sudSplQ = splitSpline.sumComponents(2, 6).GetIntegratedSpline(false);

			//cout << "RADEK gluon sudakov " << endl;
			//sudSplG.printNodes();



			//cout << "RADEK quark sudakov " << endl;
			//sudSplQ.printNodes();

			//sudInvG = sudSplG.GetInverseSpline(12*N, false);
			//sudInvQ = sudSplQ.GetInverseSpline(12*N, false);

			sudInvG = sudSplG.InverseCubic(6*N);
			sudInvQ = sudSplQ.InverseCubic(6*N);
			sudSplG.CheckInverse(sudInvG, false);
			sudSplQ.CheckInverse(sudInvQ, false);

			LnQ2min = LnQ2min_;
			LnQ2max = LnQ2max_;
			sudMinG = 0.;
			sudMaxG = sudSplG.Eval(LnQ2max);
			sudMinQ = 0.;
			sudMaxQ = sudSplQ.Eval(LnQ2max);


			//Init Splitting spline
			int Nf = alphaSpline::GetNf( (LnQ2min+LnQ2max)/2. );
			splitting = splittingUnintegrated(40000, Nf, zmin_,  max(fZmax_(LnQ2max_),fZmax_(LnQ2min)) );
			//splitting.TestPrecision(1, 113);

			//Alpha spline
			AlphaSpline = alphaSpline(Nf, N, LnQ2min, LnQ2max);
			cout << "Splines for " << Nf <<" generated." << endl;
	}
	SudakovItem()  {}
	Double getSudMaxG() const {return sudMaxG;}
	Double getSudMaxQ() const {return sudMaxQ;}


public:
	splittingGeneral splitSpline;   //Integral of splittings spline (cubic)
	Spline<Double> sudSplG, sudSplQ; //Sudakov spline (cubic)
	Spline<Double> sudInvG, sudInvQ; //Inverse Sudakov spline (linear)
	splittingUnintegrated splitting; //Splitting func spline (linear)
	alphaSpline AlphaSpline; //Splitting func aS

	Double LnQ2min, LnQ2max;
	Double sudMinG, sudMaxG;
	Double sudMinQ, sudMaxQ;


};

class MassTreshold {
	Double splitQQ, splitHQ, splitGQ;
	Double splitHG, splitGG;
	Double splitQQabs, splitHQabs, splitGQabs;
	Double splitHGabs, splitGGabs;

	Double aS2pi;

	Double sudMinG, sudMaxG;
	Double sudMinQ, sudMaxQ;

public:
	MassTreshold(Double LnQ2, int nf, Double zmin, Double zmax) {Init(LnQ2, nf, zmin,zmax);}
	MassTreshold() {}

	Double getSudMaxG() const {return sudMaxG;}
	Double getSudMaxQ() const {return sudMaxQ;}
	Double getAlphaS2pi() const {return aS2pi;}

	valarray<Double> GetSplittings()
	         {return {splitHG, splitGG, splitQQ, splitHQ, splitGQ}; }
	valarray<Double> GetSplittingsAbs()
	         {return {splitHGabs, splitGGabs, splitQQabs, splitHQabs, splitGQabs}; }

	void Init(Double LnQ2, int nf, Double zmin, Double zmax) {
		sudMinQ=sudMinG = 0;
		if(alphaSpline::iOrder != 3) {
			sudMaxQ=sudMaxG = 0;
			splitQQ=splitHQ=splitGQ = 0;
			splitGG=splitHG = 0;
			aS2pi = 0.1/(2*M_PI);
			return;
		}

		aS2pi = alphaSpline::alphaS(LnQ2, nf)/2./M_PI;

		double err = 0;
		splitQQ = aS2pi*aS2pi*IntegrateSplitting( DiscQQ, 0, zmin, zmax,  err);
		splitHQ = aS2pi*aS2pi*IntegrateSplitting( DiscHQ, 0, zmin, zmax,  err);
		splitGQ = aS2pi*aS2pi*IntegrateSplitting( DiscGQ, 0, zmin, zmax,  err);

		splitGG = aS2pi*aS2pi*IntegrateSplitting( DiscGG, 0, zmin, zmax,  err);
		splitHG = aS2pi*aS2pi*IntegrateSplitting( DiscHG, 0, zmin, zmax,  err);

		sudMaxQ = splitQQ + splitHQ + splitGQ;
		sudMaxG = splitGG + splitHG;

		splitQQabs = aS2pi*aS2pi*IntegrateSplitting( DiscQQ, 1, zmin, zmax,  err);
		splitHQabs = aS2pi*aS2pi*IntegrateSplitting( DiscHQ, 1, zmin, zmax,  err);
		splitGQabs = aS2pi*aS2pi*IntegrateSplitting( DiscGQ, 1, zmin, zmax,  err);

		splitGGabs = aS2pi*aS2pi*IntegrateSplitting( DiscGG, 1, zmin, zmax,  err);
		splitHGabs = aS2pi*aS2pi*IntegrateSplitting( DiscHG, 1, zmin, zmax,  err);




	}

};


class Sudakov {

public:
	void Init(int N, Double LnQ2min_, Double LnQ2max_, Double zmin,
	           function<Double(Double)> fZmax_,  function<Double(Double)> fZmaxD_);



	Double GetLnQ2forGluon(Double sud, int &category);
	Double GetLnQ2forQuark(Double sud, int &category);
	Double GetSudakovForGluon(Double LnQ2);
	Double GetSudakovForQuark(Double LnQ2);

	MassTreshold &GetMassTresh(int nf) {return massTresh[nf];}

	Double ConvertToGluonSud(int nf, Double sud);

	Double GetSudAboveTrQ(int nf) { return nextafter(sudakovSumQ[2*nf-1], 1e20);}

	valarray<Double> getSplittingsI(Double LnQ2);
	valarray<Double> getSplittingsDiscI();


	Double getSplitting(int type, int nf, Double aS, Double z)
	                     { return sudItems[nf].splitting.EvalOne(type, aS, z); }
	Double getAlphaS2pi(int nf, Double LnQ2)
	                     { return sudItems[nf].AlphaSpline.Eval(LnQ2, false); }

private:
	SudakovItem sudItems[7];//3,4,5,6
	MassTreshold massTresh[7];

	Double sudakovSumG[13];
	Double sudakovSumQ[13];

};

#endif
