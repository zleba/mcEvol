#ifndef _tmd_
#define _tmd_

#include "Spline.h"
#include "Sudakov.h"
#include "sudakovSpline.h"
#include <cmath>
#include <functional>
#include "TH1D.h"
#include <iostream>

#include "Vec2.h"

struct StartPDF {

	Double Int[7], IntSum[7], sum;
	function< Double(Double) > pdf[7];

	Double AD, ADbar, AU;
	
	void InitDistribution();
	Double genFirst(bool DistributeEqualy, Double &x, int &id);
};

//struct Vec2 {
	//Double px, py;
	//Double norm2() {return px*px+py*py; }
//};

struct Branching {
		Double x, LogQ2;
		Vec2 pT[3], pTmy[3];
		Double sud;
		Double weight;
		int id;

        
};



class TMD {
	
	public:
		TMD(int N, int IOrder, Double q2min, Double q2max); 
		void Init(int fl=-50);
		bool Evolve();
		void Print() const { cout <<"Branching "<< br.id<<" "<<exp(br.LogQ2) <<" "<< br.sud <<", w="<<br.weight<< endl; }
		const Branching & GetBr() const {return br; }
		const Branching & GetBrOld() const {return brOld; }

		Double CalculateZsplitting(Branching &Br);
		Double CalculateZsplittingDisc(Branching &Br, int category);

		void PrintSudakovs() const;
		void CompareSudakovs();

		static Double zmin;//, zmax;
		static function<Double(Double)> fZmax, fZmaxD;

		Double getSudakov(Double lnQ2, int fl);

	private:

		Double GetZ(Double LnQ2, Double &z) const;
		Double GetZlin(Double LnQ2, Double &z) const;

		StartPDF sPDF;

		Branching br, brOld;
		//Double x, LogQ2;
		//Double px, py;
		//Double sud;
		//Double weight;
		//int id;
		int iOrder;

		Double LogQ2min, LogQ2max;
		Double sudMaxG, sudMaxQ;


		sudakovSpline sudSplG;
		sudakovSpline sudSplQ;

		Spline<Double> sudInvG;
		Spline<Double> sudInvQ;


		Sudakov sudakovNew;
		//sudakovNew.Init(1000, LogQ2min, LogQ2max, zmin,
							 //[](Double){return 1-1e-5;}, [](Double){return 0;} );

};





#endif
