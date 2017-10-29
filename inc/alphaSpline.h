#ifndef _alphaSpline_
#define _alphaSpline_

#include "Spline.h"
#include <cmath>

class alphaSpline : public Spline<Double> {
public:
	alphaSpline(int nf, int _N, Double _minX, Double _maxX) :
	            Spline(_N, _minX, _maxX) {NF=nf; FillTogether(); }
	alphaSpline() {}

	static Double b0[7], b1[7], b2[7];

	static Double b[7];
	static Double c[7];
	static Double LnM[7];
	static Double LnLambdaSave[7];

	static int  iOrder;

	//Double func(Double x) override;
	//Double funcD(Double x) override;

	pair<Double,Double> funcBoth(Double LnQ2) override
	               {return alphaS2PiWithD(LnQ2, NF); }


	static void FixMasses(Double mc, Double mb, Double mt)
	            { LnM[4] = log(mc); LnM[5] = log(mb); LnM[6] = log(mt); }

	static void FixParameters(int Order, Double value, int iFix, Double mz = -1);
	static Double alphaS(Double LnQ2, int Nf=-1);
	static pair<Double,Double> alphaS2PiWithD(Double LnQ2, int nfUser);


	static Double GetLambda1(int nf, double a);
	static Double GetLambda2(int nf, double a);
	static Double GetLambda3(int nf, double a);


	static int GetNf(Double LnQ2) {
		Double LnQ = 0.5*LnQ2;
		if(LnQ < LnM[4])
			return 3;
		else if(LnQ < LnM[5])
			return 4;
		else if(LnQ < LnM[6])
			return 5;
		else
			return 6;
	}
	int NF;


};

#endif
