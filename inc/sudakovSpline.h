#ifndef _sudakovSpline_
#define _sudakovSpline_

#include "Spline.h"
#include "alphaSpline.h"

#include "Splittings.h"

class sudakovSpline : public Spline<Double> {
public:
	sudakovSpline(int _N, int Flavour, int IOrder, Double _minX, Double _maxX);

	Double func(Double x) override;
	Double funcD(Double x) override;

	Double InsideSud(Double x) const;

	void printNodes() const override;

	void getSplittings(Double aS, int nf,
	             Double &Pgg, Double &Pqg, Double &Pgq,
	             Double &PqIqJ, Double &PqIqbJ, Double &PqIqI, Double &PqIqbI) const;

private:
	int flavour; //0 for gluon, - for anti-quarks
	int iOrder; //1-LO, 2-NLO, 2-NNLO
	Double LogMc2, LogMb2, LogMt2;
	Double Pgg0, Pqq0, Pgq0, Pqg0;
	Double Pgg1[7], Pqg1[7], Pgq1[7], PqIqJ1[7], PqIqI1[7], PqIqbI1[7];
	Double Pgg2[7], Pqg2[7], Pgq2[7], PqIqJ2[7], PqIqbJ2[7], PqIqI2[7], PqIqbI2[7];

};



#endif
