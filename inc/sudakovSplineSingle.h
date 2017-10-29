#ifndef _sudakovSplineSingle_
#define _sudakovSplineSingle_

#include "Spline.h"
#include "alphaSpline.h"

class sudakovSplineSingle : public Spline<Double> {
public:
	sudakovSplineSingle(int _N, int Flavour, int IOrder, Double LnQ2min, Double LnQ2max);

	Double func(Double x) override;
	Double funcD(Double x) override;

	Double InsideSud(Double x) const;

	void printNodes() const override;

	void getGluonSplittings(Double aS2pi, Double *Pgg_, Double *Pqg_) const;
	void getQuarkSplittings(Double aS2pi, Double *Pgq_, Double *PqIqJ_,
	                        Double *PqIqbJ_, Double *PqIqI_, Double *PqIqbI_) const;

private:
	int flavour; //0 for gluon, - for anti-quarks
	int iOrder; //1-LO, 2-NLO, 2-NNLO
	int Nf;
	Double Pgg0, Pqq0, Pgq0, Pqg0;
	Double Pgg1, Pqg1, Pgq1, PqIqJ1, PqIqI1, PqIqbI1;
	Double Pgg2, Pqg2, Pgq2, PqIqJ2, PqIqbJ2, PqIqI2, PqIqbI2;

};

//LO splittings

Double SplittingPgg0(Double z);
Double SplittingPqg0(Double z);
Double SplittingPgq0(Double z);
Double SplittingPqq0(Double z);


Double SplittingPgg0(Double zmin, Double zmax);
Double SplittingPqg0(Double zmin, Double zmax);
Double SplittingPgq0(Double zmin, Double zmax);
Double SplittingPqq0(Double zmin, Double zmax);


//NLO splittings

Double SplittingPgg1(Double z, int nf);
Double SplittingPqg1(Double z, int nf);
Double SplittingPgq1(Double z, int nf);
Double SplittingPqIqJ1(Double z, int nf);
Double SplittingPqIqI1(Double z, int nf);
Double SplittingPqIqbI1(Double z, int nf);


Double SplittingPgg1(Double zmin, Double zmax, int nf);
Double SplittingPqg1(Double zmin, Double zmax, int nf);
Double SplittingPgq1(Double zmin, Double zmax, int nf);
Double SplittingPqIqJ1(Double zmin, Double zmax, int nf);
Double SplittingPqIqI1(Double zmin, Double zmax, int nf);
Double SplittingPqIqbI1(Double zmin, Double zmax, int nf);


//NNLO splittings

Double SplittingPgg2(Double z, int nf);
Double SplittingPqg2(Double z, int nf);
Double SplittingPgq2(Double z, int nf);
Double SplittingPqIqJ2(Double z, int nf);
Double SplittingPqIqbJ2(Double z, int nf);
Double SplittingPqIqI2(Double z, int nf);
Double SplittingPqIqbI2(Double z, int nf);


Double SplittingPgg2(Double zmin, Double zmax, int nf);
Double SplittingPqg2(Double zmin, Double zmax, int nf);
Double SplittingPgq2(Double zmin, Double zmax, int nf);
Double SplittingPqIqJ2(Double zmin, Double zmax, int nf);
Double SplittingPqIqbJ2(Double zmin, Double zmax, int nf);
Double SplittingPqIqI2(Double zmin, Double zmax, int nf);
Double SplittingPqIqbI2(Double zmin, Double zmax, int nf);




#endif
