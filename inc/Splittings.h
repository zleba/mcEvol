#ifndef Splittings_H_
#define Splittings_H_

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


double DiscHQ(double z, bool isAbs=false);
double DiscQQ(double z, bool isAbs=false);
double DiscGQ(double z, bool isAbs=false);

double DiscGG(double z, bool isAbs=false);
double DiscHG(double z, bool isAbs=false);

#endif
