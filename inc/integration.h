#ifndef _Integration_
#define _Integration_

#include <functional>

typedef double Double;

using namespace std;

Double IntegrateSplitting(function<Double(Double,int)> fun, int Nf, Double zmin, Double zmax,  Double &err);
Double Integral( function<Double(Double)> fun, Double xmin, Double xmax, Double &err );
Double Integral61( function<Double(Double)> fun, Double xmin, Double xmax, Double &err );

#endif
