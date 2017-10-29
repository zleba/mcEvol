
typedef double Double;
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "integration.h"


extern "C"  double gg1sfunc_(double *X, int *NF);
extern "C"  double pp1sfunc_(double *X, int *NF);
extern "C"  double pm1sfunc_(double *X, int *NF);
extern "C"  double gf1sfunc_(double *X, int *NF);
extern "C"  double fg1sfunc_(double *X, int *NF);
extern "C"  double ff1sfunc_(double *X, int *NF);




// g -> g
Double SplittingPgg1(Double z, int nf)
{
	//Multiplied by z
	return z * gg1sfunc_(&z,&nf);
}



Double SplittingPgg1(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPgg1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}


// g -> q
Double SplittingPqg1(Double z, int nf)
{
	//Multiplied by z, for one flavour
	//Note diffrent Furmalski-Petronzio notation
	return z * gf1sfunc_(&z, &nf) / (2*nf);
}

Double SplittingPqg1(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqg1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}



// q -> g
Double SplittingPgq1(Double z, int nf)
{
	//Note diffrent Furmalski-Petronzio notation
	return z * fg1sfunc_(&z, &nf);
}

Double SplittingPgq1(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPgq1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}


// qI -> qJ
Double SplittingPqIqJ1(Double z, int nf)
{
	//covers q -> q, q -> qb
	//covers qb -> qb, qb -> q
	//return z * ff1sfunc_(&z, &nf)/(2*nf);
	return z * ( ff1sfunc_(&z, &nf) -  pp1sfunc_(&z,&nf) )/(2*nf);
}

Double SplittingPqIqJ1(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqJ1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}


// qI -> qI
Double SplittingPqIqI1(Double z, int nf)
{
	//covers qI -> qI, qbI -> qbI
	//return z * ( ff1sfunc_(&z, &nf)/(2*nf) + ( pp1sfunc_(&z,&nf) + pm1sfunc_(&z,&nf) )/2. );
	return z * ( ff1sfunc_(&z, &nf)/(2*nf) +  pp1sfunc_(&z,&nf)*(1./2-1./(2*nf)) + pm1sfunc_(&z,&nf)* 1./2. );
}

Double SplittingPqIqI1(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqI1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}



// qI -> qbI
Double SplittingPqIqbI1(Double z, int nf)
{
	//covers qI -> qI, qbI -> qbI
	//return z * ( ff1sfunc_(&z, &nf)/(2*nf) + ( pp1sfunc_(&z,&nf) - pm1sfunc_(&z,&nf) )/2. );
	return z * ( ff1sfunc_(&z, &nf)/(2*nf) +   pp1sfunc_(&z,&nf)*(1./2-1./(2*nf)) - pm1sfunc_(&z,&nf)/2. );
}

Double SplittingPqIqbI1(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqbI1), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge error " << err << endl;
	return res;
}

