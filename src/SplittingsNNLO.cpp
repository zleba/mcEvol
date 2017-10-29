
typedef double Double;
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "integration.h"


extern "C"  double p2gga_(double *X, int *NF); // g -> g
extern "C"  double p2ggb_(double *X, int *NF); // g -> g
extern "C"  double p2ggc_(double *X, int *NF); // g -> g
extern "C"  double p2gqa_(double *X, int *NF); // q -> g
extern "C"  double p2qga_(double *X, int *NF); // g -> q

//q -> q
extern "C"  double p2psa_(double *X, int *NF); // (pure singlet)
extern "C"  double p2nspa_(double *X, int *NF); //P2_ns+
extern "C"  double p2nsma_(double *X, int *NF); //P2_ns-
extern "C"  double p2nsb_(double *X, int *NF); //P2_ns- and P2_ns+ (singular piece)
extern "C"  double p2nssa_(double *X, int *NF); //P2_nsv - P2_ns- (= Pns_S)


// qI -> qJ
Double SplittingPqIqJ2(Double z, int nf)
{
	//covers q -> q, 
	//covers qb -> qb
	return 1/8.*z * ( p2psa_(&z, &nf) + p2nssa_(&z,&nf) )/(2*nf) ;
}

Double SplittingPqIqJ2(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqJ2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}


// qI -> qbJ
Double SplittingPqIqbJ2(Double z, int nf)
{
	//covers q -> qb, 
	//covers qb -> q
	return 1/8.*z * ( p2psa_(&z, &nf) - p2nssa_(&z,&nf) )/(2*nf) ;
}

Double SplittingPqIqbJ2(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqbJ2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}




// qI -> qI
Double SplittingPqIqI2(Double z, int nf)
{
	//covers qI -> qI, qbI -> qbI
	return 1/8.*z * ( p2nsb_(&z,&nf) +  1./2*( p2nspa_(&z,&nf)+p2nsma_(&z,&nf) )   +    ( p2psa_(&z, &nf) + p2nssa_(&z,&nf) )/(2*nf)  );
}


Double SplittingPqIqI2(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqI2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}


// qI -> qbI
Double SplittingPqIqbI2(Double z, int nf)
{
	//covers qI -> qbI, qbI -> qI
	return 1/8.*z * ( 1./2*(p2nspa_(&z,&nf)-p2nsma_(&z,&nf))   +    ( p2psa_(&z, &nf) - p2nssa_(&z,&nf) )/(2*nf)   );
}


Double SplittingPqIqbI2(Double zmin, Double zmax, int nf)
{
	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqIqbI2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}





// g -> g
Double SplittingPgg2(Double z, int nf)
{
	//Multiplied by z
	return 1/8.*z * ( p2gga_(&z,&nf) + p2ggb_(&z,&nf)  );
}



Double SplittingPgg2(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPgg2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}


// g -> q
Double SplittingPqg2(Double z, int nf)
{
	//Multiplied by z, for one flavour
	//Note diffrent Furmalski-Petronzio notation
	return 1/8.*z * p2qga_(&z, &nf) / (2*nf);
}

Double SplittingPqg2(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPqg2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}


// q -> g
Double SplittingPgq2(Double z, int nf)
{
	//Note diffrent Furmalski-Petronzio notation
	return 1/8.*z * p2gqa_(&z, &nf);
}

Double SplittingPgq2(Double zmin, Double zmax, int nf)
{

	Double err = 0;
	Double res = IntegrateSplitting(static_cast<Double(*)(Double,int)>(SplittingPgq2), nf, zmin, zmax, err );
	if(err > 1e-8)
		cout << "Huge NNLO error " << err << endl;
	return res;
}


//Discontinuity functions


//Ahg (heavy quark discontinuity)
extern "C" double a2qg_(double *Z, double *FS2, double *HM2, double *OR);

//Ahq (heavy quark discontinuity)
extern "C" double a2qq_(double *Z, double *FS2, double *HM2, double *OR);
extern "C" double a2hga_(double *Y); //approximative function

//Aqq (light quark discontinuity
extern "C" double a2qqns_(double *Z, double *FS2, double *HM2,double *OR); //regular
extern "C" double softq2_(double *Z, double *FS2, double *HM2,double *OR); //singular
extern "C" double corq2_(double *Z, double *FS2, double *HM2, double *OR); //Delta piece

//Agg (to gluon discontinuity)
extern "C" double a2gg_(double *Z, double *FS2, double *HM2, double *OR); //regular
extern "C" double softg_(double *Z, double *FS2, double *HM2, double *OR); //singular
extern "C" double corg2_(double *Z, double *FS2, double *HM2, double *OR); //Delta piece


//Agq (to gluon discontinuity)
extern "C" double a2gq_(double *Z, double *FS2, double *HM2, double *OR); //regular


double DiscHQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;
	
	//double res = a2hga_(&z);  //approximative
	double res = a2qq_(&z, &dummy, &dummy, &order);//precise

	if(isAbs) res = abs(res);

	return res/4.*z;
}

double DiscQQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;

	double res = 0;
	res += a2qqns_(&z, &dummy, &dummy, &order); //regular
	res += softq2_(&z, &dummy, &dummy, &order); //singular

	if(isAbs) res = abs(res);

	return res/4.*z;
}

double DiscGQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;

	//Agq (to gluon discontinuity)
	double res = a2gq_(&z, &dummy, &dummy, &order); //regular

	if(isAbs) res = abs(res);

	return res/4.*z;
}

double DiscGG(double z, bool isAbs)
{
	double dummy=1;
	double order = 3;
	
	double res = 0;
	res += a2gg_(&z, &dummy, &dummy, &order); //regular
	res += softg_(&z, &dummy, &dummy, &order); //singular

	if(isAbs) res = abs(res);

	return res/4.*z;
}

double DiscHG(double z, bool isAbs)
{
	double dummy=1;
	double order = 3;
	
	//double res = a2hga_(&z);  //approximative
	double res = a2qg_(&z, &dummy, &dummy, &order);

	if(isAbs) res = abs(res);

	return res/4.*z;
}
