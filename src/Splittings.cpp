
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

#include <cmath>
using std::abs;


double DiscHQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;
	
	//double res = a2hga_(&z);  //approximative
	double res = a2qq_(&z, &dummy, &dummy, &order);//precise
 
	if(isAbs) res = abs(res);

	return res/4.;
}

double DiscQQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;

	double res = 0;
	res += a2qqns_(&z, &dummy, &dummy, &order); //regular
	res += softq2_(&z, &dummy, &dummy, &order); //singular

	if(isAbs) res = abs(res);

	return res/4.;
}

double DiscGQ(double z, bool isAbs)
{
	double dummy = 1;
	double order = 3;

	//Agq (to gluon discontinuity)
	double res = a2gq_(&z, &dummy, &dummy, &order); //regular

	if(isAbs) res = abs(res);

	return res/4.;
}

double DiscGG(double z, bool isAbs)
{
	double dummy=1;
	double order = 3;
	
	double res = 0;
	res += a2gg_(&z, &dummy, &dummy, &order); //regular
	res += softg_(&z, &dummy, &dummy, &order); //singular

	if(isAbs) res = abs(res);

	return res/4.;
}

double DiscHG(double z, bool isAbs)
{
	double dummy=1;
	double order = 3;
	
	//double res = a2hga_(&z);  //approximative
	double res = a2qg_(&z, &dummy, &dummy, &order);

	if(isAbs) res = abs(res);

	return res/4.;
}
