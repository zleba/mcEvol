#include <cmath>
using namespace std;

typedef double Double;

// g -> g
Double SplittingPgg0(Double z)
{
	//Multiplied by z
	return 6. * ( -2.*z + 1/(1-z)  +z*z*(1-z) );
}

Double SplittingPgg0(Double zmin, Double zmax)
{
	//Multiplied by z, integrated
	Double zmin2 = zmin*zmin;
	Double zmax2 = zmax*zmax;
	Double zmin3 = zmin2*zmin;
	Double zmax3 = zmax2*zmax;

	Double Dz2 = zmax2 - zmin2;
	Double Dz3 = zmax3 - zmin3;
	Double Dz4 = zmax*zmax3 - zmin*zmin3;

	return 6. * ( -Dz2 +  log((1-zmin)/(1-zmax))  + 1./3*Dz3 - 1./4*Dz4  );
}


// g -> q
Double SplittingPqg0(Double z)
{
	//Multiplied by z, for one flavour
	return z/2. * ( z*z + (1-z)*(1-z) );

}

Double SplittingPqg0(Double zmin, Double zmax)
{
	//Multiplied by z, for one flavour, integrated

	//Multiplied by z, integrated
	Double zmin2 = zmin*zmin;
	Double zmax2 = zmax*zmax;
	Double zmin3 = zmin2*zmin;
	Double zmax3 = zmax2*zmax;

	Double Dz2 = zmax2 - zmin2;
	Double Dz3 = zmax3 - zmin3;
	Double Dz4 = zmax*zmax3 - zmin*zmin3;

	return 1./2 * (2./4. * Dz4 - 2./3.*Dz3 + 1./2*Dz2 );
}

// q -> g
Double SplittingPgq0(Double z)
{
	//Multiplied by z
	return 4./3. * (1 + (1-z)*(1-z) );
}

Double SplittingPgq0(Double zmin, Double zmax)
{
	//Multiplied by z, integrated
	Double Dz3 = (1-zmax)*(1-zmax)*(1-zmax) - (1-zmin)*(1-zmin)*(1-zmin);
	return 4./3. * (zmax-zmin - 1./3.*Dz3 );
}


// q -> q
Double SplittingPqq0(Double z)
{
	//Multiplied by z
	return 4./3. * z*(1 + z*z) /(1-z);
}


Double SplittingPqq0(Double zmin, Double zmax)
{

	Double zmin2 = zmin*zmin;
	Double zmax2 = zmax*zmax;
	Double zmin3 = zmin2*zmin;
	Double zmax3 = zmax2*zmax;

	Double Dz2 = zmax2 - zmin2;
	Double Dz3 = zmax3 - zmin3;

	//Multiplied by z
	return 4./3. * ( -1./3.*Dz3 - 1./2*Dz2 -2*(zmax-zmin) + 2*log((1-zmin)/(1-zmax)) );

}
