#include "Sudakov.h"
#include <cassert>


void Sudakov::Init(int N, Double LnQ2min, Double LnQ2max, Double zmin,
				  function<Double(Double)> fZmax,  function<Double(Double)> fZmaxD)
{

	Double LnMc2 = 2.*alphaSpline::LnM[4];
	Double LnMb2 = 2.*alphaSpline::LnM[5];
	Double LnMt2 = 2.*alphaSpline::LnM[6];

	assert(LnMc2 < LnMb2);
	assert(LnMb2 < LnMt2);

	//Only the simplest case, TODO generalize
	int nNodes[7];
	nNodes[3] = (LnMc2 - LnQ2min) / (LnQ2max - LnQ2min) * N;
	nNodes[4] = (LnMb2 - LnMc2)   / (LnQ2max - LnQ2min) * N;
	nNodes[5] = (LnMt2 - LnMb2)   / (LnQ2max - LnQ2min) * N;
	nNodes[6] = (LnQ2max - LnMt2) / (LnQ2max - LnQ2min) * N;

	for(int i = 3; i <=6; ++i)
		assert(nNodes[i] > 4);


	sudItems[3] = SudakovItem(nNodes[3], LnQ2min, LnMc2,   zmin, fZmax, fZmaxD);
	sudItems[4] = SudakovItem(nNodes[4], LnMc2,   LnMb2,   zmin, fZmax, fZmaxD);
	sudItems[5] = SudakovItem(nNodes[5], LnMb2,   LnMt2,   zmin, fZmax, fZmaxD);
	sudItems[6] = SudakovItem(nNodes[6], LnMt2,   LnQ2max, zmin, fZmax, fZmaxD);

	massTresh[4] = MassTreshold(LnMc2, 4, zmin, fZmax(LnMc2) );
	massTresh[5] = MassTreshold(LnMb2, 5, zmin, fZmax(LnMb2) );
	massTresh[6] = MassTreshold(LnMt2, 6, zmin, fZmax(LnMt2) );


	Double sudakovsG[13]={0}, sudakovsQ[13]={0};
	for(int nf = 3; nf <= 6; ++nf) {
		sudakovsG[2*nf] = sudItems[nf].getSudMaxG();
		sudakovsQ[2*nf] = sudItems[nf].getSudMaxQ();
		if(nf >=4) {
			sudakovsG[2*nf-1] = massTresh[nf].getSudMaxG();
			sudakovsQ[2*nf-1] = massTresh[nf].getSudMaxQ();
		}
	}

	sudakovSumG[5] = 0;
	sudakovSumQ[5] = 0;
	for(int i = 6; i <=12; ++i ) {
		sudakovSumG[i] = sudakovsG[i] + sudakovSumG[i-1];
		sudakovSumQ[i] = sudakovsQ[i] + sudakovSumQ[i-1];
	}


}

Double Sudakov::GetLnQ2forGluon(Double sud, int &category) 
{
	//cout << "RADEK "<< sud <<" : "<< sudakovSumG[11]<<" "<< sudakovSumG[12] << endl;
	//cout <<"Sudakov sum: : "<< endl;
	//for(int i = 5; i <= 12; ++i)
		//cout << i << " "<< sudakovSumG[i] << endl;


	for(int i = 6; i <= 12 ; ++i)
		if(sud <= sudakovSumG[i] ) {
			category = i;
			if(i%2 == 1)
				return 2.*alphaSpline::LnM[(i+1)/2];
			else {
				Double localSud = sud - sudakovSumG[i-1];
				return sudItems[i/2].sudInvG.Eval(localSud, false);//cubic
			}
		}

	category = 12;
	return  sudItems[6].LnQ2max;
}

Double Sudakov::GetLnQ2forQuark(Double sud, int &category)
{
	for(int i = 6; i <= 12 ; ++i)
		if(sud <= sudakovSumQ[i] ) {
			category = i;
			if(i%2 == 1)
				return 2.*alphaSpline::LnM[(i+1)/2];
			else {
				Double localSud = sud - sudakovSumQ[i-1];
				return sudItems[i/2].sudInvQ.Eval(localSud, false);//cubic
			}
		}
	category = 12;
	return  sudItems[6].LnQ2max;
}

Double Sudakov::GetSudakovForGluon(Double LnQ2)
{
	int nf = alphaSpline::GetNf(LnQ2);
	Double localSud = sudItems[nf].sudSplG.Eval(LnQ2, false);//cubic
	return localSud + sudakovSumG[2*nf-1];
}

Double Sudakov::GetSudakovForQuark(Double LnQ2)
{
	int nf = alphaSpline::GetNf(LnQ2);
	Double localSud = sudItems[nf].sudSplQ.Eval(LnQ2, false);//cubic
	return localSud + sudakovSumQ[2*nf-1];
}


Double Sudakov::ConvertToGluonSud(int nf, Double sud)
{
	Double localSudQ = sud - sudakovSumQ[2*nf-2];
	Double localSudG = localSudQ /  GetMassTresh(nf).getSudMaxQ() *
	                                GetMassTresh(nf).getSudMaxG();

	return localSudG + sudakovSumG[2*nf-2];
}

valarray<Double> Sudakov::getSplittingsI(Double LnQ2)
{
	int nf = alphaSpline::GetNf(LnQ2);
	return sudItems[nf].splitSpline.Eval(LnQ2, false);
}
valarray<Double> getSplittingsDiscI()
{
	return valarray<Double>({0.,0.});
}
