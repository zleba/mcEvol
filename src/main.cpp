#include "tmd.h"
#include "sudakovSpline.h"
#include "SplittingsIntegralSpline.h"

#include "TH1D.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>


int RandI(int iMax);

class Histo {
	public:
		void Init();
        //Filling histogram with starting flavour flIn
        void Fill(const Branching &br,const Branching &brOld, int flIn);
		void Count();
	private:
		TH1D *hNev;

		TH1D *hLogXSc[8][13];
		TH1D *hLogXker[8][13][13];

		TH1D *hLogPtSc[8][3][13];
		TH1D *hLogPtNewSc[8][3][13];

		//For three x values
		TH1D *hLogPtScX001[8][3][13];
		TH1D *hLogPtNewScX001[8][3][13];
		TH1D *hLogPtScX0001[8][3][13];
		TH1D *hLogPtNewScX0001[8][3][13];
		TH1D *hLogPtScX00001[8][3][13];
		TH1D *hLogPtNewScX00001[8][3][13];

		vector<Double> Scales;
		vector<Double> LnScales;
};


void Histo::Init()
{
	#define SF TString::Format

	hNev = new TH1D("hNev", "hNev", 1, -0.5, 0.5);

	Scales = {2.00000001, 1e1, 1e3, 1e5, 1e7, 1e8 };
	LnScales.resize(Scales.size());
	for(unsigned i = 0; i < Scales.size(); ++i)
		LnScales[i] = log(Scales[i]);

	for(unsigned sc = 0; sc < Scales.size(); ++sc) {
        for(int i = 0; i < 13; ++i) 
        for(int j = 0; j < 13; ++j) 
            hLogXker[sc][i][j] = new TH1D( SF("hXKer%dSc_%d_%d",sc,i,j), SF("hXKer%dSc_%d_%d",sc,i,j), 20, -5, 0);
    }

	for(int i = 0; i < 13; ++i) 
	for(unsigned sc = 0; sc < Scales.size(); ++sc) {
		hLogXSc[sc][i]  = new TH1D( SF("hXSc%d_%d",sc,i), SF("hX%.0f_%d",Scales[sc],i), 20, -5, 0);
		for(int j = 0; j < 3; ++j) {
			hLogPtSc[sc][j][i] = new TH1D( SF("hLogPt%dSc%d_%d",j,sc,i), SF("hLogPt%dSc%d_%d",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtNewSc[sc][j][i] = new TH1D( SF("hLogPt%dNewSc%d_%d",j,sc,i), SF("hLogPt%dNewSc%d_%d",j,int(Scales[sc]),i), 120, -1, 5);

			hLogPtScX001[sc][j][i] = new TH1D( SF("hLogPt%dSc%d_%dX001",j,sc,i), SF("hLogPt%dSc%d_%dX001",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtNewScX001[sc][j][i] = new TH1D( SF("hLogPt%dNewSc%d_%dX001",j,sc,i), SF("hLogPt%dNewSc%d_%dX001",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtScX0001[sc][j][i] = new TH1D( SF("hLogPt%dSc%d_%dX0001",j,sc,i), SF("hLogPt%dSc%d_%dX0001",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtNewScX0001[sc][j][i] = new TH1D( SF("hLogPt%dNewSc%d_%dX0001",j,sc,i), SF("hLogPt%dNewSc%d_%dX0001",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtScX00001[sc][j][i] = new TH1D( SF("hLogPt%dSc%d_%dX00001",j,sc,i), SF("hLogPt%dSc%d_%dX00001",j,int(Scales[sc]),i), 120, -1, 5);
			hLogPtNewScX00001[sc][j][i] = new TH1D( SF("hLogPt%dNewSc%d_%dX00001",j,sc,i), SF("hLogPt%dNewSc%d_%dX00001",j,int(Scales[sc]),i), 120, -1, 5);
		}
	}

}


void Histo::Fill(const Branching &br,const Branching &brOld, int flIn)
{

	for(unsigned sc = 0; sc < LnScales.size(); ++sc) {
		if( brOld.LogQ2 < LnScales[sc] && br.LogQ2 > LnScales[sc] ) {
			Double logX = log10(brOld.x);
			hLogXSc[sc][brOld.id+6] ->Fill(logX, brOld.weight);
            hLogXker[sc][flIn+6][brOld.id+6]->Fill(logX, brOld.weight);
			for(int j = 0; j < 3; ++j) {
				Double pt2    = brOld.pT[j].norm2();
				Double pt2New = brOld.pTmy[j].norm2();
				hLogPtSc[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2), brOld.weight);
				hLogPtNewSc[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2New), brOld.weight);

				if(abs(logX + 2) <  0.25) {
					hLogPtScX001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2), brOld.weight);
					hLogPtNewScX001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2New), brOld.weight);
				}
				else if(abs(logX + 3) <  0.25) {
					hLogPtScX0001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2), brOld.weight);
					hLogPtNewScX0001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2New), brOld.weight);
				}
				else if(abs(logX + 4) <  0.25) {
					hLogPtScX00001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2), brOld.weight);
					hLogPtNewScX00001[sc][j][brOld.id+6]->Fill( 0.5*log10(pt2New), brOld.weight);
				}

			}
		}
	}


}

void Histo::Count()
{
	hNev->Fill(0);
}


int main()
{
	int iOrder = 1;

	alphaSpline::FixMasses( sqrt(2.8979862045712714), sqrt(24.904679585130495),	sqrt(28611.524722221096) );
	//alphaSpline::FixMasses(1e20* sqrt(2.8979862045712714),1e21* sqrt(24.904679585130495),1e22*sqrt(28611.524722221096) );

	alphaSpline::FixParameters(iOrder, 0.118, 5, sqrt(8315.25) );

	//DrawAlphaS();

	//return 0;

	//PlotSplittings();

	//return 0;
	/*
	splittingGeneral spl(30, log(2), log(2.8), 1e-6,
	             //[](Double qL){return 1. - 0.3*exp(-0.5*qL); },
					 //[](Double qL){return 0.3*0.5*exp(-0.5*qL); }  );
	             [](Double qL){return 1. - 1e-5; },
					 [](Double qL){return 0.0; }  );

	//spl.printNodes();

	Spline<Double> splG = spl.sumComponents(2,6).GetIntegratedSpline(false);

	splG.printNodes();
	*/
	/*
	cout << setprecision(10) ;
	cout << "Test point 0.845      : " << spl.Eval(0.845) <<" "<< spl.Eval(0.845, false) <<  endl;
	cout << "Test point 0.845 true : " << spl.funcBoth(0.845).first <<  endl;
	*/

	//return 0;
	//cout << SplittingPgg0(0.4) + 0.1 * SplittingPgg1(0.4, 5) << endl;
	//cout <<"P0 P1 "<<setprecision(16) <<   SplittingPqg0(0.4) << " "<< 0.1*SplittingPqIqJ1(0.4, 5) * 8  << endl;

	//return 0;

	/*
	cout <<"a100000 "<<setprecision(10)<<  alphaSpline::alphaS(log(1e5)) <<  endl;

	cout <<"a8315 "<<  alphaSpline::alphaS(log(8315.25)) <<  endl;
	cout <<"a100 "<<  alphaSpline::alphaS(log(100)) <<  endl;
	cout <<"a25Up  "<< alphaSpline::alphaS(log(24.904680)) << endl;
	cout <<"a25Dn  "<< alphaSpline::alphaS(log(24.904679)) << endl;
	cout <<"a20  "<< alphaSpline::alphaS(log(20)) << endl;
	cout <<"a2   "<< alphaSpline::alphaS(log(2)) << endl;
	cout <<"a2.69 "<< alphaSpline::alphaS(log(2.69)) << endl;
	return 0;
	*/


	/*
	myGame test;
	test.getLamba(sqrt(8315.25), 0.118);
	test.getAlpha( sqrt(8315.25) );
	test.getAlpha( sqrt(30) );
	return 0;
	*/

	//cout << alphaSpline::alphaS(log(9) ) /2./M_PI << endl;
	//cout << alphaSpline::alphaS(log(30) ) << endl;
	//cout << "Mz " << alphaSpline::alphaS(log(8315.25) ) << endl;
	//return 0;

	//TMD::fZmax  = [](Double Lq2) { return 1. - 0.3*exp(-Lq2/2); };
	//TMD::fZmaxD = [](Double Lq2) { return 0.3/2.*exp(-Lq2/2); };

	TMD::fZmax  = [](Double Lq2) { return 1-1e-6;};
	TMD::fZmaxD = [](Double Lq2) { return 0;};

	Double sCMS = 14000*14000;
	TMD tmd(3000, iOrder, 2.,  sCMS ); //Number of nodes, order, minQ2, maxQ2

	cout << "RADEK " << endl;

    int id = RandI(100000);
    cout << " My id is " << id << endl;

	TH1::SetDefaultSumw2();
	TFile *file;
	if(iOrder == 3)
		file= TFile::Open("histoNNLOh.root","recreate");
	else if(iOrder == 2)
		file= TFile::Open("histoNLOh.root","recreate");
	else
		file= TFile::Open(SF("histoLOkernel%d.root",id),"recreate");




	//DrawEvolution(tmd);

	//return 0;

	//DrawSudakov(&tmd);
	//return 0;


	//tmd.PrintSudakovs();
	//tmd.CompareSudakovs();

	//return 0;

	/*
	StartPDF sPDF;
	sPDF.InitDistribution();

	for(int i = 0; i < 1000; ++i) {
		Double x;
		int id;
		sPDF.genFirst(x, id);
	}
	return 0;
	*/
	Histo histos;
	histos.Init();

	for(int i = 0; i < 7*200000; ++i) {
        int fl = i % 7 - 3;
		tmd.Init(); //Random flavour and x accoring to PDF
		//tmd.Init(fl); //Flavour fl and x close to 1 (for KERNEL)
		histos.Count();
		//cout << "Event " << i << endl;
		do {
			//tmd.Print();
			tmd.Evolve();
			histos.Fill( tmd.GetBr(), tmd.GetBrOld(), fl);
		} while( tmd.GetBr().LogQ2 <= log(1e8) && tmd.GetBr().x > 1e-6  );
		if( (i >> 18)<<18  == i)
			cout << "New event "<< i  << endl;
	}

	file->Write();

	file->Close();

	return 0;

	sudakovSpline sudSpl(200, 0, false, log(2), log(100));

	sudSpl.printNodes();

	Spline<Double> sudInv = sudSpl.GetInverseSpline(400);

	cout << "Inverse nodes " << endl;
	sudInv.printNodes();

	return 0;
}
