#include "tmd.h"
#include "sudakovSpline.h"
#include "SplittingsIntegralSpline.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TFrame.h"

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


void PlotSplittings()
{
	//TFile *fGraphs = new TFile("Graphs.root","recreate");

	TGraph *gr0[7], *gr1[7], *gr2[7];
	for(int i = 0; i < 7; ++i) {
		gr0[i] = new TGraph();
		gr1[i] = new TGraph();
		gr2[i] = new TGraph();
	}

	int nf = 5;
	double as= alphaSpline:: alphaS(log(50*50))/2./M_PI;
	int i = 0;
	for(double th = -6; th <=6; th+=0.01) {

		Double z = (1 + tanh(th) )/2.;
	 
	 	Double C = 1;//z*(1-z);
		valarray<Double> splits0 = C*as*splittingGeneral::getSplittings0(nf, z);
		valarray<Double> splits1 = C*as*as*splittingGeneral::getSplittings1(nf, z);
		valarray<Double> splits2 = C*as*as*as*splittingGeneral::getSplittings2(nf, z);
		

		splits0[0]*=(1-z);// (1+log(1./z));
		splits1[0]*=(1-z);// (1+log(1./z));
		splits2[0]*=(1-z);// (1+log(1./z));

		//splits0[1]*=1./pow(1+log(1./(1-z)),2) *1./ (1+log(1./z));
		//splits1[1]*=1./pow(1+log(1./(1-z)),2) *1./ (1+log(1./z));
		//splits2[1]*=1./pow(1+log(1./(1-z)),2) *1./ (1+log(1./z));

		splits0[2]*=(1-z);
		splits1[2]*=(1-z);
		splits2[2]*=(1-z);

		splits0[5]*=(1-z);
		splits1[5]*=(1-z);
		splits2[5]*=(1-z);

		//splits0[6]/=(1+log(1./z));
		//splits1[6]/=(1+log(1./z));
		//splits2[6]/=(1+log(1./z));

		//splits0[4]/=(1+log(1./z));
		//splits1[4]/=(1+log(1./z));
		//splits2[4]/=(1+log(1./z));

		//splits0[3]/=(1+log(1./z));
		//splits1[3]/=(1+log(1./z));
		//splits2[3]/=(1+log(1./z));



		for(int j = 0; j < 7; ++j) {
			gr2[j]->SetLineColor(kBlack);
			gr1[j]->SetLineColor(kRed);
			gr0[j]->SetLineColor(kGreen);

			gr0[j]->SetPoint(i, th, splits2[j]);
			gr1[j]->SetPoint(i, th, splits2[j]+splits1[j]);
			gr2[j]->SetPoint(i, th, splits0[j]+splits1[j]+splits2[j]);
		}
		++i;
	}
	TCanvas *can = new TCanvas("can","canvas");
	can->Divide(2,2);

	gr2[0]->SetMinimum(-0.2);
	gr2[2]->SetMinimum(-0.05);

	gr2[2]->SetMinimum(-0.05);

	gr2[0]->SetTitle("Pgg*(1-z)");
	gr2[1]->SetTitle("Pqg");
	gr2[2]->SetTitle("Pgq*(1-z)");
	gr2[3]->SetTitle("PqIqJ");
	gr2[4]->SetTitle("PqIqbJ");
	gr2[5]->SetTitle("PqIqI*(1-z)");
	gr2[6]->SetTitle("PqIqbI");


	for(i=0; i < 4; ++i) {
		can->cd(i+1);
		double min2 = gr2[i]->GetMinimum();
		double min1 = gr1[i]->GetMinimum();
		double min0 = gr0[i]->GetMinimum();

		double max2 = gr2[i]->GetMaximum();
		double max1 = gr1[i]->GetMaximum();
		double max0 = gr0[i]->GetMaximum();
		double Min = min({min0, min1, min2});
		double Max = max({max0, max1, max2});

		//gr2[i]->SetMinimum(-1.2*abs(Min));
		//gr2[i]->SetMaximum(+1.2*abs(Max));

		gr2[i]->Draw("ac");
		gr1[i]->Draw("c same");
		gr0[i]->Draw("c same");
	}
	can->SaveAs("splittings1.eps");
	TCanvas *dan = new TCanvas("dan","canvas");
	dan->Divide(2,2);
	for(i=0; i < 3; ++i) {
		dan->cd(i+1);
		gr2[4+i]->Draw("ac");
		gr1[4+i]->Draw("c same");
		gr0[4+i]->Draw("c same");
	}
	dan->SaveAs("splittings2.eps");

	//gr0[0]->Write("ahoj00");
	//gr2[2]->Write("ahoj22");

	//fGraphs->Write();
	//fGraphs->Close();

}


void DrawSudakov(TMD *tmd)
{

	TGraph *grS;
	grS = new TGraph();

	int fl = 1;

	int i = 0;
	for(double lnQ2 = log(2); lnQ2 <= 2*log(10); lnQ2+= 0.0001) {
		grS->SetPoint(i, exp(0.5*lnQ2), exp(-tmd->getSudakov(lnQ2,fl)) );
		++i;
	}
	TCanvas *sud = new TCanvas("sud","sudakov", 1000, 1000);
	sud->SetLogx();
	sud->SetLogy();

	gStyle->SetOptStat(0);
	TH1D *h = new TH1D("hSud", "hSud", 1, sqrt(2), 10);
	h->Draw("axis");
	h->GetXaxis()->SetTitle("#mu [GeV]");
	if(fl == 0)
		h->GetYaxis()->SetTitle("#Delta_{g}");
	else
		h->GetYaxis()->SetTitle("#Delta_{q}");
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleOffset(1.2);

	double minMain = (fl==0) ? 1e-5 : 5e-3;

	h->GetYaxis()->SetRangeUser(minMain, 1.2);
	h->GetXaxis()->SetMoreLogLabels();

	grS->DrawClone("l same");

	TLine *line = new TLine;
	line->SetLineStyle(2);
	auto plotLine = [&](double m, double x) {
		double y = grS->Eval(x);
		line->DrawLine(x, m, x, y);
	};

	plotLine(minMain, sqrt(2.8979862045712714) );
	plotLine(minMain, sqrt(24.904679585130495) );

	pair<double,double> MinMax;


	double posX=0.27, posY = 0.58;
	TPad *padC = new TPad("charmS", "charmS", posX, posY, posX+0.3, posY+0.3);
	
	padC->Draw();

	padC->cd();
	padC->SetLeftMargin(0.28);
	padC->SetBottomMargin(0.2);
	padC->SetTopMargin(0.15);
	TH1D *hMy = new TH1D("hCS", "Cmass", 1, 1.65, 1.75);
	hMy->SetTitle("charm");
	hMy->GetXaxis()->SetNdivisions(303);
	hMy->GetYaxis()->SetNdivisions(303);
	hMy->GetYaxis()->SetLabelSize(0.11);
	hMy->GetXaxis()->SetLabelSize(0.11);
	hMy->GetXaxis()->SetLabelOffset(0.03);
	hMy->GetYaxis()->SetLabelOffset(0.03);

	hMy->GetXaxis()->SetTickSize(0.10);
	hMy->GetYaxis()->SetTickSize(0.10);

	hMy->Draw("axis");

	MinMax = (fl==0) ? make_pair(0.1,0.3) : make_pair(0.39, 0.6);
	hMy->GetYaxis()->SetRangeUser(MinMax.first, MinMax.second);
	grS->DrawClone("l same");
	plotLine(MinMax.first, sqrt(2.8979862045712714));

	sud->cd();
	
	posX = 0.57; posY = 0.3;
	TPad *padB = new TPad("bottomS", "bottom", posX, posY, posX+0.3, posY+0.3);
	padB->Draw();

	padB->cd();
	padB->SetLeftMargin(0.28);
	padB->SetBottomMargin(0.2);
	padB->SetTopMargin(0.15);
	TH1D *hMyB = new TH1D("hBS", "Bmass", 1, 4.95, 5.05);
	hMyB->SetTitle("bottom");
	hMyB->GetYaxis()->SetDecimals();

	hMyB->GetXaxis()->SetNdivisions(303);
	hMyB->GetYaxis()->SetNdivisions(303);
	hMyB->GetYaxis()->SetLabelSize(0.11);
	hMyB->GetXaxis()->SetLabelSize(0.11);
	hMyB->GetXaxis()->SetLabelOffset(0.03);
	hMyB->GetYaxis()->SetLabelOffset(0.03);

	hMyB->GetXaxis()->SetTickSize(0.10);
	hMyB->GetYaxis()->SetTickSize(0.10);

	hMyB->Draw("axis");
	MinMax = (fl==0) ? make_pair(1e-4,2e-4) : make_pair(1.9e-2, 2.2e-2);
	hMyB->GetYaxis()->SetRangeUser(MinMax.first, MinMax.second);
	grS->Draw("l same");
	plotLine(MinMax.first, sqrt(24.904679585130495));

	sud->cd();

	grS->DrawClone("l same");



	sud->SaveAs("sudakov.eps");

}

void DrawAlphaS()
{
	
	TGraph *grS;
	grS = new TGraph();

	int i = 0;
	double step = 0.01;
	for(double lnQ2 = log(1); lnQ2 <= 2.*log(1000); lnQ2+= step) {
		if(exp(0.5*lnQ2) > 1.65 && exp(0.5*lnQ2) < 1.75)
			step = 0.0001;
		else if(exp(0.5*lnQ2) > 4.95 && exp(0.5*lnQ2) < 5.05) 
			step = 0.0001;
		else
			step = 0.01;

		grS->SetPoint(i, exp(0.5*lnQ2), alphaSpline::alphaS(lnQ2) );
		++i;
	}
	
	TCanvas *alpha = new TCanvas("alpha","alpha", 1000, 1000);
	alpha->SetLogx();
	grS->GetXaxis()->SetTitle("#mu [GeV]");
	grS->GetYaxis()->SetTitle("#alpha_{s}");
	grS->GetXaxis()->SetTitleOffset(1.2);
	grS->GetYaxis()->SetTitleOffset(1.2);
	grS->GetYaxis()->SetDecimals();
	grS->GetYaxis()->SetNdivisions(505);
	grS->DrawClone("ac");


	TLine *line = new TLine;
	line->SetLineStyle(2);
	
	auto plotLine = [&](double m, double x) {
		double y = grS->Eval(x);
		line->DrawLine(x, m, x, y);
	};

	plotLine(0.05, sqrt(2.8979862045712714));
	plotLine(0.05, sqrt(24.904679585130495));
	plotLine(0.05, sqrt(28611.524722221096));

	gStyle->SetOptStat(0);

	double posX=0.2, posY = 0.55;
	TPad *padC = new TPad("charm", "charm", posX, posY, posX+0.3, posY+0.3);
	padC->Draw();

	padC->cd();
	padC->SetLeftMargin(0.28);
	padC->SetBottomMargin(0.2);
	TH1D *hMy = new TH1D("hC", "Cmass", 1, 1.65, 1.75);
	hMy->SetTitle("charm");
	hMy->GetXaxis()->SetNdivisions(303);
	hMy->GetYaxis()->SetNdivisions(303);
	hMy->GetYaxis()->SetLabelSize(0.11);
	hMy->GetXaxis()->SetLabelSize(0.11);
	hMy->GetXaxis()->SetLabelOffset(0.03);
	hMy->GetYaxis()->SetLabelOffset(0.03);

	hMy->GetXaxis()->SetTickSize(0.10);
	hMy->GetYaxis()->SetTickSize(0.10);

	hMy->Draw("axis");
	hMy->GetYaxis()->SetRangeUser(0.32,0.34);
	grS->DrawClone("l same");
	plotLine(0.32, sqrt(2.8979862045712714));

	alpha->cd();
	
	posX = 0.55; posY = 0.3;
	TPad *padB = new TPad("bottom", "bottom", posX, posY, posX+0.3, posY+0.3);
	padB->Draw();

	padB->cd();
	padB->SetLeftMargin(0.28);
	padB->SetBottomMargin(0.2);
	TH1D *hMyB = new TH1D("hB", "Bmass", 1, 4.95, 5.05);
	hMyB->SetTitle("bottom");
	hMyB->GetYaxis()->SetDecimals();

	hMyB->GetXaxis()->SetNdivisions(303);
	hMyB->GetYaxis()->SetNdivisions(303);
	hMyB->GetYaxis()->SetLabelSize(0.11);
	hMyB->GetXaxis()->SetLabelSize(0.11);
	hMyB->GetXaxis()->SetLabelOffset(0.03);
	hMyB->GetYaxis()->SetLabelOffset(0.03);

	hMyB->GetXaxis()->SetTickSize(0.10);
	hMyB->GetYaxis()->SetTickSize(0.10);

	hMyB->Draw("axis");
	hMyB->GetYaxis()->SetRangeUser(0.210,0.220);
	grS->Draw("l same");
	plotLine(0.210, sqrt(24.904679585130495));











	alpha->SaveAs("alpha.eps");



}

void DrawEvolution(TMD &tmd)
{
	int nEv = 100;

	vector<TGraph *> gr(nEv);

	int col = 1;
	for(int i = 0; i < nEv; ++i) {
		gr[i] = new TGraph;
		gr[i]->SetMarkerColor(col);
		gr[i]->SetMarkerStyle(20);
		gr[i]->SetMarkerSize(0.8);
		gr[i]->SetLineColor(col);
		++col;
		if(col == 10)
			++col;

	}
	cout << "Hela" << endl;
	for(int i = 0; i < nEv; ++i) {
		tmd.Init(0); //gluon at 1
		//cout << "Event " << i << endl;
		int j = 0;
		do {
			//tmd.Print();
			//cout << "Ahoj " << tmd.GetBr().LogQ2 <<" "<<  tmd.GetBr().x << " "  << endl;
			gr[i]->SetPoint(j++, log10( tmd.GetBr().x), 0.5*tmd.GetBr().LogQ2/log(10) );
			tmd.Evolve();
		} while( tmd.GetBr().LogQ2 <= log(1e8) && tmd.GetBr().x > 1e-6  );
		gr[i]->SetPoint(j++, log10( tmd.GetBr().x), 0.5*tmd.GetBr().LogQ2/log(10) );
	}
	cout << "Hela2" << endl;

	TCanvas *can = new TCanvas("can", "can", 1000, 1000);
	can->SetLeftMargin(0.14);
	//cout << "Ahoj" << endl;
	//exit(0);

	TH2D *h = new TH2D("hTemp", "", 1, -4, 0, 1, 0, 4);
	gStyle->SetOptStat(0);
	h->Draw("axis");
	//h->GetXaxis()->SetTickLe
	h->GetXaxis()->SetTickSize(0.0);
	h->GetYaxis()->SetTickSize(0.0);
	h->GetXaxis()->SetLabelSize(0.0);
	h->GetYaxis()->SetLabelSize(0.0);

	TGaxis *axisX = new TGaxis(-4.0, 0.0, 0.0, 0.0,
                            1e-4, 1,510,"G");
	axisX->SetTitleFont(42);
	axisX->SetLabelFont(42);
	axisX->SetTitle("x");
	axisX->SetTitleOffset(1.2);

	axisX->Draw();

	TGaxis *axisY = new TGaxis(-4.0, 0.0, -4.0, 4.0,
                            1, 1e4,510,"G");
	axisY->SetTitleFont(42);
	axisY->SetLabelFont(42);
	axisY->SetTitle("#mu [GeV]");
	axisY->SetTitleOffset(1.4);

	axisY->Draw();



	gr[0]->Draw("lp same");
	for(int i = 0; i < nEv; ++i)
		gr[i]->Draw("lp same");


	vector<double> allScales;
	for(int i = 0; i < nEv; ++i)
		for(int j = 0; j < gr[i]->GetN(); ++j) {
			Double x, y;
			gr[i]->GetPoint(j, x, y);
			if(y > log10(sqrt(2)))
			allScales.push_back(y);
		}

	sort(allScales.begin(), allScales.end());

	TGraph *grMean = new TGraph;

	for(int n = 0; n < allScales.size(); ++n) {
		//for every pdf determine val
		double mean = 0;

		for(int i = 0; i < nEv; ++i) {
			Double x, y;
			Double xOld=0;
			int j;
			for(j= 0; j < gr[i]->GetN(); ++j) {
				gr[i]->GetPoint(j, x, y);
				if(y > allScales[n])
					break;
				xOld = x;
			}
			double val = xOld;
			mean += val;
		}
		cout << allScales[n] <<" "<< mean/nEv << endl;
		grMean->SetPoint(n, mean/nEv, allScales[n]);
	}

	grMean->SetLineStyle(1);
	grMean->SetLineWidth(4);
	grMean->Draw("l same");


	//TFrame *frU = (TFrame *) gPad->FindObject("TFrame");
	//frU->SetFillStyle(0);
	//frU->Draw();
	//gPad->Update();

	cout << TMD::fZmax(2) << endl;
	if( TMD::fZmax(2) == 1-1e-2)
		can->SaveAs("radecek2.eps");
	else if( TMD::fZmax(2) == 1-1e-6)
		can->SaveAs("radecek6.eps");
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

	for(int i = 0; i < 7*2000000; ++i) {
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
