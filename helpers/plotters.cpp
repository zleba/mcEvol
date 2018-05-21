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
