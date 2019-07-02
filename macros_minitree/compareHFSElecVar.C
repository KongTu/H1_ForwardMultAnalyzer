#include "RiceStyle.h"

using namespace std;

void compareHFSElecVar(){

	TFile* file[3];
	file[0] = new TFile("../minitree_output/Pn_hist_data_check.root");
	file[1] = new TFile("../minitree_output/Pn_hist_django_check.root");
	file[2] = new TFile("../minitree_output/Pn_hist_rapgap_check.root");


	TH1D* h_PtBal[3];
	TH1D* h_hfsEnergy[3];
	TH1D* h_hfsPt[3];
	TH1D* h_hfsPz[3];
	TH1D* h_elecEnergy[3];
	TH1D* h_elecPt[3];
	TH1D* h_elecPz[3];
	TH1D* h_hfsEpz[3];
	TH1D* h_elecEpz[3];

	for(int ifile = 0; ifile < 3; ifile++){
		h_PtBal[ifile] = (TH1D*) file[ifile]->Get("h_PtBal");
		h_hfsEnergy[ifile] = (TH1D*) file[ifile]->Get("h_hfsEnergy");
		h_hfsPt[ifile] = (TH1D*) file[ifile]->Get("h_hfsPt");
		h_hfsPz[ifile] = (TH1D*) file[ifile]->Get("h_hfsPz");
		h_elecEnergy[ifile] = (TH1D*) file[ifile]->Get("h_elecEnergy");
		h_elecPt[ifile] = (TH1D*) file[ifile]->Get("h_elecPt");
		h_elecPz[ifile] = (TH1D*) file[ifile]->Get("h_elecPz");
		h_hfsEpz[ifile] = (TH1D*) file[ifile]->Get("h_hfsEpz");
		h_elecEpz[ifile] = (TH1D*) file[ifile]->Get("h_elecEpz");

		h_PtBal[ifile]->Rebin(10);
		h_hfsEnergy[ifile]->Rebin(10);
		h_hfsPt[ifile]->Rebin(10);
		h_hfsPz[ifile]->Rebin(10);
		h_elecEnergy[ifile]->Rebin(2);
		h_elecEnergy[ifile]->GetXaxis()->SetRangeUser(5,40);
		h_elecPt[ifile]->Rebin(2);
		h_elecPz[ifile]->Rebin(2);
		h_elecPz[ifile]->GetXaxis()->SetRangeUser(-40,5);
		h_hfsEpz[ifile]->Rebin(10);
		h_hfsEpz[ifile]->GetXaxis()->SetRangeUser(0,60);
		h_elecEpz[ifile]->Rebin(10);
		h_elecEpz[ifile]->GetXaxis()->SetRangeUser(10,70);

		if( ifile == 0 ) {
			h_PtBal[ifile]->SetMarkerStyle(20);
			h_hfsEnergy[ifile]->SetMarkerStyle(20);
			h_hfsPt[ifile]->SetMarkerStyle(20);
			h_hfsPz[ifile]->SetMarkerStyle(20);
			h_elecEnergy[ifile]->SetMarkerStyle(20);
			h_elecPt[ifile]->SetMarkerStyle(20);
			h_elecPz[ifile]->SetMarkerStyle(20);
			h_hfsEpz[ifile]->SetMarkerStyle(20);
			h_elecEpz[ifile]->SetMarkerStyle(20);
		}
		else if (ifile == 1) {
			h_PtBal[ifile]->SetMarkerStyle(24);
			h_hfsEnergy[ifile]->SetMarkerStyle(24);
			h_hfsPt[ifile]->SetMarkerStyle(24);
			h_hfsPz[ifile]->SetMarkerStyle(24);
			h_elecEnergy[ifile]->SetMarkerStyle(24);
			h_elecPt[ifile]->SetMarkerStyle(24);
			h_elecPz[ifile]->SetMarkerStyle(24);
			h_hfsEpz[ifile]->SetMarkerStyle(24);
			h_elecEpz[ifile]->SetMarkerStyle(24);
		}
		else {
			h_PtBal[ifile]->SetMarkerStyle(25);
			h_hfsEnergy[ifile]->SetMarkerStyle(25);
			h_hfsPt[ifile]->SetMarkerStyle(25);
			h_hfsPz[ifile]->SetMarkerStyle(25);
			h_elecEnergy[ifile]->SetMarkerStyle(25);
			h_elecPt[ifile]->SetMarkerStyle(25);
			h_elecPz[ifile]->SetMarkerStyle(25);
			h_hfsEpz[ifile]->SetMarkerStyle(25);
			h_elecEpz[ifile]->SetMarkerStyle(25);

			h_PtBal[ifile]->SetStats(kFALSE);
			h_hfsEnergy[ifile]->SetStats(kFALSE);
			h_hfsPt[ifile]->SetStats(kFALSE);
			h_hfsPz[ifile]->SetStats(kFALSE);
			h_elecEnergy[ifile]->SetStats(kFALSE);
			h_elecPt[ifile]->SetStats(kFALSE);
			h_elecPz[ifile]->SetStats(kFALSE);
			h_hfsEpz[ifile]->SetStats(kFALSE);
			h_elecEpz[ifile]->SetStats(kFALSE);
		}

	}


	TCanvas* c1 = new TCanvas("c1","c1",1,1,800,800);
	c1->Divide(3,3,0.001,0.001);

	c1->cd(1);
	gPad->SetLogy(1);
	h_hfsEnergy[2]->DrawNormalized("P");
	h_hfsEnergy[1]->DrawNormalized("Psame");
	h_hfsEnergy[0]->DrawNormalized("Psame");

	c1->cd(2);
	gPad->SetLogy(1);
	h_hfsPt[2]->DrawNormalized("P");
	h_hfsPt[1]->DrawNormalized("Psame");
	h_hfsPt[0]->DrawNormalized("Psame");

	c1->cd(3);
	gPad->SetLogy(1);
	h_hfsPz[2]->DrawNormalized("P");
	h_hfsPz[1]->DrawNormalized("Psame");
	h_hfsPz[0]->DrawNormalized("Psame");

	c1->cd(4);
	h_elecEnergy[2]->DrawNormalized("P");
	h_elecEnergy[1]->DrawNormalized("Psame");
	h_elecEnergy[0]->DrawNormalized("Psame");

	c1->cd(5);
	gPad->SetLogy(1);
	h_elecPt[2]->DrawNormalized("P");
	h_elecPt[1]->DrawNormalized("Psame");
	h_elecPt[0]->DrawNormalized("Psame");

	
	c1->cd(6);
	h_elecPz[2]->DrawNormalized("P");
	h_elecPz[1]->DrawNormalized("Psame");
	h_elecPz[0]->DrawNormalized("Psame");

	c1->cd(7);
	h_PtBal[2]->DrawNormalized("P");
	h_PtBal[1]->DrawNormalized("Psame");
	h_PtBal[0]->DrawNormalized("Psame");

	c1->cd(8);
	h_hfsEpz[2]->DrawNormalized("P");
	h_hfsEpz[1]->DrawNormalized("Psame");
	h_hfsEpz[0]->DrawNormalized("Psame");

	TLegend *w5 = new TLegend(0.35,0.47,0.85,0.75);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(17);
	w5->SetTextFont(45);
	w5->AddEntry(h_elecEpz[0], "data", "P");
	w5->AddEntry(h_elecEpz[1], "django ", "P");
	w5->AddEntry(h_elecEpz[2], "rapgap ", "P");
	w5->Draw("same");

	c1->cd(9);
	h_elecEpz[2]->DrawNormalized("P");
	h_elecEpz[1]->DrawNormalized("Psame");
	h_elecEpz[0]->DrawNormalized("Psame");


	











}