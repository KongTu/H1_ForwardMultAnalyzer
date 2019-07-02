#include "RiceStyle.h"

using namespace std;

void deriveRapidityWeight(){

	int total_Nevents = 11895778;
	TFile* file = new TFile("../minitree_output/Pn_hist_data.root");
	TH1D* h_eta[4];
	for(int i=0;i<4;i++){

		h_eta[i] = (TH1D*) file->Get(Form("h_eta_%d",i));
		h_eta[i]->SetMarkerStyle(24);
		h_eta[i]->SetMarkerColor(i+1);
		h_eta[i]->Scale( 100./6 * (1./h_eta[i]->GetEntries()) );
	}

	TCanvas* c1 = new TCanvas("c1","c1",1,1,500,500);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.06);
	// gPad->SetLogx(1);
	
	TH1D* base1 = makeHist("base1", "", "#eta", "Normalized dN/d#eta", 100,-3,3,kBlack);
	base1->GetYaxis()->SetRangeUser(0, 1);
	base1->GetXaxis()->SetRangeUser(-2, 2.5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.0);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.6);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.6);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base1->Draw();
	h_eta[0]->SetMarkerStyle(20);

	h_eta[0]->Draw("Psame");
	h_eta[1]->Draw("Psame");
	h_eta[2]->Draw("Psame");
	h_eta[3]->Draw("Psame");

	TLegend *w5 = new TLegend(0.27,0.69,0.65,0.85);
	w5->SetLineColor(kWhite);
	w5->SetFillColor(0);
	w5->SetTextSize(20);
	w5->SetTextFont(45);
	w5->AddEntry(h_eta[0], "total  ", "P");
	w5->AddEntry(h_eta[1], "0.075<y<0.15  ", "P");
	w5->AddEntry(h_eta[2], "0.15<y<0.3 ","P");
	w5->AddEntry(h_eta[3], "0.3<y<0.6 ","P");
	w5->Draw("same");

	double yield_low = h_eta[0]->Integral( h_eta[0]->FindBin(-1.2), h_eta[0]->FindBin(0.2) );
	double yield_mid = h_eta[0]->Integral( h_eta[0]->FindBin(-0.5), h_eta[0]->FindBin(0.9) );
	double yield_hig = h_eta[0]->Integral( h_eta[0]->FindBin(0.2), h_eta[0]->FindBin(1.6) );

	cout << "yiled low " << yield_low/1.4 << endl;
	cout << "yiled mid " << yield_mid/1.4 << endl;
	cout << "yiled hig " << yield_hig/1.4 << endl;

	c1->Print("eta_dist_h1.pdf");

}