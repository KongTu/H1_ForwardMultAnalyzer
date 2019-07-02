#include "RiceStyle.h"

using namespace std;

void plotEpz_zvertex(TString filename = "data"){

	TString histo_name;
	TString eta_bins[6]={"-1.2","0.2","-0.5","0.9","0.2","1.6"};
	TString Q2_bins[5]={"5","10","20","40","100"};
	TString y_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	TFile* file[10];
	file[0] = new TFile("../minitree_output/Pn_hist_"+filename+".root");

	//Q2,y,mult
	TH1D* h_Epz[4][4][2];
	TH1D* h_zvertex[4][4][2];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<2;k++){
				h_Epz[i][j][k] = (TH1D*) file[0]->Get(Form("h_Epz_%d_%d_%d",i,j,k));
				h_Epz[i][j][k]->Rebin(4);
				h_zvertex[i][j][k] = (TH1D*) file[0]->Get(Form("h_zvertex_%d_%d_%d",i,j,k));
			}
		}
	}


	for(int iQ2=0;iQ2<4;iQ2++){
		for(int iy=0;iy<4;iy++){

			TCanvas* c1 = new TCanvas(Form("c1_%d_%d",iQ2,iy),"c1",1,1,600,600);
			gPad->SetTicks();
			gPad->SetLeftMargin(0.15);
			gPad->SetBottomMargin(0.13);
			TH1D* base1 = makeHist("base1", "", "E-p_{z} (GeV)", "Normalized", 100,30,80,kBlack);
			base1->GetYaxis()->SetRangeUser(0.00005, 0.58);
			base1->GetXaxis()->SetTitleColor(kBlack);

			fixedFontHist1D(base1,1.2,1.4);

			base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.2);
			base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.2);
			base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.2);
			base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.2);
			base1->GetXaxis()->SetNdivisions(4,6,0);
			base1->GetYaxis()->SetNdivisions(4,6,0);

			h_Epz[iQ2][iy][0]->SetMarkerStyle(24);
			h_Epz[iQ2][iy][1]->SetMarkerStyle(25);

			base1->Draw();
			h_Epz[iQ2][iy][0]->DrawNormalized("Psame");
			h_Epz[iQ2][iy][1]->DrawNormalized("Psame");

			TLegend *w5 = new TLegend(0.65,0.67,0.85,0.75);
			w5->SetLineColor(kWhite);
			w5->SetFillColor(0);
			w5->SetTextSize(17);
			w5->SetTextFont(45);
			w5->AddEntry(h_Epz[iQ2][iy][0], "low mult", "P");
			w5->AddEntry(h_Epz[iQ2][iy][1], "high mult ", "P");
			w5->Draw("same");

			TLatex* r47 = new TLatex(0.65,0.85, y_bins[iy]+"<y<"+y_bins[iy+1]);
			r47->SetNDC();
			r47->SetTextSize(20);
			r47->SetTextFont(44);
			r47->Draw("same");

			TLatex* r50 = new TLatex(0.75,0.92, "Q^{2}("+Q2_bins[iQ2]+","+Q2_bins[iQ2+1]+")");
			c1->cd(4);
			r50->SetNDC();
			r50->SetTextSize(20);
			r50->SetTextFont(44);
			r50->Draw("same");

			TLatex* r51 = new TLatex(0.2,0.92, "ep 27.5x920 GeV");
			c1->cd(1);
			r51->SetNDC();
			r51->SetTextSize(20);
			r51->SetTextFont(44);
			r51->Draw("same");

			TLatex* r52 = new TLatex(0.2,0.85, "H1");
			c1->cd(1);
			r52->SetNDC();
			r52->SetTextSize(22);
			r52->SetTextFont(43);
			r52->Draw("same");

			// c1->Print(Form("../minitree_output/figures/Epz/"+filename+"_Epz_mult_%d_%d.pdf",iQ2,iy));

			TCanvas* c2 = new TCanvas(Form("c2_%d_%d",iQ2,iy),"c2",1,1,600,600);
			gPad->SetTicks();
			gPad->SetLeftMargin(0.15);
			gPad->SetBottomMargin(0.13);
			TH1D* base2 = makeHist("base2", "", "zvertex (cm)", "Normalized", 100,-40,40,kBlack);
			base2->GetYaxis()->SetRangeUser(0.000, 0.1);
			base2->GetXaxis()->SetTitleColor(kBlack);

			fixedFontHist1D(base2,1.2,1.4);

			base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.2);
			base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.2);
			base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.2);
			base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.2);
			base2->GetXaxis()->SetNdivisions(4,6,0);
			base2->GetYaxis()->SetNdivisions(4,6,0);

			h_zvertex[iQ2][iy][0]->SetMarkerStyle(24);
			h_zvertex[iQ2][iy][0]->Rebin(4);
			h_zvertex[iQ2][iy][1]->SetMarkerStyle(25);
			h_zvertex[iQ2][iy][1]->Rebin(4);

			base2->Draw();
			h_zvertex[iQ2][iy][0]->DrawNormalized("Psame");
			h_zvertex[iQ2][iy][1]->DrawNormalized("Psame");

			TLegend *w55 = new TLegend(0.65,0.67,0.85,0.75);
			w55->SetLineColor(kWhite);
			w55->SetFillColor(0);
			w55->SetTextSize(17);
			w55->SetTextFont(45);
			w55->AddEntry(h_zvertex[iQ2][iy][0], "low mult", "P");
			w55->AddEntry(h_zvertex[iQ2][iy][1], "high mult ", "P");
			w55->Draw("same");

			r47->Draw("same");
			r50->Draw("same");
			r51->Draw("same");
			r52->Draw("same");


			// c2->Print(Form("../minitree_output/figures/zvertex/"+filename+"_zvertex_mult_%d_%d.pdf",iQ2,iy));


		}
	}
	
}