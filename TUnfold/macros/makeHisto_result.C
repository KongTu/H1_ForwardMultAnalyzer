#include "RiceStyle.h"
#include <string>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

const double y_beam = 7.6;

void makeHisto_result(const int Q2_BIN = 0){

	gStyle->SetErrorX(0);
	
	TFile* file = new TFile("../rootfiles/results_unfoldingOutput_yetalab.root");
	TString histo_name;
	TString eta_bins[6]={"-1.2","0.2","-0.5","0.9","0.2","1.6"};
	TString Q2_bins[5]={"5","10","20","40","100"};
	TString y_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	TString x_average_Q2_5_10[]={"0.0014","0.0007","0.00033","0.00017"};
	TString x_average_Q2_10_20[]={"0.0023","0.0013","0.00066","0.00034"};
	TString x_average_Q2_20_40[]={"0.0052","0.0026","0.0013","0.00068"};
	TString x_average_Q2_40_100[]={"0.011","0.0056","0.0028","0.0013"};

	TH1D* hist_unfolded_from_django[3][4][4];
	TH1D* hist_unfolded_from_rapgap[3][4][4];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 4; k++){

				histo_name= "data_from_djangoh/hist_unfolded_("+eta_bins[2*i]+"<=eta<"+eta_bins[2*i+1]+
				")("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
				y_bins[k+1]+")";

				hist_unfolded_from_django[i][j][k] = (TH1D*) file->Get( histo_name );
				hist_unfolded_from_django[i][j][k]->SetMarkerStyle(24);
				hist_unfolded_from_django[i][j][k]->SetMarkerColor(kBlack);
				hist_unfolded_from_django[i][j][k]->SetLineColor(kBlack);

				histo_name = "data_from_rapgap/hist_unfolded_("+eta_bins[2*i]+"<=eta<"+eta_bins[2*i+1]+
				")("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
				y_bins[k+1]+")";

				hist_unfolded_from_rapgap[i][j][k] = (TH1D*) file->Get( histo_name );
				hist_unfolded_from_rapgap[i][j][k]->SetMarkerStyle(25);
				hist_unfolded_from_rapgap[i][j][k]->SetMarkerColor(kBlue);
				hist_unfolded_from_rapgap[i][j][k]->SetLineColor(kBlue);
			}
		}
	}

	TFile* file_MC_django = new TFile("../../minitree_output/Pn_hist_django.root");
	TFile* file_MC_rapgap = new TFile("../../minitree_output/Pn_hist_rapgap.root");
	TH1D* hist_mc_django[4][4][3];
	TH1D* hist_mc_rapgap[4][4][3];
	for(int i=0;i<3;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				
				hist_mc_django[j][k][i] = (TH1D*) file_MC_django->Get(Form("h_Pn_GEN_%d_%d_%d",j,k,i));
				hist_mc_rapgap[j][k][i] = (TH1D*) file_MC_rapgap->Get(Form("h_Pn_GEN_%d_%d_%d",j,k,i));
				
				hist_mc_django[j][k][i]->SetFillStyle(1001);
				hist_mc_django[j][k][i]->SetFillColorAlpha(kGreen-2,0.4);
				hist_mc_django[j][k][i]->SetMarkerColor(kGreen-2);

				hist_mc_rapgap[j][k][i]->SetFillStyle(1001);
				hist_mc_rapgap[j][k][i]->SetFillColorAlpha(kBlue-2,0.4);
				hist_mc_rapgap[j][k][i]->SetMarkerColor(kBlue-2);
			}
		}
	}

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	c1->Divide(4,3,0,0);

	TH1D* base1 = makeHist("base1", "", "N", "P(N)", 100,0,20,kBlack);
	base1->GetYaxis()->SetRangeUser(0.000005, 7);
	base1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base1,3,3.8);

	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.0);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.0);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.0);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.0);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);
	
	vector<double> S_hadron_x;

	TF1 *tf2 = new TF1("tf2","[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])",0,15);
	tf2->SetParameter(0,0.5);
	tf2->SetParameter(1,0.08);
	tf2->SetParameter(2,0.7);
	tf2->SetLineWidth(2);

	double latex_starting_eta[]={0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54};
	double latex_starting_y[]={0.48,0.49,0.525,0.56,0.48,0.49,0.525,0.56,0.46,0.47,0.51,0.55};

	int sub_panels = 1;
	for(int ieta = 0; ieta < 3; ieta++){
		for(int iy = 0; iy < 4; iy++){

			c1->cd(sub_panels);
			
			if(sub_panels==1||sub_panels==5||sub_panels==9) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=9) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(1);
			base1->Draw("");
			hist_unfolded_from_django[ieta][Q2_BIN][iy]->Draw("Psame");
			hist_unfolded_from_rapgap[ieta][Q2_BIN][iy]->Draw("Psame");

			hist_mc_django[Q2_BIN][iy][ieta]->DrawNormalized("E3 Psame");
			hist_mc_rapgap[Q2_BIN][iy][ieta]->DrawNormalized("E3 Psame");

			if( sub_panels < 5 ){
				TLatex* r46 = new TLatex(latex_starting_eta[sub_panels-1],0.76, eta_bins[2*ieta]+"<#eta<"+eta_bins[2*ieta+1]);
				r46->SetNDC();
				r46->SetTextSize(18);
				r46->SetTextFont(44);
				r46->Draw("same");
			}
			else{
				TLatex* r46 = new TLatex(latex_starting_eta[sub_panels-1],0.9, eta_bins[2*ieta]+"<#eta<"+eta_bins[2*ieta+1]);
				r46->SetNDC();
				r46->SetTextSize(18);
				r46->SetTextFont(44);
				r46->Draw("same");
			}

			if(sub_panels < 5){
				TLatex* r47 = new TLatex(latex_starting_y[sub_panels-1],0.67, y_bins[iy]+"<y<"+y_bins[iy+1]);
				r47->SetNDC();
				r47->SetTextSize(18);
				r47->SetTextFont(44);
				r47->Draw("same");
			}
			else{
				TLatex* r47 = new TLatex(latex_starting_y[sub_panels-1],0.82, y_bins[iy]+"<y<"+y_bins[iy+1]);
				r47->SetNDC();
				r47->SetTextSize(18);
				r47->SetTextFont(44);
				r47->Draw("same");
			}

			string eta1_string = (string) eta_bins[2*ieta];
			string eta2_string = (string) eta_bins[2*ieta+1];

			double eta_1 = std::stod( eta1_string );
			double eta_2 = std::stod( eta2_string );
			double x1_from_eta = TMath::Exp(eta_1-y_beam);
			double x2_from_eta = TMath::Exp(eta_2-y_beam);
			TString print_x = to_string(x1_from_eta)+"<x<"+to_string(x2_from_eta);

			TLatex* r48 = new TLatex(0.24,0.9, print_x);
			r48->SetNDC();
			r48->SetTextSize(18);
			r48->SetTextFont(44);
			// r48->Draw("same");


			TLatex* r49 = 0;
			if( sub_panels < 5 ){
				if(Q2_BIN==0) r49 = new TLatex(0.54,0.67, "#LTx_{b}#GT="+x_average_Q2_5_10[iy]);
				if(Q2_BIN==1) r49 = new TLatex(0.54,0.67, "#LTx_{b}#GT="+x_average_Q2_10_20[iy]);
				if(Q2_BIN==2) r49 = new TLatex(0.54,0.67, "#LTx_{b}#GT="+x_average_Q2_20_40[iy]);
				if(Q2_BIN==3) r49 = new TLatex(0.54,0.67, "#LTx_{b}#GT="+x_average_Q2_40_100[iy]);
			}
			else{
				if(Q2_BIN==0) r49 = new TLatex(0.54,0.82, "#LTx_{b}#GT="+x_average_Q2_5_10[iy]);
				if(Q2_BIN==1) r49 = new TLatex(0.54,0.82, "#LTx_{b}#GT="+x_average_Q2_10_20[iy]);
				if(Q2_BIN==2) r49 = new TLatex(0.54,0.82, "#LTx_{b}#GT="+x_average_Q2_20_40[iy]);
				if(Q2_BIN==3) r49 = new TLatex(0.54,0.82, "#LTx_{b}#GT="+x_average_Q2_40_100[iy]);
			}
			r49->SetNDC();
			r49->SetTextSize(18);
			r49->SetTextFont(44);
			// r49->Draw("same");

			if(Q2_BIN==0){
				if(sub_panels==3||sub_panels==6||sub_panels==9){
					TFitResultPtr r  = hist_unfolded_from_django[ieta][Q2_BIN][iy]->Fit("tf2","RMEsame");
				}
			}
			if(Q2_BIN==1){
				if(sub_panels==4||sub_panels==7||sub_panels==10){
					TFitResultPtr r  = hist_unfolded_from_django[ieta][Q2_BIN][iy]->Fit("tf2","RMEsame");
				}
			}
			if(Q2_BIN==2){
				if(sub_panels==8||sub_panels==11){
					TFitResultPtr r  = hist_unfolded_from_django[ieta][Q2_BIN][iy]->Fit("tf2","RMEsame");
				}
			}
			if(Q2_BIN==3){
				if(sub_panels==12){
					TFitResultPtr r  = hist_unfolded_from_django[ieta][Q2_BIN][iy]->Fit("tf2","RMEsame");
				}
			}
			
			double S_hadron = 0.;
			if( tf2->GetParameter(0)!=0.5 ){
				for(int n = 0; n < 15; n++){
					double pn = tf2->Eval(n);
					S_hadron += -pn*TMath::Log(pn);
				}		
				S_hadron_x.push_back(S_hadron);
				cout << "x ~ " << print_x << " and S_hadron ~ " << S_hadron << endl;
			}
			sub_panels++;


		}
		TLatex* r50 = new TLatex(0.63,0.87, "Q^{2}("+Q2_bins[Q2_BIN]+","+Q2_bins[Q2_BIN+1]+")");
		c1->cd(4);
		r50->SetNDC();
		r50->SetTextSize(18);
		r50->SetTextFont(44);
		r50->Draw("same");

		TLatex* r51 = new TLatex(0.22,0.87, "ep 27.5x920 GeV");
		c1->cd(1);
		r51->SetNDC();
		r51->SetTextSize(18);
		r51->SetTextFont(44);
		r51->Draw("same");

		TLatex* r52 = new TLatex(0.25,0.74, "H1");
		c1->cd(1);
		r52->SetNDC();
		r52->SetTextSize(20);
		r52->SetTextFont(44);
		r52->Draw("same");

		c1->cd(12);
		TLegend *w5 = new TLegend(0.05,0.22,0.5,0.5);
		w5->SetLineColor(kWhite);
		w5->SetFillColor(0);
		w5->SetTextSize(17);
		w5->SetTextFont(45);
		w5->AddEntry(tf2, "NBD fit ","L");
		w5->AddEntry(hist_mc_django[0][0][0], "django ","F");
		w5->AddEntry(hist_mc_rapgap[0][0][0], "rapgap ","F");
		w5->AddEntry(hist_unfolded_from_django[0][0][0], "unfolded data django", "P");
		w5->AddEntry(hist_unfolded_from_rapgap[0][0][0], "unfolded data rapgap ", "P");
		w5->Draw("same");
	}

	cout << "------ S of hadron -------" << endl;
	cout << "{" << S_hadron_x[0] << "," << S_hadron_x[1] << "," << S_hadron_x[2] << "}" << endl;

	c1->Print("../figures/Pn_Q2-"+Q2_bins[Q2_BIN]+"_"+Q2_bins[Q2_BIN+1]+".pdf");

}