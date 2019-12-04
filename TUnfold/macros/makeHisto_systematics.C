#include "RiceStyle.h"
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

double calColError(double Ea, double Eb, double Sa, double Sb){

	double temp = Ea/Eb;
	double temp2 = (Sa*Sa)/(Ea*Ea) + (Sb*Sb)/(Eb*Eb);
	double temp3 = (2.*Sa*Sa)/(Ea*Eb);

	return temp*(sqrt(TMath::Abs(temp2-temp3)) );
}

//remember Eplus0p5 --- hadminus1p0 does not have rapgap! 
void makeHisto_systematics(const int Q2_BIN = 0, const bool doDjango_ = true, const bool doMCcomparison_ = false,const bool doRadCorr_ = false){

	gStyle->SetErrorX(0);

	TString filename_1 = "default_hadCali_hadEplus1p0";
	TString filename_2 = "default_hadCali_hadEminus1p0";
	
	TFile* file[10];
	file[0] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_default_hadCaliNew.root");
	file[1] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_"+filename_1+".root");
	file[2] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_"+filename_2+".root");

	TString filename_allEta_1 = "hadCali_Eplus0p5";
	TString filename_allEta_2 = "hadCali_Eminus0p5";

	double yaxis_range_1 = 0.5;
	double yaxis_range_2 = 2;

	TFile* file_allEta[10];
	file_allEta[0] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_default_singleEtaBin_hadCaliNew.root");
	file_allEta[1] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_default_singleEtaBin_"+filename_allEta_1+".root");
	file_allEta[2] = new TFile("../rootfiles/results_unfoldingOutput_yetalab_default_singleEtaBin_"+filename_allEta_2+".root");

	//radiative corrections
	TFile* file_RadCorr_rapgap = new TFile("./RadCorr_RAPGAP.root");
	TH1D* hCorr_Pn_GEN[4][4][4];
	TH1D* hCorr_Pn_GEN_HCM[4][4];
	for(int iQ2=0;iQ2<4;iQ2++){
		for(int iy=0;iy<4;iy++){
			hCorr_Pn_GEN_HCM[iQ2][iy] = (TH1D*) file_RadCorr_rapgap->Get(Form("hCorr_Pn_GEN_HCM_%d_%d",iQ2,iy));
			for(int ieta=0;ieta<4;ieta++){
				hCorr_Pn_GEN[iQ2][iy][ieta] = (TH1D*) file_RadCorr_rapgap->Get(Form("hCorr_Pn_GEN_%d_%d_%d",iQ2,iy,ieta));
			}
		}
	}
	TFile* file_RadCorr_django = new TFile("./RadCorr_DJANGO.root");
	TH1D* django_hCorr_Pn_GEN[4][4][4];
	TH1D* django_hCorr_Pn_GEN_HCM[4][4];
	for(int iQ2=0;iQ2<4;iQ2++){
		for(int iy=0;iy<4;iy++){
			django_hCorr_Pn_GEN_HCM[iQ2][iy] = (TH1D*) file_RadCorr_django->Get(Form("hCorr_Pn_GEN_HCM_%d_%d",iQ2,iy));
			for(int ieta=0;ieta<4;ieta++){
				django_hCorr_Pn_GEN[iQ2][iy][ieta] = (TH1D*) file_RadCorr_django->Get(Form("hCorr_Pn_GEN_%d_%d_%d",iQ2,iy,ieta));
			}
		}
	}

	TString histo_name;
	TString eta_bins[6]={"-1.2","0.2","-0.5","0.9","0.2","1.6"};
	TString Q2_bins[5]={"5","10","20","40","100"};
	TString y_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	TString x_average_Q2_5_10[]={"0.0014","0.0007","0.00033","0.00017"};
	TString x_average_Q2_10_20[]={"0.0023","0.0013","0.00066","0.00034"};
	TString x_average_Q2_20_40[]={"0.0052","0.0026","0.0013","0.00068"};
	TString x_average_Q2_40_100[]={"0.011","0.0056","0.0028","0.0013"};

	TH1D* hist_unfolded_from_django[3][3][4][4];//file, eta, Q2, y bins
	TH1D* hist_unfolded_from_rapgap[3][3][4][4];
	for(int ifile = 0; ifile < 3; ifile++){
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					histo_name= "data_from_djangoh/hist_unfolded_("+eta_bins[2*i]+"<=eta<"+eta_bins[2*i+1]+
					")("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
					y_bins[k+1]+")";

					hist_unfolded_from_django[ifile][i][j][k] = (TH1D*) file[ifile]->Get( histo_name );
					hist_unfolded_from_django[ifile][i][j][k]->SetMarkerStyle(24);
					hist_unfolded_from_django[ifile][i][j][k]->SetMarkerColor(kBlack);
					hist_unfolded_from_django[ifile][i][j][k]->SetLineColor(kBlack);

					histo_name = "data_from_rapgap/hist_unfolded_("+eta_bins[2*i]+"<=eta<"+eta_bins[2*i+1]+
					")("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
					y_bins[k+1]+")";

					hist_unfolded_from_rapgap[ifile][i][j][k] = (TH1D*) file[ifile]->Get( histo_name );
					hist_unfolded_from_rapgap[ifile][i][j][k]->SetMarkerStyle(25);
					hist_unfolded_from_rapgap[ifile][i][j][k]->SetMarkerColor(kBlue);
					hist_unfolded_from_rapgap[ifile][i][j][k]->SetLineColor(kBlue);
				}
			}
		}
	}

	TH1D* hist_unfolded_from_rapgap_allEta[3][4][4];
	TH1D* hist_unfolded_from_django_allEta[3][4][4];
	for(int ifile = 0; ifile < 3; ifile++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 4; k++){

				histo_name= "data_from_djangoh/hist_unfolded_(-1.2<=eta<1.6)("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
				y_bins[k+1]+")";

				hist_unfolded_from_django_allEta[ifile][j][k] = (TH1D*) file_allEta[ifile]->Get( histo_name );
				hist_unfolded_from_django_allEta[ifile][j][k]->SetMarkerStyle(24);
				hist_unfolded_from_django_allEta[ifile][j][k]->SetMarkerColor(kBlack);
				hist_unfolded_from_django_allEta[ifile][j][k]->SetLineColor(kBlack);

				histo_name = "data_from_rapgap/hist_unfolded_(-1.2<=eta<1.6)("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
				y_bins[k+1]+")";

				hist_unfolded_from_rapgap_allEta[ifile][j][k] = (TH1D*) file_allEta[ifile]->Get( histo_name );
				hist_unfolded_from_rapgap_allEta[ifile][j][k]->SetMarkerStyle(25);
				hist_unfolded_from_rapgap_allEta[ifile][j][k]->SetMarkerColor(kBlue);
				hist_unfolded_from_rapgap_allEta[ifile][j][k]->SetLineColor(kBlue);

				//radiative correction to Rapgap based on Rapgap
				for(int ibin=0;ibin<hist_unfolded_from_rapgap_allEta[ifile][j][k]->GetNbinsX();ibin++){
					double value = hist_unfolded_from_rapgap_allEta[ifile][j][k]->GetBinContent(ibin+1);
					double error = hist_unfolded_from_rapgap_allEta[ifile][j][k]->GetBinError(ibin+1);
					double bincenter = hist_unfolded_from_rapgap_allEta[ifile][j][k]->GetBinCenter(ibin+1);
					double corr = hCorr_Pn_GEN[j][k][3]->GetBinContent( hCorr_Pn_GEN[j][k][3]->FindBin(bincenter) );
					if(doRadCorr_){
						hist_unfolded_from_rapgap_allEta[ifile][j][k]->SetBinContent(ibin+1, value*corr);
						hist_unfolded_from_rapgap_allEta[ifile][j][k]->SetBinError(ibin+1, error*corr);
					}
				}
				//radiative correction to Django based on Rapgap
				for(int ibin=0;ibin<hist_unfolded_from_django_allEta[ifile][j][k]->GetNbinsX();ibin++){
					double value = hist_unfolded_from_django_allEta[ifile][j][k]->GetBinContent(ibin+1);
					double error = hist_unfolded_from_django_allEta[ifile][j][k]->GetBinError(ibin+1);
					double bincenter = hist_unfolded_from_django_allEta[ifile][j][k]->GetBinCenter(ibin+1);
					double corr = hCorr_Pn_GEN[j][k][3]->GetBinContent( hCorr_Pn_GEN[j][k][3]->FindBin(bincenter) );
					if(doRadCorr_){
						hist_unfolded_from_django_allEta[ifile][j][k]->SetBinContent(ibin+1, value*corr);
						hist_unfolded_from_django_allEta[ifile][j][k]->SetBinError(ibin+1, error*corr);
					}
				}
			}
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","c1",1,1,1000,800);
	c1->Divide(4,3,0,0);

	TH1D* base1 = makeHist("base1", "", "N", "P(N) ratio", 100,0,35,kBlack);
	base1->GetYaxis()->SetRangeUser(yaxis_range_1, yaxis_range_2);
	base1->GetXaxis()->SetRangeUser(0.0, 20);
	base1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base1,3,3.8);

	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.0);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.0);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.0);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.0);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	double latex_starting_eta[]={0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54};
	double latex_starting_y[]={0.48,0.49,0.525,0.56,0.48,0.49,0.525,0.56,0.46,0.47,0.51,0.55};
	
	TH1D* temp_tight_django[16];
	TH1D* temp_loose_django[16];
	TH1D* temp_tight_rapgap[16];
	TH1D* temp_loose_rapgap[16];

	int sub_panels = 1;
	for(int ieta = 0; ieta < 3; ieta++){
		for(int iy = 0; iy < 4; iy++){

			c1->cd(sub_panels);
			
			if(sub_panels==1||sub_panels==5||sub_panels==9) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=9) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(0);
			base1->Draw("");

			temp_tight_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->Clone(Form("temp_tight_django_%d",sub_panels));
			temp_loose_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->Clone(Form("temp_loose_django_%d",sub_panels));
			temp_tight_django[sub_panels-1]->SetMarkerStyle(24);
			temp_loose_django[sub_panels-1]->SetMarkerStyle(25);

			temp_tight_rapgap[sub_panels-1] = (TH1D*) hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->Clone(Form("temp_tight_rapgap_%d",sub_panels));
			temp_loose_rapgap[sub_panels-1] = (TH1D*) hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->Clone(Form("temp_loose_rapgap_%d",sub_panels));
			temp_tight_rapgap[sub_panels-1]->SetMarkerStyle(24);
			temp_loose_rapgap[sub_panels-1]->SetMarkerStyle(25);
				
			//django
			for(int ibin=0;ibin<temp_tight_django[sub_panels-1]->GetNbinsX();ibin++){
				
				double defau = hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double defau_err = hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( defau <= 0 ) continue;
				
				double tight = hist_unfolded_from_django[1][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double tight_err = hist_unfolded_from_django[1][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( tight <= 0 ) continue;
				temp_tight_django[sub_panels-1]->SetBinContent(ibin+1, tight/defau);
				temp_tight_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,tight,defau_err,tight_err));

				double loose = hist_unfolded_from_django[2][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double loose_err = hist_unfolded_from_django[2][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( loose <= 0 ) continue;

				temp_loose_django[sub_panels-1]->SetBinContent(ibin+1, loose/defau);
				temp_loose_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,loose,defau_err,loose_err));
			}

				
			//rapgap
			for(int ibin=0;ibin<temp_tight_rapgap[sub_panels-1]->GetNbinsX();ibin++){
				
				double defau = hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double defau_err = hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( defau == 0 ) continue;
				
				double tight = hist_unfolded_from_rapgap[1][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double tight_err = hist_unfolded_from_rapgap[1][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( tight <= 0 ) continue;

				temp_tight_rapgap[sub_panels-1]->SetBinContent(ibin+1, tight/defau);
				temp_tight_rapgap[sub_panels-1]->SetBinError(ibin+1, calColError(defau,tight,defau_err,tight_err));

				double loose = hist_unfolded_from_rapgap[2][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
				double loose_err = hist_unfolded_from_rapgap[2][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
				if( loose <= 0 ) continue;

				temp_loose_rapgap[sub_panels-1]->SetBinContent(ibin+1, loose/defau);
				temp_loose_rapgap[sub_panels-1]->SetBinError(ibin+1, calColError(defau,loose,defau_err,loose_err));
			}
			if( doDjango_ ){
				temp_tight_django[sub_panels-1]->Draw("Psame");
				temp_loose_django[sub_panels-1]->Draw("Psame");

				if(sub_panels==1){
					temp_tight_django[sub_panels-1]->Fit("pol0","","",0,6);
					temp_loose_django[sub_panels-1]->Fit("pol0","","",0,6);
				}
				else{
					temp_tight_django[sub_panels-1]->Fit("pol0");
					temp_loose_django[sub_panels-1]->Fit("pol0");
				}
				
				TF1* fitFunction_tight = (TF1*) temp_tight_django[sub_panels-1]->GetFunction("pol0");
				TF1* fitFunction_loose = (TF1*) temp_loose_django[sub_panels-1]->GetFunction("pol0");
				double latex_height = 0.57;
				if( sub_panels > 4 ) latex_height = 0.72;

				double tight_fit_value = fitFunction_tight->GetParameter(0);
				double loose_fit_value = fitFunction_loose->GetParameter(0);
				double fit_value = loose_fit_value;
				if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
					fit_value = tight_fit_value;
				}

				TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.3f",fit_value) );
				rFit->SetNDC();
				rFit->SetTextSize(18);
				rFit->SetTextFont(44);
				rFit->Draw("same");


			}
			else{
				temp_tight_rapgap[sub_panels-1]->Draw("Psame");
				temp_loose_rapgap[sub_panels-1]->Draw("Psame");

				if(sub_panels==1){
					temp_tight_rapgap[sub_panels-1]->Fit("pol0","","",0,6);
					temp_loose_rapgap[sub_panels-1]->Fit("pol0","","",0,6);
				}
				else{
					temp_tight_rapgap[sub_panels-1]->Fit("pol0");
					temp_loose_rapgap[sub_panels-1]->Fit("pol0");
				}

				TF1* fitFunction_tight = (TF1*) temp_tight_rapgap[sub_panels-1]->GetFunction("pol0");
				TF1* fitFunction_loose = (TF1*) temp_loose_rapgap[sub_panels-1]->GetFunction("pol0");
				double latex_height = 0.57;
				if( sub_panels > 4 ) latex_height = 0.72;

				double tight_fit_value = fitFunction_tight->GetParameter(0);
				double loose_fit_value = fitFunction_loose->GetParameter(0);
				double fit_value = loose_fit_value;
				if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
					fit_value = tight_fit_value;
				}

				TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.3f",fit_value) );
				rFit->SetNDC();
				rFit->SetTextSize(18);
				rFit->SetTextFont(44);
				rFit->Draw("same");

			}

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
			sub_panels++;
			}
			
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
		TLegend *w5 = new TLegend(0.05,0.26,0.5,0.45);
		w5->SetLineColor(kWhite);
		w5->SetFillColor(0);
		w5->SetTextSize(17);
		w5->SetTextFont(45);
		if( doDjango_ ){
			w5->AddEntry(temp_tight_django[0], "django "+filename_1+"/default ","P");
			w5->AddEntry(temp_loose_django[0], "django "+filename_2+"/default ","P");
		}
		else{
			w5->AddEntry(temp_tight_rapgap[0], "rapgap "+filename_1+"/default ","P");
			w5->AddEntry(temp_loose_rapgap[0], "rapgap "+filename_2+"/default ","P");
		}
		
		w5->Draw("same");


		TCanvas* c3 = new TCanvas("c3","c3",1,1,1000,1000);
		c3->Divide(4,4,0,0);
		TH1D* base3 = (TH1D*) base1->Clone("base3");
		base3->GetYaxis()->SetRangeUser(yaxis_range_1, yaxis_range_2);
		base3->GetXaxis()->SetRangeUser(0,30);

		double latex_starting_eta_1[]={0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54};
		double latex_starting_y_1[]={0.48,0.49,0.48,0.48,0.49,0.53,0.53,0.53,0.51,0.525,0.56,0.53,0.52,0.51,0.5,0.55};

		TH1D* temp_allEta_tight_django[16];
		TH1D* temp_allEta_loose_django[16];
		TH1D* temp_allEta_tight_rapgap[16];
		TH1D* temp_allEta_loose_rapgap[16];

		sub_panels = 1;
		for(int iy = 0; iy < 4; iy++){
		for(int iQ2 = 0; iQ2 < 4; iQ2++){

			c3->cd(sub_panels);
			
			if(sub_panels==1||sub_panels==5||sub_panels==9||sub_panels==13) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=13) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(0);
			base3->Draw("");
			
			temp_allEta_tight_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django_allEta[1][iQ2][iy]->Clone(Form("temp_allEta_tight_django_%d",sub_panels));
			temp_allEta_loose_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django_allEta[2][iQ2][iy]->Clone(Form("temp_allEta_loose_django_%d",sub_panels));
			temp_allEta_tight_django[sub_panels-1]->SetMarkerStyle(24);
			temp_allEta_loose_django[sub_panels-1]->SetMarkerStyle(25);
			temp_allEta_tight_rapgap[sub_panels-1] = (TH1D*) hist_unfolded_from_rapgap_allEta[1][iQ2][iy]->Clone(Form("temp_allEta_tight_rapgap_%d",sub_panels));
			temp_allEta_loose_rapgap[sub_panels-1] = (TH1D*) hist_unfolded_from_rapgap_allEta[2][iQ2][iy]->Clone(Form("temp_allEta_loose_rapgap_%d",sub_panels));
			temp_allEta_tight_rapgap[sub_panels-1]->SetMarkerStyle(24);
			temp_allEta_loose_rapgap[sub_panels-1]->SetMarkerStyle(25);

		
			//django
			for(int ibin=0;ibin<temp_allEta_tight_django[sub_panels-1]->GetNbinsX();ibin++){
				
				double defau = hist_unfolded_from_django_allEta[0][iQ2][iy]->GetBinContent(ibin+1);
				double defau_err = hist_unfolded_from_django_allEta[0][iQ2][iy]->GetBinError(ibin+1);
				if( defau <= 0 ) continue;
				
				double tight = hist_unfolded_from_django_allEta[1][iQ2][iy]->GetBinContent(ibin+1);
				double tight_err = hist_unfolded_from_django_allEta[1][iQ2][iy]->GetBinError(ibin+1);
				if( tight <= 0 ) continue;
				temp_allEta_tight_django[sub_panels-1]->SetBinContent(ibin+1, tight/defau);
				temp_allEta_tight_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,tight,defau_err,tight_err));

				double loose = hist_unfolded_from_django_allEta[2][iQ2][iy]->GetBinContent(ibin+1);
				double loose_err = hist_unfolded_from_django_allEta[2][iQ2][iy]->GetBinError(ibin+1);
				if( loose <= 0 ) continue;

				temp_allEta_loose_django[sub_panels-1]->SetBinContent(ibin+1, loose/defau);
				temp_allEta_loose_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,loose,defau_err,loose_err));
			}

				
			//rapgap
			for(int ibin=0;ibin<temp_allEta_tight_rapgap[sub_panels-1]->GetNbinsX();ibin++){
				
				double defau = hist_unfolded_from_rapgap_allEta[0][iQ2][iy]->GetBinContent(ibin+1);
				double defau_err = hist_unfolded_from_rapgap_allEta[0][iQ2][iy]->GetBinError(ibin+1);
				if( defau == 0 ) continue;
				
				double tight = hist_unfolded_from_rapgap_allEta[1][iQ2][iy]->GetBinContent(ibin+1);
				double tight_err = hist_unfolded_from_rapgap_allEta[1][iQ2][iy]->GetBinError(ibin+1);
				if( tight <= 0 ) continue;

				temp_allEta_tight_rapgap[sub_panels-1]->SetBinContent(ibin+1, tight/defau);
				temp_allEta_tight_rapgap[sub_panels-1]->SetBinError(ibin+1, calColError(defau,tight,defau_err,tight_err));

				double loose = hist_unfolded_from_rapgap_allEta[2][iQ2][iy]->GetBinContent(ibin+1);
				double loose_err = hist_unfolded_from_rapgap_allEta[2][iQ2][iy]->GetBinError(ibin+1);
				if( loose <= 0 ) continue;

				temp_allEta_loose_rapgap[sub_panels-1]->SetBinContent(ibin+1, loose/defau);
				temp_allEta_loose_rapgap[sub_panels-1]->SetBinError(ibin+1, calColError(defau,loose,defau_err,loose_err));
			}
			if( doDjango_ ){
				temp_allEta_tight_django[sub_panels-1]->Draw("Psame");
				temp_allEta_loose_django[sub_panels-1]->Draw("Psame");
				if( sub_panels==8 ){
					temp_allEta_tight_django[sub_panels-1]->Fit("pol0","","",0,20.);
					temp_allEta_loose_django[sub_panels-1]->Fit("pol0","","",0,20.);
				}
				else if( sub_panels==16 ){
					temp_allEta_tight_django[sub_panels-1]->Fit("pol0","","",1,20.);
					temp_allEta_loose_django[sub_panels-1]->Fit("pol0","","",1,20.);
				}
				else{
					temp_allEta_tight_django[sub_panels-1]->Fit("pol0");
					temp_allEta_loose_django[sub_panels-1]->Fit("pol0");
				}

				TF1* fitFunction_tight = (TF1*) temp_allEta_tight_django[sub_panels-1]->GetFunction("pol0");
				TF1* fitFunction_loose = (TF1*) temp_allEta_loose_django[sub_panels-1]->GetFunction("pol0");
				double latex_height = 0.57;
				if( sub_panels > 4 ) latex_height = 0.72;

				double tight_fit_value = fitFunction_tight->GetParameter(0);
				double loose_fit_value = fitFunction_loose->GetParameter(0);
				double fit_value = loose_fit_value;
				if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
					fit_value = tight_fit_value;
				}

				TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.4f",fit_value) );
				rFit->SetNDC();
				rFit->SetTextSize(18);
				rFit->SetTextFont(44);
				rFit->Draw("same");
				
			}
			else{
				temp_allEta_tight_rapgap[sub_panels-1]->Draw("Psame");
				temp_allEta_loose_rapgap[sub_panels-1]->Draw("Psame");

				if( sub_panels==8 ){
					temp_allEta_tight_rapgap[sub_panels-1]->Fit("pol0","","",0,20.);
					temp_allEta_loose_rapgap[sub_panels-1]->Fit("pol0","","",0,20.);
				}
				else if( sub_panels==16 ){
					temp_allEta_tight_rapgap[sub_panels-1]->Fit("pol0","","",1,20.);
					temp_allEta_loose_rapgap[sub_panels-1]->Fit("pol0","","",1,20.);
				}
				else{
					temp_allEta_tight_rapgap[sub_panels-1]->Fit("pol0");
					temp_allEta_loose_rapgap[sub_panels-1]->Fit("pol0");
				}

				TF1* fitFunction_tight = (TF1*) temp_allEta_tight_rapgap[sub_panels-1]->GetFunction("pol0");
				TF1* fitFunction_loose = (TF1*) temp_allEta_loose_rapgap[sub_panels-1]->GetFunction("pol0");
				double latex_height = 0.57;
				if( sub_panels > 4 ) latex_height = 0.72;

				double tight_fit_value = fitFunction_tight->GetParameter(0);
				double loose_fit_value = fitFunction_loose->GetParameter(0);
				double fit_value = loose_fit_value;
				if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
					fit_value = tight_fit_value;
				}

				TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.4f",fit_value) );
				rFit->SetNDC();
				rFit->SetTextSize(18);
				rFit->SetTextFont(44);
				rFit->Draw("same");
			}

			if( sub_panels < 5 ){
				TLatex* r46 = new TLatex(latex_starting_y_1[sub_panels-1],0.76, y_bins[iy]+"<y<"+y_bins[iy+1]);
				r46->SetNDC();
				r46->SetTextSize(18);
				r46->SetTextFont(44);
				r46->Draw("same");
			}
			else{
				TLatex* r46 = new TLatex(latex_starting_y_1[sub_panels-1],0.9, y_bins[iy]+"<y<"+y_bins[iy+1]);
				r46->SetNDC();
				r46->SetTextSize(18);
				r46->SetTextFont(44);
				r46->Draw("same");
			}

			if(sub_panels < 5){
				TLatex* r47 = new TLatex(latex_starting_eta_1[sub_panels-1],0.65, Q2_bins[iQ2]+"<Q^{2}<"+Q2_bins[iQ2+1]);
				r47->SetNDC();
				r47->SetTextSize(18);
				r47->SetTextFont(44);
				r47->Draw("same");
			}
			else{
				TLatex* r47 = new TLatex(latex_starting_eta_1[sub_panels-1],0.80, Q2_bins[iQ2]+"<Q^{2}<"+Q2_bins[iQ2+1]);
				r47->SetNDC();
				r47->SetTextSize(18);
				r47->SetTextFont(44);
				r47->Draw("same");
			}
			
			sub_panels++;
		}

		}

		TLatex* r503 = new TLatex(0.53,0.87, "-1.2 < #eta < 1.6");
		c3->cd(4);
		r503->SetNDC();
		r503->SetTextSize(18);
		r503->SetTextFont(44);
		r503->Draw("same");

		TLatex* r513 = new TLatex(0.22,0.87, "ep 27.5x920 GeV");
		c3->cd(1);
		r513->SetNDC();
		r513->SetTextSize(18);
		r513->SetTextFont(44);
		r513->Draw("same");

		TLatex* r523 = new TLatex(0.25,0.74, "H1");
		c3->cd(1);
		r523->SetNDC();
		r523->SetTextSize(20);
		r523->SetTextFont(44);
		r523->Draw("same");

		c3->cd(16);
		TLegend *w53 = new TLegend(0.05,0.26,0.5,0.45);
		w53->SetLineColor(kWhite);
		w53->SetFillColor(0);
		w53->SetTextSize(17);
		w53->SetTextFont(45);
		if( doDjango_ ){
			w53->AddEntry(temp_allEta_tight_django[0], "django "+filename_allEta_1+"/default ","P");
			w53->AddEntry(temp_allEta_loose_django[0], "django "+filename_allEta_2+"/default ","P");
		}
		else{
			w53->AddEntry(temp_allEta_tight_rapgap[0], "rapgap "+filename_allEta_1+"/default ","P");
			w53->AddEntry(temp_allEta_loose_rapgap[0], "rapgap "+filename_allEta_2+"/default ","P");
		}
		w53->Draw("same");

		if(doMCcomparison_){

			TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,800);
			c2->Divide(4,3,0,0);
			TH1D* base2 = (TH1D*) base1->Clone("base2");
			base2->GetYaxis()->SetRangeUser(0,3);
			TH1D* temp_django[16];

			int sub_panels = 1;
			for(int ieta = 0; ieta < 3; ieta++){
				for(int iy = 0; iy < 4; iy++){

					c2->cd(sub_panels);
					
					if(sub_panels==1||sub_panels==5||sub_panels==9) gPad->SetLeftMargin(0.2);
					if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
					if(sub_panels>=9) gPad->SetBottomMargin(0.2);
					gPad->SetTicks();
					gPad->SetLogy(0);
					base2->Draw("");

					temp_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->Clone(Form("temp_django_%d",sub_panels));
					temp_django[sub_panels-1]->SetMarkerStyle(24);

					//django
					for(int ibin=0;ibin<temp_django[sub_panels-1]->GetNbinsX();ibin++){
						
						double defau = hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
						double defau_err = hist_unfolded_from_django[0][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
						if( defau <= 0 ) continue;
						
						double rapgap = hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->GetBinContent(ibin+1);
						double rapgap_err = hist_unfolded_from_rapgap[0][ieta][Q2_BIN][iy]->GetBinError(ibin+1);
						if( rapgap <= 0 ) continue;
					
						temp_django[sub_panels-1]->SetBinContent(ibin+1, rapgap/defau);
						temp_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,rapgap,defau_err,rapgap_err));
					}

					temp_django[sub_panels-1]->Draw("Psame");
					if(sub_panels==1) temp_django[sub_panels-1]->Fit("pol0","","",0,6);
					else temp_django[sub_panels-1]->Fit("pol0");

					TF1* fitFunction_tight = (TF1*) temp_django[sub_panels-1]->GetFunction("pol0");
					TF1* fitFunction_loose = (TF1*) temp_django[sub_panels-1]->GetFunction("pol0");
					double latex_height = 0.57;
					if( sub_panels > 4 ) latex_height = 0.72;

					double tight_fit_value = fitFunction_tight->GetParameter(0);
					double loose_fit_value = fitFunction_loose->GetParameter(0);
					double fit_value = loose_fit_value;
					if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
						fit_value = tight_fit_value;
					}

					TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.3f",fit_value) );
					rFit->SetNDC();
					rFit->SetTextSize(18);
					rFit->SetTextFont(44);
					rFit->Draw("same");

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
					sub_panels++;
					}
					
				}
				c2->cd(4);
				r50->Draw("same");
				c2->cd(1);
				r51->Draw("same");
				c2->cd(1);
				r52->Draw("same");
				c2->cd(12);
				TLegend *w6 = new TLegend(0.1,0.68,0.5,0.75);
				w6->SetLineColor(kWhite);
				w6->SetFillColor(0);
				w6->SetTextSize(17);
				w6->SetTextFont(45);
				w6->AddEntry(temp_django[0], "rapgap/django ","P");
				w6->Draw("same");

				// c2->Print(Form("../figures/systematics/MCmodels/Pn_rapgapDiffractive_sys_fit_%d.pdf", Q2_BIN));

				TCanvas* c4 = new TCanvas("c4","c4",1,1,1000,1000);
				c4->Divide(4,4,0,0);
				TH1D* base4 = (TH1D*) base3->Clone("base4");
				base4->GetYaxis()->SetRangeUser(0,3);
				TH1D* temp_allEta_django[16];

				sub_panels = 1;
				for(int iy = 0; iy < 4; iy++){
				for(int iQ2 = 0; iQ2 < 4; iQ2++){

					c4->cd(sub_panels);
					if(sub_panels==1||sub_panels==5||sub_panels==9||sub_panels==13) gPad->SetLeftMargin(0.2);
					if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
					if(sub_panels>=13) gPad->SetBottomMargin(0.2);
					gPad->SetTicks();
					gPad->SetLogy(0);
					base4->Draw("");

					temp_allEta_django[sub_panels-1] = (TH1D*) hist_unfolded_from_django_allEta[0][iQ2][iy]->Clone(Form("temp_allEta_django_%d",sub_panels));
					temp_allEta_django[sub_panels-1]->SetMarkerStyle(24);

					//django
					for(int ibin=0;ibin<temp_allEta_django[sub_panels-1]->GetNbinsX();ibin++){
						
						double defau = hist_unfolded_from_django_allEta[0][iQ2][iy]->GetBinContent(ibin+1);
						double defau_err = hist_unfolded_from_django_allEta[0][iQ2][iy]->GetBinError(ibin+1);
						if( defau <= 0 ) continue;
						
						double rapgap = hist_unfolded_from_rapgap_allEta[0][iQ2][iy]->GetBinContent(ibin+1);
						double rapgap_err = hist_unfolded_from_rapgap_allEta[0][iQ2][iy]->GetBinError(ibin+1);
						if( rapgap <= 0 ) continue;
					
						temp_allEta_django[sub_panels-1]->SetBinContent(ibin+1, rapgap/defau);
						temp_allEta_django[sub_panels-1]->SetBinError(ibin+1, calColError(defau,rapgap,defau_err,rapgap_err));
					}

					temp_allEta_django[sub_panels-1]->Draw("Psame");
					

					if(sub_panels==8) temp_allEta_django[sub_panels-1]->Fit("pol0","","",0,18);
					else if(sub_panels==16) temp_allEta_django[sub_panels-1]->Fit("pol0","","",1,18);
					else temp_allEta_django[sub_panels-1]->Fit("pol0","","",0,20);

					TF1* fitFunction_tight = (TF1*) temp_allEta_django[sub_panels-1]->GetFunction("pol0");
					TF1* fitFunction_loose = (TF1*) temp_allEta_django[sub_panels-1]->GetFunction("pol0");
					double latex_height = 0.57;
					if( sub_panels > 4 ) latex_height = 0.72;

					double tight_fit_value = fitFunction_tight->GetParameter(0);
					double loose_fit_value = fitFunction_loose->GetParameter(0);
					double fit_value = loose_fit_value;
					if( fabs(tight_fit_value - 1.0) > fabs(loose_fit_value - 1.0) ){
						fit_value = tight_fit_value;
					}

					TLatex* rFit = new TLatex(latex_starting_y[sub_panels-1],latex_height, Form("fit ~ %.3f",fit_value) );
					rFit->SetNDC();
					rFit->SetTextSize(18);
					rFit->SetTextFont(44);
					rFit->Draw("same");

					if( sub_panels < 5 ){
						TLatex* r46 = new TLatex(latex_starting_y_1[sub_panels-1],0.76, y_bins[iy]+"<y<"+y_bins[iy+1]);
						r46->SetNDC();
						r46->SetTextSize(18);
						r46->SetTextFont(44);
						r46->Draw("same");
					}
					else{
						TLatex* r46 = new TLatex(latex_starting_y_1[sub_panels-1],0.9, y_bins[iy]+"<y<"+y_bins[iy+1]);
						r46->SetNDC();
						r46->SetTextSize(18);
						r46->SetTextFont(44);
						r46->Draw("same");
					}

					if(sub_panels < 5){
						TLatex* r47 = new TLatex(latex_starting_eta_1[sub_panels-1],0.65, Q2_bins[iQ2]+"<Q^{2}<"+Q2_bins[iQ2+1]);
						r47->SetNDC();
						r47->SetTextSize(18);
						r47->SetTextFont(44);
						r47->Draw("same");
					}
					else{
						TLatex* r47 = new TLatex(latex_starting_eta_1[sub_panels-1],0.80, Q2_bins[iQ2]+"<Q^{2}<"+Q2_bins[iQ2+1]);
						r47->SetNDC();
						r47->SetTextSize(18);
						r47->SetTextFont(44);
						r47->Draw("same");
					}
			
					sub_panels++;
					}
					
				}
				c4->cd(4);
				r503->Draw("same");
				c4->cd(1);
				r51->Draw("same");
				c4->cd(1);
				r52->Draw("same");
				c4->cd(16);
				TLegend *w64 = new TLegend(0.1,0.68,0.5,0.75);
				w64->SetLineColor(kWhite);
				w64->SetFillColor(0);
				w64->SetTextSize(17);
				w64->SetTextFont(45);
				w64->AddEntry(temp_allEta_django[0], "rapgap/django ","P");
				w64->Draw("same");
				
				// c4->Print(Form("../figures/systematics/MCmodels/Pn_rapgapDiffractive_allEta_sys_fit.pdf"));

				delete c4;
				TCanvas* c7 = new TCanvas("c7","c7",1,1,600,600);
				gPad->SetTicks();
				gPad->SetLeftMargin(0.1);
				gPad->SetRightMargin(0.1);
				gPad->SetBottomMargin(0.15);
				gPad->SetTopMargin(0.15);

				//testing radiative effect only
				TFile* file_rapgap = new TFile("data_django_rad.root");
				TFile* file_django = new TFile("data_rapgap_rad.root");

				TGraphErrors* gr_Nch[4][2];
				TGraphErrors* gr_Var[4][2];

				for(int i=0;i<4;i++){
					gr_Nch[i][0] = (TGraphErrors*) file_rapgap->Get(Form(";%d",i+1));
					gr_Nch[i][1] = (TGraphErrors*) file_django->Get(Form(";%d",i+1));

					gr_Var[i][0] = (TGraphErrors*) file_rapgap->Get(Form(";%d",i+5));
					gr_Var[i][1] = (TGraphErrors*) file_django->Get(Form(";%d",i+5));
				}

				TGraphErrors* gr_data_ratio[4][2];
				TGraphErrors* gr_data_variance_ratio[4][2];
				for(int iQ2=0;iQ2<4;iQ2++){
					for(int j=0;j<2;j++){
						gr_data_ratio[iQ2][j]=new TGraphErrors();
						gr_data_variance_ratio[iQ2][j]=new TGraphErrors();
					}
				}
				for(int iQ2=0;iQ2<4;iQ2++){
					for(int iy=0;iy<4;iy++){
						double x,y,x1,y1;
						gr_Nch[iQ2][0]->GetPoint(iy,x,y);
						gr_Nch[iQ2][1]->GetPoint(iy,x1,y1);
						gr_data_ratio[iQ2][0]->SetPoint(iy, x, y1/y);
						
						gr_Var[iQ2][0]->GetPoint(iy,x,y);
						gr_Var[iQ2][1]->GetPoint(iy,x1,y1);
						gr_data_variance_ratio[iQ2][0]->SetPoint(iy, x, y1/y);

					}
				}

				TH1D* base7 = makeHist("base7", "", "y", "ratio #LT N_{ch} #GT", 100,0,0.6,kBlack);
				base7->GetYaxis()->SetRangeUser(0.9, 1.1);
				base7->GetXaxis()->SetTitleColor(kBlack);
				fixedFontHist1D(base7,1.2,1.0);

				base7->GetYaxis()->SetTitleSize(base7->GetYaxis()->GetTitleSize()*1.3);
				base7->GetXaxis()->SetTitleSize(base7->GetXaxis()->GetTitleSize()*1.3);
				base7->GetYaxis()->SetLabelSize(base7->GetYaxis()->GetLabelSize()*1.3);
				base7->GetXaxis()->SetLabelSize(base7->GetXaxis()->GetLabelSize()*1.3);
				base7->GetXaxis()->SetNdivisions(4,6,0);
				base7->GetYaxis()->SetNdivisions(4,6,0);
				base7->Draw();

				gr_data_ratio[0][0]->SetMarkerStyle(20);
				gr_data_ratio[0][0]->SetMarkerSize(1.4);
				gr_data_ratio[0][0]->SetMarkerColor(kBlack);

				gr_data_ratio[1][0]->SetMarkerStyle(20);
				gr_data_ratio[1][0]->SetMarkerSize(1.4);
				gr_data_ratio[1][0]->SetMarkerColor(kBlack);

				gr_data_ratio[2][0]->SetMarkerStyle(20);
				gr_data_ratio[2][0]->SetMarkerSize(1.4);
				gr_data_ratio[2][0]->SetMarkerColor(kBlack);

				gr_data_ratio[3][0]->SetMarkerStyle(24);
				gr_data_ratio[3][0]->SetMarkerSize(1.4);
				gr_data_ratio[3][0]->SetMarkerColor(kRed);
		
				gr_data_ratio[0][0]->Draw("PEsame");
				gr_data_ratio[1][0]->Draw("PEsame");
				gr_data_ratio[2][0]->Draw("PEsame");
				gr_data_ratio[3][0]->Draw("PEsame");

				// c7->Print(Form("../figures/systematics/radiation/Nch.pdf"));
				// return;
		}
	
		// c1->Print(Form("../figures/systematics/HadronEnergyScale/Pn_sys_fit_%d.pdf", Q2_BIN));
		// c3->Print(Form("../figures/systematics/HadronEnergyScale/Pn_allEta_sys_fit.pdf"));

	//Nch 1st and 2nd moment.
	TCanvas* temp[3][4][4]; 

	TF1 *tf3[3][4][4];
	sub_panels = 1;
	for(int ifile = 0;ifile<3;ifile++){
		for(int iy = 0; iy < 4; iy++){
		for(int iQ2 = 0; iQ2 < 4; iQ2++){

			// temp[ifile][iy][iQ2] = new TCanvas();
			
			if(sub_panels==1||sub_panels==5||sub_panels==9||sub_panels==13) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=13) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(1);

			//single NBD
			// tf3[ifile][iQ2][iy] = new TF1(Form("tf3_%d_%d",iQ2,iy),"[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])",0,25);
			//double NBD
			tf3[ifile][iQ2][iy] = new TF1(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])+[3]*([4]/[4]+1)*ROOT::Math::negative_binomial_pdf(x[0],[4],[5])",0,25);
			tf3[ifile][iQ2][iy]->SetParameter(0,0.5);
			tf3[ifile][iQ2][iy]->SetParameter(1,0.081);
			tf3[ifile][iQ2][iy]->SetParameter(2,0.7);
			tf3[ifile][iQ2][iy]->SetParameter(3,1.0);
			tf3[ifile][iQ2][iy]->SetParameter(4,0.3);
			tf3[ifile][iQ2][iy]->SetParameter(5,1.6);
			tf3[ifile][iQ2][iy]->SetLineWidth(2);
			if( sub_panels == 1 || sub_panels == 2 || sub_panels == 3 || sub_panels == 6 ){
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,15);
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,15);
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,15);
			}
			else{
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,25);
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,25);
				hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Fit(Form("tf3_%d_%d_%d",ifile,iQ2,iy),"RME0","",4,25);
			}
			sub_panels++;
		}
		}
	}

	TGraphErrors* gr_data[3][4];
	TGraphErrors* gr_data_variance[3][4];
	double yREC_es_ave[] = {0.05464,0.1095,0.2179,0.4121};

	for(int ifile=0;ifile<3;ifile++){
	for(int iQ2=0;iQ2<4;iQ2++){
		gr_data[ifile][iQ2] = new TGraphErrors();
		gr_data_variance[ifile][iQ2] = new TGraphErrors();

		//problem with the data is not unit bin width
		for(int iy=0;iy<4;iy++){
			//data
			TH1D* histForMultiplicity = (TH1D*) hist_unfolded_from_rapgap_allEta[ifile][iQ2][iy]->Clone(Form("histForMultiplicity_%d_%d_%d",ifile,iQ2,iy));
			double sum = 0.;
			double sum_pn = 0.;
			double err2 = 0.;
			double bin_to_stop = 25;
			//bin to stop
			for(int ibin=0;ibin<histForMultiplicity->GetNbinsX();ibin++){
				if( (histForMultiplicity->GetBinContent(ibin+1) < 1e-5 || ibin==histForMultiplicity->GetNbinsX()-1) && ibin > 0 ){
					bin_to_stop = histForMultiplicity->GetBinCenter(ibin+1);
					break;
				}
			}
			//use the first 4 bins in data and after N > 3, use NBD extrapolations.
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){
				double N = ibin;
				double Pn = 0.;
				double Pn_err = 0;
				if( ibin < 4 ){
					N = histForMultiplicity->GetBinCenter(ibin+1);
					Pn = histForMultiplicity->GetBinContent(ibin+1);
					Pn_err = histForMultiplicity->GetBinError(ibin+1);
					sum += N*Pn;
					sum_pn += Pn;
					err2 += (N*Pn_err)*(N*Pn_err);
				}
				else{
					Pn = tf3[ifile][iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += N*Pn;
					sum_pn += Pn;
					err2 += (N*Pn_err)*(N*Pn_err);
				}
			}
			gr_data[ifile][iQ2]->SetPoint(iy, yREC_es_ave[iy], sum/sum_pn );
			gr_data[ifile][iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			double first_moment = sum/sum_pn;
		
			//use mean to compute variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0; ibin<(int)bin_to_stop+1;ibin++){
				double N = ibin;
				double Pn = 0.;
				double Pn_err = 0;
				if( ibin < 4 ){
					N = histForMultiplicity->GetBinCenter(ibin+1);
					Pn = histForMultiplicity->GetBinContent(ibin+1);
					Pn_err = histForMultiplicity->GetBinError(ibin+1);
					sum += TMath::Power((N-first_moment),2)*Pn;
					sum_pn += Pn;
					err2 += (TMath::Power((N-first_moment),2)*Pn_err)*(TMath::Power((N-first_moment),2)*Pn_err);
				}
				else{
					Pn = tf3[ifile][iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += TMath::Power((N-first_moment),2)*Pn;
					sum_pn += Pn;
					err2 += (TMath::Power((N-first_moment),2)*Pn_err)*(TMath::Power((N-first_moment),2)*Pn_err);
				}
			}
			gr_data_variance[ifile][iQ2]->SetPoint(iy, yREC_es_ave[iy], sum/sum_pn );
			gr_data_variance[ifile][iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			
			}
		}
		}
		// delete temp;
		TCanvas* c5 = new TCanvas("c5","c5",1,1,600,600);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.1);
		gPad->SetRightMargin(0.1);
		gPad->SetBottomMargin(0.15);
		gPad->SetTopMargin(0.15);

		TGraphErrors* gr_data_ratio[4][2];
		TGraphErrors* gr_data_variance_ratio[4][2];
		for(int iQ2=0;iQ2<4;iQ2++){
			for(int j=0;j<2;j++){
				gr_data_ratio[iQ2][j]=new TGraphErrors();
				gr_data_variance_ratio[iQ2][j]=new TGraphErrors();
			}
		}
		for(int iQ2=0;iQ2<4;iQ2++){
			for(int iy=0;iy<4;iy++){
				double x,y,x1,y1,x2,y2;
				gr_data[0][iQ2]->GetPoint(iy,x,y);
				gr_data[1][iQ2]->GetPoint(iy,x1,y1);
				gr_data[2][iQ2]->GetPoint(iy,x2,y2);

				gr_data_ratio[iQ2][0]->SetPoint(iy, x, y1/y);
				gr_data_ratio[iQ2][1]->SetPoint(iy, x, y2/y);
				
				gr_data_variance[0][iQ2]->GetPoint(iy,x,y);
				gr_data_variance[1][iQ2]->GetPoint(iy,x1,y1);
				gr_data_variance[2][iQ2]->GetPoint(iy,x2,y2);

				gr_data_variance_ratio[iQ2][0]->SetPoint(iy, x, y1/y);
				gr_data_variance_ratio[iQ2][1]->SetPoint(iy, x, y2/y);

			}
		}

		TH1D* base5 = makeHist("base5", "", "y", "ratio #LT N_{ch} #GT", 100,0,0.6,kBlack);
		base5->GetYaxis()->SetRangeUser(0.9, 1.1);
		base5->GetXaxis()->SetTitleColor(kBlack);
		fixedFontHist1D(base5,1.2,1.0);

		base5->GetYaxis()->SetTitleSize(base5->GetYaxis()->GetTitleSize()*1.3);
		base5->GetXaxis()->SetTitleSize(base5->GetXaxis()->GetTitleSize()*1.3);
		base5->GetYaxis()->SetLabelSize(base5->GetYaxis()->GetLabelSize()*1.3);
		base5->GetXaxis()->SetLabelSize(base5->GetXaxis()->GetLabelSize()*1.3);
		base5->GetXaxis()->SetNdivisions(4,6,0);
		base5->GetYaxis()->SetNdivisions(4,6,0);
		base5->Draw();

		gr_data_ratio[0][0]->SetMarkerStyle(20);
		gr_data_ratio[0][0]->SetMarkerSize(1.4);
		gr_data_ratio[0][0]->SetMarkerColor(kBlack);

		gr_data_ratio[0][1]->SetMarkerStyle(30);
		gr_data_ratio[0][1]->SetMarkerSize(1.4);
		gr_data_ratio[0][1]->SetMarkerColor(kBlue);

		gr_data_ratio[1][0]->SetMarkerStyle(20);
		gr_data_ratio[1][0]->SetMarkerSize(1.4);
		gr_data_ratio[1][0]->SetMarkerColor(kBlack);

		gr_data_ratio[1][1]->SetMarkerStyle(30);
		gr_data_ratio[1][1]->SetMarkerSize(1.4);
		gr_data_ratio[1][1]->SetMarkerColor(kBlue);

		gr_data_ratio[2][0]->SetMarkerStyle(20);
		gr_data_ratio[2][0]->SetMarkerSize(1.4);
		gr_data_ratio[2][0]->SetMarkerColor(kBlack);

		gr_data_ratio[2][1]->SetMarkerStyle(30);
		gr_data_ratio[2][1]->SetMarkerSize(1.4);
		gr_data_ratio[2][1]->SetMarkerColor(kBlue);

		gr_data_ratio[3][0]->SetMarkerStyle(24);
		gr_data_ratio[3][0]->SetMarkerSize(1.4);
		gr_data_ratio[3][0]->SetMarkerColor(kRed);

		gr_data_ratio[3][1]->SetMarkerStyle(29);
		gr_data_ratio[3][1]->SetMarkerSize(1.8);
		gr_data_ratio[3][1]->SetMarkerColor(kGreen-2);

		gr_data_ratio[0][0]->Draw("PEsame");
		gr_data_ratio[0][1]->Draw("PEsame");
		gr_data_ratio[1][0]->Draw("PEsame");
		gr_data_ratio[1][1]->Draw("PEsame");
		gr_data_ratio[2][0]->Draw("PEsame");
		gr_data_ratio[2][1]->Draw("PEsame");
		gr_data_ratio[3][0]->Draw("PEsame");
		gr_data_ratio[3][1]->Draw("PEsame");


		TCanvas* c6 = new TCanvas("c6","c6",1,1,600,600);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.1);
		gPad->SetRightMargin(0.1);
		gPad->SetBottomMargin(0.15);
		gPad->SetTopMargin(0.15);

		TH1D* base6 = makeHist("base6", "", "y", "ratio Var(N_{ch})", 100,0,0.6,kBlack);
		base6->GetYaxis()->SetRangeUser(0.9, 1.1);
		base6->GetXaxis()->SetTitleColor(kBlack);
		fixedFontHist1D(base6,1.2,1.0);

		base6->GetYaxis()->SetTitleSize(base6->GetYaxis()->GetTitleSize()*1.3);
		base6->GetXaxis()->SetTitleSize(base6->GetXaxis()->GetTitleSize()*1.3);
		base6->GetYaxis()->SetLabelSize(base6->GetYaxis()->GetLabelSize()*1.3);
		base6->GetXaxis()->SetLabelSize(base6->GetXaxis()->GetLabelSize()*1.3);
		base6->GetXaxis()->SetNdivisions(4,6,0);
		base6->GetYaxis()->SetNdivisions(4,6,0);
		base6->Draw();

		gr_data_variance_ratio[0][0]->SetMarkerStyle(20);
		gr_data_variance_ratio[0][0]->SetMarkerSize(1.4);
		gr_data_variance_ratio[0][0]->SetMarkerColor(kBlack);

		gr_data_variance_ratio[0][1]->SetMarkerStyle(30);
		gr_data_variance_ratio[0][1]->SetMarkerSize(1.4);
		gr_data_variance_ratio[0][1]->SetMarkerColor(kBlue);

		gr_data_variance_ratio[1][0]->SetMarkerStyle(20);
		gr_data_variance_ratio[1][0]->SetMarkerSize(1.4);
		gr_data_variance_ratio[1][0]->SetMarkerColor(kBlack);

		gr_data_variance_ratio[1][1]->SetMarkerStyle(30);
		gr_data_variance_ratio[1][1]->SetMarkerSize(1.4);
		gr_data_variance_ratio[1][1]->SetMarkerColor(kBlue);

		gr_data_variance_ratio[2][0]->SetMarkerStyle(20);
		gr_data_variance_ratio[2][0]->SetMarkerSize(1.4);
		gr_data_variance_ratio[2][0]->SetMarkerColor(kBlack);

		gr_data_variance_ratio[2][1]->SetMarkerStyle(30);
		gr_data_variance_ratio[2][1]->SetMarkerSize(1.4);
		gr_data_variance_ratio[2][1]->SetMarkerColor(kBlue);

		gr_data_variance_ratio[3][0]->SetMarkerStyle(24);
		gr_data_variance_ratio[3][0]->SetMarkerSize(1.4);
		gr_data_variance_ratio[3][0]->SetMarkerColor(kRed);

		gr_data_variance_ratio[3][1]->SetMarkerStyle(29);
		gr_data_variance_ratio[3][1]->SetMarkerSize(1.8);
		gr_data_variance_ratio[3][1]->SetMarkerColor(kGreen-2);

		gr_data_variance_ratio[0][0]->Draw("PEsame");
		gr_data_variance_ratio[0][1]->Draw("PEsame");
		gr_data_variance_ratio[1][0]->Draw("PEsame");
		gr_data_variance_ratio[1][1]->Draw("PEsame");
		gr_data_variance_ratio[2][0]->Draw("PEsame");
		gr_data_variance_ratio[2][1]->Draw("PEsame");
		gr_data_variance_ratio[3][0]->Draw("PEsame");
		gr_data_variance_ratio[3][1]->Draw("PEsame");
		
		// TFile* output = new TFile("data_django_rad.root","RECREATE");
		// gr_data[0][0]->Write();
		// gr_data[0][1]->Write();
		// gr_data[0][2]->Write();
		// gr_data[0][3]->Write();

		// gr_data_variance[0][0]->Write();
		// gr_data_variance[0][1]->Write();
		// gr_data_variance[0][2]->Write();
		// gr_data_variance[0][3]->Write();

		// c5->Print(Form("../figures/systematics/ElectronEnergyScale/Nch.pdf"));
		// c6->Print(Form("../figures/systematics/ElectronEnergyScale/Var.pdf"));

	
}