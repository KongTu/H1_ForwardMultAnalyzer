#include "RiceStyle.h"
#include <string>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

TH1D* shiftedHist( TH1D* hist ){

	TString name = hist->GetName();
	TH1D* hist_new = new TH1D(name,name,30,-0.5,29.5);
	for(int i=0;i<hist->GetNbinsX();i++){

		double bincenter=hist->GetBinCenter(i+1);
		double binvalue=hist->GetBinContent(i+1);
		double binwidth=hist->GetBinWidth(i+1);
		double binerror=hist->GetBinError(i+1);

		hist_new->SetBinContent(i+1, binvalue);
		hist_new->SetBinError(i+1, binerror);
	}

	return hist_new;
}
const double y_beam = 7.6;
const double sys_Nch = 0.054; //systematics of <Nch> is 6.8%
const double sys_Var = 0.067; //systematics of Var is 7.9%

void makeHisto_result_HCM(const int ifile = 0, const bool draw_ = false, const int pickMCEG_ = 0, const bool doRadCorr_ = true){

	gStyle->SetErrorX(0);
	
	TFile* file_allEta[10];
	file_allEta[ifile] = new TFile("../rootfiles/results_unfoldingOutput_yetaStar_hadCaliNew.root");

	//radiative corrections
	TFile* file_RadCorr_rapgap = new TFile("./RadCorr_RAPGAP.root");
	TH1D* hCorr_Pn_GEN[4][4][4];
	TH1D* hCorr_Pn_GEN_HCM[4][4];
	for(int iQ2=0;iQ2<4;iQ2++){
		for(int iy=0;iy<4;iy++){
			hCorr_Pn_GEN_HCM[iQ2][iy] = (TH1D*) file_RadCorr_rapgap->Get(Form("hCorr_Pn_GEN_HCM_%d_%d",iQ2,iy));
			hCorr_Pn_GEN_HCM[iQ2][iy] = shiftedHist( hCorr_Pn_GEN_HCM[iQ2][iy] );
			for(int ieta=0;ieta<4;ieta++){
				hCorr_Pn_GEN[iQ2][iy][ieta] = (TH1D*) file_RadCorr_rapgap->Get(Form("hCorr_Pn_GEN_%d_%d_%d",iQ2,iy,ieta));
				hCorr_Pn_GEN[iQ2][iy][ieta] = shiftedHist ( hCorr_Pn_GEN[iQ2][iy][ieta] );
			}
		}
	}
	TFile* file_RadCorr_django = new TFile("./RadCorr_DJANGO.root");
	TH1D* django_hCorr_Pn_GEN[4][4][4];
	TH1D* django_hCorr_Pn_GEN_HCM[4][4];
	for(int iQ2=0;iQ2<4;iQ2++){
		for(int iy=0;iy<4;iy++){
			django_hCorr_Pn_GEN_HCM[iQ2][iy] = (TH1D*) file_RadCorr_django->Get(Form("hCorr_Pn_GEN_HCM_%d_%d",iQ2,iy));
			django_hCorr_Pn_GEN_HCM[iQ2][iy] = shiftedHist( django_hCorr_Pn_GEN_HCM[iQ2][iy] );
			for(int ieta=0;ieta<4;ieta++){
				django_hCorr_Pn_GEN[iQ2][iy][ieta] = (TH1D*) file_RadCorr_django->Get(Form("hCorr_Pn_GEN_%d_%d_%d",iQ2,iy,ieta));
				django_hCorr_Pn_GEN[iQ2][iy][ieta] = shiftedHist( django_hCorr_Pn_GEN[iQ2][iy][ieta] );
			}
		}
	}
	//end radiative MC corr files

	TString radname = "radCorr";
	if(!doRadCorr_) radname = "noRadCorr";
	TString MCEGname = "rapgap";
	if( pickMCEG_ == 1 ) MCEGname = "django";
	if( pickMCEG_ == 2 ) MCEGname = "pythia";

	TString histo_name;
	TString Q2_bins[5]={"5","10","20","40","100"};
	TString y_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	TString x_average_Q2_5_10[]={"0.0014","0.0007","0.00033","0.00017"};
	TString x_average_Q2_10_20[]={"0.0023","0.0013","0.00066","0.00034"};
	TString x_average_Q2_20_40[]={"0.0052","0.0026","0.0013","0.00068"};
	TString x_average_Q2_40_100[]={"0.011","0.0056","0.0028","0.0013"};

	TH1D* hist_unfolded_from_rapgap_allEta[4][4];
	TH1D* hist_unfolded_from_django_allEta[4][4];
	for(int j = 0; j < 4; j++){
		for(int k = 0; k < 4; k++){

			histo_name= "data_from_djangoh/hist_unfolded_(0<=etaStar<4)("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
			y_bins[k+1]+")";

			hist_unfolded_from_django_allEta[j][k] = (TH1D*) file_allEta[ifile]->Get( histo_name );
			hist_unfolded_from_django_allEta[j][k]->SetMarkerStyle(24);
			hist_unfolded_from_django_allEta[j][k]->SetMarkerColor(kBlack);
			hist_unfolded_from_django_allEta[j][k]->SetLineColor(kBlack);

			histo_name = "data_from_rapgap/hist_unfolded_(0<=etaStar<4)("+Q2_bins[j]+"<=Q2MC_es_mini<"+Q2_bins[j+1]+")("+y_bins[k]+"<=yMC_es_mini<"+
			y_bins[k+1]+")";

			hist_unfolded_from_rapgap_allEta[j][k] = (TH1D*) file_allEta[ifile]->Get( histo_name );
			hist_unfolded_from_rapgap_allEta[j][k]->SetMarkerStyle(20);
			hist_unfolded_from_rapgap_allEta[j][k]->SetMarkerColor(kBlack);
			hist_unfolded_from_rapgap_allEta[j][k]->SetLineColor(kBlack);

			//radiative correction to Rapgap based on Rapgap
			for(int ibin=0;ibin<hist_unfolded_from_rapgap_allEta[j][k]->GetNbinsX();ibin++){
				double value = hist_unfolded_from_rapgap_allEta[j][k]->GetBinContent(ibin+1);
				double error = hist_unfolded_from_rapgap_allEta[j][k]->GetBinError(ibin+1);
				double bincenter = hist_unfolded_from_rapgap_allEta[j][k]->GetBinCenter(ibin+1);
				double corr = hCorr_Pn_GEN_HCM[j][k]->GetBinContent( hCorr_Pn_GEN_HCM[j][k]->FindBin(bincenter) );
				if(doRadCorr_){
					hist_unfolded_from_rapgap_allEta[j][k]->SetBinContent(ibin+1, value*corr);
					hist_unfolded_from_rapgap_allEta[j][k]->SetBinError(ibin+1, error*corr);
				}
			}
			//radiative correction to Django based on Django
			for(int ibin=0;ibin<hist_unfolded_from_django_allEta[j][k]->GetNbinsX();ibin++){
				double value = hist_unfolded_from_django_allEta[j][k]->GetBinContent(ibin+1);
				double error = hist_unfolded_from_django_allEta[j][k]->GetBinError(ibin+1);
				double bincenter = hist_unfolded_from_django_allEta[j][k]->GetBinCenter(ibin+1);
				double corr = hCorr_Pn_GEN_HCM[j][k]->GetBinContent( hCorr_Pn_GEN_HCM[j][k]->FindBin(bincenter) );
				if(doRadCorr_){
					hist_unfolded_from_django_allEta[j][k]->SetBinContent(ibin+1, value*corr);
					hist_unfolded_from_django_allEta[j][k]->SetBinError(ibin+1, error*corr);
				}
			}
			

		}
	}
	TFile* file_MC_django = new TFile("../../minitree_output/Pn_hist_django_hadCaliNew_reweigh.root");
	TFile* file_MC_rapgap = new TFile("../../minitree_output/Pn_hist_rapgap_hadCaliNew_reweigh.root");
	TFile* file_MC_pythia = new TFile("./Pythia8235_H1paper.root");
	TH1D* hist_mc_django[4][4];
	TH1D* hist_mc_rapgap[4][4];
	TH1D* hist_mc_pythia[4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
				
			hist_mc_django[i][j] = (TH1D*) file_MC_django->Get(Form("h_Pn_GEN_HCM_%d_%d",i,j));
			hist_mc_rapgap[i][j] = (TH1D*) file_MC_rapgap->Get(Form("h_Pn_GEN_HCM_%d_%d",i,j));
			hist_mc_pythia[i][j] = (TH1D*) file_MC_pythia->Get(Form("Pn_HCM_%d_%d",i,j));
			
			hist_mc_django[i][j] = shiftedHist( hist_mc_django[i][j] );
			hist_mc_rapgap[i][j] = shiftedHist( hist_mc_rapgap[i][j] );
			hist_mc_pythia[i][j] = shiftedHist( hist_mc_pythia[i][j] );

			hist_mc_django[i][j]->SetFillStyle(1001);
			hist_mc_django[i][j]->SetFillColorAlpha(kGreen-2,0.4);
			hist_mc_django[i][j]->SetMarkerColor(kGreen-2);

			hist_mc_rapgap[i][j]->SetFillStyle(1001);
			hist_mc_rapgap[i][j]->SetFillColorAlpha(kBlue-2,0.4);
			hist_mc_rapgap[i][j]->SetMarkerColor(kBlue-2);

			hist_mc_pythia[i][j]->SetFillStyle(1001);
			hist_mc_pythia[i][j]->SetFillColorAlpha(kRed-2,0.4);
			hist_mc_pythia[i][j]->SetMarkerColor(kRed-2);

			//radiative correction to Rapgap based on Rapgap
			for(int ibin=0;ibin<hist_mc_rapgap[i][j]->GetNbinsX();ibin++){
				double value = hist_mc_rapgap[i][j]->GetBinContent(ibin+1);
				double error = hist_mc_rapgap[i][j]->GetBinError(ibin+1);
				double bincenter = hist_mc_rapgap[i][j]->GetBinCenter(ibin+1);
				double corr = hCorr_Pn_GEN_HCM[i][j]->GetBinContent( hCorr_Pn_GEN_HCM[i][j]->FindBin(bincenter) );
				if(doRadCorr_){
					hist_mc_rapgap[i][j]->SetBinContent(ibin+1, value*corr);
					hist_mc_rapgap[i][j]->SetBinError(ibin+1, error*corr);
				}
			}
			//radiative correction to Django based on Django
			for(int ibin=0;ibin<hist_mc_django[i][j]->GetNbinsX();ibin++){
				double value = hist_mc_django[i][j]->GetBinContent(ibin+1);
				double error = hist_mc_django[i][j]->GetBinError(ibin+1);
				double bincenter = hist_mc_django[i][j]->GetBinCenter(ibin+1);
				double corr = django_hCorr_Pn_GEN_HCM[i][j]->GetBinContent( django_hCorr_Pn_GEN_HCM[i][j]->FindBin(bincenter) );
				if(doRadCorr_){
					hist_mc_django[i][j]->SetBinContent(ibin+1, value*corr);
					hist_mc_django[i][j]->SetBinError(ibin+1, error*corr);
				}
			}
			
		}
	}
	

	TCanvas* c2 = new TCanvas("c2","c2",1,1,1000,1000);
	c2->Divide(4,4,0,0);

	TH1D* base2 = makeHist("base2", "", "N", "P(N)", 100,-2,33,kBlack);
	base2->GetYaxis()->SetRangeUser(0.000005, 7);
	base2->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base2,3,3.8);

	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.0);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.0);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.0);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.0);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);

	double latex_starting_eta_1[]={0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54,0.56,0.54,0.54,0.54};
	double latex_starting_y_1[]={0.48,0.49,0.48,0.48,0.49,0.53,0.53,0.53,0.51,0.525,0.56,0.53,0.52,0.51,0.5,0.55};

	std::vector<double> S_hadron_x;
	TF1 *tf3[4][4];
	int sub_panels = 1;
	for(int iy = 0; iy < 4; iy++){
		for(int iQ2 = 0; iQ2 < 4; iQ2++){

			c2->cd(sub_panels);
			
			if(sub_panels==1||sub_panels==5||sub_panels==9||sub_panels==13) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=13) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(1);
			base2->Draw("");
			// hist_unfolded_from_django_allEta[iQ2][iy]->Draw("Psame");
			hist_unfolded_from_rapgap_allEta[iQ2][iy]->Draw("PE0X0same");

			TH1D* hist_unfolded_from_rapgap_systematics = (TH1D*) hist_unfolded_from_rapgap_allEta[iQ2][iy]->Clone("hsys");
			for(int isys=0;isys<hist_unfolded_from_rapgap_systematics->GetNbinsX();isys++){
				double N = hist_unfolded_from_rapgap_systematics->GetBinCenter(isys+1);
				double Pn = hist_unfolded_from_rapgap_systematics->GetBinContent(isys+1);
				double error = 0.;
				if( N<2 ) error = 0.221*Pn;
				if( N>=2 && N<5 ) error = 0.06*Pn;
				if( N>=5 && N<15 ) error = 0.06*Pn;
				if( N>=15 ) error = 0.233*Pn;

				hist_unfolded_from_rapgap_systematics->SetBinError( isys+1, error );
			}
			gStyle->SetErrorX( 0.5 );
			hist_unfolded_from_rapgap_systematics->SetFillColorAlpha(kGray+1,0.8);
  			hist_unfolded_from_rapgap_systematics->Draw("e2same");


			hist_mc_django[iQ2][iy]->DrawNormalized("PE3 Psame");
			hist_mc_rapgap[iQ2][iy]->DrawNormalized("PE3 Psame");
			hist_mc_pythia[iQ2][iy]->DrawNormalized("PE3 Psame");

			//single NBD
			// tf3[iQ2][iy] = new TF1(Form("tf3_%d_%d",iQ2,iy),"[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])",0,30);
			//double NBD
			tf3[iQ2][iy] = new TF1(Form("tf3_%d_%d",iQ2,iy),"[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])+[3]*([4]/[4]+1)*ROOT::Math::negative_binomial_pdf(x[0],[4],[5])",0,25);
			tf3[iQ2][iy]->SetParameter(0,0.6);
			tf3[iQ2][iy]->SetParameter(1,0.08);
			tf3[iQ2][iy]->SetParameter(2,0.86);
			tf3[iQ2][iy]->SetParameter(3,1.4);
			tf3[iQ2][iy]->SetParameter(4,0.5);
			tf3[iQ2][iy]->SetParameter(5,1.9);
			tf3[iQ2][iy]->SetLineWidth(2);
			if( sub_panels == 8 ){
				// tf3[iQ2][iy]->SetParameter(0,0.6);
				// tf3[iQ2][iy]->SetParameter(1,0.08);
				// tf3[iQ2][iy]->SetParameter(2,0.8);

				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,18);
				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,18);
				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,18);

			}
			else{
				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,25);
				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,25);
				hist_unfolded_from_rapgap_allEta[iQ2][iy]->Fit(Form("tf3_%d_%d",iQ2,iy),"RME0","",4,25);
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

	TLatex* r50 = new TLatex(0.53,0.87, "0 < #eta* < 4");
	c2->cd(4);
	r50->SetNDC();
	r50->SetTextSize(18);
	r50->SetTextFont(44);
	r50->Draw("same");

	TLatex* r51 = new TLatex(0.22,0.87, "ep 27.5x920 GeV");
	c2->cd(1);
	r51->SetNDC();
	r51->SetTextSize(18);
	r51->SetTextFont(44);
	r51->Draw("same");

	TLatex* r52 = new TLatex(0.25,0.74, "H1");
	c2->cd(1);
	r52->SetNDC();
	r52->SetTextSize(20);
	r52->SetTextFont(44);
	r52->Draw("same");

	c2->cd(16);
	TLegend *w6 = new TLegend(0.05,0.22,0.5,0.5);
	w6->SetLineColor(kWhite);
	w6->SetFillColor(0);
	w6->SetTextSize(17);
	w6->SetTextFont(45);
	// w6->AddEntry(tf3[2][2],"NBD fit","L");
	w6->AddEntry(hist_mc_django[0][0], "django ","F");
	w6->AddEntry(hist_mc_rapgap[0][0], "rapgap ","F");
	w6->AddEntry(hist_mc_pythia[0][0], "pythia ","F");
	// w6->AddEntry(hist_unfolded_from_django_allEta[0][0], "unfolded data django", "P");
	w6->AddEntry(hist_unfolded_from_rapgap_allEta[0][0], "unfolded data ", "P");//rapgap
	w6->Draw("same");

	TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.1);
	gPad->SetLogx(0);
	
	double xmin_axis = 35.; double xmax_axis=250.;
	TH1D* base3 = makeHist("base3", "", "W (GeV)", "#LT N_{ch} #GT", 100,xmin_axis,xmax_axis,kBlack);
	base3->GetYaxis()->SetRangeUser(0, 8);
	base3->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base3,1.2,1.0);

	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.3);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.3);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);

	TGraphErrors* gr_data[4];
	TGraphErrors* gr_django[4];
	TGraphErrors* gr_rapgap[4];
	TGraphErrors* gr_pythia[4];

	TGraphErrors* gr_data_variance[4];
	TGraphErrors* gr_django_variance[4];
	TGraphErrors* gr_rapgap_variance[4];
	TGraphErrors* gr_pythia_variance[4];

	TGraphErrors* gr_data_entropy[4];
	TGraphErrors* gr_django_entropy[4];
	TGraphErrors* gr_rapgap_entropy[4];
	TGraphErrors* gr_pythia_entropy[4];

	double yREC_es_ave[] = {0.05464,0.1095,0.2179,0.4121};
	double W_es_ave[] = {72,101,147,205};

	for(int iQ2=0;iQ2<4;iQ2++){
		gr_data[iQ2] = new TGraphErrors();
		gr_django[iQ2] = new TGraphErrors();
		gr_rapgap[iQ2] = new TGraphErrors();
		gr_pythia[iQ2] = new TGraphErrors();

		gr_data_variance[iQ2] = new TGraphErrors();
		gr_django_variance[iQ2] = new TGraphErrors();
		gr_rapgap_variance[iQ2] = new TGraphErrors();
		gr_pythia_variance[iQ2] = new TGraphErrors();

		gr_data_entropy[iQ2] = new TGraphErrors();
		gr_django_entropy[iQ2] = new TGraphErrors();
		gr_rapgap_entropy[iQ2] = new TGraphErrors();
		gr_pythia_entropy[iQ2] = new TGraphErrors();
		
		//problem with the data is not unit bin width
		for(int iy=0;iy<4;iy++){
			//data
			TH1D* histForMultiplicity = (TH1D*) hist_unfolded_from_django_allEta[iQ2][iy]->Clone(Form("histForMultiplicity_%d_%d",iQ2,iy));
			double sum = 0.;
			double sum_pn = 0.;
			double err2 = 0.;
			double bin_to_stop = histForMultiplicity->GetBinCenter(histForMultiplicity->FindLastBinAbove(0.,1));
			
			//use the first 4 bins in data and after N > 3, use NBD extrapolations.
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){
				double N = ibin;
				double Pn = 0.;
				double Pn_err = 0;
				if( ibin == histForMultiplicity->GetBinCenter(ibin+1) ){
					N = histForMultiplicity->GetBinCenter(ibin+1);
					Pn = histForMultiplicity->GetBinContent(ibin+1);
					Pn_err = histForMultiplicity->GetBinError(ibin+1);
					sum += N*Pn;
					sum_pn += Pn;
					err2 += (N*Pn_err)*(N*Pn_err);
				}
				else{
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += N*Pn;
					sum_pn += Pn;
					err2 += (N*Pn_err)*(N*Pn_err);
				}
			}
			gr_data[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_data[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			double first_moment = sum/sum_pn;
			gr_data[iQ2]->SetPointError(iy, 0., sqrt(first_moment*first_moment*sys_Nch*sys_Nch+(sqrt(err2)/sum_pn)*(sqrt(err2)/sum_pn)) );

			//use mean to compute variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0; ibin<(int)bin_to_stop+1; ibin++){
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
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += TMath::Power((N-first_moment),2)*Pn;
					sum_pn += Pn;
					err2 += (TMath::Power((N-first_moment),2)*Pn_err)*(TMath::Power((N-first_moment),2)*Pn_err);
				}
			}
			gr_data_variance[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_data_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			double second_moment = sum/sum_pn;
			gr_data_variance[iQ2]->SetPointError(iy, 0., sqrt(second_moment*second_moment*sys_Var*sys_Var+(sqrt(err2)/sum_pn)*(sqrt(err2)/sum_pn)) );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			for(int ibin=0; ibin<(int)bin_to_stop+1; ibin++){
				double N = ibin;
				double Pn = 0.;
				double Pn_err = 0;
				if( ibin < 4 ){
					N = histForMultiplicity->GetBinCenter(ibin+1);
					Pn = histForMultiplicity->GetBinContent(ibin+1);
					Pn_err = histForMultiplicity->GetBinError(ibin+1);
					sum += -Pn*TMath::Log(Pn);
					err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
				}
				else{
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += -Pn*TMath::Log(Pn);
					err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
				}
			}
			gr_data_entropy[iQ2]->SetPoint(iy, W_es_ave[iy], sum );
			gr_data_entropy[iQ2]->SetPointError(iy, 0., sqrt((1./first_moment)*(1./first_moment)*sys_Nch*sys_Nch+err2 ));

			//django
			TH1D* histForMultiplicity_django = (TH1D*) hist_mc_django[iQ2][iy]->Clone(Form("histForMultiplicity_django_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_django->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_django[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_django[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;
			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django->GetBinWidth(ibin+1);
				double N = histForMultiplicity_django->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django->GetBinError(ibin+1);
				
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				sum_pn += Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_django_variance[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_django_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django->GetBinWidth(ibin+1);
				double N = histForMultiplicity_django->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django->GetBinError(ibin+1);
				
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_django_entropy[iQ2]->SetPoint(iy, W_es_ave[iy], sum );
			gr_django_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );


			//rapgap
			TH1D* histForMultiplicity_rapgap = (TH1D*) hist_mc_rapgap[iQ2][iy]->Clone(Form("histForMultiplicity_rapgap_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap->GetBinWidth(ibin+1);
				double N = histForMultiplicity_rapgap->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_rapgap[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_rapgap[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;

			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap->GetBinWidth(ibin+1);
				double N = histForMultiplicity_rapgap->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_rapgap_variance[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_rapgap_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap->GetBinWidth(ibin+1);
				double N = histForMultiplicity_rapgap->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap->GetBinError(ibin+1);
				
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_rapgap_entropy[iQ2]->SetPoint(iy, W_es_ave[iy], sum );
			gr_rapgap_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );

			//pythia
			TH1D* histForMultiplicity_pythia = (TH1D*) hist_mc_pythia[iQ2][iy]->Clone(Form("histForMultiplicity_pythia_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia->GetBinWidth(ibin+1);
				double N = histForMultiplicity_pythia->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_pythia[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_pythia[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;

			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia->GetBinWidth(ibin+1);
				double N = histForMultiplicity_pythia->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_pythia_variance[iQ2]->SetPoint(iy, W_es_ave[iy], sum/sum_pn );
			gr_pythia_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia->GetBinWidth(ibin+1);
				double N = histForMultiplicity_pythia->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia->GetBinError(ibin+1);
				
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_pythia_entropy[iQ2]->SetPoint(iy, W_es_ave[iy], sum );
			gr_pythia_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );
		}
	}

	base3->GetXaxis()->SetLabelOffset(999);
	base3->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis2 = new TGaxis(xmin_axis,
	                            0.0,
	                            xmax_axis,
	                            0.0,
	                            0.0,
	                            xmax_axis,
	                            510,"");
	newaxis2->SetLabelOffset(0.01);
	newaxis2->SetLabelFont(42);
	newaxis2->SetLabelSize( 0.035 );

	base3->Draw();
	newaxis2->Draw("SS");
	gPad->Update();

	TLine* l1[10];
    l1[0] = new TLine(W_es_ave[0],7.9, W_es_ave[0], 8);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(W_es_ave[1],7.9, W_es_ave[1], 8);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(W_es_ave[2],7.9, W_es_ave[2], 8);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[2] = new TLine(W_es_ave[3],7.9, W_es_ave[3], 8);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

	gr_data[0]->SetMarkerStyle(20);
	gr_data[0]->SetMarkerSize(1.4);
	gr_data[0]->SetMarkerColor(kBlack);

	gr_data[1]->SetMarkerStyle(30);
	gr_data[1]->SetMarkerSize(1.4);
	gr_data[1]->SetMarkerColor(kBlue);

	gr_data[2]->SetMarkerStyle(24);
	gr_data[2]->SetMarkerSize(1.4);
	gr_data[2]->SetMarkerColor(kRed);

	gr_data[3]->SetMarkerStyle(29);
	gr_data[3]->SetMarkerSize(1.8);
	gr_data[3]->SetMarkerColor(kGreen-2);

	gr_data[0]->Draw("Psame");
	gr_data[1]->Draw("Psame");
	gr_data[2]->Draw("Psame");
	gr_data[3]->Draw("Psame");

	gr_django[0]->SetLineStyle(2);
	gr_django[0]->SetLineWidth(4);
	gr_django[0]->SetMarkerSize(1.4);
	gr_django[0]->SetLineColor(kBlack);

	gr_django[1]->SetLineStyle(2);
	gr_django[1]->SetLineWidth(4);
	gr_django[1]->SetMarkerSize(1.4);
	gr_django[1]->SetLineColor(kBlue);

	gr_django[2]->SetLineStyle(2);
	gr_django[2]->SetLineWidth(4);
	gr_django[2]->SetMarkerSize(1.4);
	gr_django[2]->SetLineColor(kRed);

	gr_django[3]->SetLineStyle(2);
	gr_django[3]->SetLineWidth(4);
	gr_django[3]->SetMarkerSize(1.6);
	gr_django[3]->SetLineColor(kGreen-2);

	gr_rapgap[0]->SetLineStyle(2);
	gr_rapgap[0]->SetLineWidth(4);
	gr_rapgap[0]->SetMarkerSize(1.4);
	gr_rapgap[0]->SetLineColor(kBlack);

	gr_rapgap[1]->SetLineStyle(3);
	gr_rapgap[1]->SetLineWidth(4);
	gr_rapgap[1]->SetMarkerSize(1.4);
	gr_rapgap[1]->SetLineColor(kBlue);

	gr_rapgap[2]->SetLineStyle(4);
	gr_rapgap[2]->SetLineWidth(4);
	gr_rapgap[2]->SetMarkerSize(1.4);
	gr_rapgap[2]->SetLineColor(kRed);

	gr_rapgap[3]->SetLineStyle(5);
	gr_rapgap[3]->SetLineWidth(4);
	gr_rapgap[3]->SetMarkerSize(1.6);
	gr_rapgap[3]->SetLineColor(kGreen-2);

	gr_pythia[0]->SetLineStyle(2);
	gr_pythia[0]->SetLineWidth(4);
	gr_pythia[0]->SetMarkerSize(1.4);
	gr_pythia[0]->SetLineColor(kBlack);

	gr_pythia[1]->SetLineStyle(3);
	gr_pythia[1]->SetLineWidth(4);
	gr_pythia[1]->SetMarkerSize(1.4);
	gr_pythia[1]->SetLineColor(kBlue);

	gr_pythia[2]->SetLineStyle(4);
	gr_pythia[2]->SetLineWidth(4);
	gr_pythia[2]->SetMarkerSize(1.4);
	gr_pythia[2]->SetLineColor(kRed);

	gr_pythia[3]->SetLineStyle(5);
	gr_pythia[3]->SetLineWidth(4);
	gr_pythia[3]->SetMarkerSize(1.6);
	gr_pythia[3]->SetLineColor(kGreen-2);

	if( pickMCEG_ == 0 ){
		gr_rapgap[0]->Draw("L same");
		gr_rapgap[1]->Draw("L same");
		gr_rapgap[2]->Draw("L same");
		gr_rapgap[3]->Draw("L same");
	}
	else if( pickMCEG_ == 1 ){
		gr_django[0]->Draw("Lsame");
		gr_django[1]->Draw("Lsame");
		gr_django[2]->Draw("Lsame");
		gr_django[3]->Draw("Lsame");
	}
	else if( pickMCEG_ == 2 ){
		gr_pythia[0]->Draw("Lsame");
		gr_pythia[1]->Draw("Lsame");
		gr_pythia[2]->Draw("Lsame");
		gr_pythia[3]->Draw("Lsame");
	}


	TLatex* r44 = new TLatex(0.8,0.84, "H1");
	r44->SetNDC();
	r44->SetTextSize(0.04);


	TLatex* r48 = new TLatex(0.13, 0.84, "ep 27.5#times920 GeV");
	r48->SetNDC();
	r48->SetTextSize(23);
	r48->SetTextFont(43);
	r48->SetTextColor(kBlack);


	TLatex* r49 = new TLatex(0.13, 0.79, "0 < #eta* < 4.0");
	r49->SetNDC();
	r49->SetTextSize(23);
	r49->SetTextFont(43);
	r49->SetTextColor(kBlack);

	TLatex* r441 = new TLatex(0.45,0.95, "#LT y #GT");
	r441->SetNDC();
	r441->SetTextSize(20);
	r441->SetTextFont(43);
	r441->SetTextColor(kBlack);

	r44->Draw("same");
	r48->Draw("same");
	r49->Draw("same");
	r441->Draw("same");

	TLatex* r630 = new TLatex(0.48,0.34, "Q^{2} (data)");
	r630->SetNDC();
	r630->SetTextSize(0.03);
	r630->Draw("same");

	TLegend *w63 = new TLegend(0.48,0.18,0.65,0.32);
	w63->SetLineColor(kWhite);
	w63->SetFillColor(0);
	w63->SetTextSize(14);
	w63->SetTextFont(45);
	w63->SetTextColor(kBlack);
	w63->AddEntry(gr_data[0], "(5,10) ","P");
	w63->AddEntry(gr_data[1], "(10,20) ","P");
	w63->AddEntry(gr_data[2], "(20,40)", "P");
	w63->AddEntry(gr_data[3], "(40,100) ", "P");
	w63->Draw("same");

	TString mc_generator = "Q^{2} (rapgap)";
	if( pickMCEG_ == 1 ) mc_generator = "Q^{2} (django)";
	if( pickMCEG_ == 2 ) mc_generator = "Q^{2} (pythia)";
	TLatex* r631 = new TLatex(0.65,0.34, mc_generator);
	r631->SetNDC();
	r631->SetTextSize(0.03);
	r631->Draw("same");

	TLegend *w631 = new TLegend(0.65,0.18,0.87,0.32);
	w631->SetLineColor(kWhite);
	w631->SetFillColor(0);
	w631->SetTextSize(14);
	w631->SetTextFont(45);
	w631->SetTextColor(kBlack);
	w631->AddEntry(gr_rapgap[0], "(5,10) ","L");
	w631->AddEntry(gr_rapgap[1], "(10,20) ","L");
	w631->AddEntry(gr_rapgap[2], "(20,40)", "L");
	w631->AddEntry(gr_rapgap[3], "(40,100) ", "L");
	w631->Draw("same");

	// yREC_es_ave[] = {0.05464,0.1095,0.2179,0.4121};
	TLatex* W[10];
  	W[0] = new TLatex(0.22, 0.91, "0.05");
    W[1] = new TLatex(0.32, 0.91, "0.11");
    W[2] = new TLatex(0.49, 0.91, "0.22");
    W[3] = new TLatex(0.71, 0.91, "0.41");
    for(int i=0;i<4;i++){
    	W[i]->SetNDC();
		W[i]->SetTextSize(17);
		W[i]->SetTextFont(43);
	}
    W[0]->Draw("same");
    W[1]->Draw("same");
    W[2]->Draw("same");
    W[3]->Draw("same");

	TCanvas* c3_2 = new TCanvas("c3_2","c3_2",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.1);
	gPad->SetLogx(0); 
	TH1D* base3_2 = (TH1D*) base3->Clone("base3_2");
	base3_2->GetYaxis()->SetRangeUser(0,18);
	base3_2->GetYaxis()->SetTitle("Var(N_{ch})");
	base3_2->Draw();
    newaxis2->Draw("SS");
	gPad->Update();
	r44->Draw("same");
	r48->Draw("same");
	r49->Draw("same");
	r441->Draw("same");

	W[0]->Draw("same");
    W[1]->Draw("same");
    W[2]->Draw("same");
    W[3]->Draw("same");

    r630->Draw("same");
	w63->Draw("same");
    r631->Draw("same");
    w631->Draw("same");

    l1[0] = new TLine(W_es_ave[0],17.75, W_es_ave[0], 18);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(W_es_ave[1],17.75, W_es_ave[1], 18);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(W_es_ave[2],17.75, W_es_ave[2], 18);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[2] = new TLine(W_es_ave[3],17.75, W_es_ave[3], 18);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    gr_data_variance[0]->SetMarkerStyle(20);
	gr_data_variance[0]->SetMarkerSize(1.4);
	gr_data_variance[0]->SetMarkerColor(kBlack);

	gr_data_variance[1]->SetMarkerStyle(30);
	gr_data_variance[1]->SetMarkerSize(1.4);
	gr_data_variance[1]->SetMarkerColor(kBlue);

	gr_data_variance[2]->SetMarkerStyle(24);
	gr_data_variance[2]->SetMarkerSize(1.4);
	gr_data_variance[2]->SetMarkerColor(kRed);

	gr_data_variance[3]->SetMarkerStyle(29);
	gr_data_variance[3]->SetMarkerSize(1.8);
	gr_data_variance[3]->SetMarkerColor(kGreen-2);

	gr_data_variance[0]->Draw("Psame");
	gr_data_variance[1]->Draw("Psame");
	gr_data_variance[2]->Draw("Psame");
	gr_data_variance[3]->Draw("Psame");

	gr_django_variance[0]->SetLineStyle(2);
	gr_django_variance[0]->SetLineWidth(4);
	gr_django_variance[0]->SetMarkerSize(1.4);
	gr_django_variance[0]->SetLineColor(kBlack);

	gr_django_variance[1]->SetLineStyle(2);
	gr_django_variance[1]->SetLineWidth(4);
	gr_django_variance[1]->SetMarkerSize(1.4);
	gr_django_variance[1]->SetLineColor(kBlue);

	gr_django_variance[2]->SetLineStyle(2);
	gr_django_variance[2]->SetLineWidth(4);
	gr_django_variance[2]->SetMarkerSize(1.4);
	gr_django_variance[2]->SetLineColor(kRed);

	gr_django_variance[3]->SetLineStyle(2);
	gr_django_variance[3]->SetLineWidth(4);
	gr_django_variance[3]->SetMarkerSize(1.6);
	gr_django_variance[3]->SetLineColor(kGreen-2);

	gr_rapgap_variance[0]->SetLineStyle(2);
	gr_rapgap_variance[0]->SetLineWidth(4);
	gr_rapgap_variance[0]->SetMarkerSize(1.4);
	gr_rapgap_variance[0]->SetLineColor(kBlack);

	gr_rapgap_variance[1]->SetLineStyle(3);
	gr_rapgap_variance[1]->SetLineWidth(4);
	gr_rapgap_variance[1]->SetMarkerSize(1.4);
	gr_rapgap_variance[1]->SetLineColor(kBlue);

	gr_rapgap_variance[2]->SetLineStyle(4);
	gr_rapgap_variance[2]->SetLineWidth(4);
	gr_rapgap_variance[2]->SetMarkerSize(1.4);
	gr_rapgap_variance[2]->SetLineColor(kRed);

	gr_rapgap_variance[3]->SetLineStyle(5);
	gr_rapgap_variance[3]->SetLineWidth(4);
	gr_rapgap_variance[3]->SetMarkerSize(1.6);
	gr_rapgap_variance[3]->SetLineColor(kGreen-2);

	gr_pythia_variance[0]->SetLineStyle(2);
	gr_pythia_variance[0]->SetLineWidth(4);
	gr_pythia_variance[0]->SetMarkerSize(1.4);
	gr_pythia_variance[0]->SetLineColor(kBlack);

	gr_pythia_variance[1]->SetLineStyle(3);
	gr_pythia_variance[1]->SetLineWidth(4);
	gr_pythia_variance[1]->SetMarkerSize(1.4);
	gr_pythia_variance[1]->SetLineColor(kBlue);

	gr_pythia_variance[2]->SetLineStyle(4);
	gr_pythia_variance[2]->SetLineWidth(4);
	gr_pythia_variance[2]->SetMarkerSize(1.4);
	gr_pythia_variance[2]->SetLineColor(kRed);

	gr_pythia_variance[3]->SetLineStyle(5);
	gr_pythia_variance[3]->SetLineWidth(4);
	gr_pythia_variance[3]->SetMarkerSize(1.6);
	gr_pythia_variance[3]->SetLineColor(kGreen-2);

	if( pickMCEG_ == 0 ){
		gr_rapgap_variance[0]->Draw("L same");
		gr_rapgap_variance[1]->Draw("L same");
		gr_rapgap_variance[2]->Draw("L same");
		gr_rapgap_variance[3]->Draw("L same");
	}
	else if( pickMCEG_ == 1 ){
		gr_django_variance[0]->Draw("Lsame");
		gr_django_variance[1]->Draw("Lsame");
		gr_django_variance[2]->Draw("Lsame");
		gr_django_variance[3]->Draw("Lsame");
	}
	else if( pickMCEG_ == 2 ){
		gr_pythia_variance[0]->Draw("Lsame");
		gr_pythia_variance[1]->Draw("Lsame");
		gr_pythia_variance[2]->Draw("Lsame");
		gr_pythia_variance[3]->Draw("Lsame");
	}

    
    TCanvas* c4 = new TCanvas("c4","c4",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.1);
	gPad->SetLogx(1);

	TH1D* base4 = makeHist("base4", "", "#LT x #GT", "#LT N_{ch} #GT", 10000,0.00001,0.02,kBlack);
	base4->GetYaxis()->SetRangeUser(0, 8.7);
	base4->GetXaxis()->SetRangeUser(0.0001,0.02);
	base4->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base4,1.2,1.0);

	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.3);
	base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.3);
	base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.3);
	base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.3);
	base4->GetXaxis()->SetNdivisions(4,6,0);
	base4->GetYaxis()->SetNdivisions(4,6,0);

	base4->Draw();

	TGraphErrors* gr_logx_data[4];
	TGraphErrors* gr_logx_django[4];
	TGraphErrors* gr_logx_rapgap[4];
	TGraphErrors* gr_logx_pythia[4];

	TGraphErrors* gr_logx_data_variance[4];
	TGraphErrors* gr_logx_django_variance[4];
	TGraphErrors* gr_logx_rapgap_variance[4];
	TGraphErrors* gr_logx_pythia_variance[4];

	TGraphErrors* gr_logx_data_entropy[4];
	TGraphErrors* gr_logx_django_entropy[4];
	TGraphErrors* gr_logx_rapgap_entropy[4];
	TGraphErrors* gr_logx_pythia_entropy[4];

	double x_ave_Q2_5_10[]={0.0014,0.0007,0.00033,0.00017};
	double x_ave_Q2_10_20[]={0.0023,0.0013,0.00066,0.00034};
	double x_ave_Q2_20_40[]={0.0052,0.0026,0.0013,0.00068};
	double x_ave_Q2_40_100[]={0.011,0.0056,0.0028,0.0013};
	double x_ave_Q2[4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			if( i==0) x_ave_Q2[i][j] = x_ave_Q2_5_10[j];
			if( i==1) x_ave_Q2[i][j] = x_ave_Q2_10_20[j];
			if( i==2) x_ave_Q2[i][j] = x_ave_Q2_20_40[j];
			if( i==3) x_ave_Q2[i][j] = x_ave_Q2_40_100[j];
		}
	}

	for(int iQ2=0;iQ2<4;iQ2++){
		gr_logx_data[iQ2] = new TGraphErrors();
		gr_logx_django[iQ2] = new TGraphErrors();
		gr_logx_rapgap[iQ2] = new TGraphErrors();
		gr_logx_pythia[iQ2] = new TGraphErrors();

		gr_logx_data_variance[iQ2] = new TGraphErrors();
		gr_logx_django_variance[iQ2] = new TGraphErrors();
		gr_logx_rapgap_variance[iQ2] = new TGraphErrors();
		gr_logx_pythia_variance[iQ2] = new TGraphErrors();

		gr_logx_data_entropy[iQ2] = new TGraphErrors();
		gr_logx_django_entropy[iQ2] = new TGraphErrors();
		gr_logx_rapgap_entropy[iQ2] = new TGraphErrors();
		gr_logx_pythia_entropy[iQ2] = new TGraphErrors();
		
		//problem with the data is not unit bin width
		for(int iy=0;iy<4;iy++){
			//data
			TH1D* histForMultiplicity = (TH1D*) hist_unfolded_from_django_allEta[iQ2][iy]->Clone(Form("histForMultiplicity_%d_%d",iQ2,iy));
			double sum = 0.;
			double sum_pn = 0.;
			double err2 = 0.;
			double bin_to_stop = histForMultiplicity->GetBinCenter(histForMultiplicity->FindLastBinAbove(0.,1));
			
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
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += N*Pn;
					sum_pn += Pn;
					err2 += (N*Pn_err)*(N*Pn_err);
				}
			}
			gr_logx_data[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_data[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			double first_moment = sum/sum_pn;
			gr_logx_data[iQ2]->SetPointError(iy, 0., sqrt(first_moment*first_moment*sys_Nch*sys_Nch+(sqrt(err2)/sum_pn)*(sqrt(err2)/sum_pn)) );

			//use mean to compute variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0; ibin<(int)bin_to_stop+1; ibin++){
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
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					sum += TMath::Power((N-first_moment),2)*Pn;
					sum_pn += Pn;
					err2 += (TMath::Power((N-first_moment),2)*Pn_err)*(TMath::Power((N-first_moment),2)*Pn_err);
				}
			}
			gr_logx_data_variance[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_data_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			double second_moment = sum/sum_pn;
			gr_logx_data_variance[iQ2]->SetPointError(iy, 0., sqrt(second_moment*second_moment*sys_Var*sys_Var+(sqrt(err2)/sum_pn)*(sqrt(err2)/sum_pn)) );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			for(int ibin=0; ibin<(int)bin_to_stop+1; ibin++){
				double N = ibin;
				double Pn = 0.;
				double Pn_err = 0;
				if( ibin < 4 ){
					N = histForMultiplicity->GetBinCenter(ibin+1);
					Pn = histForMultiplicity->GetBinContent(ibin+1);
					Pn_err = histForMultiplicity->GetBinError(ibin+1);
					sum += -Pn*TMath::Log(Pn);
					err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
				}
				else{
					Pn = tf3[iQ2][iy]->Eval(N);
					int nbin = histForMultiplicity->FindBin(N);
					Pn_err = histForMultiplicity->GetBinError(nbin);
					if( Pn < 1e-5 ) continue;
					sum += -Pn*TMath::Log(Pn);
					err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
				}
			}
			gr_logx_data_entropy[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum );
			gr_logx_data_entropy[iQ2]->SetPointError(iy, 0., sqrt(sys_Nch*sys_Nch+err2 ));


			//django
			TH1D* histForMultiplicity_django = (TH1D*) hist_mc_django[iQ2][iy]->Clone(Form("histForMultiplicity_django_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_django->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_logx_django[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_django[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;
			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_django->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django->GetBinError(ibin+1);
				
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				sum_pn += Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_logx_django_variance[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_django_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			TH1D* histForMultiplicity_django_entropy = (TH1D*) histForMultiplicity_django->Clone("histForMultiplicity_django_entropy");
			histForMultiplicity_django_entropy->Scale(1./(histForMultiplicity_django_entropy->Integral()));
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_django_entropy->GetBinWidth(ibin+1);
				double N = histForMultiplicity_django_entropy->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_django_entropy->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_django_entropy->GetBinError(ibin+1);
				
				if(Pn < 1e-5) continue;
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_logx_django_entropy[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum );
			gr_logx_django_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );

			//rapgap
			TH1D* histForMultiplicity_rapgap = (TH1D*) hist_mc_rapgap[iQ2][iy]->Clone(Form("histForMultiplicity_rapgap_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_rapgap->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_logx_rapgap[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_rapgap[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;

			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_rapgap->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_logx_rapgap_variance[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_rapgap_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			TH1D* histForMultiplicity_rapgap_entropy = (TH1D*) histForMultiplicity_rapgap->Clone("histForMultiplicity_rapgap_entropy");
			histForMultiplicity_rapgap_entropy->Scale(1./(histForMultiplicity_rapgap_entropy->Integral()));
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_rapgap_entropy->GetBinWidth(ibin+1);
				double N = histForMultiplicity_rapgap_entropy->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_rapgap_entropy->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_rapgap_entropy->GetBinError(ibin+1);
				if( Pn < 1e-5) continue;
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_logx_rapgap_entropy[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum );
			gr_logx_rapgap_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );

			//pythia
			TH1D* histForMultiplicity_pythia = (TH1D*) hist_mc_pythia[iQ2][iy]->Clone(Form("histForMultiplicity_pythia_%d_%d",iQ2,iy));
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_pythia->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += N*Pn*weight;
				err2 += (N*Pn_err*weight)*(N*Pn_err*weight);
			}
			gr_logx_pythia[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_pythia[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );
			first_moment = sum/sum_pn;

			//from mean to variance
			sum = 0.;
			sum_pn = 0.;
			err2 = 0.;
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia->GetBinWidth(ibin+1);
				//shift bin center to integer.
				double N = histForMultiplicity_pythia->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia->GetBinError(ibin+1);

				sum_pn += Pn*weight;
				sum += TMath::Power((N-first_moment),2)*Pn*weight;
				err2 += (TMath::Power((N-first_moment),2)*Pn_err*weight)*(TMath::Power((N-first_moment),2)*Pn_err*weight);
			}
			gr_logx_pythia_variance[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum/sum_pn );
			gr_logx_pythia_variance[iQ2]->SetPointError(iy, 0., sqrt(err2)/sum_pn );

			//use Pn to compute entropy
			sum = 0.;
			err2 = 0.;
			TH1D* histForMultiplicity_pythia_entropy = (TH1D*) histForMultiplicity_pythia->Clone("histForMultiplicity_pythia_entropy");
			histForMultiplicity_pythia_entropy->Scale(1./(histForMultiplicity_pythia_entropy->Integral()));
			for(int ibin=0;ibin<(int)bin_to_stop+1;ibin++){

				double weight = histForMultiplicity_pythia_entropy->GetBinWidth(ibin+1);
				double N = histForMultiplicity_pythia_entropy->GetBinCenter(ibin+1);
				double Pn = histForMultiplicity_pythia_entropy->GetBinContent(ibin+1);
				double Pn_err = histForMultiplicity_pythia_entropy->GetBinError(ibin+1);
				
				if(Pn < 1e-5) continue;
				sum += -Pn*TMath::Log(Pn);
				err2 += (Pn_err*Pn_err)*((1.+TMath::Log(Pn))*(1.+TMath::Log(Pn)) );
			}
			gr_logx_pythia_entropy[iQ2]->SetPoint(iy, x_ave_Q2[iQ2][iy], sum );
			gr_logx_pythia_entropy[iQ2]->SetPointError(iy, 0., sqrt(err2) );


		}
	}

	gr_logx_data[0]->SetMarkerStyle(20);
	gr_logx_data[0]->SetMarkerSize(1.4);
	gr_logx_data[0]->SetMarkerColor(kBlack);

	gr_logx_data[1]->SetMarkerStyle(30);
	gr_logx_data[1]->SetMarkerSize(1.4);
	gr_logx_data[1]->SetMarkerColor(kBlue);

	gr_logx_data[2]->SetMarkerStyle(24);
	gr_logx_data[2]->SetMarkerSize(1.4);
	gr_logx_data[2]->SetMarkerColor(kRed);

	gr_logx_data[3]->SetMarkerStyle(29);
	gr_logx_data[3]->SetMarkerSize(1.8);
	gr_logx_data[3]->SetMarkerColor(kGreen-2);

	gr_logx_data[0]->Draw("Psame");
	gr_logx_data[1]->Draw("Psame");
	gr_logx_data[2]->Draw("Psame");
	gr_logx_data[3]->Draw("Psame");

	gr_logx_django[0]->SetLineStyle(2);
	gr_logx_django[0]->SetLineWidth(4);
	gr_logx_django[0]->SetMarkerSize(1.4);
	gr_logx_django[0]->SetLineColor(kBlack);

	gr_logx_django[1]->SetLineStyle(2);
	gr_logx_django[1]->SetLineWidth(4);
	gr_logx_django[1]->SetMarkerSize(1.4);
	gr_logx_django[1]->SetLineColor(kBlue);

	gr_logx_django[2]->SetLineStyle(2);
	gr_logx_django[2]->SetLineWidth(4);
	gr_logx_django[2]->SetMarkerSize(1.4);
	gr_logx_django[2]->SetLineColor(kRed);

	gr_logx_django[3]->SetLineStyle(2);
	gr_logx_django[3]->SetLineWidth(4);
	gr_logx_django[3]->SetMarkerSize(1.6);
	gr_logx_django[3]->SetLineColor(kGreen-2);

	gr_logx_rapgap[0]->SetLineStyle(2);
	gr_logx_rapgap[0]->SetLineWidth(4);
	gr_logx_rapgap[0]->SetMarkerSize(1.4);
	gr_logx_rapgap[0]->SetLineColor(kBlack);

	gr_logx_rapgap[1]->SetLineStyle(3);
	gr_logx_rapgap[1]->SetLineWidth(4);
	gr_logx_rapgap[1]->SetMarkerSize(1.4);
	gr_logx_rapgap[1]->SetLineColor(kBlue);

	gr_logx_rapgap[2]->SetLineStyle(4);
	gr_logx_rapgap[2]->SetLineWidth(4);
	gr_logx_rapgap[2]->SetMarkerSize(1.4);
	gr_logx_rapgap[2]->SetLineColor(kRed);

	gr_logx_rapgap[3]->SetLineStyle(5);
	gr_logx_rapgap[3]->SetLineWidth(4);
	gr_logx_rapgap[3]->SetMarkerSize(1.6);
	gr_logx_rapgap[3]->SetLineColor(kGreen-2);

	gr_logx_pythia[0]->SetLineStyle(2);
	gr_logx_pythia[0]->SetLineWidth(4);
	gr_logx_pythia[0]->SetMarkerSize(1.4);
	gr_logx_pythia[0]->SetLineColor(kBlack);

	gr_logx_pythia[1]->SetLineStyle(3);
	gr_logx_pythia[1]->SetLineWidth(4);
	gr_logx_pythia[1]->SetMarkerSize(1.4);
	gr_logx_pythia[1]->SetLineColor(kBlue);

	gr_logx_pythia[2]->SetLineStyle(4);
	gr_logx_pythia[2]->SetLineWidth(4);
	gr_logx_pythia[2]->SetMarkerSize(1.4);
	gr_logx_pythia[2]->SetLineColor(kRed);

	gr_logx_pythia[3]->SetLineStyle(5);
	gr_logx_pythia[3]->SetLineWidth(4);
	gr_logx_pythia[3]->SetMarkerSize(1.6);
	gr_logx_pythia[3]->SetLineColor(kGreen-2);

	if( pickMCEG_ == 0 ){
		gr_logx_rapgap[0]->Draw("L same");
		gr_logx_rapgap[1]->Draw("L same");
		gr_logx_rapgap[2]->Draw("L same");
		gr_logx_rapgap[3]->Draw("L same");
	}
	else if( pickMCEG_ == 1 ){
		gr_logx_django[0]->Draw("Lsame");
		gr_logx_django[1]->Draw("Lsame");
		gr_logx_django[2]->Draw("Lsame");
		gr_logx_django[3]->Draw("Lsame");
	}
	else if( pickMCEG_ == 2 ){
		gr_logx_pythia[0]->Draw("Lsame");
		gr_logx_pythia[1]->Draw("Lsame");
		gr_logx_pythia[2]->Draw("Lsame");
		gr_logx_pythia[3]->Draw("Lsame");
	}


	r44->Draw("same");
	r48->Draw("same");
	r49->Draw("same");
	// r441->Draw("same");

	r630->Draw("same");
	w63->Draw("same");

	r631->Draw("same");
	w631->Draw("same");

	TCanvas* c4_2 = new TCanvas("c4_2","c4_2",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.1);
	gPad->SetLogx(1);

	TH1D* base4_2 = (TH1D*) base4->Clone("base4_2");
	base4_2->GetYaxis()->SetRangeUser(-13,30);
	base4_2->GetYaxis()->SetTitle("Var(N_{ch})");
	base4_2->Draw();
	r44->Draw("same");
	r48->Draw("same");
	r49->Draw("same");
	// r441->Draw("same");

	r630->Draw("same");
	w63->Draw("same");

	r631->Draw("same");
	w631->Draw("same");

	gr_logx_data_variance[0]->SetMarkerStyle(20);
	gr_logx_data_variance[0]->SetMarkerSize(1.4);
	gr_logx_data_variance[0]->SetMarkerColor(kBlack);

	gr_logx_data_variance[1]->SetMarkerStyle(30);
	gr_logx_data_variance[1]->SetMarkerSize(1.4);
	gr_logx_data_variance[1]->SetMarkerColor(kBlue);

	gr_logx_data_variance[2]->SetMarkerStyle(24);
	gr_logx_data_variance[2]->SetMarkerSize(1.4);
	gr_logx_data_variance[2]->SetMarkerColor(kRed);

	gr_logx_data_variance[3]->SetMarkerStyle(29);
	gr_logx_data_variance[3]->SetMarkerSize(1.8);
	gr_logx_data_variance[3]->SetMarkerColor(kGreen-2);

	gr_logx_data_variance[0]->Draw("Psame");
	gr_logx_data_variance[1]->Draw("Psame");
	gr_logx_data_variance[2]->Draw("Psame");
	gr_logx_data_variance[3]->Draw("Psame");

	gr_logx_django_variance[0]->SetLineStyle(2);
	gr_logx_django_variance[0]->SetLineWidth(4);
	gr_logx_django_variance[0]->SetMarkerSize(1.4);
	gr_logx_django_variance[0]->SetLineColor(kBlack);

	gr_logx_django_variance[1]->SetLineStyle(2);
	gr_logx_django_variance[1]->SetLineWidth(4);
	gr_logx_django_variance[1]->SetMarkerSize(1.4);
	gr_logx_django_variance[1]->SetLineColor(kBlue);

	gr_logx_django_variance[2]->SetLineStyle(2);
	gr_logx_django_variance[2]->SetLineWidth(4);
	gr_logx_django_variance[2]->SetMarkerSize(1.4);
	gr_logx_django_variance[2]->SetLineColor(kRed);

	gr_logx_django_variance[3]->SetLineStyle(2);
	gr_logx_django_variance[3]->SetLineWidth(4);
	gr_logx_django_variance[3]->SetMarkerSize(1.6);
	gr_logx_django_variance[3]->SetLineColor(kGreen-2);

	gr_logx_rapgap_variance[0]->SetLineStyle(2);
	gr_logx_rapgap_variance[0]->SetLineWidth(4);
	gr_logx_rapgap_variance[0]->SetMarkerSize(1.4);
	gr_logx_rapgap_variance[0]->SetLineColor(kBlack);

	gr_logx_rapgap_variance[1]->SetLineStyle(3);
	gr_logx_rapgap_variance[1]->SetLineWidth(4);
	gr_logx_rapgap_variance[1]->SetMarkerSize(1.4);
	gr_logx_rapgap_variance[1]->SetLineColor(kBlue);

	gr_logx_rapgap_variance[2]->SetLineStyle(4);
	gr_logx_rapgap_variance[2]->SetLineWidth(4);
	gr_logx_rapgap_variance[2]->SetMarkerSize(1.4);
	gr_logx_rapgap_variance[2]->SetLineColor(kRed);

	gr_logx_rapgap_variance[3]->SetLineStyle(5);
	gr_logx_rapgap_variance[3]->SetLineWidth(4);
	gr_logx_rapgap_variance[3]->SetMarkerSize(1.6);
	gr_logx_rapgap_variance[3]->SetLineColor(kGreen-2);

	gr_logx_pythia_variance[0]->SetLineStyle(2);
	gr_logx_pythia_variance[0]->SetLineWidth(4);
	gr_logx_pythia_variance[0]->SetMarkerSize(1.4);
	gr_logx_pythia_variance[0]->SetLineColor(kBlack);

	gr_logx_pythia_variance[1]->SetLineStyle(3);
	gr_logx_pythia_variance[1]->SetLineWidth(4);
	gr_logx_pythia_variance[1]->SetMarkerSize(1.4);
	gr_logx_pythia_variance[1]->SetLineColor(kBlue);

	gr_logx_pythia_variance[2]->SetLineStyle(4);
	gr_logx_pythia_variance[2]->SetLineWidth(4);
	gr_logx_pythia_variance[2]->SetMarkerSize(1.4);
	gr_logx_pythia_variance[2]->SetLineColor(kRed);

	gr_logx_pythia_variance[3]->SetLineStyle(5);
	gr_logx_pythia_variance[3]->SetLineWidth(4);
	gr_logx_pythia_variance[3]->SetMarkerSize(1.6);
	gr_logx_pythia_variance[3]->SetLineColor(kGreen-2);

	if( pickMCEG_ == 0 ){
		gr_logx_rapgap_variance[0]->Draw("L same");
		gr_logx_rapgap_variance[1]->Draw("L same");
		gr_logx_rapgap_variance[2]->Draw("L same");
		gr_logx_rapgap_variance[3]->Draw("L same");
	}
	else if( pickMCEG_ == 1 ){
		gr_logx_django_variance[0]->Draw("Lsame");
		gr_logx_django_variance[1]->Draw("Lsame");
		gr_logx_django_variance[2]->Draw("Lsame");
		gr_logx_django_variance[3]->Draw("Lsame");
	}
	else if( pickMCEG_ == 2 ){
		gr_logx_pythia_variance[0]->Draw("Lsame");
		gr_logx_pythia_variance[1]->Draw("Lsame");
		gr_logx_pythia_variance[2]->Draw("Lsame");
		gr_logx_pythia_variance[3]->Draw("Lsame");
	}

	TCanvas* c4_3 = new TCanvas("c4_3","c4_3",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.1);
	gPad->SetLogx(1);

	TH1D* base4_3 = (TH1D*) base4->Clone("base4_3");
	base4_3->GetYaxis()->SetRangeUser(0,5.3);
	base4_3->GetXaxis()->SetRangeUser(0.00001,0.02);
	base4_3->GetYaxis()->SetTitle("S_{hadron}");
	base4_3->Draw();
	r44->Draw("same");
	r48->Draw("same");
	r49->Draw("same");
	// r441->Draw("same");

	TLatex* r6301 = new TLatex(0.17,0.34, "Q^{2} (data)");
	r6301->SetNDC();
	r6301->SetTextSize(0.03);
	r6301->Draw("same");

	TLegend *w6312 = new TLegend(0.17,0.18,0.35,0.32);
	w6312->SetLineColor(kWhite);
	w6312->SetFillColor(0);
	w6312->SetTextSize(14);
	w6312->SetTextFont(45);
	w6312->SetTextColor(kBlack);
	w6312->AddEntry(gr_data[0], "(5,10) ","P");
	w6312->AddEntry(gr_data[1], "(10,20) ","P");
	w6312->AddEntry(gr_data[2], "(20,40)", "P");
	w6312->AddEntry(gr_data[3], "(40,100) ", "P");
	w6312->Draw("same");


	TLatex* r6311 = new TLatex(0.32,0.34, mc_generator);
	r6311->SetNDC();
	r6311->SetTextSize(0.03);
	r6311->Draw("same");

	TLegend *w6311 = new TLegend(0.33,0.18,0.53,0.32);
	w6311->SetLineColor(kWhite);
	w6311->SetFillColor(0);
	w6311->SetTextSize(14);
	w6311->SetTextFont(45);
	w6311->SetTextColor(kBlack);
	w6311->AddEntry(gr_rapgap[0], "(5,10) ","L");
	w6311->AddEntry(gr_rapgap[1], "(10,20) ","L");
	w6311->AddEntry(gr_rapgap[2], "(20,40)", "L");
	w6311->AddEntry(gr_rapgap[3], "(40,100) ", "L");
	w6311->Draw("same");

	gr_logx_data_entropy[0]->SetMarkerStyle(20);
	gr_logx_data_entropy[0]->SetMarkerSize(1.4);
	gr_logx_data_entropy[0]->SetMarkerColor(kBlack);

	gr_logx_data_entropy[1]->SetMarkerStyle(30);
	gr_logx_data_entropy[1]->SetMarkerSize(1.4);
	gr_logx_data_entropy[1]->SetMarkerColor(kBlue);

	gr_logx_data_entropy[2]->SetMarkerStyle(24);
	gr_logx_data_entropy[2]->SetMarkerSize(1.4);
	gr_logx_data_entropy[2]->SetMarkerColor(kRed);

	gr_logx_data_entropy[3]->SetMarkerStyle(29);
	gr_logx_data_entropy[3]->SetMarkerSize(1.8);
	gr_logx_data_entropy[3]->SetMarkerColor(kGreen-2);

	gr_logx_data_entropy[0]->Draw("Psame");
	gr_logx_data_entropy[1]->Draw("Psame");
	gr_logx_data_entropy[2]->Draw("Psame");
	gr_logx_data_entropy[3]->Draw("Psame");

	gr_logx_django_entropy[0]->SetLineStyle(2);
	gr_logx_django_entropy[0]->SetLineWidth(4);
	gr_logx_django_entropy[0]->SetMarkerSize(1.4);
	gr_logx_django_entropy[0]->SetLineColor(kBlack);

	gr_logx_django_entropy[1]->SetLineStyle(2);
	gr_logx_django_entropy[1]->SetLineWidth(4);
	gr_logx_django_entropy[1]->SetMarkerSize(1.4);
	gr_logx_django_entropy[1]->SetLineColor(kBlue);

	gr_logx_django_entropy[2]->SetLineStyle(2);
	gr_logx_django_entropy[2]->SetLineWidth(4);
	gr_logx_django_entropy[2]->SetMarkerSize(1.4);
	gr_logx_django_entropy[2]->SetLineColor(kRed);

	gr_logx_django_entropy[3]->SetLineStyle(2);
	gr_logx_django_entropy[3]->SetLineWidth(4);
	gr_logx_django_entropy[3]->SetMarkerSize(1.6);
	gr_logx_django_entropy[3]->SetLineColor(kGreen-2);

	gr_logx_rapgap_entropy[0]->SetLineStyle(2);
	gr_logx_rapgap_entropy[0]->SetLineWidth(4);
	gr_logx_rapgap_entropy[0]->SetMarkerSize(1.4);
	gr_logx_rapgap_entropy[0]->SetLineColor(kBlack);

	gr_logx_rapgap_entropy[1]->SetLineStyle(3);
	gr_logx_rapgap_entropy[1]->SetLineWidth(4);
	gr_logx_rapgap_entropy[1]->SetMarkerSize(1.4);
	gr_logx_rapgap_entropy[1]->SetLineColor(kBlue);

	gr_logx_rapgap_entropy[2]->SetLineStyle(4);
	gr_logx_rapgap_entropy[2]->SetLineWidth(4);
	gr_logx_rapgap_entropy[2]->SetMarkerSize(1.4);
	gr_logx_rapgap_entropy[2]->SetLineColor(kRed);

	gr_logx_rapgap_entropy[3]->SetLineStyle(5);
	gr_logx_rapgap_entropy[3]->SetLineWidth(4);
	gr_logx_rapgap_entropy[3]->SetMarkerSize(1.6);
	gr_logx_rapgap_entropy[3]->SetLineColor(kGreen-2);

	gr_logx_pythia_entropy[0]->SetLineStyle(2);
	gr_logx_pythia_entropy[0]->SetLineWidth(4);
	gr_logx_pythia_entropy[0]->SetMarkerSize(1.4);
	gr_logx_pythia_entropy[0]->SetLineColor(kBlack);

	gr_logx_pythia_entropy[1]->SetLineStyle(3);
	gr_logx_pythia_entropy[1]->SetLineWidth(4);
	gr_logx_pythia_entropy[1]->SetMarkerSize(1.4);
	gr_logx_pythia_entropy[1]->SetLineColor(kBlue);

	gr_logx_pythia_entropy[2]->SetLineStyle(4);
	gr_logx_pythia_entropy[2]->SetLineWidth(4);
	gr_logx_pythia_entropy[2]->SetMarkerSize(1.4);
	gr_logx_pythia_entropy[2]->SetLineColor(kRed);

	gr_logx_pythia_entropy[3]->SetLineStyle(5);
	gr_logx_pythia_entropy[3]->SetLineWidth(4);
	gr_logx_pythia_entropy[3]->SetMarkerSize(1.6);
	gr_logx_pythia_entropy[3]->SetLineColor(kGreen-2);

	if( pickMCEG_ == 0 ){
		gr_logx_rapgap_entropy[0]->Draw("L same");
		gr_logx_rapgap_entropy[1]->Draw("L same");
		gr_logx_rapgap_entropy[2]->Draw("L same");
		gr_logx_rapgap_entropy[3]->Draw("L same");
	}
	else if( pickMCEG_ == 1 ){
		gr_logx_django_entropy[0]->Draw("Lsame");
		gr_logx_django_entropy[1]->Draw("Lsame");
		gr_logx_django_entropy[2]->Draw("Lsame");
		gr_logx_django_entropy[3]->Draw("Lsame");
	}
	else if( pickMCEG_ == 2 ){
		gr_logx_pythia_entropy[0]->Draw("Lsame");
		gr_logx_pythia_entropy[1]->Draw("Lsame");
		gr_logx_pythia_entropy[2]->Draw("Lsame");
		gr_logx_pythia_entropy[3]->Draw("Lsame");
	}

	TFile* PDF_TGraph = new TFile("PDF_TGraph.root");
	TGraphErrors* mstw_Q2_2 = (TGraphErrors*) PDF_TGraph->Get("gr1");
	TGraphErrors* mstw_Q2_10 = (TGraphErrors*) PDF_TGraph->Get("gr2");

	mstw_Q2_2->SetFillStyle(1001);
	mstw_Q2_2->SetFillColorAlpha(kGreen-2,0.4);
	mstw_Q2_2->SetMarkerStyle(24);
	mstw_Q2_2->SetLineColor(kWhite);

	mstw_Q2_10->SetFillStyle(1001);
	mstw_Q2_10->SetFillColorAlpha(kRed-2,0.4);
	mstw_Q2_10->SetMarkerStyle(25);
	mstw_Q2_10->SetLineColor(kWhite);

	mstw_Q2_2->Draw("P e3 same");
	mstw_Q2_10->Draw("P e3 same");

	TLegend *w6313 = new TLegend(0.57,0.18,0.8,0.3);
	w6313->SetLineColor(kWhite);
	w6313->SetFillColor(0);
	w6313->SetTextSize(18);
	w6313->SetTextFont(45);
	w6313->SetTextColor(kBlack);
	w6313->AddEntry(mstw_Q2_2, "ln(xG)  Q^{2} = 2 GeV^{2} ","PF");
	w6313->AddEntry(mstw_Q2_10, "ln(xG)  Q^{2} = 10 GeV^{2} ","PF");
	w6313->Draw("same");

	//KNO scaling
	TCanvas* c2_1 = new TCanvas("c2_1","c2_1",1,1,1000,1000);
	c2_1->Divide(4,4,0,0);

	TH1D* base2_1 = makeHist("base2_1", "", "z = N/#LTN#GT", "#Psi(z)", 100,-1,6.5,kBlack);
	base2_1->GetYaxis()->SetRangeUser(0.000007, 7);
	base2_1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base2_1,3,3.8);

	base2_1->GetYaxis()->SetTitleSize(base2_1->GetYaxis()->GetTitleSize()*1.0);
	base2_1->GetXaxis()->SetTitleSize(base2_1->GetXaxis()->GetTitleSize()*1.0);
	base2_1->GetYaxis()->SetLabelSize(base2_1->GetYaxis()->GetLabelSize()*1.0);
	base2_1->GetXaxis()->SetLabelSize(base2_1->GetXaxis()->GetLabelSize()*1.0);
	base2_1->GetXaxis()->SetNdivisions(4,6,0);
	base2_1->GetYaxis()->SetNdivisions(4,6,0);

	TGraphErrors* hist_KNO_scale[4][4];
	TGraphErrors* hist_KNO_scale_systematics[4][4];
	TGraphErrors* django_KNO_scale[4][4];
	TGraphErrors* pythia_KNO_scale[4][4];
	TGraphErrors* rapgap_KNO_scale[4][4];
	sub_panels = 1;
	for(int iy = 0; iy < 4; iy++){
		for(int iQ2 = 0; iQ2 < 4; iQ2++){
			hist_KNO_scale[iQ2][iy] = new TGraphErrors();
			django_KNO_scale[iQ2][iy] = new TGraphErrors();
			pythia_KNO_scale[iQ2][iy] = new TGraphErrors();
			rapgap_KNO_scale[iQ2][iy] = new TGraphErrors();

			c2_1->cd(sub_panels);
			
			if(sub_panels==1||sub_panels==5||sub_panels==9||sub_panels==13) gPad->SetLeftMargin(0.2);
			if(sub_panels==1||sub_panels==2||sub_panels==3||sub_panels==4) gPad->SetTopMargin(0.15);
			if(sub_panels>=13) gPad->SetBottomMargin(0.2);
			gPad->SetTicks();
			gPad->SetLogy(1);
			base2_1->Draw("");
			double x_tem,mean_multiplicity;
			gr_data[iQ2]->GetPoint(iy,x_tem,mean_multiplicity);
			double mean_multiplicity_error = gr_data[iQ2]->GetErrorY(iy);
			double binwidth[100];
			for(int ibin=0;ibin<hist_unfolded_from_rapgap_allEta[iQ2][iy]->GetNbinsX();ibin++){

				double bin_center = hist_unfolded_from_rapgap_allEta[iQ2][iy]->GetBinCenter(ibin+1);
				double pn = hist_unfolded_from_rapgap_allEta[iQ2][iy]->GetBinContent(ibin+1);
				double pn_error = hist_unfolded_from_rapgap_allEta[iQ2][iy]->GetBinError(ibin+1);
				binwidth[ibin] = hist_unfolded_from_rapgap_allEta[iQ2][iy]->GetBinWidth(ibin+1)/mean_multiplicity;
				double value = mean_multiplicity*pn;
				double z = bin_center/mean_multiplicity;
				double error = mean_multiplicity*pn_error;

				if( value < 1e-5) continue;
				hist_KNO_scale[iQ2][iy]->SetPoint(ibin,z,value);
				hist_KNO_scale[iQ2][iy]->SetPointError(ibin,0,error);

			}


			hist_KNO_scale_systematics[iQ2][iy] = (TGraphErrors*) hist_KNO_scale[iQ2][iy]->Clone("hsys_kno");
			for(int isys=0;isys<hist_KNO_scale_systematics[iQ2][iy]->GetN();isys++){
				double z,psi;
				hist_KNO_scale_systematics[iQ2][iy]->GetPoint(isys, z, psi);
				double N = z*mean_multiplicity;
				double error = 0.;
				if( N<2 ) error = 0.221*psi;
				if( N>=2 && N<5 ) error = 0.06*psi;
				if( N>=5 && N<15 ) error = 0.06*psi;
				if( N>=15 ) error = 0.233*psi;

				hist_KNO_scale_systematics[iQ2][iy]->SetPoint(isys,z,psi);
				hist_KNO_scale_systematics[iQ2][iy]->SetPointError(isys,binwidth[isys]/2.,error);
			}
			gStyle->SetErrorX( 0.5 );
			hist_KNO_scale_systematics[iQ2][iy]->SetFillColorAlpha(kGray+1,0.8);
  			hist_KNO_scale_systematics[iQ2][iy]->Draw("e2same");

  			hist_KNO_scale[iQ2][iy]->SetMarkerStyle(20);
			hist_KNO_scale[iQ2][iy]->Draw("PE0Zsame");

			x_tem = 0.;mean_multiplicity=0.;mean_multiplicity_error=0.;
			gr_django[iQ2]->GetPoint(iy,x_tem,mean_multiplicity);
			mean_multiplicity_error = gr_django[iQ2]->GetErrorY(iy);
			hist_mc_django[iQ2][iy]->Scale(1./hist_mc_django[iQ2][iy]->Integral());
			
			for(int ibin=0;ibin<hist_mc_django[iQ2][iy]->GetNbinsX();ibin++){
				double bin_center = hist_mc_django[iQ2][iy]->GetBinCenter(ibin+1);
				double pn = hist_mc_django[iQ2][iy]->GetBinContent(ibin+1);
				double pn_error = hist_mc_django[iQ2][iy]->GetBinError(ibin+1);
				double value = mean_multiplicity*pn;
				double z = bin_center/mean_multiplicity;
				double error = mean_multiplicity*pn_error;
				if( value < 2e-4) continue;
				django_KNO_scale[iQ2][iy]->SetPoint(ibin,z,value);
				django_KNO_scale[iQ2][iy]->SetPointError(ibin,0,error);

			}
			x_tem = 0.;mean_multiplicity=0.;mean_multiplicity_error=0.;
			gr_rapgap[iQ2]->GetPoint(iy,x_tem,mean_multiplicity);
			mean_multiplicity_error = gr_rapgap[iQ2]->GetErrorY(iy);
			hist_mc_rapgap[iQ2][iy]->Scale(1./hist_mc_rapgap[iQ2][iy]->Integral());
			
			for(int ibin=0;ibin<hist_mc_rapgap[iQ2][iy]->GetNbinsX();ibin++){
				double bin_center = hist_mc_rapgap[iQ2][iy]->GetBinCenter(ibin+1);
				double pn = hist_mc_rapgap[iQ2][iy]->GetBinContent(ibin+1);
				double pn_error = hist_mc_rapgap[iQ2][iy]->GetBinError(ibin+1);
				double value = mean_multiplicity*pn;
				double z = bin_center/mean_multiplicity;
				double error = mean_multiplicity*pn_error;
				if( value < 2e-4) continue;
				rapgap_KNO_scale[iQ2][iy]->SetPoint(ibin,z,value);
				rapgap_KNO_scale[iQ2][iy]->SetPointError(ibin,0,error);

			}
			x_tem = 0.;mean_multiplicity=0.;mean_multiplicity_error=0.;
			gr_pythia[iQ2]->GetPoint(iy,x_tem,mean_multiplicity);
			mean_multiplicity_error = gr_pythia[iQ2]->GetErrorY(iy);
			hist_mc_pythia[iQ2][iy]->Scale(1./hist_mc_pythia[iQ2][iy]->Integral());

			for(int ibin=0;ibin<hist_mc_pythia[iQ2][iy]->GetNbinsX();ibin++){
				double bin_center = hist_mc_pythia[iQ2][iy]->GetBinCenter(ibin+1);
				double pn = hist_mc_pythia[iQ2][iy]->GetBinContent(ibin+1);
				double pn_error = hist_mc_pythia[iQ2][iy]->GetBinError(ibin+1);
				double value = mean_multiplicity*pn;
				double z = bin_center/mean_multiplicity;
				double error = mean_multiplicity*pn_error;
				if( value < 2e-4) continue;
				pythia_KNO_scale[iQ2][iy]->SetPoint(ibin,z,value);
				pythia_KNO_scale[iQ2][iy]->SetPointError(ibin,0,error);

			}
			django_KNO_scale[iQ2][iy]->SetFillStyle(1001);
			django_KNO_scale[iQ2][iy]->SetFillColorAlpha(kGreen-2,0.4);
			django_KNO_scale[iQ2][iy]->SetMarkerColor(kGreen-2);

			rapgap_KNO_scale[iQ2][iy]->SetFillStyle(1001);
			rapgap_KNO_scale[iQ2][iy]->SetFillColorAlpha(kBlue-2,0.4);
			rapgap_KNO_scale[iQ2][iy]->SetMarkerColor(kBlue-2);

			pythia_KNO_scale[iQ2][iy]->SetFillStyle(1001);
			pythia_KNO_scale[iQ2][iy]->SetFillColorAlpha(kRed-2,0.4);
			pythia_KNO_scale[iQ2][iy]->SetMarkerColor(kRed-2);
			
			django_KNO_scale[iQ2][iy]->Draw("PE3 Psame");
			rapgap_KNO_scale[iQ2][iy]->Draw("PE3 Psame");
			pythia_KNO_scale[iQ2][iy]->Draw("PE3 Psame");

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
	w6->Draw("same");
	c2_1->cd(4);
	r50->Draw("same");
	c2_1->cd(1);
	r51->Draw("same");
	c2_1->cd(1);
	r52->Draw("same");

	//KNO scaling
	TCanvas* c2_2 = new TCanvas("c2_2","c2_2",1,1,600,600);
	c2_2->Divide(2,2,0,0);
	TH1D* base2_2 = (TH1D*) base2_1->Clone("base2_2");
	fixedFontHist1D(base2_2,2,2);
	base2_2->GetYaxis()->SetNdivisions(3,2,0);

	c2_2->cd(1);
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);

	base2_2->Draw();
	hist_KNO_scale[0][0]->SetMarkerStyle(20);
	hist_KNO_scale[0][1]->SetMarkerStyle(21);
	hist_KNO_scale[0][2]->SetMarkerStyle(24);
	hist_KNO_scale[0][3]->SetMarkerStyle(25);

	hist_KNO_scale[0][0]->Draw("Psame");
	hist_KNO_scale[0][1]->Draw("Psame");
	hist_KNO_scale[0][2]->Draw("Psame");
	hist_KNO_scale[0][3]->Draw("Psame");

	hist_KNO_scale_systematics[0][0]->Draw("e3same");
	hist_KNO_scale_systematics[0][1]->Draw("e3same");
	hist_KNO_scale_systematics[0][2]->Draw("e3same");
	hist_KNO_scale_systematics[0][3]->Draw("e3same");


	c2_2->cd(2);
	gPad->SetLogy(1);
	gPad->SetRightMargin(0.1);

	base2_2->Draw();
	hist_KNO_scale[1][0]->SetMarkerStyle(20);
	hist_KNO_scale[1][1]->SetMarkerStyle(21);
	hist_KNO_scale[1][2]->SetMarkerStyle(24);
	hist_KNO_scale[1][3]->SetMarkerStyle(25);

	hist_KNO_scale[1][0]->Draw("Psame");
	hist_KNO_scale[1][1]->Draw("Psame");
	hist_KNO_scale[1][2]->Draw("Psame");
	hist_KNO_scale[1][3]->Draw("Psame");

	hist_KNO_scale_systematics[1][0]->Draw("e3same");
	hist_KNO_scale_systematics[1][1]->Draw("e3same");
	hist_KNO_scale_systematics[1][2]->Draw("e3same");
	hist_KNO_scale_systematics[1][3]->Draw("e3same");

	c2_2->cd(3);
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);

	base2_2->Draw();
	hist_KNO_scale[2][0]->SetMarkerStyle(20);
	hist_KNO_scale[2][1]->SetMarkerStyle(21);
	hist_KNO_scale[2][2]->SetMarkerStyle(24);
	hist_KNO_scale[2][3]->SetMarkerStyle(25);

	hist_KNO_scale[2][0]->Draw("Psame");
	hist_KNO_scale[2][1]->Draw("Psame");
	hist_KNO_scale[2][2]->Draw("Psame");
	hist_KNO_scale[2][3]->Draw("Psame");

	hist_KNO_scale_systematics[2][0]->Draw("e3same");
	hist_KNO_scale_systematics[2][1]->Draw("e3same");
	hist_KNO_scale_systematics[2][2]->Draw("e3same");
	hist_KNO_scale_systematics[2][3]->Draw("e3same");

	c2_2->cd(4);
	gPad->SetLogy(1);
	gPad->SetRightMargin(0.1);
	gPad->SetBottomMargin(0.15);

	base2_2->Draw();
	hist_KNO_scale[3][0]->SetMarkerStyle(20);
	hist_KNO_scale[3][1]->SetMarkerStyle(21);
	hist_KNO_scale[3][2]->SetMarkerStyle(24);
	hist_KNO_scale[3][3]->SetMarkerStyle(25);

	hist_KNO_scale[3][0]->Draw("Psame");
	hist_KNO_scale[3][1]->Draw("Psame");
	hist_KNO_scale[3][2]->Draw("Psame");
	hist_KNO_scale[3][3]->Draw("Psame");

	hist_KNO_scale_systematics[3][0]->Draw("e3same");
	hist_KNO_scale_systematics[3][1]->Draw("e3same");
	hist_KNO_scale_systematics[3][2]->Draw("e3same");
	hist_KNO_scale_systematics[3][3]->Draw("e3same");

	c2_2->cd(2);
	TLatex* r444 = new TLatex(0.8,0.92, "H1");
	r444->SetNDC();
	r444->SetTextSize(20);
	r444->SetTextFont(63);
	r444->SetTextColor(kBlack);
	r444->Draw("same");
	c2_2->cd(1);
	TLatex* r488 = new TLatex(0.18, 0.92, "ep 27.5#times920 GeV");
	r488->SetNDC();
	r488->SetTextSize(18);
	r488->SetTextFont(43);
	r488->SetTextColor(kBlack);
	r488->Draw("same");
	c2_2->cd(2);
	TLatex* r499 = new TLatex(0.05, 0.92, "0 < #eta* < 4.0");
	r499->SetNDC();
	r499->SetTextSize(18);
	r499->SetTextFont(43);
	r499->SetTextColor(kBlack);
	r499->Draw("same");

	c2_2->cd(1);
	TLatex* ry1 = new TLatex(0.18,0.05, "Q^{2} (5,10) GeV^{2}");
	ry1->SetNDC();
	ry1->SetTextSize(18);
	ry1->SetTextFont(44);
	ry1->SetTextColor(kBlack);
	ry1->Draw("same");

	c2_2->cd(2);
	TLatex* ry2 = new TLatex(0.05,0.05, "Q^{2} (10,20) GeV^{2}");
	ry2->SetNDC();
	ry2->SetTextSize(18);
	ry2->SetTextFont(44);
	ry2->SetTextColor(kBlack);
	ry2->Draw("same");

	c2_2->cd(3);
	TLatex* ry3 = new TLatex(0.18,0.2, "Q^{2} (20,40) GeV^{2}");
	ry3->SetNDC();
	ry3->SetTextSize(18);
	ry3->SetTextFont(44);
	ry3->SetTextColor(kBlack);
	ry3->Draw("same");

	c2_2->cd(4);
	TLatex* ry4 = new TLatex(0.05,0.2, "Q^{2} (40,100) GeV^{2}");
	ry4->SetNDC();
	ry4->SetTextSize(18);
	ry4->SetTextFont(44);
	ry4->SetTextColor(kBlack);
	ry4->Draw("same");

	c2_2->cd(4);

	TLatex* ry5 = new TLatex(0.53,0.9, "#LT W #GT in GeV");
	ry5->SetNDC();
	ry5->SetTextSize(18);
	ry5->SetTextFont(44);
	ry5->SetTextColor(kBlack);
	ry5->Draw("same");
	TLegend *w66 = new TLegend(0.57,0.63,0.89,0.85);
	w66->SetLineColor(kWhite);
	w66->SetFillColor(0);
	w66->SetTextSize(17);
	w66->SetTextFont(45);
	w66->AddEntry(hist_KNO_scale[3][0], "72 ","P");
	w66->AddEntry(hist_KNO_scale[3][1], "101 ","P");
	w66->AddEntry(hist_KNO_scale[3][2], "141 ","P");
	w66->AddEntry(hist_KNO_scale[3][3], "205 ", "P");//rapgap
	w66->Draw("same");

	if(draw_) c2->Print(Form("../figures/results_finalv3/Pn_allEtaStar_"+radname+".pdf"));
	if(draw_) c2_1->Print(Form("../figures/results_finalv3/PsiZ_allEtaStar_"+radname+".pdf"));
	if(draw_) c2_2->Print(Form("../figures/results_finalv3/PsiZ_combine_allEtaStar_"+radname+".pdf"));
	if(draw_) c3->Print(Form("../figures/results_finalv3/Mean_Nch_allEtaStar_"+radname+"_MCEG_"+MCEGname+".pdf"));
    if(draw_) c3_2->Print(Form("../figures/results_finalv3/Var_Nch_allEtaStar_"+radname+"_MCEG_"+MCEGname+".pdf"));
	if(draw_) c4->Print(Form("../figures/results_finalv3/Mean_Nch_allEtaStar_logX_"+radname+"_MCEG_"+MCEGname+".pdf"));
	if(draw_) c4_2->Print(Form("../figures/results_finalv3/Var_Nch_allEtaStar_logX_"+radname+"_MCEG_"+MCEGname+".pdf"));
	if(draw_) c4_3->Print(Form("../figures/results_finalv3/EE_allEtaStar_logX_"+radname+"_MCEG_"+MCEGname+".pdf"));

	// TFile outfile_temp("example_unfoldedHist.root","RECREATE");
	// hist_unfolded_from_rapgap_allEta[3][3]->Write();
	// outfile_temp.Write();

}
