#include "RiceStyle.h"
#include <string>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

void testExtrapolation(){

	TF1* tf2 = new TF1(Form("tf2_%d_%d",0,0),"[0]*([1]/[1]+1)*ROOT::Math::negative_binomial_pdf(x[0],[1],[2])",0,15);
	tf2->SetParameter(0,4.76997e-01);
	tf2->SetParameter(1,4.26874e-01);
	tf2->SetParameter(2,4.75684e+00);
	tf2->SetLineWidth(2);

	double pnbins[]={0,1,2,3,4,5,6,7,8,10,14,18,22,28,32};
	int Npnbins = sizeof(pnbins)/sizeof(pnbins[0]) - 1;

	TH1D* hGen = new TH1D("hGen","hGen",30,0,30);
	TH1D* hRec = new TH1D("hRec","hRec",Npnbins,pnbins);

	const int Nevents = 1e7;

	for(int i = 0; i <Nevents; i++){
		double value = tf2->GetRandom();
		hGen->Fill( value );
		hRec->Fill( value );
	}
	for(int j=0;j<hGen->GetNbinsX();j++){
		hGen->SetBinContent(j+1, hGen->GetBinContent(j+1)/hGen->GetBinWidth(j+1));
	}
	hGen->Scale(1./Nevents);
	hGen->Draw("");
	for(int j=0;j<hRec->GetNbinsX();j++){
		hRec->SetBinContent(j+1, hRec->GetBinContent(j+1)/hRec->GetBinWidth(j+1));
	}
	hRec->Scale(1./Nevents);
	hRec->SetLineColor(kRed);
	hRec->Draw("same");

	/*
	generated truth
	*/

	// <n**j> jth moment of n 
	double sum = 0.;
	double sum_pn[4];
	for(int i=0;i<4;i++){sum_pn[i]=0.;}//initialize sums

	for(int i=0;i<hGen->GetNbinsX();i++){
		double N = hGen->GetBinCenter(i+1) - hGen->GetBinWidth(i+1)/2.;
		double Pn = hGen->GetBinContent(i+1);
		sum_pn[0] += TMath::Power(N,1)*Pn;
		sum_pn[1] += TMath::Power(N,2)*Pn;
		sum_pn[2] += TMath::Power(N,3)*Pn;
		sum_pn[3] += TMath::Power(N,4)*Pn;
		sum += Pn;
	}
	cout << "<n**1> = " << sum_pn[0]/sum << endl;
	cout << "<n**2> = " << sum_pn[1]/sum << endl;
	cout << "<n**3> = " << sum_pn[2]/sum << endl;
	cout << "<n**4> = " << sum_pn[3]/sum << endl;

	cout << "normalized moments ~ " << endl;
	cout << "<n**1>/<n>**1 = " << (sum_pn[0]/sum) / TMath::Power(sum_pn[0]/sum, 1) << endl;
	cout << "<n**2>/<n>**2 = " << (sum_pn[1]/sum) / TMath::Power(sum_pn[0]/sum, 2) << endl;
	cout << "<n**3>/<n>**3 = " << (sum_pn[2]/sum) / TMath::Power(sum_pn[0]/sum, 3) << endl;
	cout << "<n**4>/<n>**4 = " << (sum_pn[3]/sum) / TMath::Power(sum_pn[0]/sum, 4) << endl;
	
	double first_moment = sum_pn[0]/sum;
	// jth central moment of <(x-x0)**j>
	sum = 0.;
	for(int i=0;i<4;i++){sum_pn[i]=0.;}//initialize sums
	for(int i=0;i<hGen->GetNbinsX();i++){
		double N = hGen->GetBinCenter(i+1) - hGen->GetBinWidth(i+1)/2.;
		double Pn = hGen->GetBinContent(i+1);
		
		sum_pn[0] += TMath::Power((N-first_moment),1)*Pn;
		sum_pn[1] += TMath::Power((N-first_moment),2)*Pn;
		sum_pn[2] += TMath::Power((N-first_moment),3)*Pn;
		sum_pn[3] += TMath::Power((N-first_moment),4)*Pn;
		sum += Pn;
	}
	cout << "<(n-n0)**1> = " << sqrt(sum_pn[0]/sum) << endl;
	cout << "<(n-n0)**2> = " << sqrt(sum_pn[1]/sum) << endl;
	cout << "<(n-n0)**3> = " << sqrt(sum_pn[2]/sum) << endl;
	cout << "<(n-n0)**4> = " << sqrt(sum_pn[3]/sum) << endl;

	/*
	reconstructed with variable bin width.
	*/

	cout << "starting reconstruction/extrapolation " << endl;
	TH1D* hRec_extrapol = new TH1D("hRec_extrapol",";N",30,0,30);
	
	TF1* f1[10];
	int number_of_extrapo = 0;
	int stop_index = 0;
	double N_index[10];
	for(int irec=0;irec<hRec->GetNbinsX();irec++){
		double Nr = hRec->GetBinCenter(irec+1);
		double Nr_last = hRec->GetBinCenter(irec);
		double Pn = hRec->GetBinContent(irec+1);
		if( Pn > 0 ) stop_index++;
		double binwidth = hRec->GetBinWidth(irec+1);
		if( binwidth != 1. ){
			
			f1[number_of_extrapo] = new TF1(Form("f1_%d",number_of_extrapo),"([0]-[1])/([2]-[3])*(x[0]-[3]) + [1]",Nr_last,Nr);
			f1[number_of_extrapo]->SetParameter(0,hRec->GetBinContent(irec+1));
			f1[number_of_extrapo]->SetParameter(1,hRec->GetBinContent(irec));
			f1[number_of_extrapo]->SetParameter(2,Nr);
			f1[number_of_extrapo]->SetParameter(3,Nr_last);
			// f1[number_of_extrapo]->Draw("same");
			N_index[number_of_extrapo] = Nr;
			number_of_extrapo++;
		}
	}
	sum = 0.;
	for(int i=0;i<4;i++){sum_pn[i]=0.;}//initialize sums
	int ibin=0;
	int iextra=0;
	for(int jbin=0;jbin<hRec_extrapol->GetNbinsX();jbin++){
		double Nt = hRec_extrapol->GetBinCenter(jbin+1);
		double Nr = hRec->GetBinCenter(ibin+1);
		double Pn = hRec->GetBinContent(ibin+1);
		if( Nt == Nr ){
			hRec_extrapol->SetBinContent( jbin+1, Pn );
			ibin++;
		}
		else if( Nt != Nr ){
			if( Nt < N_index[iextra] ){
				Pn = f1[iextra]->Eval(Nt);
			}
			else{
				iextra++;
				ibin++;
				Pn = f1[iextra]->Eval(Nt);
			}
			
			hRec_extrapol->SetBinContent( jbin+1, Pn );
		}	
		double binwidth = hRec_extrapol->GetBinWidth(jbin+1);
		sum_pn[0] += TMath::Power(Nt-binwidth/2.,1)*Pn;
		sum_pn[1] += TMath::Power(Nt-binwidth/2.,2)*Pn;
		sum_pn[2] += TMath::Power(Nt-binwidth/2.,3)*Pn;
		sum_pn[3] += TMath::Power(Nt-binwidth/2.,4)*Pn;
		sum += Pn;
	}
	hRec_extrapol->SetLineColor(kBlue);
	hRec_extrapol->Draw("same");

	first_moment = sum_pn[0]/sum;

	cout << "<n**1> = " << sum_pn[0]/sum << endl;
	cout << "<n**2> = " << sum_pn[1]/sum << endl;
	cout << "<n**3> = " << sum_pn[2]/sum << endl;
	cout << "<n**4> = " << sum_pn[3]/sum << endl;

	cout << "normalized moments ~ " << endl;
	cout << "<n**1>/<n>**1 = " << (sum_pn[0]/sum) / TMath::Power(sum_pn[0]/sum, 1) << endl;
	cout << "<n**2>/<n>**2 = " << (sum_pn[1]/sum) / TMath::Power(sum_pn[0]/sum, 2) << endl;
	cout << "<n**3>/<n>**3 = " << (sum_pn[2]/sum) / TMath::Power(sum_pn[0]/sum, 3) << endl;
	cout << "<n**4>/<n>**4 = " << (sum_pn[3]/sum) / TMath::Power(sum_pn[0]/sum, 4) << endl;
	
	first_moment = sum_pn[0]/sum;
	// jth central moment of <(x-x0)**j>
	sum = 0.;
	for(int i=0;i<4;i++){sum_pn[i]=0.;}//initialize sums
	for(int i=0;i<hRec_extrapol->GetNbinsX();i++){
		double N = hRec_extrapol->GetBinCenter(i+1) - hRec_extrapol->GetBinWidth(i+1)/2.;
		double Pn = hRec_extrapol->GetBinContent(i+1);
		
		sum_pn[0] += TMath::Power((N-first_moment),1)*Pn;
		sum_pn[1] += TMath::Power((N-first_moment),2)*Pn;
		sum_pn[2] += TMath::Power((N-first_moment),3)*Pn;
		sum_pn[3] += TMath::Power((N-first_moment),4)*Pn;
		sum += Pn;
	}
	cout << "<(n-n0)**1> = " << sqrt(sum_pn[0]/sum) << endl;
	cout << "<(n-n0)**2> = " << sqrt(sum_pn[1]/sum) << endl;
	cout << "<(n-n0)**3> = " << sqrt(sum_pn[2]/sum) << endl;
	cout << "<(n-n0)**4> = " << sqrt(sum_pn[3]/sum) << endl;

	// TH1D* hKNO = new TH1D("hKNO",";z",30,0,3);
	// for(int i=0;i<hGen->GetNbinsX();i++){
	// 	double N = hGen->GetBinCenter(i+1) - hGen->GetBinWidth(i+1)/2.;
	// 	double z = N/first_moment;
	// 	int bin = hKNO->FindBin( z );
	// 	double Pn = hGen->GetBinContent(i+1);
	// 	double value = Pn*first_moment;
	// 	hKNO->SetBinContent(bin, value);
	// }

	// TCanvas* c2 = new TCanvas();
	// hKNO->SetMarkerStyle(20);
	// hKNO->Draw("P");




}