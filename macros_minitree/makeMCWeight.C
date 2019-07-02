#include "RiceStyle.h"

using namespace std;

void makeMCWeight(){

	TH2D* DATA_Q2vsX = 0;
	TH1D* DATA_vtxZ = 0;
	TH2D* DATA2_Q2vsX = 0;
	TH1D* DATA2_vtxZ = 0;
	
	TH2D* RAPGAP_Q2vsX = 0;
	TH1D* RAPGAP_vtxZ = 0;
	TH2D* DJANGO_Q2vsX = 0;
	TH1D* DJANGO_vtxZ = 0;

	TFile* file_data = new TFile("../minitree_output/Pn_hist_data.root");
	DATA_Q2vsX = (TH2D*) file_data->Get("h_Q2vsX_1");
	DATA_Q2vsX->Scale( 1.0/DATA_Q2vsX->Integral() );
	DATA2_Q2vsX = (TH2D*) DATA_Q2vsX->Clone("h_Q2vsX_2");

	DATA_vtxZ = (TH1D*) file_data->Get("h_vtxZ_1");
	DATA_vtxZ->Scale( 1.0/DATA_vtxZ->Integral() );
	DATA2_vtxZ = (TH1D*) DATA_vtxZ->Clone("h_vtxZ_2");

	TFile* file_mc_rapgap = 0;
	TFile* file_mc_django = 0;
	file_mc_rapgap = new TFile("../minitree_output/Pn_hist_rapgap.root");
	file_mc_django = new TFile("../minitree_output/Pn_hist_django.root");

	RAPGAP_Q2vsX = (TH2D*) file_mc_rapgap->Get("h_Q2vsX_1");
	RAPGAP_Q2vsX->Scale( 1.0/RAPGAP_Q2vsX->Integral() );

	RAPGAP_vtxZ = (TH1D*) file_mc_rapgap->Get("h_vtxZ_1");
	RAPGAP_vtxZ->Scale( 1.0/RAPGAP_vtxZ->Integral() );

	DATA_Q2vsX->Divide( RAPGAP_Q2vsX );
	DATA_vtxZ->Divide( RAPGAP_vtxZ );

	DJANGO_Q2vsX = (TH2D*) file_mc_django->Get("h_Q2vsX_1");
	DJANGO_Q2vsX->Scale( 1.0/DJANGO_Q2vsX->Integral() );

	DJANGO_vtxZ = (TH1D*) file_mc_django->Get("h_vtxZ_1");
	DJANGO_vtxZ->Scale( 1.0/DJANGO_vtxZ->Integral() );

	DATA2_Q2vsX->Divide( DJANGO_Q2vsX );
	DATA2_vtxZ->Divide( DJANGO_vtxZ );

	cout << " double Q2vsX_weight_django[]={";
	for(int i=0;i<DATA2_Q2vsX->GetNbinsX();i++){
		for(int j=0;j<DATA2_Q2vsX->GetNbinsY();j++){

			double value = DATA2_Q2vsX->GetBinContent(i+1,j+1);
			cout << value << ",";
		}
	}
	cout << endl;


	cout << " double vtxz_weight_rapgap[]={";
	for(int i = 0; i < DATA_vtxZ->GetNbinsX(); i++){
		double value = DATA_vtxZ->GetBinContent(i+1);
		double error = DATA_vtxZ->GetBinError(i+1);

		cout << value << ",";
	}
	cout << endl;

	cout << " double vtxz_weight_django[]={";
	for(int i = 0; i < DATA2_vtxZ->GetNbinsX(); i++){
		double value = DATA2_vtxZ->GetBinContent(i+1);
		double error = DATA2_vtxZ->GetBinError(i+1);

		cout << value << ",";
	}
	cout << endl;

	TFile out("MCweight.root","RECREATE");

	DATA_Q2vsX->Write();
	DATA_vtxZ->Write();
	DATA2_Q2vsX->Write();
	DATA2_vtxZ->Write();
	


	

}