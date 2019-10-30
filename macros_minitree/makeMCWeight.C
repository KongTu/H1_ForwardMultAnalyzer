#include "RiceStyle.h"

using namespace std;

void makeMCWeight(){

	TH1D* DATA_y = 0;
	TH1D* DATA_vtxZ = 0;
	TH1D* DATA2_y = 0;
	TH1D* DATA2_vtxZ = 0;
	
	TH1D* RAPGAP_y = 0;
	TH1D* RAPGAP_vtxZ = 0;
	TH1D* DJANGO_y = 0;
	TH1D* DJANGO_vtxZ = 0;

	TFile* file_data = new TFile("../minitree_output/Pn_hist_data_hadCali.root");
	DATA_y = (TH1D*) file_data->Get("h_y_1");
	DATA_y->Scale( 1.0/DATA_y->Integral() );
	DATA2_y = (TH1D*) DATA_y->Clone("h_y_2");

	DATA_vtxZ = (TH1D*) file_data->Get("h_vtxZ_1");
	DATA_vtxZ->Scale( 1.0/DATA_vtxZ->Integral() );
	DATA2_vtxZ = (TH1D*) DATA_vtxZ->Clone("h_vtxZ_2");

	TFile* file_mc_rapgap = 0;
	TFile* file_mc_django = 0;
	file_mc_rapgap = new TFile("../minitree_output/Pn_hist_rapgap_hadCali.root");
	file_mc_django = new TFile("../minitree_output/Pn_hist_django_hadCali.root");

	RAPGAP_y = (TH1D*) file_mc_rapgap->Get("h_y_1");
	RAPGAP_y->Scale( 1.0/RAPGAP_y->Integral() );

	RAPGAP_vtxZ = (TH1D*) file_mc_rapgap->Get("h_vtxZ_1");
	RAPGAP_vtxZ->Scale( 1.0/RAPGAP_vtxZ->Integral() );

	DATA_y->Divide( RAPGAP_y );
	DATA_vtxZ->Divide( RAPGAP_vtxZ );

	DJANGO_y = (TH1D*) file_mc_django->Get("h_y_1");
	DJANGO_y->Scale( 1.0/DJANGO_y->Integral() );

	DJANGO_vtxZ = (TH1D*) file_mc_django->Get("h_vtxZ_1");
	DJANGO_vtxZ->Scale( 1.0/DJANGO_vtxZ->Integral() );

	DATA2_y->Divide( DJANGO_y );
	DATA2_vtxZ->Divide( DJANGO_vtxZ );

	cout << " double y_weight_rapgap[]={";
	for(int i=0;i<DATA_y->GetNbinsX();i++){

		double value = DATA_y->GetBinContent(i+1);
		cout << value << ",";
	}
	cout << "}" << endl;

	cout << " double y_weight_django[]={";
	for(int i=0;i<DATA2_y->GetNbinsX();i++){

		double value = DATA2_y->GetBinContent(i+1);
		cout << value << ",";
		
	}
	cout << "}" << endl;


	cout << " double vtxz_weight_rapgap[]={";
	for(int i = 0; i < DATA_vtxZ->GetNbinsX(); i++){
		double value = DATA_vtxZ->GetBinContent(i+1);
		double error = DATA_vtxZ->GetBinError(i+1);

		cout << value << ",";
	}
	cout << "}" << endl;

	cout << " double vtxz_weight_django[]={";
	for(int i = 0; i < DATA2_vtxZ->GetNbinsX(); i++){
		double value = DATA2_vtxZ->GetBinContent(i+1);
		double error = DATA2_vtxZ->GetBinError(i+1);

		cout << value << ",";
	}
	cout << "}" << endl;

	TFile out("MCweight_hadCali.root","RECREATE");

	DATA_y->Write();
	DATA_vtxZ->Write();
	DATA2_y->Write();
	DATA2_vtxZ->Write();
	


	

}