#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>
#include <TChain.h>
#include <TMath.h>
#include <TRandom.h>
#include <TF1.h>
#include <string>
#include "TLorentzVector.h"

using namespace std;

void testMiniTree(){

	TFile* file = 0;
	file = new TFile("../new_output/mc_highE_DJANGOH_fullReweight_Tree.root");

	TTree* tree = (TTree*) file->Get("miniTree");
	Int_t typeChgREC_mini[200];
	Int_t nRECtrack_mini;
	tree->SetBranchAddress("typeChgREC_mini",&typeChgREC_mini);
	tree->SetBranchAddress("nRECtrack_mini",&nRECtrack_mini);

	TH1D* h_typeChgTemp = new TH1D("h_typeChgTemp", "h_typeChgTemp ",9,-4,5);
	for(int ievent = 0; ievent < tree->GetEntries(); ievent++){
		tree->GetEntry(ievent);
		for(int itrk = 0; itrk < nRECtrack_mini; itrk++){
			if( typeChgREC_mini[itrk] == 0 ) cout << "typeChgREC_mini[itrk] ~ " << typeChgREC_mini[itrk] << endl;
			h_typeChgTemp->Fill( typeChgREC_mini[itrk] );
		}
	}
	
	TFile * output = new TFile("testMiniTree.root","RECREATE");
	h_typeChgTemp->Write();
	output->Write();
}