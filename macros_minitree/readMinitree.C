#include "RiceStyle.h"

using namespace std;

TH2D* h_Q2vsX_1 = new TH2D("h_Q2vsX_1","h_Q2vsX_1",1000,0.00001,0.01,50,1,100);
TH1D* h_vtxZ_1  = new TH1D("h_vtxZ_1","h_vtxZ_1", 100,-50,50);

void readMinitree(const bool isdata_ = true, const bool isdjango_ = false){

	TFile* file = 0;
	if(isdata_){
		file = new TFile("/Users/kong/google_drive/BNL_folder/Work/DESY/H1_Collaboration/Analysis/H1_ForwardMultAnalyzer/new_output/data_highE_0607_noReweight_Tree_new_check.root");
	}
	else if(isdjango_) {
		file = new TFile("/Users/kong/google_drive/BNL_folder/Work/DESY/H1_Collaboration/Analysis/H1_ForwardMultAnalyzer/new_output/mc_highE_DJANGOH_noReweight_Tree_new_check.root");
	}
	else{
		file = new TFile("/Users/kong/google_drive/BNL_folder/Work/DESY/H1_Collaboration/Analysis/H1_ForwardMultAnalyzer/new_output/mc_highE_RAPGAP_noReweight_Tree_new_check.root");
	}
	TTree* tree = (TTree*) file->Get("miniTree");

	double Q2min = 5.;
	double Q2max = 100.;
	double ymin = 0.0375;
	double ymax = 0.6;
	double ybins[] = {0.0375,0.075,0.15,0.3,0.6};
	double eta_bins[] = {-1.2,0.2,-0.5,0.9,0.2,1.6};
	double Q2_bins[] = {5,10,20,40,100};

	TString seta_bins[6]={"-1.2","0.2","-0.5","0.9","0.2","1.6"};
	TString sQ2_bins[5]={"5","10","20","40","100"};
	TString sy_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	Float_t xMC_es_mini,yMC_es_mini,Q2MC_es_mini;
	Float_t w_mini;
	
	const int nMCtrack_MAX=400;
	// if there is no MC info, nMCtrack is set to zero
	Int_t nMCtrack_mini;
	Float_t pxMC_mini[nMCtrack_MAX];
	Float_t pyMC_mini[nMCtrack_MAX];
	Float_t pzMC_mini[nMCtrack_MAX];
	Float_t etaMC_mini[nMCtrack_MAX];
	Float_t chargeMC_mini[nMCtrack_MAX];
   	
	tree->SetBranchAddress("yMC_es_mini",&yMC_es_mini);
	tree->SetBranchAddress("Q2MC_es_mini",&Q2MC_es_mini);

	tree->SetBranchAddress("nMCtrack_mini",&nMCtrack_mini);
	tree->SetBranchAddress("pxMC_mini",&pxMC_mini);
	tree->SetBranchAddress("pyMC_mini",&pyMC_mini);
	tree->SetBranchAddress("pzMC_mini",&pzMC_mini);
	tree->SetBranchAddress("etaMC_mini",&etaMC_mini);

	Float_t xREC_es_mini,yREC_es_mini,Q2REC_es_mini;
	const int nRECtrack_MAX = 200;
	Int_t eventpass_mini;
	Int_t nRECtrack_mini;
	Int_t typeChgREC_mini[nRECtrack_MAX];
	Float_t etaAsymREC_mini;
	Float_t EpzREC_mini;
	Int_t totalMultREC_mini;
	Float_t vertex_mini[3];
	Float_t hfsEREC_mini, hfsPtREC_mini, hfsPzREC_mini;
 	Float_t elecEREC_mini, elecPtREC_mini, elecPzREC_mini;

	Float_t elecTrackMatchPhiREC_mini;
  	Float_t elecTrackMatchThetaREC_mini;

	Float_t pxREC_mini[nRECtrack_MAX];
	Float_t pyREC_mini[nRECtrack_MAX];
	Float_t pzREC_mini[nRECtrack_MAX];
	Float_t pREC_mini[nRECtrack_MAX];
	Float_t peREC_mini[nRECtrack_MAX];
	Float_t etaREC_mini[nRECtrack_MAX];
	Float_t phiREC_mini[nRECtrack_MAX];
	Int_t   bestMatchBSTrack_mini[nRECtrack_MAX]; 

	Float_t ptStarREC_mini[nRECtrack_MAX];
	Float_t etaStarREC_mini[nRECtrack_MAX];
	Float_t phiStarREC_mini[nRECtrack_MAX];

	Float_t nucliaREC_mini[nRECtrack_MAX];
	Float_t dmatchREC_mini[nRECtrack_MAX];
	Int_t imatchREC_mini[nRECtrack_MAX];
	Int_t passREC_mini[nRECtrack_MAX];
	Int_t passTightREC_mini[nRECtrack_MAX];
   	Int_t passLooseREC_mini[nRECtrack_MAX];
	
	tree->SetBranchAddress("eventpass_mini",&eventpass_mini);
	tree->SetBranchAddress("nRECtrack_mini",&nRECtrack_mini);
	tree->SetBranchAddress("w_mini",&w_mini);
	tree->SetBranchAddress("EpzREC_mini",&EpzREC_mini);
	tree->SetBranchAddress("totalMultREC_mini",&totalMultREC_mini);
	tree->SetBranchAddress("vertex_mini",&vertex_mini);
	tree->SetBranchAddress("hfsEREC_mini",&hfsEREC_mini);
	tree->SetBranchAddress("hfsPtREC_mini",&hfsPtREC_mini);
	tree->SetBranchAddress("hfsPzREC_mini",&hfsPzREC_mini);
	tree->SetBranchAddress("elecEREC_mini",&elecEREC_mini);
	tree->SetBranchAddress("elecPtREC_mini",&elecPtREC_mini);
	tree->SetBranchAddress("elecPzREC_mini",&elecPzREC_mini);

	tree->SetBranchAddress("pxREC_mini",&pxREC_mini);
	tree->SetBranchAddress("pyREC_mini",&pyREC_mini);
	tree->SetBranchAddress("pzREC_mini",&pzREC_mini);

	tree->SetBranchAddress("etaREC_mini",&etaREC_mini);
	tree->SetBranchAddress("passREC_mini",passREC_mini);
	tree->SetBranchAddress("passTightREC_mini",passTightREC_mini);
	tree->SetBranchAddress("passLooseREC_mini",passLooseREC_mini);

	tree->SetBranchAddress("xREC_es_mini",&xREC_es_mini);
	tree->SetBranchAddress("yREC_es_mini",&yREC_es_mini);
	tree->SetBranchAddress("Q2REC_es_mini",&Q2REC_es_mini);
	
	/*
	eta distribution in different x. 
	*/
	TH1D* h_eta[4];
	for(int j=0;j<4;j++){
		h_eta[j]= new TH1D(Form("h_eta_%d",j),Form("h_eta_%d",j),100,-3,3);
	}
	/*
	P(n) distribution in different Q2, y, and eta bins.
	*/
	TH1D* h_Pn[4][4][4];
	TH1D* h_Pn_GEN[4][4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				h_Pn[i][j][k] = new TH1D(Form("h_Pn_%d_%d_%d",i,j,k),Form("h_Pn_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_GEN[i][j][k] = new TH1D(Form("h_Pn_GEN_%d_%d_%d",i,j,k),Form("h_Pn_GEN_%d_%d_%d",i,j,k),30,0,30);
			}
		}
	}
	/*
	E-pz and vertex control plots for different multiplicities
	*/
	TH1D* h_Epz[4][4][2];
	TH1D* h_zvertex[4][4][2];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<2;k++){
				h_Epz[i][j][k] = new TH1D(Form("h_Epz_%d_%d_%d",i,j,k),Form("h_Epz_%d_%d_%d",i,j,k),400,0,200);
				h_zvertex[i][j][k] = new TH1D(Form("h_zvertex_%d_%d_%d",i,j,k),Form("h_zvertex_%d_%d_%d",i,j,k),200,-40,40);
			}
		}
	}

	/*
	Elec and hfs kinematic variables
	*/

	TH1D* h_PtBal = new TH1D("h_PtBal","h_PtBal",200,-5,5);
	TH1D* h_hfsEnergy = new TH1D("h_hfsEnergy","h_hfsEnergy",500,0,300);
	TH1D* h_hfsPt = new TH1D("h_hfsPt","h_hfsPt",200,0,25);
	TH1D* h_hfsPz = new TH1D("h_hfsPz","h_hfsPz",500,-50,300);
	TH1D* h_elecEnergy = new TH1D("h_elecEnergy","h_elecEnergy",200,0,100);
	TH1D* h_elecPt = new TH1D("h_elecPt","h_elecPt",200,0,25);
	TH1D* h_elecPz = new TH1D("h_elecPz","h_elecPz",200,-100,25);
	TH1D* h_elecEpz = new TH1D("h_elecEpz","h_elecEpz",500,0,100);
	TH1D* h_hfsEpz = new TH1D("h_hfsEpz","h_hfsEpz",500,0,100);

	cout << "==== Total number of events ~ " << tree->GetEntries() << " ========" << endl;
	for(int ievent = 0; ievent < tree->GetEntries(); ievent++){

		tree->GetEntry(ievent);
		
		if( !isdata_ ){
			int Q2_INDEX = -1;
			for(int iQ2 = 0; iQ2 < 4; iQ2++){
				if(Q2MC_es_mini > Q2_bins[iQ2] && Q2MC_es_mini < Q2_bins[iQ2+1] ){
					Q2_INDEX = iQ2;
				}
			}
			int y_INDEX = -1;
			for(int iy = 0; iy < 4; iy++){
				if(yMC_es_mini > ybins[iy] && yMC_es_mini < ybins[iy+1] ){
					y_INDEX = iy;
				}
			}
			int n_particle_eta[4] = {0,0,0,0};
			for(int itrk = 0; itrk < nMCtrack_mini; itrk++){
				if( TMath::Hypot(pxMC_mini[itrk],pyMC_mini[itrk]) < 0.15 ) continue;
				for(int ieta = 0; ieta < 3; ieta++){
					if( etaMC_mini[itrk] > eta_bins[2*ieta] && etaMC_mini[itrk] < eta_bins[2*ieta+1] ){
						n_particle_eta[ieta]++;
					}
				}
				if( etaMC_mini[itrk] > -1.2 && etaMC_mini[itrk] < 1.6 ){
					n_particle_eta[3]++;
				}
			}
			if(Q2_INDEX >=0 && y_INDEX >= 0){
				h_Pn_GEN[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta[0], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta[1], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta[2], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta[3], w_mini );
			}
		}

		if( Q2REC_es_mini > Q2max || Q2REC_es_mini < Q2min ) continue;
		if( yREC_es_mini > ymax || yREC_es_mini < ymin ) continue;
		if( eventpass_mini != 1 ) continue;

		h_Q2vsX_1->Fill( xREC_es_mini, Q2REC_es_mini, w_mini);
		h_vtxZ_1->Fill( vertex_mini[2], w_mini );
		
		h_PtBal->Fill( hfsPtREC_mini-elecPtREC_mini, w_mini);
		h_hfsEnergy->Fill( hfsEREC_mini, w_mini);
		h_hfsPt->Fill( hfsPtREC_mini, w_mini);
		h_hfsPz->Fill( hfsPzREC_mini, w_mini);
		h_elecEnergy->Fill( elecEREC_mini, w_mini);
		h_elecPt->Fill( elecPtREC_mini, w_mini);
		h_elecPz->Fill( elecPzREC_mini, w_mini);
		h_hfsEpz->Fill( hfsEREC_mini - hfsPzREC_mini, w_mini);
		h_elecEpz->Fill( elecEREC_mini - elecPzREC_mini, w_mini);

		int Q2_INDEX = -1;
		for(int iQ2 = 0; iQ2 < 4; iQ2++){
			if(Q2REC_es_mini > Q2_bins[iQ2] && Q2REC_es_mini < Q2_bins[iQ2+1] ){
				Q2_INDEX = iQ2;
			}
		}
		int y_INDEX = -1;
		for(int iy = 0; iy < 4; iy++){
			if(yREC_es_mini > ybins[iy] && yREC_es_mini < ybins[iy+1] ){
				y_INDEX = iy;
			}
		}
		int n_particle_eta[4] = {0,0,0,0};
		for(int itrk = 0; itrk < nRECtrack_mini; itrk++){
			
			if(passREC_mini[itrk]!=1) continue;
			for(int ieta = 0; ieta < 3; ieta++){
				if( etaREC_mini[itrk] > eta_bins[2*ieta] && etaREC_mini[itrk] < eta_bins[2*ieta+1] ){
					n_particle_eta[ieta]++;
				}
			}
			if( etaREC_mini[itrk] > -1.2 && etaREC_mini[itrk] < 1.6 ) n_particle_eta[3]++;

			h_eta[0]->Fill( etaREC_mini[itrk], w_mini );
			for(int iy=0;iy<3;iy++){
				if(yREC_es_mini>ybins[iy] && yREC_es_mini<ybins[iy+1]){
					h_eta[iy+1]->Fill( etaREC_mini[itrk], w_mini );
				}
			}
		}

		if(Q2_INDEX >=0 && y_INDEX >= 0){
			h_Pn[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta[0], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta[1], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta[2], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta[3], w_mini );
			
			if( totalMultREC_mini < 15 && totalMultREC_mini >= 0 ){
				h_Epz[Q2_INDEX][y_INDEX][0]->Fill( EpzREC_mini );
				h_zvertex[Q2_INDEX][y_INDEX][0]->Fill( vertex_mini[2] );
			}
			if( totalMultREC_mini > 15 ){
				h_Epz[Q2_INDEX][y_INDEX][1]->Fill( EpzREC_mini );
				h_zvertex[Q2_INDEX][y_INDEX][1]->Fill( vertex_mini[2] );
			}
			
		}
		
	}

	TString outname;
	if( isdata_ ){
		outname = "../minitree_output/Pn_hist_data.root";
	}
	else if( isdjango_){
		outname = "../minitree_output/Pn_hist_django.root";
	}
	else{
		outname = "../minitree_output/Pn_hist_rapgap.root";
	}
	TFile f1(outname,"RECREATE");

	for(int j=0;j<4;j++){
		h_eta[j]->Write();
	}
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				h_Pn[i][j][k]->Write();
				h_Pn_GEN[i][j][k]->Write();
			}
		}
	}
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<2;k++){
				h_Epz[i][j][k]->Write();
				h_zvertex[i][j][k]->Write();
			}
		}
	}

	h_Q2vsX_1->Write();
	h_vtxZ_1->Write();
	h_PtBal->Write();
	h_hfsEnergy->Write();
	h_hfsPt->Write();
	h_hfsPz->Write();
	h_elecEnergy->Write();
	h_elecPt->Write();
	h_elecPz->Write();
	h_hfsEpz->Write();
	h_elecEpz->Write();


}