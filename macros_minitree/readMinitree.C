#include "../header/mainAnalysis.h"

using namespace std;
void readMinitree(const int ifile_ = 0, const bool isReweigh = false, const bool isScatElec_ = false){

	TFile* file = 0;

	if( !isScatElec_ ){
		if(ifile_ == 0){
			file = new TFile("../new_output/data_highE_0607_noReweight_Tree_hadCali.root");
		}
		else if(ifile_ == 1) {

			if(!isReweigh) file = new TFile("../new_output/mc_highE_DJANGOH_noReweight_Tree_hadCali.root");
			else file = new TFile("../new_output/mc_highE_DJANGOH_fullReweight_Tree_hadCali.root");
		}
		else if(ifile_ == 2){
			if(!isReweigh) file = new TFile("../new_output/mc_highE_RAPGAP_noReweight_Tree_hadCali.root");
			else file = new TFile("../new_output/mc_highE_RAPGAP_fullReweight_Tree_hadCali_diffractive10per.root");
		}
	}
	else{
		if(ifile_ == 0){
			file = new TFile("../new_output/data_highE_0607_noReweight_Tree_scatElec.root");
		}
		else if(ifile_ == 1) {

			if(!isReweigh) file = new TFile("../new_output/mc_highE_DJANGOH_noReweight_Tree_scatElec.root");
			else file = new TFile("../new_output/mc_highE_DJANGOH_fullReweight_Tree_scatElec.root");
		}
		else if(ifile_ == 2){
			if(!isReweigh) file = new TFile("../new_output/mc_highE_RAPGAP_noReweight_Tree_scatElec.root");
			else file = new TFile("../new_output/mc_highE_RAPGAP_fullReweight_Tree_scatElec.root");
		}
	}
	TTree* tree = (TTree*) file->Get("miniTree");

	double Q2min = 5.;
	double Q2max = 100.;
	double ymin = 0.075;
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
	Float_t etaStarMC_mini[nMCtrack_MAX];
	Float_t chargeMC_mini[nMCtrack_MAX];
   	
	tree->SetBranchAddress("yMC_es_mini",&yMC_es_mini);
	tree->SetBranchAddress("Q2MC_es_mini",&Q2MC_es_mini);

	tree->SetBranchAddress("nMCtrack_mini",&nMCtrack_mini);
	tree->SetBranchAddress("pxMC_mini",&pxMC_mini);
	tree->SetBranchAddress("pyMC_mini",&pyMC_mini);
	tree->SetBranchAddress("pzMC_mini",&pzMC_mini);
	tree->SetBranchAddress("etaMC_mini",&etaMC_mini);
	tree->SetBranchAddress("etaStarMC_mini",&etaStarMC_mini);

	Float_t xREC_es_mini,yREC_es_mini,Q2REC_es_mini;
	const int nRECtrack_MAX = 200;
	Int_t eventpass_mini;
	Int_t eventpassTight_mini;
 	Int_t eventpassLoose_mini;
	Int_t nRECtrack_mini;
	Int_t typeChgREC_mini[nRECtrack_MAX];
	Float_t etaAsymREC_mini;
	Float_t EpzREC_mini;
	Int_t totalMultREC_mini;
	Float_t vertex_mini[3];
	Float_t hfsEREC_mini, hfsPtREC_mini, hfsPzREC_mini;
 	Float_t elecEREC_mini, elecPtREC_mini, elecPzREC_mini;

	Int_t elecChargeREC_mini;
	Float_t eElectronBeam_mini;

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
	tree->SetBranchAddress("eventpassTight_mini",&eventpassTight_mini);
	tree->SetBranchAddress("eventpassLoose_mini",&eventpassLoose_mini);
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
	
	tree->SetBranchAddress("eElectronBeam_mini",&eElectronBeam_mini);
	
	if( isScatElec_ ){
		tree->SetBranchAddress("elecChargeREC_mini",&elecChargeREC_mini);
	}
	
	tree->SetBranchAddress("pxREC_mini",&pxREC_mini);
	tree->SetBranchAddress("pyREC_mini",&pyREC_mini);
	tree->SetBranchAddress("pzREC_mini",&pzREC_mini);

	tree->SetBranchAddress("etaREC_mini",&etaREC_mini);
	tree->SetBranchAddress("etaStarREC_mini",&etaStarREC_mini);
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
	TH2D* h_Pn_cor[4][4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				h_Pn[i][j][k] = new TH1D(Form("h_Pn_%d_%d_%d",i,j,k),Form("h_Pn_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_GEN[i][j][k] = new TH1D(Form("h_Pn_GEN_%d_%d_%d",i,j,k),Form("h_Pn_GEN_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_cor[i][j][k] = new TH2D(Form("h_Pn_cor_%d_%d_%d",i,j,k),Form("h_Pn_cor_%d_%d_%d",i,j,k),30,0,30,30,0,30);
			}
		}
	}
	/*
	P(n) distribution in different Q2, y bins within 0<eta*<4.
	*/
	TH1D* h_Pn_HCM[4][4];
	TH1D* h_Pn_GEN_HCM[4][4];
	TH2D* h_Pn_cor_HCM[4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			h_Pn_HCM[i][j] = new TH1D(Form("h_Pn_HCM_%d_%d",i,j),Form("h_Pn_HCM_%d_%d",i,j),30,0,30);
			h_Pn_GEN_HCM[i][j] = new TH1D(Form("h_Pn_GEN_HCM_%d_%d",i,j),Form("h_Pn_GEN_HCM_%d_%d",i,j),30,0,30);
			h_Pn_cor_HCM[i][j] = new TH2D(Form("h_Pn_cor_HCM_%d_%d",i,j),Form("h_Pn_cor_HCM_%d_%d",i,j),30,0,30,30,0,30);
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

	TH2D* h_Q2vsX_1 = new TH2D("h_Q2vsX_1","h_Q2vsX_1",1000,0.00001,0.01,50,1,100);
	TH1D* h_vtxZ_1  = new TH1D("h_vtxZ_1","h_vtxZ_1", 100,-50,50);
	TH1D* h_y_1  = new TH1D("h_y_1","h_y_1", 200,0,1);

	TH2D* h_Q2vsX_2 = new TH2D("h_Q2vsX_2","h_Q2vsX_2",1000,0.00001,0.01,50,1,100);
	TH1D* h_vtxZ_2  = new TH1D("h_vtxZ_2","h_vtxZ_2", 100,-50,50);
	TH1D* h_y_2  = new TH1D("h_y_2","h_y_2", 200,0,1);

	/*
	Elec and hfs kinematic variables
	*/
	TH1D* h_etamax_all = new TH1D("h_etamax_all","h_etamax_all",200,-3,3);
	TH1D* h_PtBal_all = new TH1D("h_PtBal","h_PtBal",200,-5,5);
	TH1D* h_hfsEnergy_all = new TH1D("h_hfsEnergy","h_hfsEnergy",500,0,300);
	TH1D* h_hfsPt_all = new TH1D("h_hfsPt","h_hfsPt",200,0,25);
	TH1D* h_hfsPz_all = new TH1D("h_hfsPz","h_hfsPz",500,-50,300);
	TH1D* h_elecEnergy_all = new TH1D("h_elecEnergy","h_elecEnergy",200,0,100);
	TH1D* h_elecPt_all = new TH1D("h_elecPt","h_elecPt",200,0,25);
	TH1D* h_elecPz_all = new TH1D("h_elecPz","h_elecPz",200,-100,25);
	TH1D* h_elecEpz_all = new TH1D("h_elecEpz","h_elecEpz",500,0,100);
	TH1D* h_hfsEpz_all = new TH1D("h_hfsEpz","h_hfsEpz",500,0,100);
	TH1D* h_hfsElecEpz_all = new TH1D("h_hfsElecEpz","h_hfsElecEpz",500,0,100);

	TH1D* h_elecPlusEnergy_all = new TH1D("h_elecPlusEnergy","h_elecPlusEnergy",200,0,100);
	TH1D* h_elecPlusPt_all = new TH1D("h_elecPlusPt","h_elecPlusPt",200,0,25);
	TH1D* h_elecPlusPz_all = new TH1D("h_elecPlusPz","h_elecPlusPz",200,-100,25);
	TH1D* h_elecPlusEpz_all = new TH1D("h_elecPlusEpz","h_elecPlusEpz",500,0,100);
	TH1D* h_hfsElecPlusEpz_all = new TH1D("h_hfsElecPlusEpz","h_hfsElecPlusEpz",500,0,100);

	TH1D* h_elecMinusEnergy_all = new TH1D("h_elecMinusEnergy","h_elecMinusEnergy",200,0,100);
	TH1D* h_elecMinusPt_all = new TH1D("h_elecMinusPt","h_elecMinusPt",200,0,25);
	TH1D* h_elecMinusPz_all = new TH1D("h_elecMinusPz","h_elecMinusPz",200,-100,25);
	TH1D* h_elecMinusEpz_all = new TH1D("h_elecMinusEpz","h_elecMinusEpz",500,0,100);
	TH1D* h_hfsElecMinusEpz_all = new TH1D("h_hfsElecMinusEpz","h_hfsElecMinusEpz",500,0,100);

	TH1D* h_neutEnergy_all = new TH1D("h_neutEnergy","h_neutEnergy",200,0,100);
	TH1D* h_neutPt_all = new TH1D("h_neutPt","h_neutPt",200,0,25);
	TH1D* h_neutPz_all = new TH1D("h_neutPz","h_neutPz",200,-100,25);
	TH1D* h_neutEpz_all = new TH1D("h_neutEpz","h_neutEpz",500,0,100);
	TH1D* h_hfsNeutEpz_all = new TH1D("h_hfsNeutEpz","h_hfsNeutEpz",500,0,100);

	TH1D* h_trackEpz_all = new TH1D("h_trackEpz","h_trackEpz",500,0,100);

	TH1D* h_ptVsMult[30];
	for(int i=0;i<30;i++){
		h_ptVsMult[i] = new TH1D(Form("h_ptVsMult_%d",i),Form("h_ptVsMult_%d",i),100,0,5);
	}
	TH1D* h_etamax[4][4][2];
	TH1D* h_PtBal[4][4][2];
	TH1D* h_hfsEnergy[4][4][2];
	TH1D* h_hfsPt[4][4][2];
	TH1D* h_hfsPz[4][4][2];
	TH1D* h_elecEnergy[4][4][2];
	TH1D* h_elecPt[4][4][2];
	TH1D* h_elecPz[4][4][2];
	TH1D* h_elecEpz[4][4][2];
	TH1D* h_hfsEpz[4][4][2];
	TH1D* h_hfsElecEpz[4][4][2];

	TH1D* h_elecPlusEnergy[4][4][2];
	TH1D* h_elecPlusPt[4][4][2];
	TH1D* h_elecPlusPz[4][4][2];
	TH1D* h_elecPlusEpz[4][4][2];
	TH1D* h_hfsElecPlusEpz[4][4][2];

	TH1D* h_elecMinusEnergy[4][4][2];
	TH1D* h_elecMinusPt[4][4][2];
	TH1D* h_elecMinusPz[4][4][2];
	TH1D* h_elecMinusEpz[4][4][2];
	TH1D* h_hfsElecMinusEpz[4][4][2];

	TH1D* h_neutEnergy[4][4][2];
	TH1D* h_neutPt[4][4][2];
	TH1D* h_neutPz[4][4][2];
	TH1D* h_neutEpz[4][4][2];
	TH1D* h_hfsNeutEpz[4][4][2];

	TH1D* h_trackEpz[4][4][2];

	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 2; k++){

				h_etamax[i][j][k] = new TH1D(Form("h_etamax_%d_%d_%d",i,j,k),"h_etamax",200,-3,3);
				h_PtBal[i][j][k] = new TH1D(Form("h_PtBal_%d_%d_%d",i,j,k),"h_PtBal",200,-5,5);
				h_hfsEnergy[i][j][k] = new TH1D(Form("h_hfsEnergy_%d_%d_%d",i,j,k),"h_hfsEnergy",500,0,300);
				h_hfsPt[i][j][k] = new TH1D(Form("h_hfsPt_%d_%d_%d",i,j,k),"h_hfsPt",200,0,25);
				h_hfsPz[i][j][k] = new TH1D(Form("h_hfsPz_%d_%d_%d",i,j,k),"h_hfsPz",500,-50,300);
				h_elecEnergy[i][j][k] = new TH1D(Form("h_elecEnergy_%d_%d_%d",i,j,k),"h_elecEnergy",200,0,100);
				h_elecPt[i][j][k] = new TH1D(Form("h_elecPt_%d_%d_%d",i,j,k),"h_elecPt",200,0,25);
				h_elecPz[i][j][k] = new TH1D(Form("h_elecPz_%d_%d_%d",i,j,k),"h_elecPz",200,-100,25);
				h_elecEpz[i][j][k] = new TH1D(Form("h_elecEpz_%d_%d_%d",i,j,k),"h_elecEpz",500,0,100);
				h_hfsEpz[i][j][k] = new TH1D(Form("h_hfsEpz_%d_%d_%d",i,j,k),"h_hfsEpz",500,0,100);
				h_hfsElecEpz[i][j][k] = new TH1D(Form("h_hfsElecEpz_%d_%d_%d",i,j,k),"h_hfsElecEpz",500,0,100);

				h_elecPlusEnergy[i][j][k] = new TH1D(Form("h_elecPlusEnergy_%d_%d_%d",i,j,k),"h_elecPlusEnergy",200,0,100);
				h_elecPlusPt[i][j][k] = new TH1D(Form("h_elecPlusPt_%d_%d_%d",i,j,k),"h_elecPlusPt",200,0,25);
				h_elecPlusPz[i][j][k] = new TH1D(Form("h_elecPlusPz_%d_%d_%d",i,j,k),"h_elecPlusPz",200,-100,25);
				h_elecPlusEpz[i][j][k] = new TH1D(Form("h_elecPlusEpz_%d_%d_%d",i,j,k),"h_elecPlusEpz",500,0,100);
				h_hfsElecPlusEpz[i][j][k] = new TH1D(Form("h_hfsElecPlusEpz_%d_%d_%d",i,j,k),"h_hfsElecPlusEpz",500,0,100);

				h_elecMinusEnergy[i][j][k] = new TH1D(Form("h_elecMinusEnergy_%d_%d_%d",i,j,k),"h_elecMinusEnergy",200,0,100);
				h_elecMinusPt[i][j][k] = new TH1D(Form("h_elecMinusPt_%d_%d_%d",i,j,k),"h_elecMinusPt",200,0,25);
				h_elecMinusPz[i][j][k] = new TH1D(Form("h_elecMinusPz_%d_%d_%d",i,j,k),"h_elecMinusPz",200,-100,25);
				h_elecMinusEpz[i][j][k] = new TH1D(Form("h_elecMinusEpz_%d_%d_%d",i,j,k),"h_elecMinusEpz",500,0,100);
				h_hfsElecMinusEpz[i][j][k] = new TH1D(Form("h_hfsElecMinusEpz_%d_%d_%d",i,j,k),"h_hfsElecMinusEpz",500,0,100);

				h_neutEnergy[i][j][k] = new TH1D(Form("h_neutEnergy_%d_%d_%d",i,j,k),"h_neutEnergy",200,0,100);
				h_neutPt[i][j][k] = new TH1D(Form("h_neutPt_%d_%d_%d",i,j,k),"h_neutPt",200,0,25);
				h_neutPz[i][j][k] = new TH1D(Form("h_neutPz_%d_%d_%d",i,j,k),"h_neutPz",200,-100,25);
				h_neutEpz[i][j][k] = new TH1D(Form("h_neutEpz_%d_%d_%d",i,j,k),"h_neutEpz",500,0,100);
				h_hfsNeutEpz[i][j][k] = new TH1D(Form("h_hfsNeutEpz_%d_%d_%d",i,j,k),"h_hfsNeutEpz",500,0,100);
			
				h_trackEpz[i][j][k] = new TH1D(Form("h_trackEpz_%d_%d_%d",i,j,k),"h_trackEpz",500,0,100);

			}
		}
	}
	cout << "==== Total number of events ~ " << tree->GetEntries() << " ========" << endl;
	for(int ievent = 0; ievent < tree->GetEntries(); ievent++){

		tree->GetEntry(ievent);
		
		int n_particle_eta[4] = {0,0,0,0};
		int n_particle_HCM = 0;
		if( ifile_ != 0 ){
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
			for(int itrk = 0; itrk < nMCtrack_mini; itrk++){
				if( TMath::Hypot(pxMC_mini[itrk],pyMC_mini[itrk]) < 0.15 ) continue;
				if( etaStarMC_mini[itrk] > 0 && etaStarMC_mini[itrk] < 4.0 ){
					n_particle_HCM++;
				}
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
				h_Pn_GEN_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM, w_mini );
				
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
		h_y_1->Fill( yREC_es_mini, w_mini);
		h_vtxZ_1->Fill( vertex_mini[2], w_mini );

	   //setting event Q2 and y indices
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
		int n_particle_eta_rec[4] = {0,0,0,0};
		int n_particle_HCM_rec = 0;
		double trk_E = 0.;
		double trk_pz = 0.;
		double etamax = -99.;
		for(int itrk = 0; itrk < nRECtrack_mini; itrk++){
			
			if(passREC_mini[itrk]!=1) continue;
			if( etaREC_mini[itrk] > etamax ){
				etamax = etaREC_mini[itrk];
			}
			if( etaStarREC_mini[itrk] > 0 && etaStarREC_mini[itrk] < 4.0 ){
				n_particle_HCM_rec++;
			}
			for(int ieta = 0; ieta < 3; ieta++){
				if( etaREC_mini[itrk] > eta_bins[2*ieta] && etaREC_mini[itrk] < eta_bins[2*ieta+1] ){
					n_particle_eta_rec[ieta]++;
				}
			}
			if( etaREC_mini[itrk] > -1.2 && etaREC_mini[itrk] < 1.6 ) {
				n_particle_eta_rec[3]++;
				trk_E += sqrt(pxREC_mini[itrk]*pxREC_mini[itrk] + pyREC_mini[itrk]*pyREC_mini[itrk] + pzREC_mini[itrk]*pzREC_mini[itrk] + 0.134*0.134);
				trk_pz += pzREC_mini[itrk];
			}

			h_eta[0]->Fill( etaREC_mini[itrk], w_mini );
			for(int iy=0;iy<3;iy++){
				if(yREC_es_mini>ybins[iy] && yREC_es_mini<ybins[iy+1]){
					h_eta[iy+1]->Fill( etaREC_mini[itrk], w_mini );
				}
			}
		}

		for(int im=0;im<30;im++){
			if( im != n_particle_eta_rec[3]) continue;
			for(int itrk = 0; itrk < nRECtrack_mini; itrk++){
				if(passREC_mini[itrk]!=1) continue;
				if( etaREC_mini[itrk] > -1.2 && etaREC_mini[itrk] < 1.6 ) {
					h_ptVsMult[im]->Fill(TMath::Hypot(pxREC_mini[itrk],pyREC_mini[itrk]), w_mini);
				}
			}
		}
		

		double eBeamEnergy = eElectronBeam_mini;
		if( ifile_ != 0 ) eBeamEnergy = 27.6;

		h_etamax_all->Fill( etamax, w_mini);
		h_PtBal_all->Fill( hfsPtREC_mini/elecPtREC_mini, w_mini);
		h_hfsEnergy_all->Fill( hfsEREC_mini, w_mini);
		h_hfsPt_all->Fill( hfsPtREC_mini, w_mini);
		h_hfsPz_all->Fill( hfsPzREC_mini, w_mini);
		h_elecEnergy_all->Fill( elecEREC_mini, w_mini);
		h_elecPt_all->Fill( elecPtREC_mini, w_mini);
		h_elecPz_all->Fill( elecPzREC_mini, w_mini);
		h_hfsEpz_all->Fill( hfsEREC_mini - hfsPzREC_mini, w_mini);
		h_elecEpz_all->Fill( (27.6/eBeamEnergy)*(elecEREC_mini - elecPzREC_mini), w_mini);
		double Epz = (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini);
		Epz = (27.6/eBeamEnergy)*Epz;
		h_hfsElecEpz_all->Fill( Epz, w_mini );

		if( isScatElec_ ){
			if( elecChargeREC_mini == +1 ){
				h_elecPlusEnergy_all->Fill( elecEREC_mini, w_mini);
				h_elecPlusPt_all->Fill( elecPtREC_mini, w_mini);
				h_elecPlusPz_all->Fill( elecPzREC_mini, w_mini);
				h_elecPlusEpz_all->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
				h_hfsElecPlusEpz_all->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);
			}
			else if (elecChargeREC_mini == -1){
				h_elecMinusEnergy_all->Fill( elecEREC_mini, w_mini);
				h_elecMinusPt_all->Fill( elecPtREC_mini, w_mini);
				h_elecMinusPz_all->Fill( elecPzREC_mini, w_mini);
				h_elecMinusEpz_all->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
				h_hfsElecMinusEpz_all->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

			}
			else{
				h_neutEnergy_all->Fill( elecEREC_mini, w_mini);
				h_neutPt_all->Fill( elecPtREC_mini, w_mini);
				h_neutPz_all->Fill( elecPzREC_mini, w_mini);
				h_neutEpz_all->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
				h_hfsNeutEpz_all->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

			}
		}
		
		h_trackEpz_all->Fill(trk_E-trk_pz, w_mini);

		if(Q2_INDEX >=0 && y_INDEX >= 0){

			if( ifile_ != 0 ){

				h_Pn_cor_HCM[Q2_INDEX][y_INDEX]->Fill(n_particle_HCM_rec,n_particle_HCM,w_mini);
				h_Pn_cor[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta_rec[0], n_particle_eta[0], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta_rec[1], n_particle_eta[1], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta_rec[2], n_particle_eta[2], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta_rec[3], n_particle_eta[3], w_mini );
			}

			h_Pn_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_rec, w_mini );

			h_Pn[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta_rec[0], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta_rec[1], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta_rec[2], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta_rec[3], w_mini );
			
			if( totalMultREC_mini < 15 && totalMultREC_mini >= 0 ){
				h_Epz[Q2_INDEX][y_INDEX][0]->Fill( (hfsEREC_mini+elecEREC_mini) - (hfsPzREC_mini+elecPzREC_mini) , w_mini );
				h_zvertex[Q2_INDEX][y_INDEX][0]->Fill( vertex_mini[2],w_mini );

				h_etamax[Q2_INDEX][y_INDEX][0]->Fill( etamax, w_mini);
				h_PtBal[Q2_INDEX][y_INDEX][0]->Fill( hfsPtREC_mini/elecPtREC_mini, w_mini);
				h_hfsEnergy[Q2_INDEX][y_INDEX][0]->Fill( hfsEREC_mini, w_mini);
				h_hfsPt[Q2_INDEX][y_INDEX][0]->Fill( hfsPtREC_mini, w_mini);
				h_hfsPz[Q2_INDEX][y_INDEX][0]->Fill( hfsPzREC_mini, w_mini);
				h_elecEnergy[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini, w_mini);
				h_elecPt[Q2_INDEX][y_INDEX][0]->Fill( elecPtREC_mini, w_mini);
				h_elecPz[Q2_INDEX][y_INDEX][0]->Fill( elecPzREC_mini, w_mini);
				h_hfsEpz[Q2_INDEX][y_INDEX][0]->Fill( hfsEREC_mini - hfsPzREC_mini, w_mini);
				h_elecEpz[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
				h_hfsElecEpz[Q2_INDEX][y_INDEX][0]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini );


				if( isScatElec_ ){
					if( elecChargeREC_mini == +1 ){
						h_elecPlusEnergy[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini, w_mini);
						h_elecPlusPt[Q2_INDEX][y_INDEX][0]->Fill( elecPtREC_mini, w_mini);
						h_elecPlusPz[Q2_INDEX][y_INDEX][0]->Fill( elecPzREC_mini, w_mini);
						h_elecPlusEpz[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsElecPlusEpz[Q2_INDEX][y_INDEX][0]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);
					}
					else if (elecChargeREC_mini == -1){
						h_elecMinusEnergy[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini, w_mini);
						h_elecMinusPt[Q2_INDEX][y_INDEX][0]->Fill( elecPtREC_mini, w_mini);
						h_elecMinusPz[Q2_INDEX][y_INDEX][0]->Fill( elecPzREC_mini, w_mini);
						h_elecMinusEpz[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsElecMinusEpz[Q2_INDEX][y_INDEX][0]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

					}
					else{
						h_neutEnergy[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini, w_mini);
						h_neutPt[Q2_INDEX][y_INDEX][0]->Fill( elecPtREC_mini, w_mini);
						h_neutPz[Q2_INDEX][y_INDEX][0]->Fill( elecPzREC_mini, w_mini);
						h_neutEpz[Q2_INDEX][y_INDEX][0]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsNeutEpz[Q2_INDEX][y_INDEX][0]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

					}
				}

				h_trackEpz[Q2_INDEX][y_INDEX][0]->Fill(trk_E-trk_pz, w_mini);


			}
			if( totalMultREC_mini > 15 ){
				h_Epz[Q2_INDEX][y_INDEX][1]->Fill( (hfsEREC_mini+elecEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini  );
				h_zvertex[Q2_INDEX][y_INDEX][1]->Fill( vertex_mini[2],w_mini );

				h_etamax[Q2_INDEX][y_INDEX][1]->Fill( etamax, w_mini);
				h_PtBal[Q2_INDEX][y_INDEX][1]->Fill( hfsPtREC_mini/elecPtREC_mini, w_mini);
				h_hfsEnergy[Q2_INDEX][y_INDEX][1]->Fill( hfsEREC_mini, w_mini);
				h_hfsPt[Q2_INDEX][y_INDEX][1]->Fill( hfsPtREC_mini, w_mini);
				h_hfsPz[Q2_INDEX][y_INDEX][1]->Fill( hfsPzREC_mini, w_mini);
				h_elecEnergy[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini, w_mini);
				h_elecPt[Q2_INDEX][y_INDEX][1]->Fill( elecPtREC_mini, w_mini);
				h_elecPz[Q2_INDEX][y_INDEX][1]->Fill( elecPzREC_mini, w_mini);
				h_hfsEpz[Q2_INDEX][y_INDEX][1]->Fill( hfsEREC_mini - hfsPzREC_mini, w_mini);
				h_elecEpz[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
				h_hfsElecEpz[Q2_INDEX][y_INDEX][1]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini );


				if( isScatElec_ ){
					if( elecChargeREC_mini == +1 ){
						h_elecPlusEnergy[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini, w_mini);
						h_elecPlusPt[Q2_INDEX][y_INDEX][1]->Fill( elecPtREC_mini, w_mini);
						h_elecPlusPz[Q2_INDEX][y_INDEX][1]->Fill( elecPzREC_mini, w_mini);
						h_elecPlusEpz[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsElecPlusEpz[Q2_INDEX][y_INDEX][1]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);
					}
					else if (elecChargeREC_mini == -1){
						h_elecMinusEnergy[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini, w_mini);
						h_elecMinusPt[Q2_INDEX][y_INDEX][1]->Fill( elecPtREC_mini, w_mini);
						h_elecMinusPz[Q2_INDEX][y_INDEX][1]->Fill( elecPzREC_mini, w_mini);
						h_elecMinusEpz[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsElecMinusEpz[Q2_INDEX][y_INDEX][1]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

					}
					else{
						h_neutEnergy[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini, w_mini);
						h_neutPt[Q2_INDEX][y_INDEX][1]->Fill( elecPtREC_mini, w_mini);
						h_neutPz[Q2_INDEX][y_INDEX][1]->Fill( elecPzREC_mini, w_mini);
						h_neutEpz[Q2_INDEX][y_INDEX][1]->Fill( elecEREC_mini - elecPzREC_mini, w_mini);
						h_hfsNeutEpz[Q2_INDEX][y_INDEX][1]->Fill( (elecEREC_mini+hfsEREC_mini) - (hfsPzREC_mini+elecPzREC_mini), w_mini);

					}
				}

				h_trackEpz[Q2_INDEX][y_INDEX][1]->Fill(trk_E-trk_pz, w_mini);
			}
			
		}
		
	}

	TString outname;
	if( !isScatElec_ ){
		if( ifile_ == 0 ){
			outname = "../minitree_output/Pn_hist_data_hadCali.root";
		}
		else if( ifile_ == 1 ){
			if(!isReweigh) outname = "../minitree_output/Pn_hist_django_hadCali.root";
			else outname = "../minitree_output/Pn_hist_django_hadCali_reweigh.root";
		}
		else if( ifile_ == 2 ){
			if(!isReweigh) outname = "../minitree_output/Pn_hist_rapgap_hadCali.root";
			else outname = "../minitree_output/Pn_hist_rapgap_hadCali_reweigh_diffractive10per.root";
		}
	}
	else{
		if( ifile_ == 0 ){
			outname = "../minitree_output/Pn_hist_data_scatElec.root";
		}
		else if( ifile_ == 1 ){
			if(!isReweigh) outname = "../minitree_output/Pn_hist_django_scatElec.root";
			else outname = "../minitree_output/Pn_hist_django_scatElec_reweigh.root";
		}
		else if( ifile_ == 2 ){
			if(!isReweigh) outname = "../minitree_output/Pn_hist_rapgap_scatElec.root";
			else outname = "../minitree_output/Pn_hist_rapgap_scatElec_reweigh.root";
		}		
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
				h_Pn_cor[i][j][k]->Write();
			}
		}
	}
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			h_Pn_HCM[i][j]->Write();
			h_Pn_GEN_HCM[i][j]->Write();
			h_Pn_cor_HCM[i][j]->Write();
			
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
	h_y_1->Write();
	h_vtxZ_1->Write();
	h_etamax_all->Write();

	h_PtBal_all->Write();
	h_hfsEnergy_all->Write();
	h_hfsPt_all->Write();
	h_hfsPz_all->Write();
	h_elecEnergy_all->Write();
	h_elecPt_all->Write();
	h_elecPz_all->Write();
	h_hfsEpz_all->Write();
	h_elecEpz_all->Write();
	h_hfsElecEpz_all->Write();

	if( isScatElec_ ){
		h_elecPlusEnergy_all->Write();
		h_elecPlusPt_all->Write();
		h_elecPlusPz_all->Write();
		h_elecPlusEpz_all->Write();
		h_hfsElecPlusEpz_all->Write();

		h_elecMinusEnergy_all->Write();
		h_elecMinusPt_all->Write();
		h_elecMinusPz_all->Write();
		h_elecMinusEpz_all->Write();
		h_hfsElecMinusEpz_all->Write();

		h_neutEnergy_all->Write();
		h_neutPt_all->Write();
		h_neutPz_all->Write();
		h_neutEpz_all->Write();
		h_hfsNeutEpz_all->Write();
	}

	h_trackEpz_all->Write();
	

	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<2;k++){
				h_etamax[i][j][k]->Write();
				h_PtBal[i][j][k]->Write();
				h_hfsEnergy[i][j][k]->Write();
				h_hfsPt[i][j][k]->Write();
				h_hfsPz[i][j][k]->Write();
				h_elecEnergy[i][j][k]->Write();
				h_elecPt[i][j][k]->Write();
				h_elecPz[i][j][k]->Write();
				h_hfsEpz[i][j][k]->Write();
				h_elecEpz[i][j][k]->Write();
				h_hfsElecEpz[i][j][k]->Write();

				if( isScatElec_ ){
					h_elecPlusEnergy[i][j][k]->Write();
					h_elecPlusPt[i][j][k]->Write();
					h_elecPlusPz[i][j][k]->Write();
					h_elecPlusEpz[i][j][k]->Write();
					h_hfsElecPlusEpz[i][j][k]->Write();

					h_elecMinusEnergy[i][j][k]->Write();
					h_elecMinusPt[i][j][k]->Write();
					h_elecMinusPz[i][j][k]->Write();
					h_elecMinusEpz[i][j][k]->Write();
					h_hfsElecMinusEpz[i][j][k]->Write();

					h_neutEnergy[i][j][k]->Write();
					h_neutPt[i][j][k]->Write();
					h_neutPz[i][j][k]->Write();
					h_neutEpz[i][j][k]->Write();
					h_hfsNeutEpz[i][j][k]->Write();
				}

				h_trackEpz[i][j][k]->Write();
			}
		}
	}
	
	for(int i=0;i<30;i++){
		h_ptVsMult[i]->Write();
	}


}
