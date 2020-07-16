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

#define PIMASS 0.13957
#define ELECTRON_MASS 0.00051
//binning before was 30,0,30;
double Pn_binning[]={-0.5,0.5,1.5,2.5,3.5,5.5,7.5,11.5,15.5,19.5,25.5,31.5};

using namespace std;
void readMinitree(const int ifile_ = 0, const bool isReweigh = false){

	//input files
	TFile* file = 0;
	if(ifile_ == 0){
		file = new TFile("../new_output/data_highE_0607_noReweight_Tree_hadCaliNewKine_June16.root");
	}else if(ifile_ == 1) {
		if(!isReweigh) file = new TFile("../new_output/mc_highE_DJANGOH_noReweight_Tree_hadCaliNew.root");
		else file = new TFile("../new_output/mc_highE_DJANGOH_fullReweight_Tree_hadCaliNewKine_forPhotonicElectron.root");
	}else if(ifile_ == 2){
		if(!isReweigh) file = new TFile("../new_output/mc_highE_RAPGAP_noReweight_Tree_hadCaliNew.root");
		else file = new TFile("../new_output/mc_highE_RAPGAP_fullReweight_Tree_hadCaliNewKine_forPhotonicElectron.root");
	}
	else if(ifile_ == 3){
		file = new TFile("../new_output/mc_highE_PYTHIA6_noReweight_Tree_hadCaliNewKine_photoproduction.root");
	}
	//output files
	TString outname;
	if( ifile_ == 0 ){
		outname = "../minitree_output/Pn_hist_data_hadCaliNewKine_June16.root";
	}else if( ifile_ == 1 ){
		if(!isReweigh) outname = "../minitree_output/Pn_hist_django_extendEtalabLooseTrack.root";
		else outname = "../minitree_output/Pn_hist_django_hadCaliNewKine_reweigh_forPhotonicElectron.root";
	}else if( ifile_ == 2 ){
		if(!isReweigh) outname = "../minitree_output/Pn_hist_rapgap_extendEtalabLooseTrack.root";
		else outname = "../minitree_output/Pn_hist_rapgap_hadCaliNewKine_reweigh_forPhotonicElectron.root";
	}
	else if( ifile_ == 3 ){
		outname = "../minitree_output/Pn_hist_pythia_hadCaliNewKine_photoproduction.root";
	}
	
	TTree* tree = (TTree*) file->Get("miniTree");

	double Q2min = 5.;
	double Q2max = 100.;
	double ymin = 0.0375;
	double ymax = 0.6;
	double ybins[] = {0.0375,0.075,0.15,0.3,0.6};
	double eta_bins[] = {-1.2,0.2,-0.5,0.9,0.2,1.6};
	double Q2_bins[] = {5,10,20,40,100};
	double electron_likelihood = 0.01;

	TString seta_bins[6]={"-1.2","0.2","-0.5","0.9","0.2","1.6"};
	TString sQ2_bins[5]={"5","10","20","40","100"};
	TString sy_bins[5]={"0.0375","0.075","0.15","0.3","0.6"};

	Float_t xMC_es_mini,yMC_es_mini,Q2MC_es_mini;
	Float_t yMC_mini,Q2MC_mini;
	Float_t w_mini;
	
	const int nMCtrack_MAX=400;
	// if there is no MC info, nMCtrack is set to zero
	Int_t nMCtrack_mini;
	Int_t isQEDcMC_mini;
	Int_t isQEDbkg_mini;
	Int_t isDIFFbkg_mini;
	Float_t dRRadPhot_mini;
	Float_t dPhiRadPhot_mini;
	Float_t elecPxMC_mini;
	Float_t elecPyMC_mini;
	Float_t elecPzMC_mini;
	Float_t elecEMC_mini;
	Float_t phoPxMC_mini[3];
	Float_t phoPyMC_mini[3];
	Float_t phoPzMC_mini[3];
	Float_t phoEMC_mini[3];

	Float_t pxMC_mini[nMCtrack_MAX];
	Float_t pyMC_mini[nMCtrack_MAX];
	Float_t pzMC_mini[nMCtrack_MAX];
	Float_t etaMC_mini[nMCtrack_MAX];
	Float_t etaStarMC_mini[nMCtrack_MAX];
	Float_t ptStarMC_mini[nMCtrack_MAX];
	Float_t chargeMC_mini[nMCtrack_MAX];
	Int_t isDaughtersMC_mini[nMCtrack_MAX];
   	
	tree->SetBranchAddress("yMC_es_mini",&yMC_es_mini);
	tree->SetBranchAddress("Q2MC_es_mini",&Q2MC_es_mini);
	tree->SetBranchAddress("yMC_mini",&yMC_mini);
	tree->SetBranchAddress("Q2MC_mini",&Q2MC_mini);
	tree->SetBranchAddress("isQEDcMC_mini",&isQEDcMC_mini);
	tree->SetBranchAddress("isQEDbkg_mini",&isQEDbkg_mini);
	tree->SetBranchAddress("isDIFFbkg_mini",&isDIFFbkg_mini);
	tree->SetBranchAddress("elecPxMC_mini",&elecPxMC_mini);
	tree->SetBranchAddress("elecPyMC_mini",&elecPyMC_mini);
	tree->SetBranchAddress("elecPzMC_mini",&elecPzMC_mini);
	tree->SetBranchAddress("elecEMC_mini",&elecEMC_mini);
	tree->SetBranchAddress("phoPxMC_mini",&phoPxMC_mini);
	tree->SetBranchAddress("phoPyMC_mini",&phoPyMC_mini);
	tree->SetBranchAddress("phoPzMC_mini",&phoPzMC_mini);
	tree->SetBranchAddress("phoEMC_mini",&phoEMC_mini);

	tree->SetBranchAddress("dRRadPhot_mini",&dRRadPhot_mini);
	tree->SetBranchAddress("dPhiRadPhot_mini",&dPhiRadPhot_mini);

	tree->SetBranchAddress("nMCtrack_mini",&nMCtrack_mini);
	tree->SetBranchAddress("pxMC_mini",&pxMC_mini);
	tree->SetBranchAddress("pyMC_mini",&pyMC_mini);
	tree->SetBranchAddress("pzMC_mini",&pzMC_mini);
	tree->SetBranchAddress("etaMC_mini",&etaMC_mini);
	tree->SetBranchAddress("etaStarMC_mini",&etaStarMC_mini);
	tree->SetBranchAddress("ptStarMC_mini",&ptStarMC_mini);
	tree->SetBranchAddress("isDaughtersMC_mini",&isDaughtersMC_mini);

	Float_t xREC_es_mini,yREC_es_mini,Q2REC_es_mini;
	const int nRECtrack_MAX = 200;
	Int_t eventpass_mini;
	Int_t nRECtrack_mini;
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
   	Int_t typeChgREC_mini[nRECtrack_MAX];

	Float_t dcaPrimeREC_mini[nRECtrack_MAX];
	Float_t startHitsRadiusREC_mini[nRECtrack_MAX];
	Float_t dedxProtonREC_mini[nRECtrack_MAX];
	Float_t dedxLikelihoodProtonREC_mini[nRECtrack_MAX];
	Float_t dedxElectronREC_mini[nRECtrack_MAX];
	Float_t dedxLikelihoodElectronREC_mini[nRECtrack_MAX];
	
	tree->SetBranchAddress("eventpass_mini",&eventpass_mini);
	tree->SetBranchAddress("nRECtrack_mini",&nRECtrack_mini);
	tree->SetBranchAddress("w_pdg_mini",&w_mini);
	tree->SetBranchAddress("EpzREC_mini",&EpzREC_mini);
	tree->SetBranchAddress("totalMultREC_mini",&totalMultREC_mini);
	tree->SetBranchAddress("vertex_mini",&vertex_mini);
	tree->SetBranchAddress("eElectronBeam_mini",&eElectronBeam_mini);
	
	tree->SetBranchAddress("pxREC_mini",&pxREC_mini);
	tree->SetBranchAddress("pyREC_mini",&pyREC_mini);
	tree->SetBranchAddress("pzREC_mini",&pzREC_mini);

	tree->SetBranchAddress("etaREC_mini",&etaREC_mini);
	tree->SetBranchAddress("etaStarREC_mini",&etaStarREC_mini);
	tree->SetBranchAddress("passREC_mini",passREC_mini);
	// tree->SetBranchAddress("passTightREC_mini",passTightREC_mini);
	// tree->SetBranchAddress("passLooseREC_mini",passLooseREC_mini);

	tree->SetBranchAddress("xREC_es_mini",&xREC_es_mini);
	tree->SetBranchAddress("yREC_es_mini",&yREC_es_mini);
	tree->SetBranchAddress("Q2REC_es_mini",&Q2REC_es_mini);
	
	tree->SetBranchAddress("typeChgREC_mini",&typeChgREC_mini);
	tree->SetBranchAddress("dcaPrimeREC_mini",&dcaPrimeREC_mini);
	tree->SetBranchAddress("startHitsRadiusREC_mini",&startHitsRadiusREC_mini);
	tree->SetBranchAddress("dedxProtonREC_mini",&dedxProtonREC_mini);
	tree->SetBranchAddress("dedxLikelihoodProtonREC_mini",&dedxLikelihoodProtonREC_mini);
	tree->SetBranchAddress("dedxElectronREC_mini",&dedxElectronREC_mini);
	tree->SetBranchAddress("dedxLikelihoodElectronREC_mini",&dedxLikelihoodElectronREC_mini);

	//define output file
	TFile *outputfile = new TFile(outname,"RECREATE");
	/*
	eta distribution in different x. 
	*/
	TH1D* h_mult = new TH1D("h_mult",";N_{trk}",50,0,50);
	TH1D* h_eta = new TH1D("h_eta",";#eta",100,-3,3);
	TH1D* h_eta_pos[4];
	TH1D* h_eta_neg[4];
	for(int j=0;j<4;j++){
		h_eta_pos[j]= new TH1D(Form("h_eta_pos_%d",j),Form("h_eta_pos_%d",j),100,-3,3);
		h_eta_neg[j]= new TH1D(Form("h_eta_neg_%d",j),Form("h_eta_neg_%d",j),100,-3,3);
	}
	/*
	P(n) distribution in different Q2, y, and eta bins.
	*/
	TH1D* h_Pn[4][4][4];
	TH1D* h_Pn_QEDc[4][4][4];
	TH1D* h_Pn_GEN[4][4][4];
	TH2D* h_Pn_cor[4][4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				h_Pn[i][j][k] = new TH1D(Form("h_Pn_%d_%d_%d",i,j,k),Form("h_Pn_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_QEDc[i][j][k] = new TH1D(Form("h_Pn_QEDc_%d_%d_%d",i,j,k),Form("h_Pn_QEDc_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_GEN[i][j][k] = new TH1D(Form("h_Pn_GEN_%d_%d_%d",i,j,k),Form("h_Pn_GEN_%d_%d_%d",i,j,k),30,0,30);
				h_Pn_cor[i][j][k] = new TH2D(Form("h_Pn_cor_%d_%d_%d",i,j,k),Form("h_Pn_cor_%d_%d_%d",i,j,k),30,0,30,30,0,30);
			}
		}
	}
	/*
	P(n) distribution in different Q2, y bins within 0<eta*<4.
	*/

	//QED Compton event
	TH1D* h_eGammaPhiMC[4][4][2][2];
	TH1D* h_eGammaPhiMC_allQ2y[2][2];
	TH1D* h_sumPtMC[4][4][2][2];
	TH1D* h_sumPtMC_allQ2y[2][2];
	TH1D* h_EpzElecPhotMC[4][4][2][2];
	TH1D* h_EpzElecPhotMC_allQ2y[2][2];
	TH2D* h_deltaPhiVsThetaMC[2][2];
	
	TH2D* h_gammaMCEnergyVsTheta_allQ2y[2][2];
	TH2D* h_elecMCEnergyVsTheta_allQ2y[2][2];

	TH2D* h_dRRadPhotVsMult[4][4]; 
	TH2D* h_dPhiRadPhotVsMult[4][4];

	TH1D* h_Pn_HCM[4][4];
	TH1D* h_Pn_genREC_HCM[4][4];	
	TH1D* h_Pn_genNextREC_HCM[4][4];	
	TH1D* h_Pn_genNotREC_HCM[4][4];	
	TH1D* h_Pn_GEN_HCM[4][4];
	TH1D* h_Pn_GEN_HCM_noSel[4][4];
	TH2D* h_Pn_cor_HCM[4][4];

	TH1D* h_K0sMass[4][4];
	TH1D* h_PhotMass[4][4][4];
	TH2D* h_AP[4][4][4];

	for(int k=0;k<2;k++){
		for(int l=0;l<2;l++){
			h_eGammaPhiMC_allQ2y[k][l] = new TH1D(Form("h_eGammaPhiMC_allQ2y_%d_%d",k,l),";#Delta#phi",100,-3.15,3.15);
			h_sumPtMC_allQ2y[k][l] = new TH1D(Form("h_sumPtMC_allQ2y_%d_%d",k,l),";sum pt",100,0,5);
			h_EpzElecPhotMC_allQ2y[k][l] = new TH1D(Form("h_EpzElecPhotMC_allQ2y_%d_%d",k,l),";E-pz",100,0,70);
			h_deltaPhiVsThetaMC[k][l] = new TH2D(Form("h_deltaPhiVsThetaMC_%d_%d",k,l),";theta;#Delta#phi",100,-3.15,3.15,100,-3.15,3.15);
			h_gammaMCEnergyVsTheta_allQ2y[k][l] = new TH2D(Form("h_gammaMCEnergyVsTheta_allQ2y_%d_%d",k,l),";theta;energy",100,-3.15,3.15,500,0,50);
			h_elecMCEnergyVsTheta_allQ2y[k][l] = new TH2D(Form("h_elecMCEnergyVsTheta_allQ2y_%d_%d",k,l),";theta;energy",100,-3.15,3.15,500,0,50);
		}
	}
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<2;k++){
				for(int l=0;l<2;l++){
					h_eGammaPhiMC[i][j][k][l] = new TH1D(Form("h_eGammaPhiMC_%d_%d_%d_%d",i,j,k,l),";#Delta#phi",100,-3.15,3.15);
					h_sumPtMC[i][j][k][l] = new TH1D(Form("h_sumPtMC_%d_%d_%d_%d",i,j,k,l),";sum pt",100,0,5);
					h_EpzElecPhotMC[i][j][k][l] = new TH1D(Form("h_EpzElecPhotMC_%d_%d_%d_%d",i,j,k,l),";E-pz",100,0,70);
				}
			}
			h_dRRadPhotVsMult[i][j]= new TH2D(Form("h_dRRadPhotVsMult_%d_%d",i,j),";mult;dR",30,0,30,300,0,15);
			h_dPhiRadPhotVsMult[i][j]= new TH2D(Form("h_dPhiRadPhotVsMult_%d_%d",i,j),";mult;dPhi",30,0,30,300,-6.28,6.28);

			h_Pn_HCM[i][j] = new TH1D(Form("h_Pn_HCM_%d_%d",i,j),Form("h_Pn_HCM_%d_%d",i,j),30,0,30);
			h_Pn_genREC_HCM[i][j] = new TH1D(Form("h_Pn_genREC_HCM_%d_%d",i,j),Form("h_Pn_genREC_HCM_%d_%d",i,j),30,0,30);
			h_Pn_genNextREC_HCM[i][j] = new TH1D(Form("h_Pn_genNextREC_HCM_%d_%d",i,j),Form("h_Pn_genNextREC_HCM_%d_%d",i,j),30,0,30);
			h_Pn_genNotREC_HCM[i][j] = new TH1D(Form("h_Pn_genNotREC_HCM_%d_%d",i,j),Form("h_Pn_genNotREC_HCM_%d_%d",i,j),30,0,30);
			h_Pn_GEN_HCM[i][j] = new TH1D(Form("h_Pn_GEN_HCM_%d_%d",i,j),Form("h_Pn_GEN_HCM_%d_%d",i,j),30,0,30);
			h_Pn_GEN_HCM_noSel[i][j] = new TH1D(Form("h_Pn_GEN_HCM_noSel_%d_%d",i,j),Form("h_Pn_GEN_HCM_noSel_%d_%d",i,j),30,0,30);
			h_Pn_cor_HCM[i][j] = new TH2D(Form("h_Pn_cor_HCM_%d_%d",i,j),Form("h_Pn_cor_HCM_%d_%d",i,j),30,0,30,30,0,30);
			
			h_K0sMass[i][j] = new TH1D(Form("h_K0sMass_%d_%d",i,j),Form("h_K0sMass_%d_%d",i,j),200,0.25,0.54);

			for(int m=0;m<4;m++){
				h_PhotMass[i][j][m] = new TH1D(Form("h_PhotMass_%d_%d_%d",i,j,m),Form("h_PhotMass_%d_%d_%d",i,j,m),80,0,0.25);
				h_AP[i][j][m] = new TH2D(Form("h_AP_%d_%d_%d",i,j,m),";#alpha;p'_{T}",100,-1,1,100,0,0.3);
			}
			
		}
	}

	TH2D* h_Q2vsX = new TH2D("h_Q2vsX","h_Q2vsX",1000,0.00001,0.01,50,1,100);
	TH1D* h_vtxZ  = new TH1D("h_vtxZ","h_vtxZ", 100,-50,50);
	TH1D* h_y  = new TH1D("h_y","h_y", 200,0,1);
	TH1D* h_y_QEDc  = new TH1D("h_y_QEDc","h_y_QEDc", 200,0,1);
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
	TH1D* h_trackEpz_all = new TH1D("h_trackEpz","h_trackEpz",500,0,100);
	TH1D* h_noSel_pt = new TH1D("h_noSel_pt","pt",100,0,10);
	TH1D* h_noSel_eta = new TH1D("h_noSel_eta","eta",100,-1.7,1.7);

	TH1D* h_dedxElectronLikehood = new TH1D("h_dedxElectronLikehood","h_dedxElectronLikehood",100,-3,3);
	TH2D* h_dedxProtonVsp = new TH2D("h_dedxProtonVsp",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxProtonVspCut = new TH2D("h_dedxProtonVspCut",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxElectronVsp = new TH2D("h_dedxElectronVsp",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxElectronVspCut = new TH2D("h_dedxElectronVspCut",";p(GeV);dE/dx",100,0,5,300,0,100);
	
	TH1D* h_TestMass = new TH1D("h_TestMass","mass",100,0,0.2);
	TH1D* h_dedxElectronThetaCut[4];
	TH1D* h_dedxElectronPtCut[4];
	for(int m=0;m<4;m++){
	 h_dedxElectronThetaCut[m] = new TH1D(Form("h_dedxElectronThetaCut_%d",m),";#theta (rad)",100,0,3.14);
	 h_dedxElectronPtCut[m] = new TH1D(Form("h_dedxElectronPtCut_%d",m),";#p_{T} (GeV)",100,0,5);
	}
	
	TH1D* h_chargedDcaPrime = new TH1D("h_chargedDcaPrime",";charge*DCA'",100,-10,10);
	TH1D* h_chargedDcaPrimeProton = new TH1D("h_chargedDcaPrimeProton",";charge*DCA'",100,-10,10);
	TH1D* h_chargeRstart = new TH1D("h_chargeRstart",";charge*R_{start}",200,-50,50);
	TH1D* h_chargeRstartProton = new TH1D("h_chargeRstartProton",";charge*R_{start}",200,-50,50);
	TH1D* h_chargeRstartNoCut = new TH1D("h_chargeRstartNoCut",";charge*R_{start}",200,-100,100);
	TH1D* h_chargeRstartProtonNoCut = new TH1D("h_chargeRstartProtonNoCut",";charge*R_{start}",200,-50,50);
	TH1D* h_chargeRstartSplit = new TH1D("h_chargeRstartSplit",";charge*R_{start}",200,-50,50);
	TH1D* h_chargeRstartProtonSplit = new TH1D("h_chargeRstartProtonSplit",";charge*R_{start}",200,-50,50);
	TH2D* h_deltaPtDeltaEta = new TH2D("h_deltaPtDeltaEta",";p_{T};#eta",200,-1,1,1000,-3.2,3.2);

	cout << "==== Total number of events ~ " << tree->GetEntries() << " ========" << endl;
	for(int ievent = 0; ievent < tree->GetEntries(); ievent++){

		tree->GetEntry(ievent);
		if( ievent%100000 == 0 )cout << "Events ~ " << ievent << endl;
		
		double n_particle_eta[4] = {0.,0.,0.,0.};
		double n_particle_HCM = 0.;
		double n_particle_HCM_noSel = 0.;
		//MC generator
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
				if( isDaughtersMC_mini[itrk]!=0 ) continue;
				if( TMath::Hypot(pxMC_mini[itrk],pyMC_mini[itrk]) > 0.15 ){
					if( fabs(etaMC_mini[itrk]) < 1.6 ){
						if( etaStarMC_mini[itrk] > 0 && etaStarMC_mini[itrk] < 4.0 ){
								n_particle_HCM++;	
						}
						else{
							n_particle_HCM_noSel++;
						}
					}
				}
				if( TMath::Hypot(pxMC_mini[itrk],pyMC_mini[itrk]) < 0.15 ) continue;
				for(int ieta = 0; ieta < 3; ieta++){
					if( etaMC_mini[itrk] > eta_bins[2*ieta] && etaMC_mini[itrk] < eta_bins[2*ieta+1] ){
						n_particle_eta[ieta]++;
					}
				}
				if( etaMC_mini[itrk] > -1.6 && etaMC_mini[itrk] < 1.6 ){
					n_particle_eta[3]++;
				}
			}
			int generator_index = -1;
			if( ifile_==1 ){
				if(ievent<3e7) generator_index=0;
				else generator_index=1;
			}
			if( ifile_==2 ){
				//no mixing diffractive
				// if( ievent>=39999984 ) continue;
				if(ievent<44999982) generator_index=0;
				else generator_index=1;
			}
			else{
				generator_index=0;
			}

			//filling QED Compton delta phi
			if(Q2_INDEX>=0 && y_INDEX >= 0 && isQEDcMC_mini == 0 && isQEDbkg_mini==0 ) {
				TLorentzVector eMC;eMC.SetPxPyPzE(elecPxMC_mini,elecPyMC_mini,elecPzMC_mini,elecEMC_mini);
				TLorentzVector eGamma, sumEgamma;
				for(int itype=0;itype<3;itype++){
					eGamma.SetPxPyPzE(phoPxMC_mini[itype],phoPyMC_mini[itype],phoPzMC_mini[itype],phoEMC_mini[itype]);
					if( eGamma.E() <= 0 ) continue;
					sumEgamma = eMC+eGamma;
					if( n_particle_eta[3] < 2 ){
						if(eGamma.E()>0.1) h_deltaPhiVsThetaMC[0][generator_index]->Fill( eGamma.Theta(),eMC.DeltaPhi(eGamma),w_mini );
						h_eGammaPhiMC_allQ2y[0][generator_index]->Fill( eMC.DeltaPhi(eGamma), w_mini );
						if( fabs(eMC.DeltaPhi(eGamma)) > 2.7 ){
							h_elecMCEnergyVsTheta_allQ2y[0][0]->Fill(eMC.Theta(),eMC.E());
							h_gammaMCEnergyVsTheta_allQ2y[0][0]->Fill(eGamma.Theta(),eGamma.E());
						}
						if( fabs(eMC.DeltaPhi(eGamma)) < 2.7 && fabs(eMC.DeltaPhi(eGamma))>0.4 ){
							h_elecMCEnergyVsTheta_allQ2y[0][1]->Fill(eMC.Theta(),eMC.E());
							h_gammaMCEnergyVsTheta_allQ2y[0][1]->Fill(eGamma.Theta(),eGamma.E());
						}
						h_sumPtMC_allQ2y[0][generator_index]->Fill( sumEgamma.Pt(), w_mini );
						h_EpzElecPhotMC_allQ2y[0][generator_index]->Fill( sumEgamma.E()-sumEgamma.Pz(), w_mini );
						h_eGammaPhiMC[Q2_INDEX][y_INDEX][0][generator_index]->Fill( eMC.DeltaPhi(eGamma), w_mini );
						h_sumPtMC[Q2_INDEX][y_INDEX][0][generator_index]->Fill( sumEgamma.Pt(), w_mini );
						h_EpzElecPhotMC[Q2_INDEX][y_INDEX][0][generator_index]->Fill( sumEgamma.E()-sumEgamma.Pz(), w_mini);
					}
					else{
						if(eGamma.E()>0.1) h_deltaPhiVsThetaMC[1][generator_index]->Fill( eGamma.Theta(),eMC.DeltaPhi(eGamma),w_mini );
						h_eGammaPhiMC_allQ2y[1][generator_index]->Fill( eMC.DeltaPhi(eGamma), w_mini );
						if( fabs(eMC.DeltaPhi(eGamma)) > 2.7 ){
							h_elecMCEnergyVsTheta_allQ2y[1][0]->Fill(eMC.Theta(),eMC.E());
							h_gammaMCEnergyVsTheta_allQ2y[1][0]->Fill(eGamma.Theta(),eGamma.E());
						}
						if( fabs(eMC.DeltaPhi(eGamma)) < 2.7 && fabs(eMC.DeltaPhi(eGamma))>0.4 ){
							h_elecMCEnergyVsTheta_allQ2y[1][1]->Fill(eMC.Theta(),eMC.E());
							h_gammaMCEnergyVsTheta_allQ2y[1][1]->Fill(eGamma.Theta(),eGamma.E());
						}
						h_sumPtMC_allQ2y[1][generator_index]->Fill( sumEgamma.Pt(), w_mini );
						h_EpzElecPhotMC_allQ2y[1][generator_index]->Fill( sumEgamma.E()-sumEgamma.Pz(), w_mini );
						h_eGammaPhiMC[Q2_INDEX][y_INDEX][1][generator_index]->Fill( eMC.DeltaPhi(eGamma), w_mini );
						h_sumPtMC[Q2_INDEX][y_INDEX][1][generator_index]->Fill( sumEgamma.Pt(), w_mini );
						h_EpzElecPhotMC[Q2_INDEX][y_INDEX][1][generator_index]->Fill( sumEgamma.E()-sumEgamma.Pz(), w_mini);
					}
				}
				
			}
			if(Q2_INDEX >=0 && y_INDEX >= 0 && isQEDcMC_mini == 0 && isQEDbkg_mini==0 ){//no QEDc event counted as radiative Gen
				//HCM frame
				h_Pn_GEN_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM, w_mini );
				h_Pn_GEN_HCM_noSel[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_noSel, w_mini );
				//lab frame
				h_Pn_GEN[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta[0], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta[1], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta[2], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta[3], w_mini );
			
				h_dPhiRadPhotVsMult[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM, dPhiRadPhot_mini);
				h_dRRadPhotVsMult[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM, dRRadPhot_mini);
			}
		}

		if( Q2REC_es_mini > Q2max || Q2REC_es_mini < Q2min ) continue;
		if( yREC_es_mini > ymax || yREC_es_mini < ymin ) continue;
		if( eventpass_mini != 1 ) continue;
		
		h_Q2vsX->Fill( xREC_es_mini, Q2REC_es_mini, w_mini);
		h_y->Fill( yREC_es_mini, w_mini);
		if( ifile_!=0 && isQEDcMC_mini==1 ) h_y_QEDc->Fill(yREC_es_mini, w_mini);
		h_vtxZ->Fill( vertex_mini[2], w_mini );

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

		vector< TLorentzVector> elecp_vect,elecm_vect;
		for(int itrk=0;itrk<nRECtrack_mini;itrk++){
			if( passREC_mini[itrk] != 1 ) continue;
			if( fabs(etaREC_mini[itrk]) > 1.6 ) continue;
			h_dedxElectronLikehood->Fill( 
				dedxLikelihoodElectronREC_mini[itrk], w_mini );
			if(dedxLikelihoodElectronREC_mini[itrk] < electron_likelihood) continue;
			TLorentzVector electron_candidate;
			double E_cand = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+pyREC_mini[itrk]*pyREC_mini[itrk]+pzREC_mini[itrk]*pzREC_mini[itrk]+ELECTRON_MASS*ELECTRON_MASS);
			electron_candidate.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_cand);
			int track_charge = typeChgREC_mini[itrk];
			if( typeChgREC_mini[itrk] > 0 ) track_charge = 1;
			else if( typeChgREC_mini[itrk] < 0 ) track_charge = -1;
			else continue;
			if( track_charge == -1 ) elecm_vect.push_back(electron_candidate);
			if( track_charge == +1 ) elecp_vect.push_back(electron_candidate);
		}
		if(Q2_INDEX>-1 && y_INDEX>-1){
			//unlike-sign pairs
			for(int icand=0;icand<elecm_vect.size();icand++){
				for(int jcand=0;jcand<elecp_vect.size();jcand++){
					TLorentzVector photon_candidate = elecm_vect[icand]+elecp_vect[jcand];
					h_PhotMass[Q2_INDEX][y_INDEX][0]->Fill( photon_candidate.M(), w_mini );
				}
			}
			//like-sign pairs
			for(int icand=0;icand<elecm_vect.size();icand++){
				for(int jcand=icand+1;jcand<elecm_vect.size();jcand++){
					TLorentzVector photon_candidate = elecm_vect[icand]+elecm_vect[jcand];
					h_PhotMass[Q2_INDEX][y_INDEX][2]->Fill( photon_candidate.M(), w_mini );
				}
			}
			for(int icand=0;icand<elecp_vect.size();icand++){
				for(int jcand=icand+1;jcand<elecp_vect.size();jcand++){
					TLorentzVector photon_candidate = elecp_vect[icand]+elecp_vect[jcand];
					h_PhotMass[Q2_INDEX][y_INDEX][2]->Fill( photon_candidate.M(), w_mini );
				}
			}
		}



		for(int itrk = 0; itrk < nRECtrack_mini; itrk++){

			//filling eta_lab
			if( passREC_mini[itrk] != 1 ) continue; 
				h_eta->Fill(etaREC_mini[itrk], w_mini);
			if( fabs(etaREC_mini[itrk]) > 1.6 ) continue;

			int chargetrack_1 = typeChgREC_mini[itrk];
			if( typeChgREC_mini[itrk] > 0 ) chargetrack_1 = 1;
			if( typeChgREC_mini[itrk] < 0 ) chargetrack_1 = -1;
	
			//double nested loops
			// for(int jtrk = itrk+1; jtrk < nRECtrack_mini; jtrk++){
				
				// if( itrk==jtrk ) continue;
				// if( passREC_mini[jtrk] != 1 ) continue;
				// if( abs(etaREC_mini[jtrk]) > 1.6 ) continue;
				
				// int chargetrack_2 = typeChgREC_mini[jtrk];
				// if( typeChgREC_mini[jtrk] > 0 ) chargetrack_2 = 1;
				// if( typeChgREC_mini[jtrk] < 0 ) chargetrack_2 = -1;

				// double E_pip = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]
				// 	+pyREC_mini[itrk]*pyREC_mini[itrk]
				// 	+pzREC_mini[itrk]*pzREC_mini[itrk]
				// 	+PIMASS*PIMASS);
				// pip.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_pip);
				// double E_pim = sqrt( pxREC_mini[jtrk]*pxREC_mini[jtrk]
				// 	+pyREC_mini[jtrk]*pyREC_mini[jtrk]
				// 	+pzREC_mini[jtrk]*pzREC_mini[jtrk]
				// 	+PIMASS*PIMASS );
				// pim.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_pim);
				// if( chargetrack_1 != chargetrack_2 ){
				// 	k0s_candidate = pip+pim;
				// 	//AP plot:
				// 	TVector3 k0s_candidate_3Vect = k0s_candidate.Vect();
				// 	double p_angle = pip.Angle(k0s_candidate_3Vect);
				// 	double m_angle = pim.Angle(k0s_candidate_3Vect);
				// 	double pt_p = pip.P()*TMath::Sin(p_angle);
				// 	double pL_p = pip.P()*TMath::Cos(p_angle);
				// 	double pt_m = pim.P()*TMath::Sin(m_angle);
				// 	double pL_m = pim.P()*TMath::Cos(m_angle);
				// 	double alpha = (pL_p-pL_m)/(pL_p+pL_m);
				// 	if(Q2_INDEX>-1 && y_INDEX>-1){
				// 		h_AP[Q2_INDEX][y_INDEX][2]->Fill(alpha, pt_p, w_mini);
				// 		h_AP[Q2_INDEX][y_INDEX][2]->Fill(alpha, pt_m, w_mini);
				// 	}
				// 	//end AP
				// }
				// else{k0s_candidate.SetPxPyPzE(-99,-99,-99,-99);}
				
				// if(dedxLikelihoodElectronREC_mini[jtrk] < electron_likelihood) continue;

				// double E_elecp = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+
				// 	pyREC_mini[itrk]*pyREC_mini[itrk]+
				// 	pzREC_mini[itrk]*pzREC_mini[itrk]+
				// 	ELECTRON_MASS*ELECTRON_MASS);
				// elecp.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_elecp);
				// double E_elecm = sqrt(pxREC_mini[jtrk]*pxREC_mini[jtrk]+
				// 	pyREC_mini[jtrk]*pyREC_mini[jtrk]+
				// 	pzREC_mini[jtrk]*pzREC_mini[jtrk]+
				// 	ELECTRON_MASS*ELECTRON_MASS);
				// elecm.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_elecm);
				// photon_candidate = elecm+elecp;
				// if( photon_candidate.M() < min_mass ) {
				// 	photon_candidate_min = photon_candidate;
				// 	min_mass = photon_candidate.M();
				// 	min_track2_charge = chargetrack_2;
				// 	elecp_min = elecp;
				// 	elecm_min = elecm;
				// }
				// if(Q2_INDEX>-1 && y_INDEX>-1){
				// 	//unlike-sign pairs
				// 	if( chargetrack_1 != chargetrack_2  ){
				// 		h_PhotMass[Q2_INDEX][y_INDEX][0]->Fill( photon_candidate.M(), w_mini );
				// 		if( photon_candidate.M() < 0.1 && photon_candidate.M() > 0. ) {
				// 			h_dedxElectronThetaCut[0]->Fill(elecp.Theta(), w_mini);
				// 			h_dedxElectronThetaCut[0]->Fill(elecm.Theta(), w_mini);

				// 			h_dedxElectronPtCut[0]->Fill(elecp.Pt(), w_mini);
				// 			h_dedxElectronPtCut[0]->Fill(elecm.Pt(), w_mini);
				// 		}
				// 	}
				// 	//like-sign pairs
				// 	if( chargetrack_1 == chargetrack_2 ){
				// 		double deltaR = elecp.DeltaR(elecm);
				// 		if( deltaR < 0.05 ) continue;
				// 		h_PhotMass[Q2_INDEX][y_INDEX][2]->Fill( photon_candidate.M(), w_mini );
				// 		if( photon_candidate.M() < 0.1 && photon_candidate.M() > 0. ) {
				// 			h_dedxElectronThetaCut[2]->Fill(elecp.Theta(), w_mini);
				// 			h_dedxElectronThetaCut[2]->Fill(elecm.Theta(), w_mini);

				// 			h_dedxElectronPtCut[2]->Fill(elecp.Pt(), w_mini);
				// 			h_dedxElectronPtCut[2]->Fill(elecm.Pt(), w_mini);
				// 		}
				// 	}
		
				// }
							
				// if(Q2_INDEX>-1 && y_INDEX>-1){
				// 	if(k0s_candidate.E()!=-99) h_K0sMass[Q2_INDEX][y_INDEX]->Fill( k0s_candidate.M(), w_mini );
				// 	if( elecp.E()!=-99 && elecm.E() != -99 && (chargetrack_1!=chargetrack_2) ){
				// 		//AP plot:
				// 		TVector3 photon_candidate_3Vect = photon_candidate.Vect();
				// 		double p_angle = elecp.Angle(photon_candidate_3Vect);
				// 		double m_angle = elecm.Angle(photon_candidate_3Vect);
				// 		double pt_p = elecp.P()*TMath::Sin(p_angle);
				// 		double pL_p = elecp.P()*TMath::Cos(p_angle);
				// 		double pt_m = elecm.P()*TMath::Sin(m_angle);
				// 		double pL_m = elecm.P()*TMath::Cos(m_angle);
				// 		double alpha = (pL_p-pL_m)/(pL_p+pL_m);
				// 		if( photon_candidate.M() < 0.1 ){
				// 			h_AP[Q2_INDEX][y_INDEX][0]->Fill(alpha, pt_p, w_mini);
				// 			h_AP[Q2_INDEX][y_INDEX][0]->Fill(alpha, pt_m, w_mini);
				// 		}
				// 		h_AP[Q2_INDEX][y_INDEX][1]->Fill(alpha, pt_p, w_mini);
				// 		h_AP[Q2_INDEX][y_INDEX][1]->Fill(alpha, pt_m, w_mini);
				// 		if(fabs(alpha)>0.8) h_TestMass->Fill(photon_candidate.M(), w_mini);
				// 		//end AP
				// 	}
				// }
		
			// }//end double loop

			// if(Q2_INDEX>-1 && y_INDEX>-1){
			// //unlike-sign pairs
			// 	if( chargetrack_1 != min_track2_charge && photon_candidate_min.E() != -99  ){
			// 		h_PhotMass[Q2_INDEX][y_INDEX][0]->Fill( photon_candidate_min.M(), w_mini );
			// 		if( photon_candidate_min.M() < 0.1 && photon_candidate_min.M() > 0. ) {
			// 			h_dedxElectronThetaCut[0]->Fill(elecp_min.Theta(), w_mini);
			// 			h_dedxElectronThetaCut[0]->Fill(elecm_min.Theta(), w_mini);

			// 			h_dedxElectronPtCut[0]->Fill(elecp_min.Pt(), w_mini);
			// 			h_dedxElectronPtCut[0]->Fill(elecm_min.Pt(), w_mini);
			// 		}
			// 	}
			// //like-sign pairs
			// 	if( chargetrack_1 == min_track2_charge && photon_candidate_min.E() != -99 ){
			// 		double deltaR = elecp_min.DeltaR(elecm_min);
			// 		if( deltaR < 0.05 ) continue;
			// 		h_PhotMass[Q2_INDEX][y_INDEX][2]->Fill( photon_candidate_min.M(), w_mini );
			// 		if( photon_candidate_min.M() < 0.1 && photon_candidate_min.M() > 0. ) {
			// 			h_dedxElectronThetaCut[2]->Fill(elecp_min.Theta(), w_mini);
			// 			h_dedxElectronThetaCut[2]->Fill(elecm_min.Theta(), w_mini);

			// 			h_dedxElectronPtCut[2]->Fill(elecp_min.Pt(), w_mini);
			// 			h_dedxElectronPtCut[2]->Fill(elecm_min.Pt(), w_mini);
			// 		}
			// 	}
		
			// }
			// photon_candidate_min.SetPxPyPzE(-99,-99,-99,-99);
			// elecp_min.SetPxPyPzE(-99,-99,-99,-99);
			// elecm_min.SetPxPyPzE(-99,-99,-99,-99);

			//Rstart without cut
			int chargetrack = typeChgREC_mini[itrk];
			if( typeChgREC_mini[itrk] >= 1 ) chargetrack = 1;
			if( typeChgREC_mini[itrk] < 0 ) chargetrack = -1;
			
			h_chargeRstartNoCut->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
			if(dedxLikelihoodProtonREC_mini[itrk] > electron_likelihood) h_chargeRstartProtonNoCut->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini);
			
			if( passREC_mini[itrk] != 1 ) continue;
			if( etaREC_mini[itrk] > etamax ){
				etamax = etaREC_mini[itrk];
			}
			if( etaStarREC_mini[itrk] < 0 || etaStarREC_mini[itrk] > 4.0 ){
				if( fabs(etaREC_mini[itrk]) < 1.6 ){
					n_particle_HCM_rec++;
					h_noSel_pt->Fill(TMath::Hypot(pxREC_mini[itrk],pyREC_mini[itrk]),w_mini);
					h_noSel_eta->Fill(etaREC_mini[itrk],w_mini);
				}
			}
			for(int ieta = 0; ieta < 3; ieta++){
				if( etaREC_mini[itrk] > eta_bins[2*ieta] && etaREC_mini[itrk] < eta_bins[2*ieta+1] ){
					n_particle_eta_rec[ieta]++;
				}
			}
			if( etaREC_mini[itrk] > -1.6 && etaREC_mini[itrk] < 1.6 ) {
				n_particle_eta_rec[3]++;
				trk_E += sqrt(pxREC_mini[itrk]*pxREC_mini[itrk] + pyREC_mini[itrk]*pyREC_mini[itrk] + pzREC_mini[itrk]*pzREC_mini[itrk] + PIMASS*PIMASS);
				trk_pz += pzREC_mini[itrk];
				// dE/dx
				double pREC = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+pyREC_mini[itrk]*pyREC_mini[itrk]+pzREC_mini[itrk]*pzREC_mini[itrk]);
				h_dedxProtonVsp->Fill( pREC, dedxProtonREC_mini[itrk], w_mini);
				h_dedxElectronVsp->Fill( pREC, dedxElectronREC_mini[itrk], w_mini);
				h_chargedDcaPrime->Fill( chargetrack*dcaPrimeREC_mini[itrk], w_mini );
				h_chargeRstart->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
				//PROTON
				if(dedxLikelihoodProtonREC_mini[itrk] > 0.003){
					h_dedxProtonVspCut->Fill( pREC, dedxProtonREC_mini[itrk], w_mini);
					h_chargedDcaPrimeProton->Fill( chargetrack*dcaPrimeREC_mini[itrk], w_mini );
					h_chargeRstartProton->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
				}
				//ELECTRON
				if(dedxLikelihoodElectronREC_mini[itrk] > electron_likelihood){
					h_dedxElectronVspCut->Fill( pREC, dedxElectronREC_mini[itrk], w_mini);
				}
			}
			for(int iy=0;iy<4;iy++){
				if(yREC_es_mini>ybins[iy] && yREC_es_mini<ybins[iy+1]){
					if(chargetrack > 0) h_eta_pos[iy]->Fill( etaREC_mini[itrk], w_mini );
					if(chargetrack < 0) h_eta_neg[iy]->Fill( etaREC_mini[itrk], w_mini );
				}
			}
		}
		
		// event level distributions.
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
		h_trackEpz_all->Fill(trk_E-trk_pz, w_mini);

		if(Q2_INDEX >=0 && y_INDEX >= 0){
			if( ifile_ != 0 ){
				//2D correlation between gen and rec mult
				h_Pn_cor_HCM[Q2_INDEX][y_INDEX]->Fill(n_particle_HCM_rec,n_particle_HCM,w_mini);
				h_Pn_cor[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta_rec[0], n_particle_eta[0], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta_rec[1], n_particle_eta[1], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta_rec[2], n_particle_eta[2], w_mini );
				h_Pn_cor[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta_rec[3], n_particle_eta[3], w_mini );
				
				//filling with signal rec==gen.
				if( n_particle_HCM_rec == n_particle_HCM ){
					h_Pn_genREC_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_rec, w_mini );
				}
				if( n_particle_HCM_rec != n_particle_HCM ){
					h_Pn_genNotREC_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_rec, w_mini );
				}
				if( n_particle_HCM_rec == n_particle_HCM+1 || n_particle_HCM_rec == n_particle_HCM-1 ){
					h_Pn_genNextREC_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_rec, w_mini );
				}
			}

			h_Pn_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM_rec, w_mini );
			h_Pn[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta_rec[0], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta_rec[1], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta_rec[2], w_mini );
			h_Pn[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta_rec[3], w_mini );
			
			h_mult->Fill( n_particle_eta_rec[3], w_mini );
		}
	}

	outputfile->Write();
	outputfile->Close();

}
