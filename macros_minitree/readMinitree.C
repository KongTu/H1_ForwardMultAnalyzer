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
		file = new TFile("../new_output/data_highE_0607_noReweight_Tree_hadCaliNewKine.root");
	}else if(ifile_ == 1) {
		if(!isReweigh) file = new TFile("../new_output/mc_highE_DJANGOH_noReweight_Tree_hadCaliNew.root");
		else file = new TFile("../new_output/mc_highE_DJANGOH_fullReweight_Tree.root");
	}else if(ifile_ == 2){
		if(!isReweigh) file = new TFile("../new_output/mc_highE_RAPGAP_noReweight_Tree_hadCaliNew.root");
		else file = new TFile("../new_output/mc_highE_RAPGAP_fullReweight_Tree_hadCaliNewKine.root");
	}
	//output files
	TString outname;
	if( ifile_ == 0 ){
		outname = "../minitree_output/Pn_hist_data_hadCaliNewKine.root";
	}else if( ifile_ == 1 ){
		if(!isReweigh) outname = "../minitree_output/Pn_hist_django_extendEtalabLooseTrack.root";
		else outname = "../minitree_output/Pn_hist_django_hadCaliNewKine_reweigh.root";
	}else if( ifile_ == 2 ){
		if(!isReweigh) outname = "../minitree_output/Pn_hist_rapgap_extendEtalabLooseTrack.root";
		else outname = "../minitree_output/Pn_hist_rapgap_hadCaliNewKine_reweigh.root";
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
	Float_t eGammaPhiMC_mini;
	Float_t sumPtMC_mini;
	Int_t isQEDcMC_mini;
	Float_t elecPxMC_mini;
	Float_t elecPyMC_mini;
	Float_t elecPzMC_mini;
	Float_t elecEMC_mini;
	Float_t phoPxMC_mini;
	Float_t phoPyMC_mini;
	Float_t phoPzMC_mini;
	Float_t phoEMC_mini;

	Float_t pxMC_mini[nMCtrack_MAX];
	Float_t pyMC_mini[nMCtrack_MAX];
	Float_t pzMC_mini[nMCtrack_MAX];
	Float_t etaMC_mini[nMCtrack_MAX];
	Float_t etaStarMC_mini[nMCtrack_MAX];
	Float_t ptStarMC_mini[nMCtrack_MAX];
	Float_t chargeMC_mini[nMCtrack_MAX];
   	
	tree->SetBranchAddress("yMC_es_mini",&yMC_es_mini);
	tree->SetBranchAddress("Q2MC_es_mini",&Q2MC_es_mini);
	tree->SetBranchAddress("eGammaPhiMC_mini",&eGammaPhiMC_mini);
	tree->SetBranchAddress("sumPtMC_mini",&sumPtMC_mini);
	tree->SetBranchAddress("isQEDcMC_mini",&isQEDcMC_mini);
	tree->SetBranchAddress("elecPxMC_mini",&elecPxMC_mini);
	tree->SetBranchAddress("elecPyMC_mini",&elecPyMC_mini);
	tree->SetBranchAddress("elecPzMC_mini",&elecPzMC_mini);
	tree->SetBranchAddress("elecEMC_mini",&elecEMC_mini);
	tree->SetBranchAddress("phoPxMC_mini",&phoPxMC_mini);
	tree->SetBranchAddress("phoPyMC_mini",&phoPyMC_mini);
	tree->SetBranchAddress("phoPzMC_mini",&phoPzMC_mini);
	tree->SetBranchAddress("phoEMC_mini",&phoEMC_mini);

	tree->SetBranchAddress("nMCtrack_mini",&nMCtrack_mini);
	tree->SetBranchAddress("pxMC_mini",&pxMC_mini);
	tree->SetBranchAddress("pyMC_mini",&pyMC_mini);
	tree->SetBranchAddress("pzMC_mini",&pzMC_mini);
	tree->SetBranchAddress("etaMC_mini",&etaMC_mini);
	tree->SetBranchAddress("etaStarMC_mini",&etaStarMC_mini);
	tree->SetBranchAddress("ptStarMC_mini",&ptStarMC_mini);

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
	tree->SetBranchAddress("elecChargeREC_mini",&elecChargeREC_mini);
	
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
				h_Pn[i][j][k] = new TH1D(Form("h_Pn_%d_%d_%d",i,j,k),Form("h_Pn_%d_%d_%d",i,j,k),11,Pn_binning);
				h_Pn_QEDc[i][j][k] = new TH1D(Form("h_Pn_QEDc_%d_%d_%d",i,j,k),Form("h_Pn_QEDc_%d_%d_%d",i,j,k),11,Pn_binning);
				h_Pn_GEN[i][j][k] = new TH1D(Form("h_Pn_GEN_%d_%d_%d",i,j,k),Form("h_Pn_GEN_%d_%d_%d",i,j,k),11,Pn_binning);
				h_Pn_cor[i][j][k] = new TH2D(Form("h_Pn_cor_%d_%d_%d",i,j,k),Form("h_Pn_cor_%d_%d_%d",i,j,k),11,Pn_binning,11,Pn_binning);
			}
		}
	}
	/*
	P(n) distribution in different Q2, y bins within 0<eta*<4.
	*/

	TH1D* h_Pn_HCM[4][4];
	TH1D* h_Pn_genREC_HCM[4][4];	
	TH1D* h_Pn_genNextREC_HCM[4][4];	
	TH1D* h_Pn_genNotREC_HCM[4][4];	
	TH1D* h_Pn_GEN_HCM[4][4];
	TH2D* h_Pn_cor_HCM[4][4];

	TH1D* h_K0sMass[4][4];
	TH1D* h_K0sMassTight[4][4];
	TH1D* h_K0sMassLoose[4][4];
	TH1D* h_PhotMass[4][4];
	TH1D* h_PhotMassTight[4][4];
	TH1D* h_PhotMassLoose[4][4];
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			
			h_Pn_HCM[i][j] = new TH1D(Form("h_Pn_HCM_%d_%d",i,j),Form("h_Pn_HCM_%d_%d",i,j),11,Pn_binning);
			h_Pn_genREC_HCM[i][j] = new TH1D(Form("h_Pn_genREC_HCM_%d_%d",i,j),Form("h_Pn_genREC_HCM_%d_%d",i,j),11,Pn_binning);
			h_Pn_genNextREC_HCM[i][j] = new TH1D(Form("h_Pn_genNextREC_HCM_%d_%d",i,j),Form("h_Pn_genNextREC_HCM_%d_%d",i,j),11,Pn_binning);
			h_Pn_genNotREC_HCM[i][j] = new TH1D(Form("h_Pn_genNotREC_HCM_%d_%d",i,j),Form("h_Pn_genNotREC_HCM_%d_%d",i,j),11,Pn_binning);
			h_Pn_GEN_HCM[i][j] = new TH1D(Form("h_Pn_GEN_HCM_%d_%d",i,j),Form("h_Pn_GEN_HCM_%d_%d",i,j),11,Pn_binning);
			h_Pn_cor_HCM[i][j] = new TH2D(Form("h_Pn_cor_HCM_%d_%d",i,j),Form("h_Pn_cor_HCM_%d_%d",i,j),11,Pn_binning,11,Pn_binning);
			
			h_K0sMass[i][j] = new TH1D(Form("h_K0sMass_%d_%d",i,j),Form("h_K0sMass_%d_%d",i,j),200,0.25,0.54);
			h_K0sMassTight[i][j] = new TH1D(Form("h_K0sMassTight_%d_%d",i,j),Form("h_K0sMassTight_%d_%d",i,j),200,0.25,0.54);
			h_K0sMassLoose[i][j] = new TH1D(Form("h_K0sMassLoose_%d_%d",i,j),Form("h_K0sMassLoose_%d_%d",i,j),200,0.25,0.54);

			h_PhotMass[i][j] = new TH1D(Form("h_PhotMass_%d_%d",i,j),Form("h_PhotMass_%d_%d",i,j),200,0,0.25);
			h_PhotMassTight[i][j] = new TH1D(Form("h_PhotMassTight_%d_%d",i,j),Form("h_PhotMassTight_%d_%d",i,j),200,0,0.25);
			h_PhotMassLoose[i][j] = new TH1D(Form("h_PhotMassLoose_%d_%d",i,j),Form("h_PhotMassLoose_%d_%d",i,j),200,0,0.25);
		}
	}

	TH2D* h_Q2vsX = new TH2D("h_Q2vsX","h_Q2vsX",1000,0.00001,0.01,50,1,100);
	TH1D* h_vtxZ  = new TH1D("h_vtxZ","h_vtxZ", 100,-50,50);
	TH1D* h_y  = new TH1D("h_y","h_y", 200,0,1);

	//QED Compton event
	TH1D* h_eGammaPhiMC = new TH1D("h_eGammaPhiMC",";#Delta#phi",100,-3.15,3.15);
	TH1D* h_sumPtMC = new TH1D("h_sumPtMC",";#Delta#phi",100,0,5);

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

	TH2D* h_dedxProtonVsp = new TH2D("h_dedxProtonVsp",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxProtonVspCut = new TH2D("h_dedxProtonVspCut",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxElectronVsp = new TH2D("h_dedxElectronVsp",";p(GeV);dE/dx",100,0,5,300,0,100);
	TH2D* h_dedxElectronVspCut = new TH2D("h_dedxElectronVspCut",";p(GeV);dE/dx",100,0,5,300,0,100);
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
		
		double n_particle_eta[4] = {0.,0.,0.,0.};
		double n_particle_HCM = 0.;
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
				
				if( TMath::Hypot(pxMC_mini[itrk],pyMC_mini[itrk]) > 0.15 ){
					if( fabs(etaMC_mini[itrk]) < 1.6 ){
						if( etaStarMC_mini[itrk] > 0 && etaStarMC_mini[itrk] < 4.0 ){
							n_particle_HCM++;
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
			//filling QED Compton delta phi
			// if(Q2_INDEX >=0 && y_INDEX >= 0 && isQEDcMC_mini == 1 ) {
			// 	h_eGammaPhiMC->Fill( eGammaPhiMC_mini );
			// 	h_sumPtMC->Fill( sumPtMC_mini );
			// }
			// if(Q2_INDEX >=0 && y_INDEX >= 0 && isQEDcMC_mini == 0 ){//no QEDc event counted as radiative Gen
			//no QEDc event counted as radiative Gen
			if(Q2_INDEX >=0 && y_INDEX >= 0 ){
				//HCM frame
				h_Pn_GEN_HCM[Q2_INDEX][y_INDEX]->Fill( n_particle_HCM, w_mini );
				//lab frame
				h_Pn_GEN[Q2_INDEX][y_INDEX][0]->Fill( n_particle_eta[0], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][1]->Fill( n_particle_eta[1], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][2]->Fill( n_particle_eta[2], w_mini );
				h_Pn_GEN[Q2_INDEX][y_INDEX][3]->Fill( n_particle_eta[3], w_mini );
			}
		}

		if( Q2REC_es_mini > Q2max || Q2REC_es_mini < Q2min ) continue;
		if( yREC_es_mini > ymax || yREC_es_mini < ymin ) continue;
		if( eventpass_mini != 1 ) continue; 
	
		h_Q2vsX->Fill( xREC_es_mini, Q2REC_es_mini, w_mini);
		h_y->Fill( yREC_es_mini, w_mini);
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

		TLorentzVector k0s_candidate(0,0,0,0);
		TLorentzVector photon_candidate(-99,-99,-99,-99);
		TLorentzVector pip(0,0,0,0),pim(0,0,0,0);
		TLorentzVector elecp(0,0,0,0),elecm(0,0,0,0);


	//double nested loops
		for(int itrk = 0; itrk < nRECtrack_mini; itrk++){

			//filling eta_lab
			if( passREC_mini[itrk] == 1 ) h_eta->Fill(etaREC_mini[itrk], w_mini);
			
			for(int jtrk = itrk+1; jtrk < nRECtrack_mini; jtrk++){
				if( itrk==jtrk ) continue;
				if( passREC_mini[itrk] == 1 && passREC_mini[jtrk] == 1  ){

					double E_pip = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]
						+pyREC_mini[itrk]*pyREC_mini[itrk]
						+pzREC_mini[itrk]*pzREC_mini[itrk]
						+PIMASS*PIMASS);
					pip.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_pip);
					double E_pim = sqrt( pxREC_mini[jtrk]*pxREC_mini[jtrk]
						+pyREC_mini[jtrk]*pyREC_mini[jtrk]
						+pzREC_mini[jtrk]*pzREC_mini[jtrk]
						+PIMASS*PIMASS );
					pim.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_pim);

					k0s_candidate = pip+pim;

					if(dedxLikelihoodElectronREC_mini[itrk] > 0.003 && dedxLikelihoodElectronREC_mini[jtrk] > 0.003){
						double E_elecp = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+
							pyREC_mini[itrk]*pyREC_mini[itrk]+
							pzREC_mini[itrk]*pzREC_mini[itrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecp.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_elecp);
						double E_elecm = sqrt(pxREC_mini[jtrk]*pxREC_mini[jtrk]+
							pyREC_mini[jtrk]*pyREC_mini[jtrk]+
							pzREC_mini[jtrk]*pzREC_mini[jtrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecm.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_elecm);
						
						photon_candidate = elecp+elecm;
					}
					
					if(Q2_INDEX>-1 && y_INDEX>-1){
						h_K0sMass[Q2_INDEX][y_INDEX]->Fill( k0s_candidate.M(), w_mini );
						h_PhotMass[Q2_INDEX][y_INDEX]->Fill( photon_candidate.M(), w_mini );
					}
				}
				if( passTightREC_mini[itrk] == 1 && passTightREC_mini[jtrk] == 1  ){

					double E_pip = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]
						+pyREC_mini[itrk]*pyREC_mini[itrk]
						+pzREC_mini[itrk]*pzREC_mini[itrk]
						+PIMASS*PIMASS);
					pip.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_pip);
					double E_pim = sqrt( pxREC_mini[jtrk]*pxREC_mini[jtrk]
						+pyREC_mini[jtrk]*pyREC_mini[jtrk]
						+pzREC_mini[jtrk]*pzREC_mini[jtrk]
						+PIMASS*PIMASS );
					pim.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_pim);

					k0s_candidate = pip+pim;

					if(dedxLikelihoodElectronREC_mini[itrk] > 0.003 && dedxLikelihoodElectronREC_mini[jtrk] > 0.003){
						double E_elecp = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+
							pyREC_mini[itrk]*pyREC_mini[itrk]+
							pzREC_mini[itrk]*pzREC_mini[itrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecp.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_elecp);
						double E_elecm = sqrt(pxREC_mini[jtrk]*pxREC_mini[jtrk]+
							pyREC_mini[jtrk]*pyREC_mini[jtrk]+
							pzREC_mini[jtrk]*pzREC_mini[jtrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecm.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_elecm);
						
						photon_candidate = elecp+elecm;
					}
					
					if(Q2_INDEX>-1 && y_INDEX>-1){
						h_K0sMassTight[Q2_INDEX][y_INDEX]->Fill( k0s_candidate.M(), w_mini );
						h_PhotMassTight[Q2_INDEX][y_INDEX]->Fill( photon_candidate.M(), w_mini );
					}
				}
				if( passLooseREC_mini[itrk] == 1 && passLooseREC_mini[jtrk] == 1  ){

					double E_pip = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]
						+pyREC_mini[itrk]*pyREC_mini[itrk]
						+pzREC_mini[itrk]*pzREC_mini[itrk]
						+PIMASS*PIMASS);
					pip.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_pip);
					double E_pim = sqrt( pxREC_mini[jtrk]*pxREC_mini[jtrk]
						+pyREC_mini[jtrk]*pyREC_mini[jtrk]
						+pzREC_mini[jtrk]*pzREC_mini[jtrk]
						+PIMASS*PIMASS );
					pim.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_pim);

					k0s_candidate = pip+pim;

					if(dedxLikelihoodElectronREC_mini[itrk] > 0.003 && dedxLikelihoodElectronREC_mini[jtrk] > 0.003){
						double E_elecp = sqrt(pxREC_mini[itrk]*pxREC_mini[itrk]+
							pyREC_mini[itrk]*pyREC_mini[itrk]+
							pzREC_mini[itrk]*pzREC_mini[itrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecp.SetPxPyPzE(pxREC_mini[itrk],pyREC_mini[itrk],pzREC_mini[itrk],E_elecp);
						double E_elecm = sqrt(pxREC_mini[jtrk]*pxREC_mini[jtrk]+
							pyREC_mini[jtrk]*pyREC_mini[jtrk]+
							pzREC_mini[jtrk]*pzREC_mini[jtrk]+
							ELECTRON_MASS*ELECTRON_MASS);
						elecm.SetPxPyPzE(pxREC_mini[jtrk],pyREC_mini[jtrk],pzREC_mini[jtrk],E_elecm);
						
						photon_candidate = elecp+elecm;
					}
					if(Q2_INDEX>-1 && y_INDEX>-1){
						h_K0sMassLoose[Q2_INDEX][y_INDEX]->Fill( k0s_candidate.M(), w_mini );
						h_PhotMassLoose[Q2_INDEX][y_INDEX]->Fill( photon_candidate.M(), w_mini );
					}
				}
				
				//fill delta pt and delta eta to see split tracks
				if( passREC_mini[itrk] == 1 && passREC_mini[jtrk] == 1  ){
					int chargetrack = typeChgREC_mini[itrk];
					if( typeChgREC_mini[itrk] > 0 ) chargetrack = 1;
					if( typeChgREC_mini[itrk] < 0 ) chargetrack = -1;
					double deltaPt = TMath::Hypot(pxREC_mini[itrk],pyREC_mini[itrk]) - TMath::Hypot(pxREC_mini[jtrk],pyREC_mini[jtrk]);
					double deltaEta = etaREC_mini[itrk] - etaREC_mini[jtrk];
					h_deltaPtDeltaEta->Fill( deltaPt, deltaEta, w_mini );
					if( fabs(deltaPt) < 0.03 && fabs(deltaEta) < 0.05 ){
						h_chargeRstartSplit->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
						if( dedxLikelihoodProtonREC_mini[itrk] > 0.003 ) h_chargeRstartProtonSplit->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
					}
				}
			}
	//end double loop

			//Rstart without cut

			int chargetrack = typeChgREC_mini[itrk];
			if( typeChgREC_mini[itrk] >= 1 ) chargetrack = 1;
			if( typeChgREC_mini[itrk] < 0 ) chargetrack = -1;
			
			h_chargeRstartNoCut->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini );
			if(dedxLikelihoodProtonREC_mini[itrk] > 0.003) h_chargeRstartProtonNoCut->Fill( chargetrack*startHitsRadiusREC_mini[itrk], w_mini);
			
			if( passREC_mini[itrk] != 1 ) continue;
			if( etaREC_mini[itrk] > etamax ){
				etamax = etaREC_mini[itrk];
			}
			if( etaStarREC_mini[itrk] > 0 && etaStarREC_mini[itrk] < 4.0 ){
				if( fabs(etaREC_mini[itrk]) < 1.6 ){
					n_particle_HCM_rec++;
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
				if(dedxLikelihoodElectronREC_mini[itrk] > 0.003){
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
