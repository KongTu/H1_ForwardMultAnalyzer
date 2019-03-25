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

int nEta=100;
int etaMin=-5;
int etaMax=5;
double etabins_low[] = {0,1.0,1.6,2.1,2.6,3.1,3.7,5.0};
double etabins_high[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.7,5.0};

//event level variables
TH2D* h_Q2vsX[2];
TH1D* h_yREC[2];
TH1D* h_vtxZ[2];
TH2D* h_vtxXY[2];
TH1D* h_EminusPz[2];
TH1D* h_weights[2];
TH1D* h_weights_gen[2];

//scattered electrons
TH1D* h_eEnergy[2];
TH1D* h_eTheta[2];
TH1D* h_hadESpacal[2];
TH1D* h_emE[2];
TH1D* h_clusterRadius[2];
TH2D* h_xyCluster[2];

//background finder
TH1D* h_bkg_vtxZ[10];
TH1D* h_bkg_Epz[10];
TH1D* h_nobkg_vtxZ[10];
TH1D* h_nobkg_Epz[10];

//track level variables
TH1D *h_eta_gen[2];
TH1D *h_eta[2];
TH1D *h_etaStar_gen_low[2];
TH1D *h_etaStar_gen_high[2];
TH1D *h_etaStar_low[2];
TH1D *h_etaStar_high[2];
TH1D *h_pt[2];
TH1D *h_ptStar[2];
TH1D *h_phi[2];
TH2D *h_etaPt[2];

//track level cut variables
TH1D *h_chi2vtx[2];
TH1D *h_chi2trk[2];
TH1D *h_Nhits[2];
TH1D *h_dcaPrime[2];
TH1D *h_dz0Prime[2];
TH1D *h_trackLength[2];
TH1D *h_rStartHits[2];
TH1D *h_trkTheta[2];
TH1D *h_zLengthHit[2];
TH1D *h_chi2Link[2];
TH1D *h_rZero[2];
TH1D *h_dcaPrimeSinTheta[2];
TH1D *h_momRes[2];

//P(N) analysis level histograms
TH1D *h_Q2REC[4][5];
TH1D *h_xREC[4][5];
TH1D *h_Pn[4][5];
TH1D *h_Pn_gen[4][5];
TH2D *h_response_Pn[4][5];
TH1D *h_Pn_all = new TH1D("h_Pn_all",";N",50,0,50);


// for(int jbin=0;jbin<4;jbin++){
//   for(int ibin=0;ibin<4;ibin++){
//      response[jbin][ibin].Setup(50,0,50);
//   }
// }


// for(int iq2 = 0; iq2 < 4; iq2++){
//   for(int ix = 0; ix < 4; ix++){
     
//      h_Q2REC[iq2][ix] = new TH1D(Form("h_Q2REC_%d_%d",iq2,ix),";Q^{2} (GeV^{2})",100,0.1,100);
//      h_xREC[iq2][ix] = new TH1D(Form("h_xREC_%d_%d",iq2,ix),";x",1000,0.00001,0.001);

//      h_Pn[iq2][ix] = new TH1D(Form("h_Pn_%d_%d",iq2,ix),";N;P(n)",50,0,50);
//      h_Pn_gen[iq2][ix] = new TH1D(Form("h_Pn_gen_%d_%d",iq2,ix),";N;P(n)",50,0,50);
     
//      h_response_Pn[iq2][ix] = new TH2D(Form("h_response_Pn_%d_%d",iq2,ix),";P(n) measured;P(n) truth",50,0,50,50,0,50);
//   }
// }