#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>

#include "RiceStyle.h"

void makeCentCorrTable() {
   
   TFile* file = new TFile("../run/ForwardMultAnalyzer_newVars.root");
   TTree *tree = (TTree*) file->Get("properties");
   //gFile->GetObject("properties",tree);
   int nEta=40;
   double etaMin=-2.5;
   double etaMax=7.;
   int nPt=100;
   double ptMin=0;
   double ptMax=10.0;

   TH2D *h_genEtaStar=new TH2D("h_genEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_genMatchEtaStar=new TH2D("h_genMatchEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_genMatchCentEtaStar=new TH2D("h_genMatchCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_genMatchCombCentEtaStar=new TH2D("h_genMatchCombCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_genMatchFwdCombCentEtaStar=new TH2D("h_genMatchFwdCombCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   
   TH2D *h_recEtaStar=new TH2D("h_recEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_recMatchEtaStar=new TH2D("h_recMatchEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_recMatchCentEtaStar=new TH2D("h_recMatchCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_recMatchCombCentEtaStar=new TH2D("h_recMatchCombCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   TH2D *h_recMatchFwdCombCentEtaStar=new TH2D("h_recMatchFwdCombCentEtaStar",";#eta*;p_{T} (GeV/c)",nEta,etaMin,etaMax,nPt,ptMin,ptMax);
   

   float Q2min=5.;
   float Q2max=100.;
   float ymin=0.1;
   float ymax=0.6;

   float etaStarMin=0;
   float etaStarMax=4;

   double zvtxOffset=0.;
   
   double ptcut=0.2;
   double chi2nTrackCut=10.0;
   double chi2nVtxCut=7.0;
   double dcaCut=5.0;
   int NhitsCut=10;
   double trackLengthCut=2.0;

   int print=100;

   if(tree) {
      Int_t vertexType;
      Float_t w;
      Float_t elecPxREC,elecPyREC;
      Float_t vertex[3];
      Float_t xREC,yREC,Q2REC;
      Int_t nMCtrack;
      Float_t ptStarMC[400];
      Float_t etaStarMC[400];
      Float_t pxMC[400];
      Float_t pyMC[400];
      Int_t imatchMC[400];

      Int_t nRECtrack;
      Int_t typeChgREC[200];
      Float_t pxREC[200];
      Float_t pyREC[200];
      Float_t ptStarREC[200];
      Float_t phiStarREC[200];
      Float_t etaStarREC[200];
      Float_t chi2vtxREC[400];
      Float_t chi2nvREC[400];
      Int_t vtxNdfREC[200];
      Int_t nvNdfREC[200];
      Int_t vtxNHitsREC[200];
      Int_t nvNHitsREC[200];
      Float_t vtxTrackLengthREC[400];
      Float_t nvTrackLengthREC[400];
      Float_t dcaPrimeREC[400];
      Float_t dz0PrimeREC[400];

      Int_t imatchREC[400];
      Float_t ptResFstSelected[400];

      tree->SetBranchAddress("w",&w);
      tree->SetBranchAddress("vertexType",&vertexType);
      tree->SetBranchAddress("vertex",vertex);
      tree->SetBranchAddress("xREC",&xREC);
      tree->SetBranchAddress("yREC",&yREC);
      tree->SetBranchAddress("Q2REC",&Q2REC);
      tree->SetBranchAddress("elecPxREC",&elecPxREC);
      tree->SetBranchAddress("elecPyREC",&elecPyREC);

      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStarMC",etaStarMC);
      tree->SetBranchAddress("ptStarMC",ptStarMC);
      tree->SetBranchAddress("imatchMC",imatchMC);

      tree->SetBranchAddress("pxMC",&pxMC);
      tree->SetBranchAddress("pyMC",&pyMC);

      tree->SetBranchAddress("nRECtrack",&nRECtrack);
      tree->SetBranchAddress("typeChgREC",typeChgREC);
      tree->SetBranchAddress("pxREC",pxREC);
      tree->SetBranchAddress("pyREC",pyREC);
      tree->SetBranchAddress("etaStarREC",etaStarREC);
      tree->SetBranchAddress("phiStarREC",phiStarREC);
      tree->SetBranchAddress("ptStarREC",ptStarREC);
      tree->SetBranchAddress("imatchREC",imatchREC);
      tree->SetBranchAddress("chi2vtxREC",chi2vtxREC);
      tree->SetBranchAddress("chi2nvREC",chi2nvREC);
      tree->SetBranchAddress("vtxNdfREC",vtxNdfREC);
      tree->SetBranchAddress("nvNdfREC",nvNdfREC);
      tree->SetBranchAddress("vtxNHitsREC",vtxNHitsREC);
      tree->SetBranchAddress("nvNHitsREC",nvNHitsREC);
      tree->SetBranchAddress("vtxTrackLengthREC",vtxTrackLengthREC);
      tree->SetBranchAddress("nvTrackLengthREC",nvTrackLengthREC);
      tree->SetBranchAddress("dcaPrimeREC",dcaPrimeREC);
      tree->SetBranchAddress("dz0PrimeREC",dz0PrimeREC);
      //tree->SetBranchAddress("ptResFstSelected",ptResFstSelected);
      
      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);
         // very simple quality cuts
         if(TMath::Abs(vertex[2]+zvtxOffset)>30.) continue;
         if(Q2REC<Q2min || Q2REC>Q2max) continue;
         if(yREC<ymin) continue;
         if(yREC>ymax) continue;

         //loop over MC tracks and match to REC: efficiency
         for(int k=0;k<nMCtrack;k++) {
            double pt = TMath::Hypot(pxMC[k],pyMC[k]);
            if(pt<ptcut) continue;
            double eta=etaStarMC[k];
            if(eta < etaStarMin || eta > etaStarMax ) continue;
            h_genEtaStar->Fill(eta,pt,w);
            int irec=imatchMC[k];
            if(irec>=0) {
               //track quality cut:
               int type=typeChgREC[irec];
               if(type<0) type= -type;
               if( chi2nvREC[irec]/nvNdfREC[irec] > chi2nTrackCut ) continue;
               if(type!=4) {
                  if( chi2vtxREC[irec]/vtxNdfREC[irec] > chi2nVtxCut ) continue;
                  if( fabs(dcaPrimeREC[irec]) > dcaCut || fabs(dz0PrimeREC[irec]) > dcaCut ) continue;
                  if( vtxNHitsREC[irec] < NhitsCut || nvNHitsREC[irec] < NhitsCut ) continue;
                  if( vtxTrackLengthREC[irec] < trackLengthCut ) continue;
               }
               if( TMath::Hypot(pxREC[irec],pyREC[irec]) < ptcut ) continue;
               if( etaStarREC[irec] > etaStarMax || etaStarREC[irec] < etaStarMin ) continue;
               //if(type==4) h_genMatchEtaStar->Fill(eta,pt,w);
               //if(type==2||type==3) h_genMatchFwdCombCentEtaStar->Fill(eta,pt,w);
               if(type==1) h_genMatchCentEtaStar->Fill(eta,pt,w);
              
            }
         }
         
         //loop over REC tracks and does not match to MC: fake
         for(int j = 0; j<nRECtrack; j++){
            int type=typeChgREC[j];
            if(type<0) type= -type;
            //track quality cut:
            if( chi2nvREC[j]/nvNdfREC[j] > chi2nTrackCut ) continue;
            if(type!=4) {
               if( chi2vtxREC[j]/vtxNdfREC[j] > chi2nVtxCut ) continue;
               if( fabs(dcaPrimeREC[j]) > dcaCut || fabs(dz0PrimeREC[j]) > dcaCut ) continue;
               if( vtxNHitsREC[j] < NhitsCut || nvNHitsREC[j] < NhitsCut ) continue;
               if( vtxTrackLengthREC[j] < trackLengthCut ) continue;
            }
            double pt = TMath::Hypot(pxREC[j],pyREC[j]);
            if( pt<ptcut ) continue;
            double eta = etaStarREC[j];
            if( eta > etaStarMax || eta < etaStarMin ) continue;
            h_recEtaStar->Fill(eta,pt,w);
            int jrec=imatchREC[j];
            if(jrec<0) {
               // if(type==4)h_recMatchEtaStar->Fill(eta,w);
               // if(type==2||type==3) h_recMatchFwdCombCentEtaStar->Fill(eta,w);
               if(type==1) h_recMatchCentEtaStar->Fill(eta,pt,w);
            }
         }

         if(print) print--;
      }
   }

   TH2D* eff = (TH2D*) h_genMatchCentEtaStar->Clone("eff");
   eff->GetYaxis()->SetTitleOffset(1.6);
   eff->GetXaxis()->SetTitleOffset(1.6);
   eff->Divide(h_genEtaStar);
   //eff->Draw("lego2");

   TH2D* fak = (TH2D*) h_recMatchCentEtaStar->Clone("fak");
   fak->GetYaxis()->SetTitleOffset(1.6);
   fak->GetXaxis()->SetTitleOffset(1.6);
   fak->Divide(h_recEtaStar);
   //fak->Draw("lego2");

   TH2D* corr = (TH2D*) eff->Clone("corr");
   for(int i = 0; i < h_genMatchCentEtaStar->GetNbinsX(); i++){
      for(int j = 0; j < h_genMatchCentEtaStar->GetNbinsY(); j++){
         
         double efficiency = eff->GetBinContent(i+1,j+1);
         double fake = fak->GetBinContent(i+1,j+1);
         double correction = 1./(efficiency*(1.-fake));
         if(efficiency*(1.-fake) <= 0) {correction = 0.;};
         
         cout << "corr " << correction << endl;
         corr->SetBinContent(i+1,j+1, correction);
      }
   }

   TCanvas* c1 = new TCanvas("c1","c1",600,600);
   gPad->SetTicks();
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.13);
   gPad->SetLogz(1);
   corr->GetYaxis()->SetTitleOffset(1.8);
   corr->GetYaxis()->SetRangeUser(0,8.0);

   corr->GetXaxis()->SetTitleOffset(1.8);
   corr->GetXaxis()->SetRangeUser(-0.2,4);
   
   corr->GetZaxis()->SetTitleOffset(1.6);
   corr->GetZaxis()->SetTitle("1/(eff*(1-fak))");
   corr->Draw("lego2");

   c1->Print("../figures/Final/CentCorrectionTable.pdf");
}
