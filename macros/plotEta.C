#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>

void plotEta(void) {
   TTree *tree;
   gFile->GetObject("properties",tree);
   int nEta=70;
   int etaMin=-7;
   int etaMax=7.;
   TH1D *h_genEtaStar=new TH1D("h_genEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchEtaStar=new TH1D("h_genMatchEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchCentEtaStar=new TH1D("h_genMatchCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchCombCentEtaStar=new TH1D("h_genMatchCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchFwdCombCentEtaStar=new TH1D("h_genMatchFwdCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_effEtaStar=new TH1D("h_effEtaStar",";#eta*;#epsilon",nEta,etaMin,etaMax);
   TH1D *h_effMatchEtaStar=new TH1D("h_effMatchEtaStar",";#eta*;#epsilon",
                                    nEta,etaMin,etaMax);
   TH1D *h_effMatchCentEtaStar=new TH1D("h_effMatchCentEtaStar",
                                        ";#eta*;#epsilon",nEta,etaMin,etaMax);
   TH1D *h_effMatchCombCentEtaStar=new TH1D("h_effMatchCombCentEtaStar",
                                            ";#eta*;#epsilon",nEta,etaMin,etaMax);
   TH1D *h_effMatchFwdCombCentEtaStar=new TH1D("h_effMatchFwdCombCentEtaStar",
                                               ";#eta*;#epsilon",nEta,etaMin,etaMax);

   double Q2min=30.;
   double ymin=0.1;
   double ymax=0.6;

   double zvtxOffset=0;

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
      Int_t imatchMC[400];
      Int_t nRECtrack;
      Int_t typeChgREC[200];
      //Float_t ptStarREC[200];
      //Float_t phiStarREC[200];
      //Float_t etaStarREC[200];
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

      tree->SetBranchAddress("nRECtrack",&nRECtrack);
      tree->SetBranchAddress("typeChgREC",typeChgREC);
      //tree->SetBranchAddress("etaStarREC",etaStarREC);
      //tree->SetBranchAddress("phiStarREC",phiStarREC);
      //tree->SetBranchAddress("ptStarREC",ptStarREC);
      
      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);
         // very simple quality cuts
         if(TMath::Abs(vertex[2]+zvtxOffset)>30.) continue;
         if(Q2REC<Q2min) continue;
         if(yREC<ymin) continue;
         if(yREC>ymax) continue;

         // discard events with certain orientation of
         //  FST acceptance  wrt electron direction
         //double phi=atan2(elecPyREC,elecPxREC);
         //if(print) {
         //   cout<<phi<<" "<<elecPyREC<<" "<<elecPxREC<<"\n";
         //}
         //if(fabs(phi)<1.) continue;
         
         // generated eta*
         //  possibly matched with a reconstructed track
         //  for different reco track types
         for(int k=0;k<nMCtrack;k++) {
            if(ptStarMC[k]<0.1) continue;
            double eta=etaStarMC[k];
            h_genEtaStar->Fill(eta,w);
            int irec=imatchMC[k];
            if(irec>=0) {
               int type=typeChgREC[irec];
               if(type<0) type= -type;
               h_genMatchEtaStar->Fill(eta,w);
               if(type<=3) h_genMatchFwdCombCentEtaStar->Fill(eta,w);
               if(type<=2) h_genMatchCombCentEtaStar->Fill(eta,w);
               if(type==1) h_genMatchCentEtaStar->Fill(eta,w);
            }
         }
         if(print) print--;
      }
   }
   TCanvas *c=new TCanvas("H1 tracks","H1 tracks",600,600);
   //c->Divide(1,2);
   //c->cd(1);

   h_effEtaStar->Divide(h_genEtaStar,h_genEtaStar);
   h_effMatchEtaStar->Divide(h_genMatchEtaStar,h_genEtaStar);
   h_effMatchFwdCombCentEtaStar->Divide(h_genMatchFwdCombCentEtaStar,h_genEtaStar);
   h_effMatchCombCentEtaStar->Divide(h_genMatchCombCentEtaStar,h_genEtaStar);
   h_effMatchCentEtaStar->Divide(h_genMatchCentEtaStar,h_genEtaStar);

   h_effEtaStar->SetLineWidth(2);
   h_effEtaStar->DrawClone();
   h_effMatchEtaStar->SetFillStyle(1001);
   h_effMatchEtaStar->SetLineColor(kBlue);
   h_effMatchEtaStar->SetFillColor(kBlue);
   h_effMatchEtaStar->Draw("SAME");
   h_effMatchFwdCombCentEtaStar->SetFillStyle(1001);
   h_effMatchFwdCombCentEtaStar->SetLineColor(kBlue-7);
   h_effMatchFwdCombCentEtaStar->SetFillColor(kBlue-7);
   h_effMatchFwdCombCentEtaStar->Draw("SAME");
   h_effMatchCombCentEtaStar->SetFillStyle(1001);
   h_effMatchCombCentEtaStar->SetLineColor(kCyan);
   h_effMatchCombCentEtaStar->SetFillColor(kCyan);
   h_effMatchCombCentEtaStar->Draw("SAME");
   h_effMatchCentEtaStar->SetFillStyle(1001);
   h_effMatchCentEtaStar->SetLineColor(kMagenta);
   h_effMatchCentEtaStar->SetFillColor(kMagenta);
   h_effMatchCentEtaStar->Draw("SAME");
   TString zvtxShift;
   if(zvtxOffset!=0.0) {
      zvtxShift=TString::Format(" |zvtx+%g|<30 cm", zvtxOffset);
   }
   TLegend *legend=new TLegend
      (0.28,0.65,0.6,0.85,TString::Format("H1 tracks Q2>%g %g<y<%g"+zvtxShift,
                                         Q2min,ymin,ymax));
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetTextSize(0.025);
   legend->AddEntry(h_effEtaStar,"gen tracks pt*>100 MeV","f");
   legend->AddEntry(h_effMatchCentEtaStar,"cent matched","f");
   legend->AddEntry(h_effMatchCombCentEtaStar,"comb matched","f");
   legend->AddEntry(h_effMatchFwdCombCentEtaStar,"fwd matched","f");
   legend->AddEntry(h_effMatchEtaStar,"FST matched","f");
   legend->Draw();
}
