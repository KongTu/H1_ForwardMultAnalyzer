#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotEta(void) {
   TTree *tree;
   gFile->GetObject("properties",tree);
   TH1D *h_genEtaStar=new TH1D("h_genEtaStar",";eta*",100,-10.,10.);
   TH1D *h_recEtaStar=new TH1D("h_recEtaStar",";eta*",100,-10.,10.);
   TH1D *h_recCentEtaStar=new TH1D("h_recCentEtaStar",";eta*",100,-10.,10.);
   TH1D *h_recCombCentEtaStar=new TH1D("h_recCombCentEtaStar",";eta*",100,-10.,10.);
   TH1D *h_recFwdCombCentEtaStar=new TH1D("h_recFwdCombCentEtaStar",";eta*",100,-10.,10.);

   double Q2min=5.;
   double ymin=0.1;
   double ymax=0.6;

   if(tree) {
      Int_t vertexType;
      Float_t w;
      Float_t vertex[3];
      Float_t xREC,yREC,Q2REC;
      Int_t nMCtrack;
      Float_t ptStarMC[200];
      Float_t etaStarMC[200];
      Int_t nRECtrack;
      Int_t typeREC[200];
      Float_t ptStarREC[200];
      Float_t etaStarREC[200];
      tree->SetBranchAddress("w",&w);
      tree->SetBranchAddress("vertexType",&vertexType);
      tree->SetBranchAddress("vertex",vertex);
      tree->SetBranchAddress("xREC",&xREC);
      tree->SetBranchAddress("yREC",&yREC);
      tree->SetBranchAddress("Q2REC",&Q2REC);

      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStarMC",etaStarMC);
      tree->SetBranchAddress("ptStarMC",ptStarMC);

      tree->SetBranchAddress("nRECtrack",&nRECtrack);
      tree->SetBranchAddress("typeREC",typeREC);
      tree->SetBranchAddress("etaStarREC",etaStarREC);
      tree->SetBranchAddress("ptStarREC",ptStarREC);
      
      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);
         // very simple quality cuts
         if(TMath::Abs(vertex[2])>30.) continue;
         if(Q2REC<Q2min) continue;
         if(yREC<ymin) continue;
         if(yREC>ymax) continue;

         // generated eta* for different track types
         for(int k=0;k<nMCtrack;k++) {
            if(ptStarMC[k]<0.1) continue;
            h_genEtaStar->Fill(etaStarMC[k],w);
         }
         // reconstructed eta* for different track types
         for(int k=0;k<nRECtrack;k++) {
            if(ptStarREC[k]<0.1) continue;
            h_recEtaStar->Fill(etaStarREC[k],w);
            if(typeREC[k]>0) {
               if(typeREC[k]<=3) {
                  h_recFwdCombCentEtaStar->Fill(etaStarREC[k],w);
                  if(typeREC[k]<=2) {
                     h_recCombCentEtaStar->Fill(etaStarREC[k],w);
                     if(typeREC[k]==1) {
                        h_recCentEtaStar->Fill(etaStarREC[k],w);
                     }
                  }
               }
            }
         }
      }
   }
   TCanvas *c=new TCanvas("H1 tracks","H1 tracks",600,600);
   h_genEtaStar->SetLineWidth(2);
   h_genEtaStar->DrawClone();
   /*h_recEtaStar->SetFillStyle(1001);
   h_recEtaStar->SetLineColor(kBlue);
   h_recEtaStar->SetFillColor(kBlue);
   h_recEtaStar->Draw("SAME"); */
   h_recFwdCombCentEtaStar->SetFillStyle(1001);
   h_recFwdCombCentEtaStar->SetLineColor(kBlue-7);
   h_recFwdCombCentEtaStar->SetFillColor(kBlue-7);
   h_recFwdCombCentEtaStar->Draw("SAME");
   h_recCombCentEtaStar->SetFillStyle(1001);
   h_recCombCentEtaStar->SetLineColor(kCyan);
   h_recCombCentEtaStar->SetFillColor(kCyan);
   h_recCombCentEtaStar->Draw("SAME");
   h_recCentEtaStar->SetFillStyle(1001);
   h_recCentEtaStar->SetLineColor(kMagenta);
   h_recCentEtaStar->SetFillColor(kMagenta);
   h_recCentEtaStar->Draw("SAME");
   h_genEtaStar->Draw("SAME");
   TLegend *legend=new TLegend
      (0.3,0.65,0.6,0.85,TString::Format("H1 tracks Q2>%f %f<y<%f",
                                         Q2min,ymin,ymax));
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->AddEntry(h_genEtaStar,"gen tracks pt*>100 MeV","f");
   legend->AddEntry(h_recCentEtaStar,"central tracks","f");
   legend->AddEntry(h_recCombCentEtaStar,"combined tracks","f");
   legend->AddEntry(h_recFwdCombCentEtaStar,"forward tracks","f");
   legend->Draw();
}
