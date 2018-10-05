#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>

#include "RiceStyle.h"

void plotEta() {
   
   TFile* file = new TFile("../run/ForwardMultAnalyzer_LowE_8833.root");
   TTree *tree = (TTree*) file->Get("properties");
   //gFile->GetObject("properties",tree);
   int nEta=70;
   int etaMin=-7;
   int etaMax=7.;
   TH1D *h_genEtaStar=new TH1D("h_genEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchEtaStar=new TH1D("h_genMatchEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchCentEtaStar=new TH1D("h_genMatchCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchCombCentEtaStar=new TH1D("h_genMatchCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_genMatchFwdCombCentEtaStar=new TH1D("h_genMatchFwdCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   
   TH1D *h_recEtaStar=new TH1D("h_recEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_recMatchEtaStar=new TH1D("h_recMatchEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_recMatchCentEtaStar=new TH1D("h_recMatchCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_recMatchCombCentEtaStar=new TH1D("h_recMatchCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   TH1D *h_recMatchFwdCombCentEtaStar=new TH1D("h_recMatchFwdCombCentEtaStar",";eta*",nEta,etaMin,etaMax);
   

   float Q2min=30.;
   float ymin=0.1;
   float ymax=0.6;

   double zvtxOffset=70.;

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
      Float_t ptStarREC[200];
      Float_t phiStarREC[200];
      Float_t etaStarREC[200];
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
      tree->SetBranchAddress("etaStarREC",etaStarREC);
      tree->SetBranchAddress("phiStarREC",phiStarREC);
      tree->SetBranchAddress("ptStarREC",ptStarREC);
      tree->SetBranchAddress("imatchREC",imatchREC);
      tree->SetBranchAddress("ptResFstSelected",ptResFstSelected);
      
      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);
         // very simple quality cuts
         if(TMath::Abs(vertex[2]+zvtxOffset)>30.) continue;
         if(Q2REC<Q2min) continue;
         if(yREC<ymin) continue;
         if(yREC>ymax) continue;

         //loop over MC tracks and match to REC: efficiency
         for(int k=0;k<nMCtrack;k++) {
            //double pt = sqrt(pxMC[k]*pxMC[k]+pyMC[k]*pyMC[k]);
            if(ptStarMC[k]< 0.1) continue;
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
         
         //loop over REC tracks and does not match to MC: fake
         for(int j = 0; j<nRECtrack; j++){
            
            if( ptStarREC[j] < 0.1) continue;
            double eta = etaStarREC[j];
            h_recEtaStar->Fill(eta, w);
            int jrec=imatchREC[j];
            if(jrec<0) {
               int type=typeChgREC[jrec];
               if(type<0) type= -type;
               h_recMatchEtaStar->Fill(eta,w);
               if(type<=3) h_recMatchFwdCombCentEtaStar->Fill(eta,w);
               if(type<=2) h_recMatchCombCentEtaStar->Fill(eta,w);
               if(type==1) h_recMatchCentEtaStar->Fill(eta,w);
            }
         }

         if(print) print--;
      }
   }

   TGraphAsymmErrors * eff[5];
   eff[0] = new TGraphAsymmErrors();
   eff[1] = new TGraphAsymmErrors();
   eff[2] = new TGraphAsymmErrors();
   eff[3] = new TGraphAsymmErrors();
   eff[4] = new TGraphAsymmErrors();

   TGraphAsymmErrors * fak[5];
   fak[0] = new TGraphAsymmErrors();
   fak[1] = new TGraphAsymmErrors();
   fak[2] = new TGraphAsymmErrors();
   fak[3] = new TGraphAsymmErrors();
   fak[4] = new TGraphAsymmErrors();

   TCanvas* c1 = new TCanvas("c1","c1",600,600);
   gPad->SetTicks();
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.13);
   gPad->SetLogy(1);
   //gStyle->SetPadBorderMode(0.1);
   //gStyle->SetOptTitle(0);

   TH1D* base2 = makeHist("base2", "", "#eta*", "Tracking Efficiency", 100,-2,1.3,kBlack);
   base2->GetYaxis()->SetRangeUser(0.001, 10);
   base2->GetXaxis()->SetTitleColor(kBlack);
   
   fixedFontHist1D(base2,1.1,1.25);

   base2->GetYaxis()->SetTitleOffset(1.3);
   base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
   base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
   base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
   base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
   base2->GetXaxis()->SetNdivisions(4,6,0);
   base2->GetYaxis()->SetNdivisions(4,6,0);

   base2->Draw();

   eff[0]->Divide(h_genEtaStar,h_genEtaStar,"cp");
   eff[1]->Divide(h_genMatchEtaStar,h_genEtaStar,"cp");
   eff[2]->Divide(h_genMatchFwdCombCentEtaStar,h_genEtaStar,"cp");
   eff[3]->Divide(h_genMatchCombCentEtaStar,h_genEtaStar,"cp");
   eff[4]->Divide(h_genMatchCentEtaStar,h_genEtaStar,"cp");

   //eff[0]->SetMarkerStyle(20);
   eff[0]->SetLineWidth(2);
   eff[0]->DrawClone("SAME");
   eff[1]->SetFillStyle(1001);
   eff[1]->SetLineColor(kBlue);
   eff[1]->SetMarkerColor(kBlue);
   eff[1]->SetFillColor(kBlue);
   eff[1]->SetMarkerStyle(20);
   eff[1]->Draw("PSAME");
   eff[2]->SetFillStyle(1001);
   eff[2]->SetLineColor(kGreen-2);
   eff[2]->SetMarkerColor(kGreen-2);
   eff[2]->SetFillColor(kGreen-2);
   eff[2]->SetMarkerStyle(20);
   eff[2]->Draw("PSAME");
   eff[3]->SetFillStyle(1001);
   eff[3]->SetLineColor(kRed);
   eff[3]->SetMarkerColor(kRed);
   eff[3]->SetFillColor(kRed);
   eff[3]->SetMarkerStyle(20);
   eff[3]->Draw("PSAME");
   eff[4]->SetFillStyle(1001);
   eff[4]->SetLineColor(kBlack);
   eff[4]->SetMarkerColor(kBlack);
   eff[4]->SetFillColor(kBlack);
   eff[4]->SetMarkerStyle(20);
   eff[4]->Draw("PSAME");
   
   TString zvtxShift;
   zvtxShift=TString::Format("|v_{z}+%.f|<30 cm", zvtxOffset);
   
   TLegend *legend1=new TLegend(0.50,0.73,0.8,0.86,Form("H1 tracks Q2>%.2f, %.2f<y<%.2f",Q2min,ymin,ymax));
   legend1->SetBorderSize(0);
   legend1->SetFillStyle(0);
   legend1->SetTextSize(0.025);
   legend1->AddEntry(eff[0],"gen tracks","L");
   legend1->AddEntry(eff[4],"cent matched","P");
   legend1->AddEntry(eff[3],"comb matched","P");
   legend1->AddEntry(eff[2],"fwd matched","P");
   legend1->AddEntry(eff[1],"FST matched","P");
   legend1->Draw("same");

   TLatex* r44 = new TLatex(0.84,0.91, "H1");
   r44->SetNDC();
   r44->SetTextSize(0.04);

   TLatex* r45 = new TLatex(0.86,0.91, "experiment");
   r45->SetNDC();
   r45->SetTextSize(21);
   r45->SetTextFont(53);

   TLatex* r46 = new TLatex(0.18,0.80, zvtxShift);
   r46->SetNDC();
   r46->SetTextSize(21);
   r46->SetTextFont(44);

   TLatex* r47 = new TLatex(0.18,0.85, "p_{T} > 0.1 GeV/c");
   r47->SetNDC();
   r47->SetTextSize(21);
   r47->SetTextFont(44);

   TLatex* r48 = new TLatex(0.16, 0.91, "ep 27.5#times460 GeV");
   r48->SetNDC();
   r48->SetTextSize(23);
   r48->SetTextFont(43);
   r48->SetTextColor(kBlack);

   r44->Draw("same");
   r46->Draw("same");
   r47->Draw("same");
   r48->Draw("same");

   TCanvas* c2 = new TCanvas("c2","c2",600,600);
   gPad->SetTicks();
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.13);
   gPad->SetLogy(1);
   //gStyle->SetPadBorderMode(0.1);
   //gStyle->SetOptTitle(0);

   TH1D* base3 = makeHist("base3", "", "#eta*", "Fake rate", 100,-2,1.3,kBlack);
   base3->GetYaxis()->SetRangeUser(0.001, 10);
   base3->GetXaxis()->SetTitleColor(kBlack);
   
   fixedFontHist1D(base3,1.1,1.25);

   base3->GetYaxis()->SetTitleOffset(1.3);
   base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
   base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
   base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
   base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
   base3->GetXaxis()->SetNdivisions(4,6,0);
   base3->GetYaxis()->SetNdivisions(4,6,0);

   base3->Draw();

   fak[0]->Divide(h_recEtaStar,h_recEtaStar,"cp");
   fak[1]->Divide(h_recMatchEtaStar,h_recEtaStar,"cp");
   fak[2]->Divide(h_recMatchFwdCombCentEtaStar,h_recEtaStar,"cp");
   fak[3]->Divide(h_recMatchCombCentEtaStar,h_recEtaStar,"cp");
   fak[4]->Divide(h_recMatchCentEtaStar,h_recEtaStar,"cp");

   //fak[0]->SetMarkerStyle(20);
   fak[0]->SetLineWidth(2);
   //fak[0]->DrawClone("SAME");
   fak[1]->SetFillStyle(1001);
   fak[1]->SetLineColor(kBlue);
   fak[1]->SetMarkerColor(kBlue);
   fak[1]->SetFillColor(kBlue);
   fak[1]->SetMarkerStyle(20);
   fak[1]->Draw("PSAME");
   fak[2]->SetFillStyle(1001);
   fak[2]->SetLineColor(kGreen-2);
   fak[2]->SetMarkerColor(kGreen-2);
   fak[2]->SetFillColor(kGreen-2);
   fak[2]->SetMarkerStyle(20);
   fak[2]->Draw("PSAME");
   fak[3]->SetFillStyle(1001);
   fak[3]->SetLineColor(kRed);
   fak[3]->SetMarkerColor(kRed);
   fak[3]->SetFillColor(kRed);
   fak[3]->SetMarkerStyle(20);
   fak[3]->Draw("PSAME");
   fak[4]->SetFillStyle(1001);
   fak[4]->SetLineColor(kBlack);
   fak[4]->SetMarkerColor(kBlack);
   fak[4]->SetFillColor(kBlack);
   fak[4]->SetMarkerStyle(20);
   fak[4]->Draw("PSAME");
   
   TLegend *legend=new TLegend(0.50,0.73,0.8,0.86,Form("H1 tracks Q2>%.2f, %.2f<y<%.2f",Q2min,ymin,ymax));
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetTextSize(0.025);
   //legend->AddEntry(fak[0],"gen tracks","L");
   legend->AddEntry(fak[4],"cent not matched","P");
   legend->AddEntry(fak[3],"comb not matched","P");
   legend->AddEntry(fak[2],"fwd not matched","P");
   legend->AddEntry(fak[1],"FST not matched","P");
   legend->Draw("same");

   r44->Draw("same");
   r46->Draw("same");
   r47->Draw("same");
   r48->Draw("same");
  
   // c1->Print("../figures/etaStar_LowE_vtx70_eff.pdf");
   // c2->Print("../figures/etaStar_LowE_vtx70_fak.pdf");
}
