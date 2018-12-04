#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>

#include "RiceStyle.h"

void plotPhiFST() {
   
   //TFile* file = new TFile("../run/ForwardMultAnalyzer_LowE_8833.root");
   TFile* file = new TFile("../run/ForwardMultAnalyzer_newVars.root");
   TTree *tree = (TTree*) file->Get("properties");
   //gFile->GetObject("properties",tree);
   int nEta=20;
   int etaMin=-4;
   int etaMax=4.;
   TH1D *h_genEtaStar=new TH1D("h_genEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_genMatchEtaStar=new TH1D("h_genMatchEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_genMatchCentEtaStar=new TH1D("h_genMatchCentEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_genMatchCombCentEtaStar=new TH1D("h_genMatchCombCentEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_genMatchFwdCombCentEtaStar=new TH1D("h_genMatchFwdCombCentEtaStar",";#phi",nEta,etaMin,etaMax);
   
   TH1D *h_recEtaStar=new TH1D("h_recEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_recMatchEtaStar=new TH1D("h_recMatchEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_recMatchCentEtaStar=new TH1D("h_recMatchCentEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_recMatchCombCentEtaStar=new TH1D("h_recMatchCombCentEtaStar",";#phi",nEta,etaMin,etaMax);
   TH1D *h_recMatchFwdCombCentEtaStar=new TH1D("h_recMatchFwdCombCentEtaStar",";#phi",nEta,etaMin,etaMax);
   

   float Q2min=5.;
   float Q2max=100.;
   float ymin=0.1;
   float ymax=0.6;
   float etaStarMin=0;
   float etaStarMax=1.5;

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
            //double pt = sqrt(pxMC[k]*pxMC[k]+pyMC[k]*pyMC[k]);
            if(TMath::Hypot(pxMC[k],pyMC[k])<ptcut) continue;
            double phi = TMath::ATan(pyMC[k]/pxMC[k]);
            if( pxMC[k] < 0 && pyMC[k] > 0 ) phi = phi+3.14;
            if( pxMC[k] < 0 && pyMC[k] < 0 ) phi = phi-3.14;
            double eta=etaStarMC[k];
            if( eta > etaStarMax || eta < etaStarMin ) continue;
            h_genEtaStar->Fill(phi,w);
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
               if(type==4) h_genMatchEtaStar->Fill(phi,w);
               if(type==2||type==3) h_genMatchFwdCombCentEtaStar->Fill(phi,w);
               if(type==1) h_genMatchCentEtaStar->Fill(phi,w);
              
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
            if( TMath::Hypot(pxREC[j],pyREC[j])<ptcut ) continue;
            double phi = TMath::ATan(pyREC[j]/pxREC[j]);
            if( pxREC[j] < 0 && pyREC[j] > 0 ) phi = phi+3.14;
            if( pxREC[j] < 0 && pyREC[j] < 0 ) phi = phi-3.14;
            double eta = etaStarREC[j];
            if( eta > etaStarMax || eta < etaStarMin ) continue;
            h_recEtaStar->Fill(phi,w);
            int jrec=imatchREC[j];
            if(jrec<0) {
               
               if(type==4)h_recMatchEtaStar->Fill(phi,w);
               if(type==2||type==3) h_recMatchFwdCombCentEtaStar->Fill(phi,w);
               if(type==1) h_recMatchCentEtaStar->Fill(phi,w);
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
   //gPad->SetLogy(1);
   //gStyle->SetPadBorderMode(0.1);
   //gStyle->SetOptTitle(0);

   TH1D* base2 = makeHist("base2", "", "#phi", "Tracking Efficiency", 100,-4,4,kBlack);
   base2->GetYaxis()->SetRangeUser(0.001, 1.4);
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
   eff[2]->Divide(h_genMatchCentEtaStar,h_genEtaStar,"cp");
   eff[3]->Divide(h_genMatchFwdCombCentEtaStar,h_genEtaStar,"cp");

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
   eff[2]->SetLineColor(kRed);
   eff[2]->SetMarkerColor(kRed);
   eff[2]->SetFillColor(kRed);
   eff[2]->SetMarkerStyle(20);
   eff[2]->Draw("PSAME");
   eff[3]->SetFillStyle(1001);
   eff[3]->SetLineColor(kGreen-3);
   eff[3]->SetMarkerColor(kGreen-3);
   eff[3]->SetFillColor(kGreen-3);
   eff[3]->SetMarkerStyle(20);
   //eff[3]->Draw("PSAME");

   
   TString zvtxShift;
   zvtxShift=TString::Format("|v_{z}+%.f|<30 cm", zvtxOffset);
   
   TLegend *legend1=new TLegend(0.50,0.73,0.8,0.86,Form("H1 tracks Q2>%.2f, %.2f<y<%.2f",Q2min,ymin,ymax));
   legend1->SetBorderSize(0);
   legend1->SetFillStyle(0);
   legend1->SetTextSize(0.025);
   legend1->AddEntry(eff[0],"gen tracks","L");
   // legend1->AddEntry(eff[4],"cent matched","P");
   //legend1->AddEntry(eff[3],"fwd(type=2,3) matched","P");
   legend1->AddEntry(eff[2],"cent(type=1) matched","P");
   legend1->AddEntry(eff[1],"FST(type=4) matched","P");
   legend1->Draw("same");

   TLatex* r44 = new TLatex(0.84,0.91, "H1");
   r44->SetNDC();
   r44->SetTextSize(0.04);

   TLatex* r45 = new TLatex(0.86,0.91, "experiment");
   r45->SetNDC();
   r45->SetTextSize(21);
   r45->SetTextFont(53);

   TLatex* r46 = new TLatex(0.18,0.75, zvtxShift);
   r46->SetNDC();
   r46->SetTextSize(21);
   r46->SetTextFont(44);

   TLatex* r47 = new TLatex(0.18,0.85, "p_{T} > 0.2 GeV/c");
   r47->SetNDC();
   r47->SetTextSize(21);
   r47->SetTextFont(44);

   TLatex* r47_1 = new TLatex(0.18,0.80, "0 < #eta* < 1.5");
   r47_1->SetNDC();
   r47_1->SetTextSize(21);
   r47_1->SetTextFont(44);

   TLatex* r48 = new TLatex(0.16, 0.91, "ep 27.5#times460 GeV");
   r48->SetNDC();
   r48->SetTextSize(23);
   r48->SetTextFont(43);
   r48->SetTextColor(kBlack);

   r44->Draw("same");
   r46->Draw("same");
   r47->Draw("same");
   r47_1->Draw("same");
   r48->Draw("same");

   TCanvas* c2 = new TCanvas("c2","c2",600,600);
   gPad->SetTicks();
   gPad->SetLeftMargin(0.15);
   gPad->SetBottomMargin(0.13);
   gPad->SetLogy(1);
   //gStyle->SetPadBorderMode(0.1);
   //gStyle->SetOptTitle(0);

   TH1D* base3 = makeHist("base3", "", "#phi", "Fake rate", 100,-4,4,kBlack);
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
   fak[2]->Divide(h_recMatchCentEtaStar,h_recEtaStar,"cp");
   fak[3]->Divide(h_recMatchFwdCombCentEtaStar,h_recEtaStar,"cp");

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
   fak[2]->SetLineColor(kRed);
   fak[2]->SetMarkerColor(kRed);
   fak[2]->SetFillColor(kRed);
   fak[2]->SetMarkerStyle(20);
   fak[2]->Draw("PSAME");
   fak[3]->SetFillStyle(1001);
   fak[3]->SetLineColor(kGreen-3);
   fak[3]->SetMarkerColor(kGreen-3);
   fak[3]->SetFillColor(kGreen-3);
   fak[3]->SetMarkerStyle(20);
   //fak[3]->Draw("PSAME");

   
   TLegend *legend=new TLegend(0.50,0.73,0.8,0.86,Form("H1 tracks Q2>%.2f, %.2f<y<%.2f",Q2min,ymin,ymax));
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetTextSize(0.025);
   //legend->AddEntry(fak[0],"gen tracks","L");
   // legend->AddEntry(fak[4],"cent not matched","P");
   //legend->AddEntry(fak[3],"fwd(type=2,3) not matched","P");
   legend->AddEntry(fak[2],"cent(type=1) not matched","P");
   legend->AddEntry(fak[1],"FST(type=4) not matched","P");
   legend->Draw("same");

   r44->Draw("same");
   r46->Draw("same");
   r47->Draw("same");
   r47_1->Draw("same");
   r48->Draw("same");

   // eff[2]->Fit("pol0","","",-0.5,1.2);
   // TF1* myFunc = (TF1*) eff[2]->GetFunction("pol0");
   // double eff_value = myFunc->GetParameter(0);

   // fak[2]->Fit("pol0","","",-0.5,1.2);
   // TF1* myFunc1 = (TF1*) fak[2]->GetFunction("pol0");
   // double fak_value = myFunc1->GetParameter(0);

   // double ff = eff_value*(1-fak_value);

   // cout << "factor " << ff << endl;
  
   c1->Print("../figures/Final/phiLab_eff_vtx70_2.pdf");
   c2->Print("../figures/Final/phiLab_fak_vtx70_2.pdf");
}
