#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <iostream>

#include "RiceStyle.h"

void plotEtaNoFST(const int iq2=0, const int iy=0, const bool doComb_=true, const bool doFwd_= false) {
   
   TChain* tree = new TChain("properties");
   //tree->Add("../run/mc_9xxx_rapgap31/*.root");
   // tree->Add("../run/test_boost/ForwardMultAnalyzer_ESigma_cleanEenergy_rapgap.root");
   tree->Add("../run/test_boost/ForwardMultAnalyzer_ESigma_cleanEenergy.root");

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
   
   double Q2_bins[5]={5,10,20,40,100};
   double y_bins[5]={0.0375,0.075,0.15,0.3,0.6};

   TString Q2BINS[4]={"(5-10)","(10,20)","(20,40)","(40,100)"};
   TString YBINS[4]={"(0.0375-0.075)","(0.075,0.15)","(0.15,0.3)","(0.3,0.6)"};

   float Q2min=Q2_bins[iq2];
   float Q2max=Q2_bins[iq2+1];
   float ymin=y_bins[iy];
   float ymax=y_bins[iy+1];

   double zvtxOffset=0.;
   double ptcut=0.15;

   int print=100;

   if(tree) {
      
      Int_t vertexType;
      Float_t w;
      Float_t trigWeightAC;
      Float_t trigWeightRW;
      Float_t elecPxREC,elecPyREC;
      Float_t hfsEREC,hfsPzREC;
      Float_t elecEREC,elecPzREC;
      Float_t elecXclusREC,elecYclusREC, elecThetaREC,elecEnergyREC,elecEfracREC,elecHfracREC;
      Float_t elecEradREC,elecEcraREC;

      Float_t vertex[3];
      Float_t xREC,yREC,Q2REC;
      Int_t nMCtrack;
      Float_t ptStarMC[400];
      Float_t etaStarMC[400];
      Float_t pxMC[400];
      Float_t pyMC[400];
      Float_t etaMC[400];
      Int_t imatchMC[400];

      Int_t nRECtrack;
      Int_t typeChgREC[200];
      Float_t pxREC[200];
      Float_t pyREC[200];
      Float_t pREC[200];
      Float_t peREC[200];
      Float_t ptStarREC[200];
      Float_t phiStarREC[200];
      Float_t etaStarREC[200];
      Float_t etaREC[200];
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
      Float_t nucliaREC[400];

      Int_t imatchREC[400];

      Float_t startHitsRadiusREC[400];
      Float_t endHitsRadiusREC[400];
      Float_t trkThetaREC[400];
      Float_t chi2TrkREC[400];
      Int_t   ndfTrkREC[400];
      Float_t zLengthHitREC[400];
      Float_t chi2LinkREC[400];
      Int_t ndfLinkREC[400];
      Float_t rZeroREC[400];

      tree->SetBranchAddress("w",&w);
      tree->SetBranchAddress("vertexType",&vertexType);
      tree->SetBranchAddress("vertex",vertex);
      tree->SetBranchAddress("trigWeightAC",&trigWeightAC);
      tree->SetBranchAddress("trigWeightRW",&trigWeightRW);
      tree->SetBranchAddress("xREC",&xREC);
      tree->SetBranchAddress("yREC",&yREC);
      tree->SetBranchAddress("Q2REC",&Q2REC);
      tree->SetBranchAddress("elecPxREC",&elecPxREC);
      tree->SetBranchAddress("elecPyREC",&elecPyREC);
      tree->SetBranchAddress("elecPzREC",&elecPzREC);
      tree->SetBranchAddress("elecEREC",&elecEREC);
      tree->SetBranchAddress("elecEradREC",&elecEradREC);
      tree->SetBranchAddress("elecEcraREC",&elecEcraREC);

      tree->SetBranchAddress("hfsPzREC",&hfsPzREC);
      tree->SetBranchAddress("hfsEREC",&hfsEREC);

      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStarMC",etaStarMC);
      tree->SetBranchAddress("ptStarMC",ptStarMC);
      tree->SetBranchAddress("imatchMC",imatchMC);

      tree->SetBranchAddress("pxMC",&pxMC);
      tree->SetBranchAddress("pyMC",&pyMC);

      tree->SetBranchAddress("imatchREC",imatchREC);

      tree->SetBranchAddress("elecXclusREC",&elecXclusREC);
      tree->SetBranchAddress("elecYclusREC",&elecYclusREC);
      tree->SetBranchAddress("elecThetaREC",&elecThetaREC);
      tree->SetBranchAddress("elecEnergyREC",&elecEnergyREC);
      tree->SetBranchAddress("elecEfracREC",&elecEfracREC);
      tree->SetBranchAddress("elecHfracREC",&elecHfracREC);
      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStarMC",etaStarMC);
      tree->SetBranchAddress("ptStarMC",ptStarMC);
      tree->SetBranchAddress("pxMC",&pxMC);
      tree->SetBranchAddress("pyMC",&pyMC);
      tree->SetBranchAddress("etaMC",&etaMC);

      tree->SetBranchAddress("nRECtrack",&nRECtrack);
      tree->SetBranchAddress("typeChgREC",typeChgREC);
      tree->SetBranchAddress("pxREC",pxREC);
      tree->SetBranchAddress("pyREC",pyREC);
      tree->SetBranchAddress("pREC",pREC);
      tree->SetBranchAddress("peREC",peREC);
      tree->SetBranchAddress("etaStarREC",etaStarREC);
      tree->SetBranchAddress("etaREC",etaREC);
      tree->SetBranchAddress("phiStarREC",phiStarREC);
      tree->SetBranchAddress("ptStarREC",ptStarREC);
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
      tree->SetBranchAddress("nucliaREC",nucliaREC);

      tree->SetBranchAddress("startHitsRadiusREC",startHitsRadiusREC);
      tree->SetBranchAddress("endHitsRadiusREC",endHitsRadiusREC);
      tree->SetBranchAddress("trkThetaREC",trkThetaREC);
      tree->SetBranchAddress("chi2TrkREC",chi2TrkREC);
      tree->SetBranchAddress("ndfTrkREC",ndfTrkREC);
      tree->SetBranchAddress("zLengthHitREC",zLengthHitREC);
      tree->SetBranchAddress("chi2LinkREC",chi2LinkREC);
      tree->SetBranchAddress("ndfLinkREC",ndfLinkREC);
      tree->SetBranchAddress("rZeroREC",rZeroREC);
      
      //Generate a random number:
      TF1* rand = new TF1("rand","1",0,1);

      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);
         // very simple quality cuts
         if(vertexType != 1 ) continue;
         if(TMath::Abs(vertex[2]+zvtxOffset)>35.0) continue;
         if(Q2REC<Q2min || Q2REC>Q2max) continue;
         if(yREC<ymin) continue;
         if(yREC>ymax) continue;

         double Epz = hfsEREC+elecEREC - (hfsPzREC+elecPzREC);
         if( Epz > 70 || Epz < 35 ) continue;

         //scattered electron cuts
         if( sqrt(elecXclusREC*elecXclusREC+elecYclusREC*elecYclusREC) > 76. ) continue;

         //loop over MC tracks and match to REC: efficiency
         for(int j=0;j<nMCtrack;j++) {
            double eta=etaMC[j];
            if( TMath::Hypot(pxMC[j],pyMC[j])<ptcut ) continue;
            h_genEtaStar->Fill(eta,w);
            int irec=imatchMC[j];
            if(irec>=0) {
               //track quality cut:
               int type=typeChgREC[irec];
               if(type<0) type= -type;
               
               /*track quality cut
               1. central tracks
               2. combined tracks
               3. forward tracks
               */
               
               if( type == 1 ){
                  if( TMath::Hypot(pxREC[irec],pyREC[irec])<ptcut ) continue;
                  if( fabs( dcaPrimeREC[irec]*TMath::Sin(trkThetaREC[irec]) ) > 2.0 ) continue;
                  if( startHitsRadiusREC[irec] > 50.0 ) continue;
                  if( fabs(startHitsRadiusREC[irec] - endHitsRadiusREC[irec]) < 10. ) continue;
                  if( vtxNHitsREC[irec] < 0 ) continue;
                  if( fabs(trkThetaREC[irec] - elecThetaREC) < 0.2 ) continue;
                  if( rand->GetRandom() > nucliaREC[irec] && rand->GetRandom() < 1.0 ) continue; 
               }      
               else if( doComb_ && type == 2 ){
                  if( TMath::Hypot(pxREC[irec],pyREC[irec])<ptcut ) continue;
                  if( pREC[irec] < 0.5 ) continue;
                  if( TMath::RadToDeg()*trkThetaREC[irec] < 10. || TMath::RadToDeg()*trkThetaREC[irec] > 30. ) continue;
                  if( fabs(dcaPrimeREC[irec]) > 5.0 ) continue;
                  if( startHitsRadiusREC[irec] > 50.0 ) continue;
                  if( vtxNHitsREC[irec] < 0 ) continue;
                  if( peREC[irec]/pREC[irec] > 99999.9 ) continue;
                  if( chi2vtxREC[irec] > 50. ) continue;
                  if( chi2LinkREC[irec] > 50. ) continue;
                  if( rand->GetRandom() > nucliaREC[irec] && rand->GetRandom() < 1.0 ) continue; 
               }  
               else if( doFwd_ && type == 3 ){
                  if( TMath::Hypot(pxREC[irec],pyREC[irec])<ptcut ) continue;
                  if( pREC[irec] < 0.5 ) continue;
                  if( TMath::RadToDeg()*trkThetaREC[irec] < 6. || TMath::RadToDeg()*trkThetaREC[irec] > 25. ) continue;
                  if( startHitsRadiusREC[irec] > 25.0 ) continue;
                  if( fabs(zLengthHitREC[irec]) < 10. ) continue;
                  if( rZeroREC[irec] > 20. ) continue;
                  if( peREC[irec]/pREC[irec] > 99999.9 ) continue;
                  if( chi2vtxREC[irec] > 25. ) continue;
                  if( chi2TrkREC[irec] > 10. ) continue;
                  if( rand->GetRandom() > nucliaREC[irec] && rand->GetRandom() < 1.0 ) continue; 
               }  
               else{
                  continue;
               }  
               
               if(type==3) h_genMatchEtaStar->Fill(eta,w);
               if(type==2) h_genMatchFwdCombCentEtaStar->Fill(eta,w);
               if(type==1) h_genMatchCentEtaStar->Fill(eta,w);
              
            }
         }
         
         //loop over REC tracks and does not match to MC: fake
         for(int j = 0; j<nRECtrack; j++){
            int type=typeChgREC[j];
            if(type<0) type= -type;
            //track quality cut:
            
            /*track quality cut
            1. central tracks
            2. combined tracks
            3. forward tracks
            */
   
            if( type == 1 ){
               if( TMath::Hypot(pxREC[j],pyREC[j])<ptcut ) continue;
               if( fabs( dcaPrimeREC[j]*TMath::Sin(trkThetaREC[j]) ) > 2.0 ) continue;
               if( startHitsRadiusREC[j] > 50.0 ) continue;
               if( fabs(startHitsRadiusREC[j] - endHitsRadiusREC[j]) < 10. ) continue;
               if( vtxNHitsREC[j] < 0 ) continue;
               if( fabs(trkThetaREC[j] - elecThetaREC) < 0.2 ) continue;
               if( rand->GetRandom() > nucliaREC[j] && rand->GetRandom() < 1.0 ) continue; 
            }      
            else if( doComb_ && type == 2 ){
               if( TMath::Hypot(pxREC[j],pyREC[j])<ptcut ) continue;
               if( pREC[j] < 0.5 ) continue;
               if( TMath::RadToDeg()*trkThetaREC[j] < 10. || TMath::RadToDeg()*trkThetaREC[j] > 30. ) continue;
               if( fabs(dcaPrimeREC[j]) > 5.0 ) continue;
               if( startHitsRadiusREC[j] > 50.0 ) continue;
               if( vtxNHitsREC[j] < 0 ) continue;
               if( peREC[j]/pREC[j] > 99999.9 ) continue;
               if( chi2vtxREC[j] > 50. ) continue;
               if( chi2LinkREC[j] > 50. ) continue;
               if( rand->GetRandom() > nucliaREC[j] && rand->GetRandom() < 1.0 ) continue; 
            }  
            else if( doFwd_ && type == 3 ){
               if( TMath::Hypot(pxREC[j],pyREC[j])<ptcut ) continue;
               if( pREC[j] < 0.5 ) continue;
               if( TMath::RadToDeg()*trkThetaREC[j] < 6. || TMath::RadToDeg()*trkThetaREC[j] > 25. ) continue;
               if( startHitsRadiusREC[j] > 25.0 ) continue;
               if( fabs(zLengthHitREC[j]) < 10. ) continue;
               if( rZeroREC[j] > 20. ) continue;
               if( peREC[j]/pREC[j] > 99999.9 ) continue;
               if( chi2vtxREC[j] > 25. ) continue;
               if( chi2TrkREC[j] > 10. ) continue;
               if( rand->GetRandom() > nucliaREC[j] && rand->GetRandom() < 1.0 ) continue; 
            }  
            else{
               continue;
            }  

            double eta = etaREC[j];
            h_recEtaStar->Fill(eta, w);
            int jrec=imatchREC[j];
            if(jrec<0) {
               if(type==3) h_recMatchEtaStar->Fill(eta,w);
               if(type==2) h_recMatchFwdCombCentEtaStar->Fill(eta,w);
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
   //gPad->SetLogy(1);
   //gStyle->SetPadBorderMode(0.1);
   //gStyle->SetOptTitle(0);

   TH1D* base2 = makeHist("base2", "", "#eta", "Tracking Efficiency #times Acceptance", 100,-2,3.,kBlack);
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
   eff[3]->Draw("PSAME");

   
   TString zvtxShift;
   zvtxShift=TString::Format("|v_{z}+%.f|<35 cm", zvtxOffset);
   
   TLegend *legend1=new TLegend(0.45,0.73,0.8,0.86,Form("H1 tracks %.1f>Q2>%.1f, %.4f<y<%.4f",Q2max,Q2min,ymin,ymax));
   legend1->SetBorderSize(0);
   legend1->SetFillStyle(0);
   legend1->SetTextSize(0.025);
   legend1->AddEntry(eff[0],"gen tracks","L");
   // legend1->AddEntry(eff[4],"cent matched","P");
   legend1->AddEntry(eff[3],"comb(type=2) matched","P");
   legend1->AddEntry(eff[2],"cent(type=1) matched","P");
   legend1->AddEntry(eff[1],"fwd(type=3) matched","P");
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

   TLatex* r47 = new TLatex(0.18,0.85, "p_{T} > 0.15 GeV/c");
   r47->SetNDC();
   r47->SetTextSize(21);
   r47->SetTextFont(44);

   TLatex* r48 = new TLatex(0.16, 0.91, "ep 27.5#times920 GeV");
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

   TH1D* base3 = makeHist("base3", "", "#eta", "Fake rate", 100,-2,3,kBlack);
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
   fak[3]->Draw("PSAME");

   
   TLegend *legend=new TLegend(0.45,0.73,0.8,0.86,Form("H1 tracks %.1f>Q2>%.1f, %.4f<y<%.4f",Q2max,Q2min,ymin,ymax));
   legend->SetBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetTextSize(0.025);
   //legend->AddEntry(fak[0],"gen tracks","L");
   // legend->AddEntry(fak[4],"cent not matched","P");
   legend->AddEntry(fak[3],"comb(type=2) not matched","P");
   legend->AddEntry(fak[2],"cent(type=1) not matched","P");
   legend->AddEntry(fak[1],"fwd(type=3) not matched","P");
   legend->Draw("same");

   r44->Draw("same");
   r46->Draw("same");
   r47->Draw("same");
   r48->Draw("same");

   // eff[2]->Fit("pol0","","",-0.5,1.2);
   // TF1* myFunc = (TF1*) eff[2]->GetFunction("pol0");
   // double eff_value = myFunc->GetParameter(0);

   // fak[2]->Fit("pol0","","",-0.5,1.2);
   // TF1* myFunc1 = (TF1*) fak[2]->GetFunction("pol0");
   // double fak_value = myFunc1->GetParameter(0);

   // double ff = eff_value*(1-fak_value);

   // cout << "factor " << ff << endl;
  
   c1->Print("./eff_final_Q2_"+Q2BINS[iq2]+"_y_"+YBINS[iy]+".pdf");
   c2->Print("./fak_final_Q2_"+Q2BINS[iq2]+"_y_"+YBINS[iy]+".pdf");

}
