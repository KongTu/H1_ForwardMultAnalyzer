#include "../header/mainAnalysis.h"

using namespace std;

void mainAnalysis_all(const bool doGenOnly_ = false, const bool doGen_ = true, const bool doRapgap_ = true, const bool doReweight_ = false, const bool doImprove_ = false, const bool doComb_=true, const bool doFwd_= false) {
   
   TChain* tree = new TChain("properties");
   
   if( doRapgap_ && doGen_ && !doGenOnly_ ){
      tree->Add("../batch/output/mc_9299/*.root");
      tree->Add("../batch/output/mc_9300/*.root");
      tree->Add("../batch/output/mc_9301/*.root");
      tree->Add("../batch/output/mc_9302/*.root");
      tree->Add("../batch/output/mc_9303/*.root");
      tree->Add("../batch/output/mc_9304/*.root");
      tree->Add("../batch/output/mc_9305/*.root");
      tree->Add("../batch/output/mc_9306/*.root");
      // tree->Add("../run/test_boost/ForwardMultAnalyzer_ESigma_cleanEenergy_rapgap.root");
   }
   else if( !doRapgap_ && doGen_ && !doGenOnly_){
      tree->Add("../batch/output/mc_8926_4/*.root");
      tree->Add("../batch/output/mc_8927_4/*.root");
      // tree->Add("../run/test_boost/ForwardMultAnalyzer_ESigma_cleanEenergy.root");
   }
   else if( !doGen_ && !doGenOnly_ ){
      tree->Add("../batch/output/data_highE_06_resubmit_baseline/*.root");
      tree->Add("../batch/output/data_highE_07_resubmit/*.root");
      // tree->Add("../data/test_boost/ForwardMultAnalyzer_confirm.root");
   }
   else if( doGen_ && doGenOnly_ && doRapgap_ ){
      tree->Add("../run/test_NONRAD/ForwardMultAnalyzer_NONRAD_RAPGAP.root");
   }
   else if( doGen_ && doGenOnly_ && !doRapgap_ ){
      tree->Add("../run/test_NONRAD/ForwardMultAnalyzer_NONRAD_DJANGOH.root");
   }

   for( int itr = 0; itr < 2; itr++){
      //event level variables
      h_Q2vsX[itr] = new TH2D(Form("h_Q2vsX_%d",itr),";xREC;Q2REC",1000,0.00001,0.01,50,1,100);
      h_yREC[itr] = new TH1D(Form("h_yREC_%d",itr),";xREC",1000,0,1);
      h_vtxZ[itr] = new TH1D(Form("h_vtxZ_%d",itr),";vtxZ (cm)",400,-200,200);
      h_vtxXY[itr] = new TH2D(Form("h_vtxXY_%d",itr),";x;y",400,-1,1,400,-1,1);
      h_EminusPz[itr] = new TH1D(Form("h_EminusPz_%d",itr),"E-pz",500,0,500);
      h_weights[itr] = new TH1D(Form("h_weights_%d",itr),"weights",1,0,1e10);
      h_weights_gen[itr] = new TH1D(Form("h_weights_gen_%d",itr),"weights",1,0,1e10);

      h_eEnergy[itr] = new TH1D(Form("h_eEnergy_%d",itr),";electron E",200,0,100);
      h_eTheta[itr] = new TH1D(Form("h_eTheta_%d",itr),";electron theta",200,0,3.2);
      h_hadESpacal[itr] = new TH1D(Form("h_hadESpacal_%d",itr),";electron had energy",200,0,5);
      h_emE[itr] = new TH1D(Form("h_emE_%d",itr),";electron ecal E",200,0,5);
      h_clusterRadius[itr] = new TH1D(Form("h_clusterRadius_%d",itr),";electron cluster radius",100,0,10);
      h_xyCluster[itr] = new TH2D(Form("h_xyCluster_%d",itr),";x (cm) ;y (cm)",40,-80,80,40,-80,80);

      //track level variables
      h_eta_gen[itr]=new TH1D(Form("h_eta_gen_%d",itr),";#eta generated",nEta,etaMin,etaMax);
      h_eta[itr]=new TH1D(Form("h_eta_%d",itr),";#eta",nEta,etaMin,etaMax);
      h_etaStar_gen_low[itr]=new TH1D(Form("h_etaStar_gen_low_%d",itr),";#eta* generated",7,etabins_low);
      h_etaStar_gen_high[itr]=new TH1D(Form("h_etaStar_gen_high_%d",itr),";#eta* generated",8,etabins_high);
      h_etaStar_low[itr]=new TH1D(Form("h_etaStar_low_%d",itr),";#eta*",7,etabins_low);
      h_etaStar_high[itr]=new TH1D(Form("h_etaStar_high_%d",itr),";#eta*",8,etabins_high);
      h_pt[itr]=new TH1D(Form("h_pt_%d",itr),";p_{T} (GeV/c)",120,0,12);
      h_ptStar[itr]=new TH1D(Form("h_ptStar_%d",itr),";p_{T}* (GeV/c)",120,0,12);
      h_phi[itr]=new TH1D(Form("h_phi_%d",itr),";#phi",100,-4,4);
      h_etaPt[itr]=new TH2D(Form("h_etaPt_%d",itr),";#eta;p_{T}",nEta,etaMin,etaMax,120,0,12);

      //track level cut variables
      h_chi2vtx[itr]=new TH1D(Form("h_chi2vtx_%d",itr),";chi2vtx",100,0,50);
      h_chi2trk[itr]=new TH1D(Form("h_chi2trk_%d",itr),";chi2trk",100,0,50);
      h_Nhits[itr]=new TH1D(Form("h_Nhits_%d",itr),";Nhits",10000,0,10000);
      h_dcaPrime[itr]=new TH1D(Form("h_dcaPrime_%d",itr),";dcaPrime",400,-20,20);
      h_dz0Prime[itr]=new TH1D(Form("h_dz0Prime_%d",itr),";dz0Prime",200,-100,100);
      h_trackLength[itr]=new TH1D(Form("h_trackLength_%d",itr),";trackLength",100,0,100);
      h_rStartHits[itr]=new TH1D(Form("h_rStartHits_%d",itr),";rStartHits",100,0,100);
      h_trkTheta[itr]=new TH1D(Form("h_trkTheta_%d",itr),";trkTheta",100,0,3.2);
      h_zLengthHit[itr]=new TH1D(Form("h_zLengthHit_%d",itr),";zLengthHit",200,-100,100);
      h_chi2Link[itr]=new TH1D(Form("h_chi2Link_%d",itr),";chi2Link",150,0,30);
      h_rZero[itr]=new TH1D(Form("h_rZero_%d",itr),";rZero",500,0,50);
      h_dcaPrimeSinTheta[itr]=new TH1D(Form("h_dcaPrimeSinTheta_%d",itr),";dcaPrimeSinTheta",200,-40,40);
      h_momRes[itr]=new TH1D(Form("h_momRes_%d",itr),";dp/p (GeV/c)",500,0,100);

   }

   //doReweighting on Q2,x,zvertex,pt,eta,Nch
   TH2D* DATA_Q2vsX = 0;
   TH2D* DATA_etaPt = 0;
   TH1D* DATA_vtxZ = 0;
   TH1D* h_Pn_all_data = 0;

   TH2D* MC_Q2vsX = 0;
   TH2D* MC_etaPt = 0;
   TH1D* MC_vtxZ = 0;
   TH1D* h_Pn_all_mc = 0;

   if( doReweight_ ){
      TFile* file_data = new TFile("../new_output/data_highE_0607_noReweight.root");
      DATA_Q2vsX = (TH2D*) file_data->Get("h_Q2vsX_1");
      DATA_Q2vsX->Scale( 1.0/DATA_Q2vsX->Integral() );
      DATA_etaPt = (TH2D*) file_data->Get("h_etaPt_1");
      DATA_etaPt->Scale( 1.0/DATA_etaPt->Integral() );
      DATA_vtxZ = (TH1D*) file_data->Get("h_vtxZ_1");
      DATA_vtxZ->Scale( 1.0/DATA_vtxZ->Integral() );

      TFile* file_mc = 0;
      if( doRapgap_ ){
         file_mc = new TFile("../new_output/mc_highE_RAPGAP_noReweight.root");
      }
      else{ 
         file_mc = new TFile("../new_output/mc_highE_DJANGOH_noReweight.root");
      }
      MC_Q2vsX = (TH2D*) file_mc->Get("h_Q2vsX_1");
      MC_Q2vsX->Scale( 1.0/MC_Q2vsX->Integral() );
      MC_etaPt = (TH2D*) file_mc->Get("h_etaPt_1");
      MC_etaPt->Scale( 1.0/MC_etaPt->Integral() );
      MC_vtxZ = (TH1D*) file_mc->Get("h_vtxZ_1");
      MC_vtxZ->Scale( 1.0/MC_vtxZ->Integral() );

      DATA_Q2vsX->Divide( MC_Q2vsX );
      DATA_etaPt->Divide( MC_etaPt );
      DATA_vtxZ->Divide( MC_vtxZ );
   }

   //initialize background finder histograms
   for(int bit = 0; bit < 10; bit++){
     
      int n = TMath::Power(2,bit);

      h_bkg_vtxZ[bit] = new TH1D(Form("h_bkg_vtxZ_bit_%d",n),";vtxZ",400,-200,200);
      h_bkg_Epz[bit] = new TH1D(Form("h_bkg_Epz_bit_%d",n),";Epz",400,0,200);     
      h_nobkg_vtxZ[bit] = new TH1D(Form("h_nobkg_vtxZ_bit_%d",n),";vtxZ",400,-200,200);
      h_nobkg_Epz[bit] = new TH1D(Form("h_nobkg_Epz_bit_%d",n),";Epz",400,0,200);   
   }

   float Q2min=5.;
   float Q2max=100.;
   float ymin=0.05;
   float ymax=0.6;
   float xmin=0.0001;
   float xmax=0.01;

   double zvtxOffset=0.;
   double ptcut=0.15;

   double q2_range[]={5.0,10.0,20.0,50.0,100.0};
   double y_range[]={-1.5,-0.5,0.5,1.5,2.5};

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
     
      Int_t ibgREC;

      Float_t vertex[3];
      Float_t xREC,yREC,Q2REC;
      Float_t xREC_es,yREC_es,Q2REC_es;
      
      Float_t simvertex[3];
      Float_t elecEMC,elecEradMC;
      Float_t xMC,yMC,Q2MC;
      Float_t xGKI,yGKI,Q2GKI;
      Float_t xMC_es,yMC_es,Q2MC_es;
      Int_t nMCtrack;
      Float_t ptStarMC[400];
      Float_t etaStarMC[400];
      Float_t pxMC[400];
      Float_t pyMC[400];
      Float_t etaMC[400];
 
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

      tree->SetBranchAddress("xGKI",&xGKI);
      tree->SetBranchAddress("yGKI",&yGKI);
      tree->SetBranchAddress("Q2GKI",&Q2GKI);
      tree->SetBranchAddress("simvertex",&simvertex);
      
      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStar2MC",etaStarMC);
      tree->SetBranchAddress("ptStar2MC",ptStarMC);
      tree->SetBranchAddress("pxMC",&pxMC);
      tree->SetBranchAddress("pyMC",&pyMC);
      tree->SetBranchAddress("etaMC",&etaMC);
      tree->SetBranchAddress("xMC",&xMC);
      tree->SetBranchAddress("yMC",&yMC);
      tree->SetBranchAddress("Q2MC",&Q2MC);
      tree->SetBranchAddress("xMC_es",&xMC_es);
      tree->SetBranchAddress("yMC_es",&yMC_es);
      tree->SetBranchAddress("Q2MC_es",&Q2MC_es);
      tree->SetBranchAddress("elecEMC",&elecEMC);
      tree->SetBranchAddress("elecEradMC",&elecEradMC);

      tree->SetBranchAddress("ibg",&ibgREC);
      tree->SetBranchAddress("elecXclusREC",&elecXclusREC);
      tree->SetBranchAddress("elecYclusREC",&elecYclusREC);
      tree->SetBranchAddress("elecThetaREC",&elecThetaREC);
      tree->SetBranchAddress("elecEnergyREC",&elecEnergyREC);
      tree->SetBranchAddress("elecEfracREC",&elecEfracREC);
      tree->SetBranchAddress("elecHfracREC",&elecHfracREC);
      
      tree->SetBranchAddress("w",&w);
      tree->SetBranchAddress("vertexType",&vertexType);
      tree->SetBranchAddress("vertex",vertex);
      tree->SetBranchAddress("trigWeightAC",&trigWeightAC);
      tree->SetBranchAddress("trigWeightRW",&trigWeightRW);
      tree->SetBranchAddress("xREC",&xREC);
      tree->SetBranchAddress("yREC",&yREC);
      tree->SetBranchAddress("Q2REC",&Q2REC);
      tree->SetBranchAddress("xREC_es",&xREC_es);
      tree->SetBranchAddress("yREC_es",&yREC_es);
      tree->SetBranchAddress("Q2REC_es",&Q2REC_es);
      tree->SetBranchAddress("elecPxREC",&elecPxREC);
      tree->SetBranchAddress("elecPyREC",&elecPyREC);
      tree->SetBranchAddress("elecPzREC",&elecPzREC);
      tree->SetBranchAddress("elecEREC",&elecEREC);
      tree->SetBranchAddress("elecEradREC",&elecEradREC);
      tree->SetBranchAddress("elecEcraREC",&elecEcraREC);

      tree->SetBranchAddress("hfsPzREC",&hfsPzREC);
      tree->SetBranchAddress("hfsEREC",&hfsEREC);
      
      tree->SetBranchAddress("nRECtrack",&nRECtrack);
      tree->SetBranchAddress("typeChgREC",typeChgREC);
      tree->SetBranchAddress("pxREC",pxREC);
      tree->SetBranchAddress("pyREC",pyREC);
      tree->SetBranchAddress("pREC",pREC);
      tree->SetBranchAddress("peREC",peREC);
      tree->SetBranchAddress("etaStar2REC",etaStarREC);
      tree->SetBranchAddress("etaREC",etaREC);
      tree->SetBranchAddress("phiStar2REC",phiStarREC);
      tree->SetBranchAddress("ptStar2REC",ptStarREC);
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
      
      tree->SetBranchAddress("imatchREC",imatchREC);

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

      double sum_weights_gen = 0.;
      double sum_weights_rec_gen = 0.;
      double sum_weights_rec_mis = 0.;
      int truth_index = -1;

      cout << "total. number of events = " << tree->GetEntries() << endl;
      for(int i=0;i<tree->GetEntries();i++) {
         tree->GetEntry(i);

         //assigning all reweights:
         double evt_weight = 1.;
         if( !doGen_ && !doReweight_ ){
            //this is data weighting for triggers
            evt_weight = w*trigWeightRW;
         }
         else if( doGen_ && doReweight_ ) {
            //reweighting MC to DATA distribution in phase space variables
            double Q2x_weight = DATA_Q2vsX->GetBinContent( DATA_Q2vsX->FindBin(xREC_es,Q2REC_es) );
            if( Q2x_weight == 0. ) Q2x_weight = 1.0;
            double vtxZ_weight = DATA_vtxZ->GetBinContent( DATA_vtxZ->FindBin( vertex[2]) );
            if( vtxZ_weight == 0. ) vtxZ_weight = 1.0;
            
            evt_weight = w*Q2x_weight*vtxZ_weight;
         }
         else if( doGen_ && !doReweight_ ) {
            //this is MC without anyreweighting
            evt_weight = w;
         }
         else{
            cout << "Error in config! " << endl;
         }

         if( doGen_ ){
            //generator level event selections:
            if(yMC_es>ymin && yMC_es<ymax){
               if(Q2MC_es>Q2min && Q2MC_es<Q2max){
            
                     //gen weights
                     sum_weights_gen += evt_weight;
                     //gen loop:          
                     for(int j=0;j<nMCtrack;j++) {
                        double eta=etaMC[j];
                        double etaStar=etaStarMC[j];
                        //do something to gen particle here 
                        if( eta > 2.5 || eta < -2.0 ) continue;
                        if(TMath::Hypot(pxMC[j],pyMC[j])<ptcut) continue;
                        // if( etaStar > 5.0 || etaStar < 0.0 ) continue;
                        // if( ptStarMC[j] > 10. || ptStarMC[j] < 0 ) continue;

                        h_eta_gen[1]->Fill( eta, evt_weight);

                        if( ptStarMC[j] > 0. && ptStarMC[j] < 1.0 ) h_etaStar_gen_low[1]->Fill( etaStar, evt_weight);
                        if( ptStarMC[j] > 1. && ptStarMC[j] < 10. ) h_etaStar_gen_high[1]->Fill( etaStar, evt_weight);      
                     }       
               }
            }
          }//end doGen_

         if( doGenOnly_ ) continue;
        
         /**RECO level starts here both MC and DATA**/
         double Epz = hfsEREC+elecEREC - (hfsPzREC+elecPzREC);

         //fill event level histogram before cuts:
         h_EminusPz[0]->Fill( Epz, evt_weight );
         h_Q2vsX[0]->Fill(xREC_es,Q2REC_es, evt_weight);
         h_yREC[0]->Fill(yREC_es, evt_weight);
         h_vtxZ[0]->Fill( vertex[2], evt_weight);
         h_vtxXY[0]->Fill(vertex[0],vertex[1],evt_weight);
   
         //trigger cuts:
         if( (doGen_ && trigWeightRW <= 0) || (!doGen_ && trigWeightAC <= 0) ) continue; //require trigger fired
         //vertex type cuts:
         if( vertexType != 1 ) continue;
         //bkg finder bit cuts:
         if( (ibgREC & 1) != 0 ) continue;
         if( (ibgREC & 2) != 0 ) continue;
         if( (ibgREC & 4) != 0 ) continue;
         if( (ibgREC & 8) != 0 ) continue;
         if( (ibgREC & 16) != 0 ) continue;
         if( (ibgREC & 32) != 0 ) continue;
         if( (ibgREC & 64) != 0 ) continue;
         //E-pz cuts:
         if( Epz > 70 || Epz < 35 ) continue;
         //vertex cuts:
         if(TMath::Abs(vertex[2]+zvtxOffset)>35.) continue;
         //phase space cuts:
         if(yREC_es<ymin) continue;
         if(yREC_es>ymax) continue;
         if(Q2REC_es<Q2min || Q2REC_es>Q2max) continue;
         //Cut electron spatial :
         if( TMath::Hypot(elecXclusREC,elecYclusREC) > 70. || TMath::Hypot(elecXclusREC,elecYclusREC) < 15. ) continue;
         if( elecEnergyREC < 12. ) continue;

         //filling all event histograms after cuts:
         for(int bit = 0; bit < 10; bit++){
               int n = TMath::Power(2,bit);
               if( (ibgREC & n) != 0 ) {
                  h_bkg_vtxZ[bit]->Fill( vertex[2] );
                  h_bkg_Epz[bit]->Fill( Epz );
               } 
               if( (ibgREC & n) == 0 ) {
                  h_nobkg_vtxZ[bit]->Fill( vertex[2] );
                  h_nobkg_Epz[bit]->Fill( Epz );
               } 
         }
         h_EminusPz[1]->Fill( Epz, evt_weight );
         h_Q2vsX[1]->Fill(xREC_es,Q2REC_es, evt_weight);
         h_yREC[1]->Fill(yREC_es, evt_weight);
         h_vtxZ[1]->Fill( vertex[2], evt_weight);
         h_vtxXY[1]->Fill(vertex[0],vertex[1],evt_weight);
   
         /*
         separate the events into signal group and background
         group. Signal group is within the phase space cut
         Background group is outside of reco phase space but within
         gen level phase space.
         */
         if( yMC_es > ymax || yMC_es < ymin || Q2MC_es > Q2max || Q2MC_es < Q2min ){
            truth_index = 0;
         }
         else{
            truth_index = 1;
         }

         //data is always 1.
         if( !doGen_ ) truth_index = 1;
         //no improve bin-by-bin correction applied.
         if( !doImprove_ ) truth_index = 1;

         //fill reco weights
         if( truth_index == 0 ){
            sum_weights_rec_mis += evt_weight;
         }
         else if (truth_index == 1){
            sum_weights_rec_gen += evt_weight;
         }
         else{ 
            cout << "something is wrong!" << endl;
         }
         

         //fill scattered electron distributions:
         h_xyCluster[truth_index]->Fill( elecXclusREC, elecYclusREC, evt_weight*pow(Q2REC_es,2) );
         h_emE[truth_index]->Fill( elecEfracREC, evt_weight );
         h_hadESpacal[truth_index]->Fill( elecHfracREC, evt_weight );
         h_eEnergy[truth_index]->Fill( elecEnergyREC, evt_weight );
         h_eTheta[truth_index]->Fill( elecThetaREC, evt_weight );
         h_clusterRadius[truth_index]->Fill( elecEcraREC, evt_weight );
         
         //loop over reconstructed tracks:
         for(int j = 0; j<nRECtrack; j++){
            int type=typeChgREC[j];
            if(type<0) type= -type;

            double etaStar = etaStarREC[j];
            double eta = etaREC[j];
            double pt = TMath::Hypot(pxREC[j],pyREC[j]);
            double phi = TMath::ATan(pyREC[j]/pxREC[j]);
            if( pxREC[j] < 0 && pyREC[j] > 0 ) phi = phi+3.14;
            if( pxREC[j] < 0 && pyREC[j] < 0 ) phi = phi-3.14;

            if( eta > 2.5 || eta < -2.0 ) continue;
            //if( etaStar > 5.0 || etaStar < 0.0 ) continue;
            //if( ptStarREC[j] > 10. || ptStarREC[j] < 0 ) continue;
            
            //filling the track distributions:
            if( type == 1 || (type == 2 && doComb_) || (type == 3 && doFwd_) ){
               h_chi2vtx[truth_index]->Fill(chi2vtxREC[j],evt_weight);
               h_chi2trk[truth_index]->Fill(chi2TrkREC[j],evt_weight);
               h_dcaPrime[truth_index]->Fill(dcaPrimeREC[j],evt_weight);
               h_dz0Prime[truth_index]->Fill(dz0PrimeREC[j],evt_weight);
               h_Nhits[truth_index]->Fill(vtxNHitsREC[j],evt_weight);
               h_trackLength[truth_index]->Fill(vtxTrackLengthREC[j],evt_weight);
              
               h_rStartHits[truth_index]->Fill( startHitsRadiusREC[j], evt_weight);
               h_trkTheta[truth_index]->Fill( trkThetaREC[j], evt_weight);
               h_zLengthHit[truth_index]->Fill( zLengthHitREC[j], evt_weight);
               h_chi2Link[truth_index]->Fill( chi2LinkREC[j], evt_weight);
               h_rZero[truth_index]->Fill( rZeroREC[j], evt_weight);
               h_dcaPrimeSinTheta[truth_index]->Fill(dcaPrimeREC[j]*TMath::Sin(trkThetaREC[j]), evt_weight);
               h_momRes[truth_index]->Fill(peREC[j]/pREC[j], evt_weight);
            }
            
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

            h_pt[truth_index]->Fill( pt, evt_weight );
            h_ptStar[truth_index]->Fill( ptStarREC[j], evt_weight );
            h_eta[truth_index]->Fill( etaREC[j], evt_weight );
            h_phi[truth_index]->Fill( phi, evt_weight );
            h_etaPt[truth_index]->Fill( etaREC[j], pt, evt_weight );

            if( ptStarREC[j] > 0. && ptStarREC[j] < 1.0 ) h_etaStar_low[truth_index]->Fill( etaStar, evt_weight );
            if( ptStarREC[j] > 1. && ptStarREC[j] < 10. ) h_etaStar_high[truth_index]->Fill( etaStar, evt_weight );
            
         }
      }  

      h_weights_gen[1]->Fill( sum_weights_gen );
      h_weights[1]->Fill( sum_weights_rec_gen );
      h_weights[0]->Fill( sum_weights_rec_mis );

      if( !doImprove_ ){
         cout << "The sum of weights in gen = " << sum_weights_gen << endl;
         cout << "The sum of weights in rec = " << sum_weights_rec_gen << endl;
      }
      else{
         cout << "The sum of weights in gen = " << sum_weights_gen << endl;
         cout << "The sum of weights in rec and gen = " << h_weights[1]->GetMean() << endl;
         cout << "The sum of weights in rec but not gen = " << h_weights[0]->GetMean() << endl;
      }
   }

   TString outfile_name;
   if( !doReweight_ ){
      if( doRapgap_ && doGen_ && !doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_noReweight.root";
      }
      else if( !doRapgap_ && doGen_ && !doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_noReweight.root";
      }
      else if( doRapgap_ && doGen_ && doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_NONRAD_noReweight.root";
      }
      else if( !doRapgap_ && doGen_ && doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_NONRAD_noReweight.root";
      }
      else{
            outfile_name = "../new_output/data_highE_0607_noReweight.root";
      }

   }
   else{
      if( doRapgap_ && doGen_ && !doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_fullReweight.root";
      }
      else if( !doRapgap_ && doGen_ && !doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_fullReweight.root";
      }
      else if( doRapgap_ && doGen_ && doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_NONRAD_fullReweight.root";
      }
      else if( !doRapgap_ && doGen_ && doGenOnly_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_NONRAD_fullReweight.root";
      }
      else{
            cout << "don't reweight data! " << endl;
            outfile_name = "../new_output/data_highE_0607_noReweight.root";
      }
   }
   
   TFile outfile( outfile_name, "RECREATE");
   
   for(int itr = 0; itr < 2; itr++){
      h_vtxZ[itr]->Write();
      h_vtxXY[itr]->Write();
      h_Q2vsX[itr]->Write();
      h_yREC[itr]->Write();
      h_EminusPz[itr]->Write();
      h_weights[itr]->Write();
      h_weights_gen[itr]->Write();

      h_eEnergy[itr]->Write();
      h_eTheta[itr]->Write();
      h_xyCluster[itr]->Write();
      h_emE[itr]->Write();
      h_hadESpacal[itr]->Write();
      h_clusterRadius[itr]->Write();

      h_pt[itr]->Write();
      h_ptStar[itr]->Write();
      h_eta[itr]->Write();
      h_etaStar_low[itr]->Write();
      h_etaStar_high[itr]->Write();
      h_phi[itr]->Write();
      h_etaPt[itr]->Write();
      h_eta_gen[itr]->Write();
      h_etaStar_gen_low[itr]->Write();
      h_etaStar_gen_high[itr]->Write();

      h_chi2vtx[itr]->Write();
      h_chi2trk[itr]->Write();
      h_Nhits[itr]->Write();
      h_dcaPrime[itr]->Write();
      h_dz0Prime[itr]->Write();
      h_trackLength[itr]->Write();
      h_rStartHits[itr]->Write();
      h_trkTheta[itr]->Write();
      h_zLengthHit[itr]->Write();
      h_chi2Link[itr]->Write();
      h_rZero[itr]->Write();
      h_dcaPrimeSinTheta[itr]->Write();
      h_momRes[itr]->Write();
   }
   for(int bit = 0; bit < 10; bit++){
      h_bkg_vtxZ[bit]->Write();
      h_bkg_Epz[bit]->Write();
      h_nobkg_vtxZ[bit]->Write();
      h_nobkg_Epz[bit]->Write();
   }
   
}
