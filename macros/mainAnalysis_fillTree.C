#include "../header/mainAnalysis.h"

using namespace std;

int getPassFlag(vector<int> trackType, vector<double> cuts, int trackQuality){

   //Generate a random number:
      TF1* rand = new TF1("rand","1",0,1);

   //trackType   
   int type = trackType[0];   
   int doComb_ = trackType[1];
   int doFwd_ = trackType[2];

   //define cut variables
   double pt = cuts[0];
   double dcaPrime = cuts[1];
   double trkTheta = cuts[2];
   double startHitsRadius = cuts[3];
   double endHitsRadius = cuts[4];
   double vtxNHits = cuts[5];
   double elecTheta = cuts[6];
   double nuclia = cuts[7];
   double p = cuts[8];
   double pe = cuts[9];
   double chi2vtx = cuts[10];
   double chi2Link = cuts[11];
   double zLengthHit = cuts[12];
   double rZero = cuts[13];
   double chi2Trk = cuts[14];

   //define pass flag
   int pass = 0;

   //define track quality, trackQuality = 0 (tight), = 1 (default), = 2 (loose)
   double ptcut = 0.15;
   double cuts_1_value[3]={1.0,2.0,3.0};
   double cuts_2_value[3]={40.,50.,60.};
   double cuts_3_value[3]={15.,10.,7.};
   double cuts_4_value[3]={3.,5.,8.};

   if( type == 1 ){
      if( pt<ptcut ) pass = 0;
      if( fabs( dcaPrime*TMath::Sin(trkTheta) ) > cuts_1_value[trackQuality] ) pass = 0;
      if( startHitsRadius > cuts_2_value[trackQuality] ) pass = 0;
      if( fabs(startHitsRadius - endHitsRadius) < cuts_3_value[trackQuality] ) pass = 0;
      if( vtxNHits < 0 ) pass = 0;
      if( fabs(trkTheta - elecTheta) < 0.2 ) pass = 0;
      if( rand->GetRandom() > nuclia && rand->GetRandom() < 1.0 ) pass = 0; 
      //pass everything
      pass = 1;
   }      
   else if( doComb_ && type == 2 ){
      if( pt<ptcut ) pass = 0;
      if( p < 0.5 ) pass = 0;
      if( TMath::RadToDeg()*trkTheta < 10. || TMath::RadToDeg()*trkTheta > 30. ) pass = 0;
      if( fabs(dcaPrime) > cuts_4_value[trackQuality] ) pass = 0;
      if( startHitsRadius > cuts_2_value[trackQuality] ) pass = 0;
      if( vtxNHits < 0 ) pass = 0;
      if( pe/p > 99999.9 ) pass = 0;
      if( chi2vtx > 50. ) pass = 0;
      if( chi2Link > 50. ) pass = 0;
      if( rand->GetRandom() > nuclia && rand->GetRandom() < 1.0 ) pass = 0; 
      //pass everything
      pass = 1;
   }  
   else if( doFwd_ && type == 3 ){
      if( pt<ptcut ) pass = 0;
      if( p < 0.5 ) pass = 0;
      if( TMath::RadToDeg()*trkTheta < 6. || TMath::RadToDeg()*trkTheta > 25. ) pass = 0;
      if( startHitsRadius > 25.0 ) pass = 0;
      if( fabs(zLengthHit) < 10. ) pass = 0;
      if( rZero > 20. ) pass = 0;
      if( pe/p > 99999.9 ) pass = 0;
      if( chi2vtx > 25. ) pass = 0;
      if( chi2Trk > 10. ) pass = 0;
      if( rand->GetRandom() > nuclia && rand->GetRandom() < 1.0 ) pass = 0; 
      //pass everything
      pass = 1;
   }  
   else{
      pass = 0;
   }

   return pass;

}

struct MyEvent {
   // general information
   Int_t run_mini,evno_mini; // run and event number
   Float_t w_mini; // event weight

   // event quality .. not complete yet
   Float_t vertex_mini[3];
   // Monte Carlo information
   Float_t xGKI_mini,yGKI_mini,Q2GKI_mini;
   Float_t xMC_es_mini,yMC_es_mini,Q2MC_es_mini;

   enum {
      nMCtrack_MAX=400
   };
   // if there is no MC info, nMCtrack is set to zero
   Int_t nMCtrack_mini;
   Float_t pxMC_mini[nMCtrack_MAX];
   Float_t pyMC_mini[nMCtrack_MAX];
   Float_t pzMC_mini[nMCtrack_MAX];
   Float_t etaMC_mini[nMCtrack_MAX];
   Float_t chargeMC_mini[nMCtrack_MAX];

   Float_t ptStarMC_mini[nMCtrack_MAX];
   Float_t etaStarMC_mini[nMCtrack_MAX];
   Float_t phiStarMC_mini[nMCtrack_MAX];

   Int_t imatchMC_mini[nMCtrack_MAX];

   // reconstructed quantities
   Float_t xREC_es_mini,yREC_es_mini,Q2REC_es_mini;
   
   enum {
      nRECtrack_MAX=200
   };
 
   Int_t eventpass_mini;
   Int_t nRECtrack_mini;
   Int_t typeChgREC_mini[nRECtrack_MAX];

   Float_t pxREC_mini[nRECtrack_MAX];
   Float_t pyREC_mini[nRECtrack_MAX];
   Float_t pzREC_mini[nRECtrack_MAX];
   Float_t pREC_mini[nRECtrack_MAX];
   Float_t peREC_mini[nRECtrack_MAX];
   Float_t etaREC_mini[nRECtrack_MAX];
   Float_t phiREC_mini[nRECtrack_MAX];

   Float_t ptStarREC_mini[nRECtrack_MAX];
   Float_t etaStarREC_mini[nRECtrack_MAX];
   Float_t phiStarREC_mini[nRECtrack_MAX];

   Float_t nucliaREC_mini[nRECtrack_MAX];
   Float_t dmatchREC_mini[nRECtrack_MAX];
   Int_t imatchREC_mini[nRECtrack_MAX];
   Int_t passREC_mini[nRECtrack_MAX];
   Int_t passTightREC_mini[nRECtrack_MAX];
   Int_t passLooseREC_mini[nRECtrack_MAX];
};

void mainAnalysis_fillTree(const bool doGen_ = true, const bool doRapgap_ = true, const bool doReweight_ = false, const bool doComb_=true, const bool doFwd_= false) {
   
   TChain* tree = new TChain("properties");
   
   if( doRapgap_ && doGen_ ){
      tree->Add("../batch/output/mc_9299/*.root");
      tree->Add("../batch/output/mc_9300/*.root");
      tree->Add("../batch/output/mc_9301/*.root");
      tree->Add("../batch/output/mc_9302/*.root");
      tree->Add("../batch/output/mc_9303/*.root");
      tree->Add("../batch/output/mc_9304/*.root");
      tree->Add("../batch/output/mc_9305/*.root");
      tree->Add("../batch/output/mc_9306/*.root");
   }
   else if( !doRapgap_ && doGen_){
      tree->Add("../batch/output/mc_8926_4/*.root");
      tree->Add("../batch/output/mc_8927_4/*.root");
   }
   else if( !doGen_ ){
      tree->Add("../batch/output/data_highE_06_resubmit/*.root");
      tree->Add("../batch/output/data_highE_07_resubmit/*.root");
   }
   else{ cout << "no files" << endl;}
  
   //doReweighting on Q2,x,zvertex,pt,eta,Nch
   TH2D* DATA_Q2vsX = 0;
   TH2D* DATA_etaPt = 0;
   TH1D* DATA_vtxZ = 0;

   TH2D* MC_Q2vsX = 0;
   TH2D* MC_etaPt = 0;
   TH1D* MC_vtxZ = 0;

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

   //start to define new miniTree:
   TTree *outtree =new TTree("miniTree","miniTree");
   MyEvent myEvent;
   outtree->Branch("run_mini",&myEvent.run_mini,"run_mini/I");
   outtree->Branch("evno_mini",&myEvent.evno_mini,"evno_mini/I");
   outtree->Branch("w_mini",&myEvent.w_mini,"w_mini/F");
   outtree->Branch("vertex_mini",myEvent.vertex_mini,"vertex_mini[3]/F");

   outtree->Branch("xGKI_mini",&myEvent.xGKI_mini,"xGKI_mini/F");
   outtree->Branch("yGKI_mini",&myEvent.yGKI_mini,"yGKI_mini/F");
   outtree->Branch("Q2GKI_mini",&myEvent.Q2GKI_mini,"Q2GKI_mini/F");
   outtree->Branch("xMC_es_mini",&myEvent.xMC_es_mini,"xMC_es_mini/F");
   outtree->Branch("yMC_es_mini",&myEvent.yMC_es_mini,"yMC_es_mini/F");
   outtree->Branch("Q2MC_es_mini",&myEvent.Q2MC_es_mini,"Q2MC_es_mini/F");

   outtree->Branch("nMCtrack_mini",&myEvent.nMCtrack_mini,"nMCtrack_mini/I");

   outtree->Branch("pxMC_mini",myEvent.pxMC_mini,"pxMC_mini[nMCtrack_mini]/F");
   outtree->Branch("pyMC_mini",myEvent.pyMC_mini,"pyMC_mini[nMCtrack_mini]/F");
   outtree->Branch("pzMC_mini",myEvent.pzMC_mini,"pzMC_mini[nMCtrack_mini]/F");
   outtree->Branch("etaMC_mini",myEvent.etaMC_mini,"etaMC_mini[nMCtrack_mini]/F");
   outtree->Branch("chargeMC_mini",myEvent.chargeMC_mini,"chargeMC_mini[nMCtrack_mini]/F");
   outtree->Branch("ptStarMC_mini",myEvent.ptStarMC_mini,"ptStarMC_mini[nMCtrack_mini]/F");
   outtree->Branch("etaStarMC_mini",myEvent.etaStarMC_mini,"etaStarMC_mini[nMCtrack_mini]/F");
   outtree->Branch("phiStarMC_mini",myEvent.phiStarMC_mini,"phiStarMC_mini[nMCtrack_mini]/F");
   outtree->Branch("imatchMC_mini",myEvent.imatchMC_mini,"imatchMC_mini[nMCtrack_mini]/I");

   outtree->Branch("xREC_es_mini",&myEvent.xREC_es_mini,"xREC_es_mini/F");
   outtree->Branch("yREC_es_mini",&myEvent.yREC_es_mini,"yREC_es_mini/F");
   outtree->Branch("Q2REC_es_mini",&myEvent.Q2REC_es_mini,"Q2REC_es_mini/F");

   outtree->Branch("eventpass_mini",&myEvent.eventpass_mini,"eventpass_mini/I");
   
   outtree->Branch("nRECtrack_mini",&myEvent.nRECtrack_mini,"nRECtrack_mini/I");
   outtree->Branch("typeChgREC_mini",myEvent.typeChgREC_mini,"typeChgREC_mini[nRECtrack_mini]/I");
   
   outtree->Branch("pxREC_mini",myEvent.pxREC_mini,"pxREC_mini[nRECtrack_mini]/F");
   outtree->Branch("pyREC_mini",myEvent.pyREC_mini,"pyREC_mini[nRECtrack_mini]/F");
   outtree->Branch("pzREC_mini",myEvent.pzREC_mini,"pzREC_mini[nRECtrack_mini]/F");
   outtree->Branch("pREC_mini",myEvent.pREC_mini,"pREC_mini[nRECtrack_mini]/F");
   outtree->Branch("peREC_mini",myEvent.peREC_mini,"peREC_mini[nRECtrack_mini]/F");
   outtree->Branch("etaREC_mini",myEvent.etaREC_mini,"etaREC_mini[nRECtrack_mini]/F");
   outtree->Branch("phiREC_mini",myEvent.phiREC_mini,"phiREC_mini[nRECtrack_mini]/F");
   outtree->Branch("ptStarREC_mini",myEvent.ptStarREC_mini,"ptStarREC_mini[nRECtrack_mini]/F");
   outtree->Branch("etaStarREC_mini",myEvent.etaStarREC_mini,"etaStarREC_mini[nRECtrack_mini]/F");
   outtree->Branch("phiStarREC_mini",myEvent.phiStarREC_mini,"phiStarREC_mini[nRECtrack_mini]/F");

   outtree->Branch("nucliaREC_mini",myEvent.nucliaREC_mini,"nucliaREC_mini[nRECtrack_mini]/F");
   outtree->Branch("dmatchREC_mini",myEvent.dmatchREC_mini,"dmatchREC_mini[nRECtrack_mini]/F");
   outtree->Branch("imatchREC_mini",myEvent.imatchREC_mini,"imatchREC_mini[nRECtrack_mini]/I");
   outtree->Branch("passREC_mini",myEvent.passREC_mini,"passREC_mini[nRECtrack_mini]/I");
   outtree->Branch("passTightREC_mini",myEvent.passTightREC_mini,"passTightREC_mini[nRECtrack_mini]/I");
   outtree->Branch("passLooseREC_mini",myEvent.passLooseREC_mini,"passLooseREC_mini[nRECtrack_mini]/I");


   float Q2min=5.;
   float Q2max=100.;
   float ymin=0.05;
   float ymax=0.6;


   double zvtxOffset=0.;
   double ptcut=0.15;

   if(tree) {
      Int_t run, evno;
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
      Float_t phiStarMC[400];
      Float_t pxMC[400];
      Float_t pyMC[400];
      Float_t pzMC[400];
      Float_t etaMC[400];
      Float_t chargeMC[400];
 
      Int_t nRECtrack;
      Int_t typeChgREC[200];
      Float_t pxREC[200];
      Float_t pyREC[200];
      Float_t pzREC[200];
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
      Int_t dmatchREC[400];
      Int_t imatchMC[400];

      Float_t startHitsRadiusREC[400];
      Float_t endHitsRadiusREC[400];
      Float_t trkThetaREC[400];
      Float_t chi2TrkREC[400];
      Int_t   ndfTrkREC[400];
      Float_t zLengthHitREC[400];
      Float_t chi2LinkREC[400];
      Int_t ndfLinkREC[400];
      Float_t rZeroREC[400];

      tree->SetBranchAddress("run",&run);
      tree->SetBranchAddress("evno",&evno);
      tree->SetBranchAddress("xGKI",&xGKI);

      tree->SetBranchAddress("xGKI",&xGKI);
      tree->SetBranchAddress("yGKI",&yGKI);
      tree->SetBranchAddress("Q2GKI",&Q2GKI);
      tree->SetBranchAddress("simvertex",&simvertex);
      
      tree->SetBranchAddress("nMCtrack",&nMCtrack);
      tree->SetBranchAddress("etaStar2MC",etaStarMC);
      tree->SetBranchAddress("ptStar2MC",ptStarMC);
      tree->SetBranchAddress("phiStar2MC",phiStarMC);
      tree->SetBranchAddress("pxMC",&pxMC);
      tree->SetBranchAddress("pyMC",&pyMC);
      tree->SetBranchAddress("pzMC",&pzMC);
      tree->SetBranchAddress("etaMC",&etaMC);
      tree->SetBranchAddress("chargeMC",&chargeMC);
      tree->SetBranchAddress("xMC",&xMC);
      tree->SetBranchAddress("yMC",&yMC);
      tree->SetBranchAddress("Q2MC",&Q2MC);
      tree->SetBranchAddress("xMC_es",&xMC_es);
      tree->SetBranchAddress("yMC_es",&yMC_es);
      tree->SetBranchAddress("Q2MC_es",&Q2MC_es);
      tree->SetBranchAddress("elecEMC",&elecEMC);
      tree->SetBranchAddress("elecEradMC",&elecEradMC);
      tree->SetBranchAddress("imatchMC",&imatchMC);

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
      tree->SetBranchAddress("pzREC",pzREC);
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
      tree->SetBranchAddress("dmatchREC",dmatchREC);

      tree->SetBranchAddress("startHitsRadiusREC",startHitsRadiusREC);
      tree->SetBranchAddress("endHitsRadiusREC",endHitsRadiusREC);
      tree->SetBranchAddress("trkThetaREC",trkThetaREC);
      tree->SetBranchAddress("chi2TrkREC",chi2TrkREC);
      tree->SetBranchAddress("ndfTrkREC",ndfTrkREC);
      tree->SetBranchAddress("zLengthHitREC",zLengthHitREC);
      tree->SetBranchAddress("chi2LinkREC",chi2LinkREC);
      tree->SetBranchAddress("ndfLinkREC",ndfLinkREC);
      tree->SetBranchAddress("rZeroREC",rZeroREC);

      cout << "total. number of events = " << tree->GetEntries() << endl;
      int Nevents = 1000;
      for(int i=0;i<Nevents;i++) {
         tree->GetEntry(i);

         //assigning all reweights:
         double evt_weight = 1.;
         if( !doGen_ && !doReweight_ ){
            //this is data weighting for triggers
            evt_weight = w*trigWeightRW;
         }
         else if( doGen_ && doReweight_ ) {
            //reweighting MC to DATA distribution in phase space variables
            double Q2x_weight = DATA_Q2vsX->GetBinContent( DATA_Q2vsX->FindBin(xMC_es,Q2MC_es) );
            if( Q2x_weight == 0. ) Q2x_weight = 1.0;
            double vtxZ_weight = DATA_vtxZ->GetBinContent( DATA_vtxZ->FindBin( simvertex[2]) );
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


         myEvent.run_mini = run;
         myEvent.evno_mini = evno;
         myEvent.w_mini = evt_weight;

         myEvent.xGKI_mini = -999.;
         myEvent.yGKI_mini = -999.;
         myEvent.Q2GKI_mini = -999.;

         myEvent.xMC_es_mini = -999.;
         myEvent.yMC_es_mini = -999.;
         myEvent.Q2MC_es_mini = -999.;

         myEvent.nMCtrack_mini = -999;

         for(int j=0;j<400;j++) {
               myEvent.pxMC_mini[j] = -999.;
               myEvent.pyMC_mini[j] = -999.;
               myEvent.pzMC_mini[j] = -999.;
               myEvent.etaMC_mini[j] = -999.;
               myEvent.chargeMC_mini[j] = -999.;
               myEvent.ptStarMC_mini[j] = -999.;
               myEvent.etaStarMC_mini[j] = -999.;
               myEvent.phiStarMC_mini[j] = -999.;
               myEvent.imatchMC_mini[j] = -999;
         }
         
         if( doGen_ ){
            //generator level event selections:
            myEvent.xGKI_mini = xGKI;
            myEvent.yGKI_mini = yGKI;
            myEvent.Q2GKI_mini = Q2GKI;

            myEvent.xMC_es_mini = xMC_es;
            myEvent.yMC_es_mini = yMC_es;
            myEvent.Q2MC_es_mini = Q2MC_es;

            myEvent.nMCtrack_mini = nMCtrack;

            for(int j=0;j<nMCtrack;j++) {
               myEvent.pxMC_mini[j] = pxMC[j];
               myEvent.pyMC_mini[j] = pyMC[j];
               myEvent.pzMC_mini[j] = pzMC[j];
               myEvent.etaMC_mini[j] = etaMC[j];
               myEvent.chargeMC_mini[j] = chargeMC[j];
               myEvent.ptStarMC_mini[j] = ptStarMC[j];
               myEvent.etaStarMC_mini[j] = etaStarMC[j];
               myEvent.phiStarMC_mini[j] = phiStarMC[j];
               myEvent.imatchMC_mini[j] = imatchMC[j];
            }
          }//end doGen_

         int event_pass = 1;

         /**RECO level starts here both MC and DATA**/
         if( (doGen_ && trigWeightRW <= 0) || (!doGen_ && trigWeightAC <= 0) ) event_pass = 0; //require trigger fired
         if( vertexType != 1 ) event_pass = 0;

         double Epz = hfsEREC+elecEREC - (hfsPzREC+elecPzREC);

         //bkg finder bit cuts:
         if( (ibgREC & 1) != 0 ) event_pass = 0;
         if( (ibgREC & 2) != 0 ) event_pass = 0;
         if( (ibgREC & 4) != 0 ) event_pass = 0;
         if( (ibgREC & 8) != 0 ) event_pass = 0;
         if( (ibgREC & 16) != 0 ) event_pass = 0;
         if( (ibgREC & 32) != 0 ) event_pass = 0;
         if( (ibgREC & 64) != 0 ) event_pass = 0;   
         //E-pz cuts
         if( Epz > 70 || Epz < 35 ) event_pass = 0;
         //vertex cuts
         if(TMath::Abs(vertex[2]+zvtxOffset)>35.) event_pass = 0;
         //kinematic cuts are not included   
         //Cut electron spatial 
         if( TMath::Hypot(elecXclusREC,elecYclusREC) > 70. || TMath::Hypot(elecXclusREC,elecYclusREC) < 15. ) event_pass = 0;
         if( elecEREC < 12. ) event_pass = 0;
   
         //after all event selection cuts:
         myEvent.xREC_es_mini = xREC_es;
         myEvent.yREC_es_mini = yREC_es;
         myEvent.Q2REC_es_mini = Q2REC_es;

         myEvent.vertex_mini[0] = vertex[0];
         myEvent.vertex_mini[1] = vertex[1];
         myEvent.vertex_mini[2] = vertex[2];

         myEvent.eventpass_mini = event_pass;
         myEvent.nRECtrack_mini = nRECtrack;

         //loop over reconstructed tracks:
         for(int j = 0; j<nRECtrack; j++){
            int type=typeChgREC[j];            
            if(type<0) type= -type;
           
            double phi = TMath::ATan(pyREC[j]/pxREC[j]);
            if( pxREC[j] < 0 && pyREC[j] > 0 ) phi = phi+3.14;
            if( pxREC[j] < 0 && pyREC[j] < 0 ) phi = phi-3.14;

            int pass_default = 0;
            int pass_tight = 0;
            int pass_loose = 0; 
            /*track quality cut
            1. central tracks
            2. combined tracks
            3. forward tracks
            */

            double ptREC = TMath::Hypot(pxREC[j],pyREC[j]);
            vector<double> cutVar{ptREC,dcaPrimeREC[j],trkThetaREC[j],startHitsRadiusREC[j],endHitsRadiusREC[j],
               (double)vtxNHitsREC[j],elecThetaREC,nucliaREC[j],pREC[j],peREC[j],chi2vtxREC[j],chi2LinkREC[j],zLengthHitREC[j],
               rZeroREC[j],chi2TrkREC[j]}; 
            
            vector<int> trackType{type, (int)doComb_, (int)doFwd_ };

            pass_tight = getPassFlag(trackType, cutVar, 0);
            pass_default = getPassFlag(trackType, cutVar, 1);
            pass_loose = getPassFlag(trackType, cutVar, 2);

            //assign values to each branch on track levels:
            myEvent.typeChgREC_mini[j] = typeChgREC[j];

            myEvent.pxREC_mini[j] = pxREC[j];
            myEvent.pyREC_mini[j] = pyREC[j];
            myEvent.pzREC_mini[j] = pzREC[j];

            myEvent.pREC_mini[j] = pREC[j];
            myEvent.peREC_mini[j] = peREC[j];
            myEvent.etaREC_mini[j] = etaREC[j];
            myEvent.phiREC_mini[j] = phi;
            myEvent.etaStarREC_mini[j] = etaStarREC[j];
            myEvent.ptStarREC_mini[j] = ptStarREC[j];
            myEvent.phiStarREC_mini[j] = phiStarREC[j];
            myEvent.nucliaREC_mini[j] = nucliaREC[j];
            myEvent.dmatchREC_mini[j] = dmatchREC[j];
            myEvent.imatchREC_mini[j] = imatchREC[j];       
            myEvent.passREC_mini[j] = pass_default; 
            myEvent.passTightREC_mini[j] = pass_tight;
            myEvent.passLooseREC_mini[j] = pass_loose;      
         }
         outtree->Fill();
      }  
   }

   TString outfile_name;
   if( !doReweight_ ){
      if( doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_noReweight_Tree.root";
      }
      else if( !doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_noReweight_Tree.root";
      }
      else{
            outfile_name = "../new_output/data_highE_0607_noReweight_Tree.root";
      }

   }
   else{
      if( doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_fullReweight_Tree.root";
      }
      else if( !doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_fullReweight_Tree.root";
      }
      else{
            cout << "don't reweight data! " << endl;
            outfile_name = "../new_output/data_highE_0607_noReweight_Tree.root";
      }
   }
   
   TFile outfile( outfile_name, "RECREATE");
   
   outtree->Write();
   
}
