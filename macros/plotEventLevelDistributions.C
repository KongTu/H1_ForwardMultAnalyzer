#include "RiceStyle.h"
#include "mainAnalysis_fillTree.C"
#include "TLorentzVector.h"

using namespace std;

TH2D* h_SpaCalXY = new TH2D("h_SpaCalXY",";x;y",40,-80,80,40,-80,80);
TH1D* h_scatElecTheta = new TH1D("h_scatElecTheta", "#theta_{e}",100,2.5,3.2);
TH1D* h_scatElecPhi = new TH1D("h_scatElecPhi",";#phi",100,-3.2,3.2);
void plotEventLevelDistributions(const int start = 0, int end = -1, const bool doGen_ = false, const bool doRapgap_ = false, const bool doReweight_ = false){


//dirty hack of the problem of using TFile with TTree here. 
	TH2D* DATA_Q2vsX = new TH2D("h_Q2vsX_1","h_Q2vsX_1",1000,0.00001,0.01,50,1,100);
	TH1D* DATA_vtxZ = new TH1D("h_vtxZ_1","h_vtxZ_1", 100,-50,50);
	if( doRapgap_ ){
		for(int i = 0; i < DATA_vtxZ->GetNbinsX(); i++){
		 	DATA_vtxZ->SetBinContent(i+1,vtxz_weight_rapgap[i]);
		}
		int count = 0;
		for(int i = 0; i < DATA_Q2vsX->GetNbinsX(); i++){
		 for(int j = 0; j < DATA_Q2vsX->GetNbinsY(); j++){
		    DATA_Q2vsX->SetBinContent(i+1,j+1,Q2vsX_weight_rapgap[count]);
		    count++;
		 }
		}
	}
	else{
		for(int i = 0; i < DATA_vtxZ->GetNbinsX(); i++){
		 DATA_vtxZ->SetBinContent(i+1,vtxz_weight_django[i]);
		}
		int count = 0;
		for(int i = 0; i < DATA_Q2vsX->GetNbinsX(); i++){
		 for(int j = 0; j < DATA_Q2vsX->GetNbinsY(); j++){
		    DATA_Q2vsX->SetBinContent(i+1,j+1,Q2vsX_weight_django[count]);
		    count++;
		 }
		}
	}
//end dirty trick

	//starting TChain;
	TChain* tree = new TChain("properties");

	if( doRapgap_ && doGen_ ){
		tree->Add("../batch/output/mc_9299_hadCali/*.root");
		tree->Add("../batch/output/mc_9300_hadCali/*.root");
		tree->Add("../batch/output/mc_9301_hadCali/*.root");
		tree->Add("../batch/output/mc_9302_hadCali/*.root");
		tree->Add("../batch/output/mc_9303_hadCali/*.root");
		tree->Add("../batch/output/mc_9304_hadCali/*.root");
		tree->Add("../batch/output/mc_9305_hadCali/*.root");
		tree->Add("../batch/output/mc_9306_hadCali/*.root");
	}
	else if( !doRapgap_ && doGen_){
		tree->Add("../batch/output/mc_8926_hadCali/*.root");
		tree->Add("../batch/output/mc_8927_hadCali/*.root");
	}
	else if( !doGen_ ){
		tree->Add("../batch/output/data_highE_06_hadCali/*.root");
		tree->Add("../batch/output/data_highE_07_hadCali/*.root");
	}
	else{ cout << "no files" << endl;}
	
	if(tree) {
    //define variables in trees 
      Int_t vertexType;
      Float_t w;
      Float_t trigWeightAC;
      Float_t trigWeightRW;
      Float_t elecPxREC,elecPyREC,elecEREC,elecPzREC;
      Float_t hfsEREC,hfsPxREC, hfsPyREC, hfsPzREC;
      Float_t elecXclusREC,elecYclusREC, elecThetaREC,elecEnergyREC,elecEfracREC,elecHfracREC;
      Float_t elecEradREC,elecEcraREC;
      Int_t ibgREC;

      Float_t vertex[3];
      Float_t xREC_es,yREC_es,Q2REC_es;
      
      Float_t simvertex[3];
      Float_t elecEMC,elecEradMC;
      Float_t xMC_es,yMC_es,Q2MC_es;

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

      Float_t startHitsRadiusREC[400];
      Float_t endHitsRadiusREC[400];
      Float_t trkThetaREC[400];
      Float_t chi2TrkREC[400];
      Int_t   ndfTrkREC[400];
      Float_t zLengthHitREC[400];
      Float_t chi2LinkREC[400];
      Int_t ndfLinkREC[400];
      Float_t rZeroREC[400];

      tree->SetBranchAddress("simvertex",&simvertex);
      tree->SetBranchAddress("xMC_es",&xMC_es);
      tree->SetBranchAddress("yMC_es",&yMC_es);
      tree->SetBranchAddress("Q2MC_es",&Q2MC_es);

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
      tree->SetBranchAddress("xREC_es",&xREC_es);
      tree->SetBranchAddress("yREC_es",&yREC_es);
      tree->SetBranchAddress("Q2REC_es",&Q2REC_es);
      tree->SetBranchAddress("elecPxREC",&elecPxREC);
      tree->SetBranchAddress("elecPyREC",&elecPyREC);
      tree->SetBranchAddress("elecPzREC",&elecPzREC);
      tree->SetBranchAddress("elecEREC",&elecEREC);
      tree->SetBranchAddress("elecEradREC",&elecEradREC);
      tree->SetBranchAddress("elecEcraREC",&elecEcraREC);

      tree->SetBranchAddress("hfsPxREC",&hfsPxREC);
      tree->SetBranchAddress("hfsPyREC",&hfsPyREC);
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
    //end   

      
      cout << "total number of events ~ " << tree->GetEntries() << endl;
      int totalEvents = 0;
      if( end == -1 ) end = tree->GetEntries();
      for(int i=start;i<end;i++) {
         tree->GetEntry(i);
         //assigning all reweights:
         double evt_weight = 1.;
        //reweighting  
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
        //end reweighting
         /**RECO level starts here both MC and DATA**/
         if( (doGen_ && trigWeightRW <= 0) || (!doGen_ && trigWeightAC <= 0) ) continue; //require trigger fired
         if( vertexType != 1 ) continue;
         //bkg finder bit cuts:
         if( (ibgREC & 1) != 0 ) continue;
         if( (ibgREC & 2) != 0 ) continue;
         if( (ibgREC & 4) != 0 ) continue;
         if( (ibgREC & 8) != 0 ) continue;
         if( (ibgREC & 16) != 0 ) continue;
         if( (ibgREC & 32) != 0 ) continue;
         if( (ibgREC & 64) != 0 ) continue;  
         //kinematic cuts are not included   
         //Cut electron spatial 
         if( TMath::Hypot(elecXclusREC,elecYclusREC) > 70. || TMath::Hypot(elecXclusREC,elecYclusREC) < 15. ) continue;
         if( elecEREC < 12. ) continue;    
         //E-pz cuts
         double Epz = hfsEREC+elecEREC - (hfsPzREC+elecPzREC);
         if( Epz > 70 || Epz < 35 ) continue;
         //vertex cuts
         if(TMath::Abs(vertex[2])>35.) continue;
         //phase space cuts
         if( Q2REC_es < 5 || Q2REC_es > 100. ) continue;
         if( yREC_es < 0.075 || yREC_es > 0.6 ) continue;

         TLorentzVector elec;
         elec.SetPxPyPzE(elecPxREC,elecPyREC,elecPzREC,elecEREC);

         //after all event selection cuts:
         h_SpaCalXY->Fill( elecXclusREC, elecYclusREC, evt_weight*pow(Q2REC_es,2) );
         h_scatElecTheta->Fill( elecThetaREC, evt_weight);
         h_scatElecPhi->Fill(elec.Phi(), evt_weight);

         double Ntracks_eta_m = 0;
         double Ntracks_eta_p = 0;
         vector<double> cutVar;
         vector<int> trackType;
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
            /*
            track quality variables are stored using arrays below
            */
            double ptREC = TMath::Hypot(pxREC[j],pyREC[j]);

            cutVar.push_back(ptREC); cutVar.push_back(dcaPrimeREC[j]); cutVar.push_back(trkThetaREC[j]);
            cutVar.push_back(startHitsRadiusREC[j]); cutVar.push_back(endHitsRadiusREC[j]); cutVar.push_back((double)vtxNHitsREC[j]);
            cutVar.push_back(elecThetaREC); cutVar.push_back(nucliaREC[j]); cutVar.push_back(pREC[j]);
            cutVar.push_back(peREC[j]); cutVar.push_back(chi2vtxREC[j]); cutVar.push_back(chi2LinkREC[j]);
            cutVar.push_back(zLengthHitREC[j]); cutVar.push_back(rZeroREC[j]); cutVar.push_back(chi2TrkREC[j]);

            trackType.push_back(type); trackType.push_back( 1 );
            trackType.push_back( 0 );

            pass_tight = getPassFlag(trackType, cutVar, 0);
            pass_default = getPassFlag(trackType, cutVar, 1);
            pass_loose = getPassFlag(trackType, cutVar, 2);

            trackType.clear();
            cutVar.clear();
  
         }

         totalEvents++;
         if( totalEvents%100000 == 0 )cout << "Events ~ " << totalEvents << endl;

      }
      cout << "Number of events processed ~ " << totalEvents << endl;  
   }

   TString outfile_name;

   if( !doReweight_ ){
      if( doRapgap_ && doGen_ ){
            outfile_name = Form("../new_output/mc_highE_RAPGAP_noReweight_Hist.root");
      }
      else if( !doRapgap_ && doGen_ ){
            outfile_name = Form("../new_output/mc_highE_DJANGOH_noReweight_Hist.root");
      }
      else{
            outfile_name = Form("../new_output/data_highE_0607_noReweight_Hist.root");
      }
   }
   else{
      if( doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_RAPGAP_fullReweight_Hist.root";
      }
      else if( !doRapgap_ && doGen_ ){
            outfile_name = "../new_output/mc_highE_DJANGOH_fullReweight_Hist.root";
      }
      else{
            cout << "don't reweight data! " << endl;
            outfile_name = "../new_output/data_highE_0607_noReweight_Hist.root";
      }
   }
   
   TFile outfile( outfile_name, "RECREATE");

   h_SpaCalXY->Write();
   h_scatElecTheta->Write();
   h_scatElecPhi->Write();



}