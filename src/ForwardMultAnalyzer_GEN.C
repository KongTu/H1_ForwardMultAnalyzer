///////////////////////////////////////////////////////
// Search for J/Psi on MODS
//
// Author     : Ingo Strauch (strauch@mail.desy.de)
// Created    : 19.12.2000
// Redone by  : Gero Flucke (gero.flucke@desy.de)
//          and Bengt Wessling (bengt.wessling@desy.de)
// Last update: $Date: 2010/10/08 13:23:57 $
//          by: $Author: msteder $
//
///////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <set>

// ROOT includes
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

// for testing
#include "H1Skeleton/H1EventFiller.h"

// H1 OO includes
#include "H1Skeleton/H1Tree.h"
#include "H1Steering/H1StdCmdLine.h"
#include "H1Arrays/H1ArrayF.h"
#include "H1Pointers/H1IntPtr.h"
#include "H1Pointers/H1BytePtr.h"
#include "H1Pointers/H1ShortPtr.h"
#include "H1Pointers/H1FloatPtr.h"
#include "H1Mods/H1PartMCArrayPtr.h"
#include "H1Mods/H1PartMC.h"
#include "H1Mods/H1GetPartMCId.h"
#include "H1Mods/H1SelVertexArrayPtr.h"
#include "H1Mods/H1SelVertex.h"
#include "H1Mods/H1PartCandArrayPtr.h"
#include "H1Mods/H1PartCand.h"
#include "H1Mods/H1PartEm.h"
#include "H1Mods/H1PartSelTrack.h"
#include "H1Mods/H1InclHfsIterator.h"
#include "H1PhysUtils/H1NuclIACor.h"
#include "H1Mods/H1PartEmArrayPtr.h"

#include "H1Geom/H1DetectorStatus.h"
#include "H1Geom/H1DBManager.h"

#include "H1Tracks/H1FSTFittedTrack.h"
#include "H1Tracks/H1FSTFittedTrackArrayPtr.h"

#include "H1Tracks/H1CombinedFittedTrack.h"
#include "H1Tracks/H1CombinedFittedTrackArrayPtr.h"

#include "H1Tracks/H1ForwardFittedTrack.h"
#include "H1Tracks/H1ForwardFittedTrackArrayPtr.h"

//#include "H1Tracks/H1FSTTrackArrayPtr.h"


//#include "H1Mods/H1GkiInfoArrayPtr.h"
//#include "H1Mods/H1GkiInfo.h"
#include <TLorentzRotation.h>
#include "H1HadronicCalibration/H1HadronicCalibration.h"
#include "H1PhysUtils/H1MakeKine.h"

#include "elecCut.C"
#include "elecCut.h"

using namespace std;

static double const ME=0.0005109989461;
static double const M_CHARGED_PION=0.13957061;

static double const ELEC_ISOLATION_CONE=0.1;

bool floatEqual(double a,double b) {
   double d1=fabs(a-b);
   double d2=fabs(a+b);
   return (d1<=1.E-5*d2);
}

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());
   return boost;
}
//boost to HCM frame with kinematics from e-sigma method
TLorentzRotation BoostToHCM_es(TLorentzVector const &eBeam_lab,
                               TLorentzVector const &pBeam_lab,
                               TLorentzVector const &eScat_lab, 
                               double Q2_es, 
                               double y_es) {

   double escat_lab_es_E = (Q2_es)/(4.*eBeam_lab.E()) + eBeam_lab.E()*(1.-y_es);
   double b_par = 4.*eBeam_lab.E()*eBeam_lab.E()*(1.-y_es)/(Q2_es);
   double escat_lab_es_theta = TMath::ACos((1.-b_par)/(1.+b_par));
   
   double escat_lab_es_pz = sqrt(escat_lab_es_E*escat_lab_es_E - ME*ME)*TMath::Cos(escat_lab_es_theta);
   double escat_lab_es_pt = sqrt(escat_lab_es_E*escat_lab_es_E - ME*ME - escat_lab_es_pz*escat_lab_es_pz);
   double escat_lab_es_eta = -TMath::Log(TMath::Tan(escat_lab_es_theta/2.));
   double phi_elec = eScat_lab.Phi();

   TLorentzVector eScat_lab_ES;
   eScat_lab_ES.SetPtEtaPhiE(escat_lab_es_pt, escat_lab_es_eta, phi_elec, escat_lab_es_E);

   //same as before
   TLorentzVector q_lab=eBeam_lab - eScat_lab_ES;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;

   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-phi_elec);
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());

   return boost;

}

void GetKinematics(TLorentzVector const &ebeam,TLorentzVector const &pbeam,
                   TLorentzVector const &escat,
                   Float_t *x,Float_t *y,Float_t *Q2) {
   TLorentzVector q=ebeam-escat;
   *Q2= -q.Mag2();
   double pq=pbeam.Dot(q);
   *y= pq/pbeam.Dot(ebeam);
   *x= *Q2/(2.*pq);
}


struct MyEvent {
   // general information
   Int_t run,evno; // run and event number
   Float_t w; // event weight

      // trigger information
   UInt_t l1l2l3ac[4];
   UInt_t l1l2l3rw[4];
   UInt_t hasActualST;
   UInt_t hasRawST;
   Float_t trigWeightRW;
   Float_t trigWeightAC;

   // background finders
   UInt_t ibg;
   UInt_t ibgfm;
   UInt_t ibgam;
   UInt_t iqn;

   // event quality .. not complete yet
   Int_t vertexType;
   Float_t vertex[3];
   Float_t beamSpot[2];
   Float_t beamTilt[2];
   Float_t eProtonBeam,eElectronBeam;

   // Monte Carlo information
   Float_t eProtonBeamMC;
   Float_t eElectronBeamMC;
   Float_t simvertex[3];
   Float_t xGKI,yGKI,Q2GKI;

   Float_t elecPxMC,elecPyMC,elecPzMC,elecEMC,elecEradMC; // scattered electron
   Float_t elecEcraREC;
   Float_t xMC,yMC,Q2MC;
   Float_t xMC_es,yMC_es,Q2MC_es;

   enum {
      nMCtrack_MAX=400
   };
   // if there is no MC info, nMCtrack is set to zero
   Int_t nMCtrackAll;
   Int_t nMCtrack;
   Int_t idMC[nMCtrack_MAX];
   Int_t idxRad;

   Float_t pxMC[nMCtrack_MAX];
   Float_t pyMC[nMCtrack_MAX];
   Float_t pzMC[nMCtrack_MAX];
   Float_t etaMC[nMCtrack_MAX];
   Float_t chargeMC[nMCtrack_MAX];

   Float_t ptStarMC[nMCtrack_MAX];
   Float_t etaStarMC[nMCtrack_MAX];
   Float_t phiStarMC[nMCtrack_MAX];
   Float_t ptStar2MC[nMCtrack_MAX];
   Float_t etaStar2MC[nMCtrack_MAX];
   Float_t phiStar2MC[nMCtrack_MAX];

   Float_t log10zMC[nMCtrack_MAX];
   Int_t imatchMC[nMCtrack_MAX];

   // reconstructed quantities
   Float_t elecPxREC,elecPyREC,elecPzREC,elecEREC,elecEradREC; // scattered electron
   Float_t elecXclusREC,elecYclusREC, elecThetaREC,elecEnergyREC,elecEfracREC,elecHfracREC;
   Int_t elecTypeREC;

   Float_t xREC,yREC,Q2REC;
   Float_t xREC_es,yREC_es,Q2REC_es;
   Float_t hfsPxREC,hfsPyREC,hfsPzREC,hfsEREC; // hadronic final state
   enum {
      nRECtrack_MAX=200
   };
   // if there is no scattered electron, no tracks are saved either
   //  (nRECtrack=0)
   Int_t nRECtrackAll;
   Int_t nRECtrack;
   Int_t typeChgREC[nRECtrack_MAX];

   Float_t pxREC[nRECtrack_MAX];
   Float_t pyREC[nRECtrack_MAX];
   Float_t pzREC[nRECtrack_MAX];
   Float_t pREC[nRECtrack_MAX];
   Float_t peREC[nRECtrack_MAX];
   Float_t etaREC[nRECtrack_MAX];

   Float_t ptStarREC[nRECtrack_MAX];
   Float_t etaStarREC[nRECtrack_MAX];
   Float_t phiStarREC[nRECtrack_MAX];
   Float_t ptStar2REC[nRECtrack_MAX];
   Float_t etaStar2REC[nRECtrack_MAX];
   Float_t phiStar2REC[nRECtrack_MAX];
   Float_t chi2vtxREC[nRECtrack_MAX];
   Int_t   vtxNdfREC[nRECtrack_MAX];
   Int_t   vtxNHitsREC[nRECtrack_MAX];
   Float_t vtxTrackLengthREC[nRECtrack_MAX];
   Float_t dcaPrimeREC[nRECtrack_MAX];
   Float_t dz0PrimeREC[nRECtrack_MAX];

   //non vertex fitted parameter not used
   Float_t chi2nvREC[nRECtrack_MAX]; 
   Int_t   nvNdfREC[nRECtrack_MAX];
   Int_t   nvNHitsREC[nRECtrack_MAX];
   Float_t nvTrackLengthREC[nRECtrack_MAX];
   
   Float_t startHitsRadiusREC[nRECtrack_MAX];
   Float_t endHitsRadiusREC[nRECtrack_MAX];
   Float_t trkThetaREC[nRECtrack_MAX];
   Float_t chi2TrkREC[nRECtrack_MAX];
   Int_t   ndfTrkREC[nRECtrack_MAX];
   Float_t zLengthHitREC[nRECtrack_MAX];
   Float_t chi2LinkREC[nRECtrack_MAX];
   Int_t ndfLinkREC[nRECtrack_MAX];
   Float_t rZeroREC[nRECtrack_MAX];

   Float_t log10zREC[nRECtrack_MAX];
   Float_t nucliaREC[nRECtrack_MAX];
   Float_t dmatchREC[nRECtrack_MAX];
   Int_t imatchREC[nRECtrack_MAX];

   // FST tracks
   Int_t nRECfstFitted;

   // auxillary information, not to be saved
   H1PartMC const *partMC[nMCtrack_MAX];
   TVector3 momREC[nRECtrack_MAX];
   TMatrix covREC[nRECtrack_MAX];
};

class DummyFiller : public H1EventFiller {
public:
   virtual void   Fill(H1Event* event) { }
};

int main(int argc, char* argv[]) {
   // parse the command line
   H1StdCmdLine opts;
   opts.Parse(&argc, argv);

   // open run selection and detector status file
   // TString goodRunFileName("SelectedRuns_HighE0607_e+p_920.root");
   // TFile goodRunFile(goodRunFileName);
   // if(!goodRunFile.IsOpen()) {
   //    cerr<<"Error: could not open file "<<goodRunFileName<<"\n";
   //    return 2;
   // }
   // H1RunList* goodRunList
   //    = (H1RunList*) goodRunFile.Get("H1RunList");
   // if(!goodRunList) {
   //    cerr<<"Error: no runlist in file - return!\n";
   //    return 2;
   // }
   // H1DetectorStatus *detectorStatus
   //    = (H1DetectorStatus*)goodRunFile.Get("MyDetectorStatus");
   // if(!detectorStatus) {
   //    cerr<<"Error: no detector status in file - return!\n";
   //    return 3;
   // }

   // Load mODS/HAT files
   H1Tree::Instance()->Open();            // this statement must be there!

   //cout << H1Tree::Instance()->SelectHat("NumJPsi>0")
   //     << " events selected " << endl;

   TFile *file=new TFile(opts.GetOutput(), "RECREATE");

   /*
   Finding out more about ISR and FSR
   */
   TH2D* h_dPhi_theta_ISR=new TH2D("h_dPhi_theta_ISR",";#theta;#Delta#phi",150,0,7,300,-7,7);
   TH2D* h_dPhi_theta_FSR=new TH2D("h_dPhi_theta_FSR",";#theta;#Delta#phi",150,0,7,300,-7,7);
   TH2D* h_dPhi_theta_noR=new TH2D("h_dPhi_theta_noR",";#theta;#Delta#phi",150,0,7,300,-7,7);

   TH1D* h_Q2diff = new TH1D("h_Q2diff",";#DeltaQ2",1000,-100,100);
   TH1D* h_Xdiff = new TH1D("h_Xdiff",";#DeltaX",1000,-1,1);
   TH1D* h_Ydiff = new TH1D("h_Ydiff",";#DeltaY",1000,-1,1);

   TH1D* h_ISR_Q2diff = new TH1D("h_ISR_Q2diff",";#DeltaQ2",1000,-100,100);
   TH1D* h_ISR_Xdiff = new TH1D("h_ISR_Xdiff",";#DeltaX",1000,-1,1);
   TH1D* h_ISR_Ydiff = new TH1D("h_ISR_Ydiff",";#DeltaY",1000,-1,1);

   TH1D* h_FSR_Q2diff = new TH1D("h_FSR_Q2diff",";#DeltaQ2",1000,-100,100);
   TH1D* h_FSR_Xdiff = new TH1D("h_FSR_Xdiff",";#DeltaX",1000,-1,1);
   TH1D* h_FSR_Ydiff = new TH1D("h_FSR_Ydiff",";#DeltaY",1000,-1,1);

   TH1D* h_noR_Q2diff = new TH1D("h_noR_Q2diff",";#DeltaQ2",1000,-100,100);
   TH1D* h_noR_Xdiff = new TH1D("h_noR_Xdiff",";#DeltaX",1000,-1,1);
   TH1D* h_noR_Ydiff = new TH1D("h_noR_Ydiff",";#DeltaY",1000,-1,1);

   TTree *output=new TTree("properties","properties");
   MyEvent myEvent;
   output->Branch("run",&myEvent.run,"run/I");
   output->Branch("evno",&myEvent.evno,"evno/I");
   output->Branch("w",&myEvent.w,"w/F");
   output->Branch("vertexType",&myEvent.vertexType,"vertexType/I");
   output->Branch("vertex",myEvent.vertex,"vertex[3]/F");
   output->Branch("beamSpot",myEvent.beamSpot,"beamSpot[2]/F");
   output->Branch("beamTilt",myEvent.beamTilt,"beamTilt[2]/F");
   output->Branch("eProtonBeam",&myEvent.eProtonBeam,"eProtonBeam/F");
   output->Branch("eElectronBeam",&myEvent.eElectronBeam,"eElectronBeam/F");
   output->Branch("eProtonBeamMC",&myEvent.eProtonBeamMC,"eProtonBeamMC/F");
   output->Branch("eElectronBeamMC",&myEvent.eElectronBeamMC,"eElectronBeamMC/F");
   
   output->Branch("l1l2l3ac",myEvent.l1l2l3ac,"l1l2l3ac[4]/i");
   output->Branch("l1l2l3rw",myEvent.l1l2l3ac,"l1l2l3rw[4]/i");
   output->Branch("trigWeightAC",&myEvent.trigWeightAC,"trigWeightAC/F");
   output->Branch("trigWeightRW",&myEvent.trigWeightRW,"trigWeightRW/F");

   output->Branch("ibg",&myEvent.ibg,"ibg/I");
   output->Branch("ibgfm",&myEvent.ibgfm,"ibgfm/I");
   output->Branch("ibgam",&myEvent.ibgam,"ibgam/I");
   output->Branch("iqn",&myEvent.iqn,"iqn/I");

   output->Branch("simvertex",myEvent.simvertex,"simvertex[3]/F");
   output->Branch("elecEradMC",&myEvent.elecEradMC,"elecEradMC/F");
   output->Branch("elecPxMC",&myEvent.elecPxMC,"elecPxMC/F");
   output->Branch("elecPyMC",&myEvent.elecPyMC,"elecPyMC/F");
   output->Branch("elecPzMC",&myEvent.elecPzMC,"elecPzMC/F");
   output->Branch("elecEMC",&myEvent.elecEMC,"elecEMC/F");
   output->Branch("xGKI",&myEvent.xGKI,"xGKI/F");
   output->Branch("yGKI",&myEvent.yGKI,"yGKI/F");
   output->Branch("Q2GKI",&myEvent.Q2GKI,"Q2GKI/F");
   output->Branch("xMC",&myEvent.xMC,"xMC/F");
   output->Branch("yMC",&myEvent.yMC,"yMC/F");
   output->Branch("Q2MC",&myEvent.Q2MC,"Q2MC/F");
   output->Branch("xMC_es",&myEvent.xMC_es,"xMC_es/F");
   output->Branch("yMC_es",&myEvent.yMC_es,"yMC_es/F");
   output->Branch("Q2MC_es",&myEvent.Q2MC_es,"Q2MC_es/F");

   output->Branch("nMCtrackAll",&myEvent.nMCtrackAll,"nMCtrackAll/I");
   output->Branch("nMCtrack",&myEvent.nMCtrack,"nMCtrack/I");
   output->Branch("idMC",myEvent.idMC,"idMC[nMCtrack]/I");
   output->Branch("idxRad",&myEvent.idxRad,"idxRad/I");

   output->Branch("pxMC",myEvent.pxMC,"pxMC[nMCtrack]/F");
   output->Branch("pyMC",myEvent.pyMC,"pyMC[nMCtrack]/F");
   output->Branch("pzMC",myEvent.pzMC,"pzMC[nMCtrack]/F");
   output->Branch("etaMC",myEvent.etaMC,"etaMC[nMCtrack]/F");
   output->Branch("chargeMC",myEvent.chargeMC,"chargeMC[nMCtrack]/F");
   output->Branch("ptStarMC",myEvent.ptStarMC,"ptStarMC[nMCtrack]/F");
   output->Branch("etaStarMC",myEvent.etaStarMC,"etaStarMC[nMCtrack]/F");
   output->Branch("phiStarMC",myEvent.phiStarMC,"phiStarMC[nMCtrack]/F");
   output->Branch("ptStar2MC",myEvent.ptStar2MC,"ptStar2MC[nMCtrack]/F");
   output->Branch("etaStar2MC",myEvent.etaStar2MC,"etaStar2MC[nMCtrack]/F");
   output->Branch("phiStar2MC",myEvent.phiStar2MC,"phiStar2MC[nMCtrack]/F");
   
   output->Branch("log10zMC",myEvent.log10zMC,"log10zMC[nMCtrack]/F");
   output->Branch("imatchMC",myEvent.imatchMC,"imatchMC[nMCtrack]/I");

   output->Branch("elecEradREC",&myEvent.elecEradREC,"elecEradREC/F");
   output->Branch("elecPxREC",&myEvent.elecPxREC,"elecPxREC/F");
   output->Branch("elecPyREC",&myEvent.elecPyREC,"elecPyREC/F");
   output->Branch("elecPzREC",&myEvent.elecPzREC,"elecPzREC/F");
   output->Branch("elecEREC",&myEvent.elecEREC,"elecEREC/F");
   output->Branch("elecEcraREC",&myEvent.elecEcraREC,"elecEcraREC/F");
   output->Branch("elecXclusREC",&myEvent.elecXclusREC,"elecXclusREC/F");
   output->Branch("elecYclusREC",&myEvent.elecYclusREC,"elecYclusREC/F");
   output->Branch("elecThetaREC",&myEvent.elecThetaREC,"elecThetaREC/F");
   output->Branch("elecTypeREC",&myEvent.elecTypeREC,"elecTypeREC/I");
   output->Branch("elecEnergyREC",&myEvent.elecEnergyREC,"elecEnergyREC/F");
   output->Branch("elecEfracREC",&myEvent.elecEfracREC,"elecEfracREC/F");
   output->Branch("elecHfracREC",&myEvent.elecHfracREC,"elecHfracREC/F");

   output->Branch("xREC",&myEvent.xREC,"xREC/F");
   output->Branch("yREC",&myEvent.yREC,"yREC/F");
   output->Branch("Q2REC",&myEvent.Q2REC,"Q2REC/F");
   output->Branch("xREC_es",&myEvent.xREC_es,"xREC_es/F");
   output->Branch("yREC_es",&myEvent.yREC_es,"yREC_es/F");
   output->Branch("Q2REC_es",&myEvent.Q2REC_es,"Q2REC_es/F");
   output->Branch("hfsPxREC",&myEvent.hfsPxREC,"hfsPxREC/F");
   output->Branch("hfsPyREC",&myEvent.hfsPyREC,"hfsPyREC/F");
   output->Branch("hfsPzREC",&myEvent.hfsPzREC,"hfsPzREC/F");
   output->Branch("hfsEREC",&myEvent.hfsEREC,"hfsEREC/F");

   output->Branch("nRECtrackAll",&myEvent.nRECtrackAll,"nRECtrackAll/I");
   output->Branch("nRECtrack",&myEvent.nRECtrack,"nRECtrack/I");
   output->Branch("typeChgREC",myEvent.typeChgREC,"typeChgREC[nRECtrack]/I");
   
   output->Branch("pxREC",myEvent.pxREC,"pxREC[nRECtrack]/F");
   output->Branch("pyREC",myEvent.pyREC,"pyREC[nRECtrack]/F");
   output->Branch("pzREC",myEvent.pzREC,"pzREC[nRECtrack]/F");
   output->Branch("pREC",myEvent.pREC,"pREC[nRECtrack]/F");
   output->Branch("peREC",myEvent.peREC,"peREC[nRECtrack]/F");
   output->Branch("etaREC",myEvent.etaREC,"etaREC[nRECtrack]/F");

   output->Branch("ptStarREC",myEvent.ptStarREC,"ptStarREC[nRECtrack]/F");
   output->Branch("etaStarREC",myEvent.etaStarREC,"etaStarREC[nRECtrack]/F");
   output->Branch("phiStarREC",myEvent.phiStarREC,"phiStarREC[nRECtrack]/F");
   output->Branch("ptStar2REC",myEvent.ptStar2REC,"ptStar2REC[nRECtrack]/F");
   output->Branch("etaStar2REC",myEvent.etaStar2REC,"etaStar2REC[nRECtrack]/F");
   output->Branch("phiStar2REC",myEvent.phiStar2REC,"phiStar2REC[nRECtrack]/F");

   output->Branch("log10zREC",myEvent.log10zREC,"log10zREC[nRECtrack]/F");
   output->Branch("chi2vtxREC",myEvent.chi2vtxREC,"chi2vtxREC[nRECtrack]/F");
   output->Branch("chi2nvREC",myEvent.chi2nvREC,"chi2nvREC[nRECtrack]/F");
   output->Branch("vtxNdfREC",myEvent.vtxNdfREC,"vtxNdfREC[nRECtrack]/I");
   output->Branch("nvNdfREC",myEvent.nvNdfREC,"nvNdfREC[nRECtrack]/I");
   output->Branch("vtxNHitsREC",myEvent.vtxNHitsREC,"vtxNHitsREC[nRECtrack]/I");
   output->Branch("nvNHitsREC",myEvent.nvNHitsREC,"nvNHitsREC[nRECtrack]/I");
   output->Branch("vtxTrackLengthREC",myEvent.vtxTrackLengthREC,"vtxTrackLengthREC[nRECtrack]/F");
   output->Branch("nvTrackLengthREC",myEvent.nvTrackLengthREC,"nvTrackLengthREC[nRECtrack]/F");
   output->Branch("dcaPrimeREC",myEvent.dcaPrimeREC,"dcaPrimeREC[nRECtrack]/F");
   output->Branch("dz0PrimeREC",myEvent.dz0PrimeREC,"dz0PrimeREC[nRECtrack]/F");
   
   output->Branch("startHitsRadiusREC",myEvent.startHitsRadiusREC,"startHitsRadiusREC[nRECtrack]/F");
   output->Branch("endHitsRadiusREC",myEvent.endHitsRadiusREC,"endHitsRadiusREC[nRECtrack]/F");
   output->Branch("trkThetaREC",myEvent.trkThetaREC,"trkThetaREC[nRECtrack]/F");
   output->Branch("chi2TrkREC",myEvent.chi2TrkREC,"chi2TrkREC[nRECtrack]/F");
   output->Branch("ndfTrkREC",myEvent.ndfTrkREC,"ndfTrkREC[nRECtrack]/I");
   output->Branch("zLengthHitREC",myEvent.zLengthHitREC,"zLengthHitREC[nRECtrack]/F");
   output->Branch("chi2LinkREC",myEvent.chi2LinkREC,"chi2LinkREC[nRECtrack]/F");
   output->Branch("ndfLinkREC",myEvent.ndfLinkREC,"ndfLinkREC[nRECtrack]/I");
   output->Branch("rZeroREC",myEvent.rZeroREC,"rZeroREC[nRECtrack]/F");

   output->Branch("nucliaREC",myEvent.nucliaREC,"nucliaREC[nRECtrack]/F");
   output->Branch("dmatchREC",myEvent.dmatchREC,"dmatchREC[nRECtrack]/F");
   output->Branch("imatchREC",myEvent.imatchREC,"imatchREC[nRECtrack]/I");

   output->Branch("nRECfstFitted",&myEvent.nRECfstFitted,"nRECfstFitted/I");

   H1ShortPtr runtype("RunType"); // 0=data, 1=MC, 2=CERN test, 3=CERN MC test
   H1FloatPtr beamx0("BeamX0");          // x position of beam spot (at z=0)
   H1FloatPtr beamy0("BeamY0");          // y position of beam spot (at z=0)

   H1FloatPtr eBeamP("EBeamP"); // proton beam energy from DMIS
   H1FloatPtr eBeamE("EBeamE"); // electron beam energy from DMIS
   H1FloatPtr beamtiltx0("BeamTiltX0");      // x slope of beam axis
   H1FloatPtr beamtilty0("BeamTiltY0");      // y slope of beam axis
   H1IntPtr run("RunNumber");
   H1IntPtr evno("EventNumber");
  
   H1BytePtr l1l2l3ac("Il1l2l3ac");
   H1BytePtr l1l2rw("Il1l2rw");
   H1BytePtr l1l3rw("Il1l3rw");

   H1IntPtr   ibg("Ibg");  // standard array of background finders from OM group (bit packed)
   H1IntPtr  ibgfm("Ibgfm");  // extra finders from OM  (bit packed)
   H1IntPtr  ibgam("Ibgam");  // background finders from DUK group (bit packed)
   H1IntPtr  iqn("Iqn");           // The Lar coherent noise flag

   H1FloatPtr Q2Gki("Q2Gki");
   H1FloatPtr xGki("XGki");
   H1FloatPtr yGki("YGki");
   H1FloatPtr weight1("Weight1");
   H1FloatPtr weight2("Weight2");

   H1FloatPtr genEnElec("GenEnElec"); //   Electron energy, combined with photon for FSR
   H1FloatPtr genPhElec("GenPhElec"); //   Electron phi, combined with photon for FSR
   H1FloatPtr genThElec("GenThElec"); //   Electron theta, combined with photon for FSR

   H1ShortPtr ivtyp("Ivtyp"); // vertex type
   H1SelVertexArrayPtr vertex; // all good vertices
   //H1PartCandArrayPtr partCand; // all good tracks
   H1PartCandArrayPtr partCandArray; // all good tracks

   H1PartMCArrayPtr mcpart;

   H1FSTFittedTrackArrayPtr fstFittedTrack;
   //H1FSTTrackArrayPtr fstNonFittedTrack;

   Int_t eventCounter = 0;

   H1HadronicCalibration *hadronicCalibration=H1HadronicCalibration::Instance();
   hadronicCalibration->ApplyHadronicCalibration(H1HadronicCalibration::eHighPtJet);
   hadronicCalibration->ApplyHadronicCalibration(kTRUE);

   // Loop over events
   static int print=10;
   while (gH1Tree->Next() && !opts.IsMaxEvent(eventCounter)) {
      ++eventCounter;

      double w=*weight1 * *weight2;
      if(print || ((eventCounter %10000)==0))  { 
         cout<<eventCounter
             <<" event "<<*run<<" "<<*evno<<" type="<<*runtype<<" weight="<<w<<"\n";
         if(!print) print=1; //print this event
      }
      myEvent.run=*run;
      myEvent.evno=*evno;
      myEvent.w=w;

      if(*runtype==1) {
         // handle MC information
         if(print) {
            // kinematic variables from GKI bank
            cout<<" xGKI="<<*xGki<<" Q2GKI="<<*Q2Gki
                <<" y="<<*yGki
                <<" w="<<w<<"\n";
         }
         myEvent.xGKI = *xGki;
         myEvent.yGKI = *yGki;
         myEvent.Q2GKI = *Q2Gki;

         double _Q2GKI = myEvent.Q2GKI;
         double _yGKI = myEvent.yGKI;
         double _xGKI = myEvent.xGKI;

         H1GetPartMCId mcPartId(&*mcpart);
         mcPartId.Fill();

         TLorentzVector ebeam_MC_lab
            (mcpart[mcPartId.GetIdxBeamElectron()]->GetFourVector());
         TLorentzVector pbeam_MC_lab
            (mcpart[mcPartId.GetIdxBeamProton()]->GetFourVector());

         myEvent.eProtonBeamMC=pbeam_MC_lab.E();
         myEvent.eElectronBeamMC=ebeam_MC_lab.E();

         TLorentzVector escat0_MC_lab
            (mcpart[mcPartId.GetIdxScatElectron()]->GetFourVector());

         /*begin test*/
         TLorentzVector radPhot_MC_lab;
         if( mcPartId.GetIdxRadPhoton() >= 0 ){
            radPhot_MC_lab = mcpart[mcPartId.GetIdxRadPhoton()]->GetFourVector();
      
            double delta_phi = escat0_MC_lab.Phi() - radPhot_MC_lab.Phi();
         
            if( mcPartId.GetRadType() == 0 ){
               h_dPhi_theta_noR->Fill(escat0_MC_lab.Theta(), delta_phi );
            }
            else if( mcPartId.GetRadType() == 1 ){
               h_dPhi_theta_ISR->Fill(escat0_MC_lab.Theta(), delta_phi );
            }
            else if( mcPartId.GetRadType() == 2 ){
               h_dPhi_theta_FSR->Fill(escat0_MC_lab.Theta(), delta_phi );
            }
            else{
               cout << "something is wrong!" << endl;
            }
         }

         //HFS 4-vectors
         //TLorentzVector hfs_MC_lab = ebeam_MC_lab+pbeam_MC_lab-escat0_MC_lab;
         double hfs_MC_E_lab = 0.;
         double hfs_MC_pz_lab = 0.;
         for(int i=0;i<mcpart.GetEntries();i++) {
            H1PartMC *part=mcpart[i];
            int pdgid = part->GetPDG();
            int status=part->GetStatus();
            float charge=part->GetCharge();
            int elec_id = mcPartId.GetIdxScatElectron();
            if( status != 0 || i == elec_id ) continue;
            //if( TMath::RadToDeg()*part->GetTheta() > 177.5 || TMath::RadToDeg()*part->GetTheta() < 6. ) continue;

            hfs_MC_E_lab += part->GetE();
            hfs_MC_pz_lab += part->GetPz();
         }

         double sigma = hfs_MC_E_lab - hfs_MC_pz_lab;

         H1MakeKine makeKin_es;
         makeKin_es.MakeESig(escat0_MC_lab.E(), escat0_MC_lab.Theta(),sigma, ebeam_MC_lab.E(), pbeam_MC_lab.E());
         
         double Q2_esigma = makeKin_es.GetQ2es();
         double y_esigma = makeKin_es.GetYes();
         double x_esigma = makeKin_es.GetXes();

         myEvent.Q2MC_es = Q2_esigma;
         myEvent.yMC_es = y_esigma;
         myEvent.xMC_es = x_esigma;

         H1MakeKine makeKin_ISR;
         H1MakeKine makeKin_FSR;
         H1MakeKine makeKin_noR;

         double Q2_ISR=Q2_esigma;
         double y_ISR=y_esigma;
         double x_ISR=x_esigma;

         double Q2_FSR=Q2_esigma;
         double y_FSR=y_esigma;
         double x_FSR=x_esigma;

         double Q2_noR=Q2_esigma;
         double y_noR=y_esigma;
         double x_noR=x_esigma;

         myEvent.idxRad = mcPartId.GetRadType();

         if( mcPartId.GetIdxRadPhoton() >= 0 ){
            
            radPhot_MC_lab = mcpart[mcPartId.GetIdxRadPhoton()]->GetFourVector();
            
            if( mcPartId.GetRadType() == 1 ){
               makeKin_ISR.MakeESig(escat0_MC_lab.E(), escat0_MC_lab.Theta(), sigma, (ebeam_MC_lab).E(), pbeam_MC_lab.E());
               Q2_ISR=makeKin_ISR.GetQ2es();
               y_ISR=makeKin_ISR.GetYes();
               x_ISR=makeKin_ISR.GetXes();

               h_ISR_Q2diff->Fill( Q2_ISR - _Q2GKI );
               h_ISR_Ydiff->Fill( y_ISR - _yGKI );
               h_ISR_Xdiff->Fill( x_ISR - _xGKI );

            }
            else if( mcPartId.GetRadType() == 2 ){
               makeKin_FSR.MakeESig((escat0_MC_lab).E(), (escat0_MC_lab).Theta(), sigma, ebeam_MC_lab.E(), pbeam_MC_lab.E());
               Q2_FSR=makeKin_FSR.GetQ2es();
               y_FSR=makeKin_FSR.GetYes();
               x_FSR=makeKin_FSR.GetXes();

               h_FSR_Q2diff->Fill( Q2_FSR - _Q2GKI );
               h_FSR_Ydiff->Fill( y_FSR - _yGKI );
               h_FSR_Xdiff->Fill( x_FSR - _xGKI );
            }
         }
         else{
            makeKin_noR.MakeESig(escat0_MC_lab.E(), escat0_MC_lab.Theta(),sigma, ebeam_MC_lab.E(), pbeam_MC_lab.E());
            Q2_noR=makeKin_noR.GetQ2es();
            y_noR=makeKin_noR.GetYes();
            x_noR=makeKin_noR.GetXes();

            h_noR_Q2diff->Fill( Q2_noR - _Q2GKI );
            h_noR_Ydiff->Fill( y_noR - _yGKI );
            h_noR_Xdiff->Fill( x_noR - _xGKI );         
         }
         //end test

         // add radiative photon(s) in a cone
         TLorentzVector escatPhot_MC_lab(escat0_MC_lab);
         set<int> isElectron;
         isElectron.insert(mcPartId.GetIdxScatElectron());
         for(int i=0;i<mcpart.GetEntries();i++) {
            H1PartMC *part=mcpart[i];
            int status=part->GetStatus();
            if((status==0 )&&(part->GetPDG()==22)) {
               TLorentzVector p(part->GetFourVector());
               if(p.DeltaR(escat0_MC_lab)<ELEC_ISOLATION_CONE) {
                  // this photon counts with the electron
                  isElectron.insert(i);
                  escatPhot_MC_lab += p;
               }
            }
         }
         myEvent.elecEradMC=escatPhot_MC_lab.E()-escat0_MC_lab.E();
         myEvent.elecPxMC=escatPhot_MC_lab.X();
         myEvent.elecPyMC=escatPhot_MC_lab.Y();
         myEvent.elecPzMC=escatPhot_MC_lab.Z();
         myEvent.elecEMC=escatPhot_MC_lab.E();

         if(print) {
            cout<<"MC scattered electron is made of "<<isElectron.size()<<" particle(s)\n";
         }

         //H1MakeKine maybe helpful
         GetKinematics(ebeam_MC_lab,pbeam_MC_lab,escatPhot_MC_lab,
                       &myEvent.xMC,&myEvent.yMC,&myEvent.Q2MC);
         TLorentzRotation boost_MC_HCM = BoostToHCM(ebeam_MC_lab,pbeam_MC_lab,escatPhot_MC_lab);
         TLorentzVector q_MC_lab(ebeam_MC_lab-escatPhot_MC_lab);
         //New boost using the e-Sigma method
         TLorentzRotation boost_MC_HCM_es = BoostToHCM_es(ebeam_MC_lab,pbeam_MC_lab,escat0_MC_lab,Q2_esigma,y_esigma);

         //difference with respect to GKI values:
         h_Xdiff->Fill( x_esigma - _xGKI );
         h_Q2diff->Fill( Q2_esigma - _Q2GKI );
         h_Ydiff->Fill( y_esigma - _yGKI );

         // final state particles
         //bool haveElectron=false;
         myEvent.nMCtrackAll=0;
         myEvent.nMCtrack=0;
         for(int i=0;i<mcpart.GetEntries();i++) {
            
            H1PartMC *part=mcpart[i];
            if(print) {
               //cout << i << " " ; part->Print();
            }
            // skip particles counted as electron
            if(isElectron.find(i)!=isElectron.end()) continue;

            int status=part->GetStatus();
            if(status==0) {
               // generator "stable" particles
               // if((!haveElectron)&&
               //    ((part->GetPDG()==11)||(part->GetPDG()== -11))) {
               //    haveElectron=true;               
               // } else 
               if(part->GetCharge()!=0.) {
                  // other charged particles
                  TLorentzVector h=part->GetFourVector();
                  double log10z=TMath::Log10((h*pbeam_MC_lab)/(q_MC_lab*pbeam_MC_lab));
                  // boost to hadronic-centre-of-mass frame
                  TLorentzVector hStar = boost_MC_HCM*h;
                  TLorentzVector hStar2 = boost_MC_HCM_es*h;
                  double etaStar=hStar.Eta();
                  double ptStar=hStar.Pt();
                  double phiStar=hStar.Phi();

                  double etaStar2=hStar2.Eta();
                  double ptStar2=hStar2.Pt();
                  double phiStar2=hStar2.Phi();

                  if(print && etaStar2 < -20) {
                     //cout << i << " " ; part->Print();
                     cout<<"MCpart "<<myEvent.nMCtrackAll
                         <<" "<<part->GetPDG()
                         <<" etaLab="<<h.Eta()
                         <<" ptLab="<<h.Pt()
                         <<" phiLab="<<h.Phi()
                         <<" ptStar="<<ptStar
                         <<" etaStar="<<etaStar
                         <<" phiStar="<<phiStar
                         <<" ptStar2="<<ptStar2
                         <<" etaStar2="<<etaStar2
                         <<" phiStar2="<<phiStar2
                         <<" Boost px "<<hStar2.Px()
                         <<" Boost py "<<hStar2.Py()
                         <<" Boost pz "<<hStar2.Pz()
                         <<" log10(z)="<<log10z<<"\n";
                  }
                  myEvent.nMCtrackAll++;
                  if(myEvent.nMCtrack<MyEvent::nMCtrack_MAX) {
                     int k=myEvent.nMCtrack;
                     myEvent.idMC[k]=part->GetPDG();
                     myEvent.pxMC[k]=h.X();
                     myEvent.pyMC[k]=h.Y();
                     myEvent.pzMC[k]=h.Z();
                     myEvent.etaMC[k]=h.Eta();
                     myEvent.chargeMC[k]=part->GetCharge();

                     myEvent.ptStarMC[k]=hStar.Pt();
                     myEvent.etaStarMC[k]=hStar.Eta();
                     myEvent.phiStarMC[k]=hStar.Phi();

                     myEvent.ptStar2MC[k]=hStar2.Pt();
                     myEvent.etaStar2MC[k]=hStar2.Eta();
                     myEvent.phiStar2MC[k]=hStar2.Phi();

                     myEvent.log10zMC[k]=log10z;
                     myEvent.imatchMC[k]=-1;
                     myEvent.partMC[k]=part;
                     myEvent.nMCtrack=k+1;
                  }
               }
            } // end loop over stable particles
         }
      }//end of MC particles


      output->Fill();
   }

    // Summary
    cout << "\nProcessed " << eventCounter
        << " events\n\n";
    cerr << "\nProcessed " << eventCounter
       << " events\n\n";

    // Write histogram to file
    output->Write();

    h_dPhi_theta_noR->Write();
    h_dPhi_theta_ISR->Write();
    h_dPhi_theta_FSR->Write();

    h_Q2diff->Write();
    h_Xdiff->Write();
    h_Ydiff->Write();

    h_ISR_Q2diff->Write();
    h_ISR_Xdiff->Write();
    h_ISR_Ydiff->Write();

    h_FSR_Q2diff->Write();
    h_FSR_Xdiff->Write();
    h_FSR_Ydiff->Write();

    h_noR_Q2diff->Write();
    h_noR_Xdiff->Write();
    h_noR_Ydiff->Write();

    //file.Close();
    delete file;

    cout << "Histograms written to " << opts.GetOutput() << endl;
    cerr << "Histograms written to " << opts.GetOutput() << endl;

    return 0;
}
