#include "VarList.h"
#include "ClassifierBinning.h"

#include <TError.h>
#include <TFile.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TKey.h>
#include <TList.h>

#ifdef PRIVATE_TUNFOLD
#include "../TUnfold_V17.7/TUnfoldDensity.h"
#else
#include <TUnfoldDensity.h>
#endif

#define NSTEP 40
//#define NSTEP 10
//#define STOP_EARLY

using namespace std;

int main(int argc, char const *argv[]) {
  // switch on histogram errors
  TH1::SetDefaultSumw2();
  gErrorIgnoreLevel=kInfo;
  gErrorAbortLevel=kError;

  TString xmlBinning;
  if(argc>=2) {
    xmlBinning=argv[1];
  } else {
    cout<<"usage: "<<argv[0]<<" binning.xml\n";
    return 1;
  }
  int toyMin=0;
  int toyMax=0;
  if(argc>=4) {
     toyMin=TString(argv[2]).Atoi();
     toyMax=TString(argv[3]).Atoi();
  }

  // read steering in XML format                                                
  TDOMParser binningParser;
  Int_t XMLerror=binningParser.ParseFile(xmlBinning);
  if(XMLerror) {
    cout<<"error="<<XMLerror<<" from TDOMParser (while parsing \""
        <<xmlBinning<<"\")\n";
    return 2;
  }
  TXMLDocument const *document=binningParser.GetXMLDocument();

  TString inputFileName;
  TString outputFileName;
  TXMLNode *root=document->GetRootNode();

  if(root && (!TString(root->GetNodeName()).CompareTo("TUnfoldBinning")) &&
     (root->GetNodeType()==TXMLNode::kXMLElementNode)) {
    for(TXMLNode *node=root->GetChildren();node;node=node->GetNextNode()) {
      if(node->GetNodeType()==TXMLNode::kXMLElementNode &&
         !TString(node->GetNodeName()).CompareTo("UnfoldingInput")) {
	inputFileName=node->GetChildren()->GetContent();
	int l;
        do {
          l=inputFileName.Length();
          inputFileName=inputFileName.Strip(TString::kLeading,'\n');
          inputFileName=inputFileName.Strip(TString::kLeading,' ');
          inputFileName=inputFileName.Strip(TString::kTrailing,'\n');
          inputFileName=inputFileName.Strip(TString::kTrailing,' ');
        } while(l!=inputFileName.Length());
	cout<<"main: inputFileName: "<<inputFileName<<"\n";
      }
      if(node->GetNodeType()==TXMLNode::kXMLElementNode &&
         !TString(node->GetNodeName()).CompareTo("UnfoldingOutput")) {
	outputFileName=node->GetChildren()->GetContent();
        int l;
        do {
          l=outputFileName.Length();
          outputFileName=outputFileName.Strip(TString::kLeading,'\n');
          outputFileName=outputFileName.Strip(TString::kLeading,' ');
          outputFileName=outputFileName.Strip(TString::kTrailing,'\n');
          outputFileName=outputFileName.Strip(TString::kTrailing,' ');
        } while(l!=outputFileName.Length());
	cout<<"main: outputFileName: "<<outputFileName<<"\n";
      }
    }
  }

  // open output file
  if(toyMin<toyMax) {
     outputFileName+=TString::Format("_%d-%d",toyMin,toyMax-1);
  }
  TFile *outputFile=new TFile(outputFileName,"recreate");
     
  //==================================================
  // read binning schemes and input histograms
  TFile *inputFile=new TFile(inputFileName);

  if(!inputFile->IsOpen()) {
    cout<<"failed to open: "<<inputFileName<<"\n";
    exit(2);
  }

  outputFile->cd();

  TUnfoldBinning *recBinning,*genBinning,*covBinning;

  inputFile->GetObject("recQ2",recBinning);
  inputFile->GetObject("genQ2",genBinning);
  inputFile->GetObject("covBinning",covBinning);

  TObjString *binningXML_object;
  inputFile->GetObject("binningXML",binningXML_object);
  if(!binningXML_object) {
     cout<<"Error: can not read binning from unfolding input\n";
  }

  if(!(recBinning && genBinning && covBinning)) {

     TString binningXML=binningXML_object->GetString();
     binningParser.ParseBuffer(binningXML,binningXML.Length());

     document=binningParser.GetXMLDocument();

     if(!genBinning) {
        genBinning=TUnfoldBinningXML::ImportXML(document,"genQ2");
     }
     if(!recBinning) {
        recBinning=TUnfoldBinningXML::ImportXML(document,"recQ2");
     }
     if(!covBinning) {
        covBinning=TUnfoldBinningXML::ImportXML(document,"covBinning");
     }

     if((!recBinning)||(!genBinning)||(!covBinning)) {
        cout<<"problem to read binning schemes\n";
     }
  }

  // save binning schemes to output file
  binningXML_object->Write("binningXML");
  recBinning->Write();
  genBinning->Write();
  covBinning->Write();

  VarList dummyVars;
  ClassifierBinning covClassifier(covBinning,dummyVars);
  ClassifierBinning recClassifier(recBinning,dummyVars);
  ClassifierBinning genClassifier(genBinning,dummyVars);

  vector<ClassifierBinning const *> etaClasses=covClassifier.Enumerate();
  vector<ClassifierBinning const *> recClasses=recClassifier.Enumerate();

  // read directories in output file
  TList *keys=inputFile->GetListOfKeys();
  TIter next(keys);
  TKey *key;
  vector<TDirectory *> source;
  while ((key = (TKey *)next())) {
     TDirectory *dir=key->ReadObject<TDirectory>();
     if(dir) {
        //cout<<"directory: "<<key->GetName()<<"\n";
        source.push_back(dir);
     }
  }

  // unfolding
  for(size_t src1=0;src1<source.size();src1++) {
     vector<TUnfoldDensity *> unfolding;
     cout<<"read histograms from "<<source[src1]->GetName()<<" neta="<<etaClasses.size()<<"\n";
     for(size_t iEta=0;iEta<etaClasses.size();iEta++) {
        TString etaName=etaClasses[iEta]->GetName();
        TH2 * hist_genRec1=0;
        //TH1 * hist_genMC1=0;
        TH1 * hist_fake1=0;
        TH1 * hist_QEDc1=0;
        source[src1]->GetObject("hist_genRec_"+etaName,hist_genRec1);
        source[src1]->GetObject("hist_fake_"+etaName,hist_fake1);
        //source[src1]->GetObject("hist_gen_"+etaName,hist_genMC1);
        source[src1]->GetObject("hist_QEDc_"+etaName,hist_QEDc1);

        if(hist_genRec1) {
           TUnfoldDensity *unfold1=
              new TUnfoldDensity
              (hist_genRec1,TUnfoldDensity::kHistMapOutputHoriz,
               TUnfoldDensity::kRegModeSize,
               TUnfoldDensity::kEConstraintNone,
               TUnfoldDensity::kDensityModeNone,
               genBinning,recBinning);
           unfolding.push_back(unfold1);
           if(hist_fake1) unfold1->SubtractBackground(hist_fake1,"fakes",1.0,0.2);
            // if(hist_QEDc1) unfold1->SubtractBackground(hist_QEDc1, "QEDc",0.33,0.2);
        }
     }

     if(!unfolding.size()) continue;
     
     for(size_t src2=0;src2<source.size();src2++) {
        if(src1==src2) {
           continue;
           // could generate toys here...
        }
        cout<<"unfolding "<<source[src2]->GetName()
            <<" using MC="<<source[src1]->GetName()<<"\n";
        TString name=TString(source[src2]->GetName())+
           "_from_"+source[src1]->GetName();
        outputFile->cd();
        TDirectoryFile *outputDir=new TDirectoryFile(name,name);
        outputFile->Add(outputDir);
        outputDir->cd();

        for(size_t recBin=0;recBin<recClasses.size();recBin++) {
           TH2 *hist_covar=0;
           source[src2]->GetObject
              ("hist_recCovar_"+recClasses[recBin]->GetName(),hist_covar);
           if(hist_covar) hist_covar->Write();
        }

        for(size_t iEta=0;iEta<etaClasses.size();iEta++) {
           TString etaName=etaClasses[iEta]->GetName();
           TH1 *hist_truth=0,*hist_rec=0;
           source[src2]->GetObject
              ("hist_gen_"+etaName,hist_truth);
           source[src2]->GetObject
              ("hist_rec_"+etaName,hist_rec);
           if(hist_rec) {
              unfolding[iEta]->SetInput(hist_rec,1.,0.);
              if(hist_truth) hist_truth->Write();
              hist_rec->Write();

              TGraph *lCurve;
              TSpline *logTauX,*logTauY,*logTauCurvature;
              double tauMin=1.E-9;
              double tauMax=1.E-3;
              int iBest=unfolding[iEta]->ScanLcurve
                 (NSTEP,tauMin,tauMax,&lCurve,&logTauX,&logTauY,
                  &logTauCurvature);
              Double_t t[1],x[1],y[1];
              logTauX->GetKnot(iBest,t[0],x[0]);
              logTauY->GetKnot(iBest,t[0],y[0]);
              TGraph *bestLcurve=new TGraph(1,x,y);
              TGraph *bestLogTauX=new TGraph(1,t,x);

              lCurve->Write("Lcurve_"+etaName);
              bestLcurve->Write("bestLcurve_"+etaName);
              logTauX->Write("logTauX_"+etaName);
              bestLogTauX->Write("bestLogTauX_"+etaName);
              logTauY->Write("LogTauY_"+etaName);
              logTauCurvature->Write("LogTauCurvature_"+etaName);


              TH2 *hist_genRec=unfolding[iEta]->GetProbabilityMatrix
                 ("hist_genRec_"+etaName);
              hist_genRec->Write();

              TH1 *hist_fakes=unfolding[iEta]->GetBackground
                 ("hist_fakes_"+etaName,"fakes");
              if(hist_fakes) hist_fakes->Write();

              TH1 *hist_genMC=unfolding[iEta]->GetBias
                 ("hist_bias_"+etaName);

              TH1 *unfolded=unfolding[iEta]->GetOutput
                 ("hist_unfolded_"+etaName);
              TH2 *ematrixTotal=unfolding[iEta]->GetEmatrixTotal
                 ("hist_ematrixTotal_"+etaName);
              TH2 *rhoIJ=unfolding[iEta]->GetRhoIJtotal
                 ("hist_rhoij_"+etaName);
              unfolded->Write();
              ematrixTotal->Write();
              rhoIJ->Write();
              hist_genMC->Write();

              // map to bin numbers
              TH2 *hist_dxdy=unfolding[iEta]->GetTransferMatrix
                 ("hist_dxdy_"+etaName);
              hist_dxdy->Write();

           }
        }
#ifdef STOP_EARLY
        break;
#endif
     }
#ifdef STOP_EARLY
     break;
#endif
  }

  delete outputFile;
  delete inputFile;

  return 0;
}
