#include <set>

#include <TError.h>
#include <TFile.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TKey.h>
#include <TList.h>

#ifdef PRIVATE_TUNFOLD
#include "../TUnfold_V17.7/TUnfoldBinningXML.h"
#else
#include <TUnfoldBinningXML.h>
#endif

#include "UnfoldingResult.h"
#include "UnfoldingTUnfoldDensity.h"
#include "UnfoldingIterative.h"

#define NSTEP 40
//#define STOP_EARLY

using namespace std;

static int const NSaveDetails=1;

void extractMultiplicity(int eventBin,int i0,int i1,Int_t const *binMap,
                         TH1 const *src,TH1 *dest,
                         TH2 const *src_e=0,TH2 *dest_e=0);


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

  TUnfoldBinning *recBinning,*genBinning;

  inputFile->GetObject("recBinning",recBinning);
  inputFile->GetObject("genBinning",genBinning);

  TObjString *binningXML_object;
  inputFile->GetObject("binningXML",binningXML_object);
  if(!binningXML_object) {
     cout<<"Error: can not read binning from unfolding input\n";
  }

  if(!(recBinning && genBinning)) {

     TString binningXML=binningXML_object->GetString();
     binningParser.ParseBuffer(binningXML,binningXML.Length());

     document=binningParser.GetXMLDocument();

     if(!genBinning) {
        genBinning=TUnfoldBinningXML::ImportXML(document,"genBinning");
     }
     if(!recBinning) {
        recBinning=TUnfoldBinningXML::ImportXML(document,"recBinning");
     }

     if((!recBinning)||(!genBinning)) {
        cout<<"problem to read binning schemes\n";
     }
  }

  // save binning schemes to output file
  binningXML_object->Write("binningXML");
  recBinning->Write();
  genBinning->Write();

  // read directories in outout file
  TList *keys=inputFile->GetListOfKeys();
  TIter next(keys);
  TKey *key;
  vector<TDirectory *> source;
  while ((key = (TKey *)next())) {
     TDirectory *dir=key->ReadObject<TDirectory>();
     if(dir) {
        source.push_back(dir);
     }
  }

  // unfolding
  for(size_t src1=0;src1<source.size();src1++) {
     TH2 *hist_genRec=0;
     TH1 *hist_genMC=0;
     TH1 *hist_fake=0;
     source[src1]->GetObject("hist_genRec",hist_genRec);
     if(!hist_genRec) continue;
     source[src1]->GetObject("hist_fake",hist_fake);
     source[src1]->GetObject("hist_gen",hist_genMC);

     TUnfoldDensity tunfoldDensity
        (hist_genRec,TUnfoldDensity::kHistMapOutputHoriz,
         TUnfoldDensity::kRegModeSize,
         TUnfoldDensity::kEConstraintNone,
         TUnfoldDensity::kDensityModeNone,
         genBinning,recBinning);

     if(hist_fake) tunfoldDensity.SubtractBackground(hist_fake,"fakes",1.0,0.2);

     for(size_t src2=0;src2<source.size();src2++) {
        if(src1==src2) {
           continue;
           // could generate toys here...
        }
        TH1 *hist_gen;
        TH1 *hist_rec;
        TH2D *hist_recCovar;
        source[src2]->GetObject("hist_gen",hist_gen);
        source[src2]->GetObject("hist_rec",hist_rec);
        source[src2]->GetObject("hist_recCovar",hist_recCovar);
        if(!(hist_rec&&hist_recCovar)) continue;
        cout<<"unfolding "<<source[src2]->GetName()
            <<" using MC="<<source[src1]->GetName()<<"\n";

        // loop over unfolding algorithms
        TString name=TString(source[src2]->GetName())
                      +"_from_"+source[src1]->GetName();
        outputFile->cd();
        TDirectoryFile *outputDir=new TDirectoryFile(name,name);

        outputFile->Add(outputDir);
        outputDir->cd();

        if(hist_gen) hist_gen->Write();
        hist_genRec->Write();
        if(hist_fake) hist_fake->Write();
        hist_rec->Write();
        hist_recCovar->Write();

        tunfoldDensity.SetInput(hist_rec,0.,1.,hist_recCovar);

        TGraph *lCurve;
        TSpline *logTauX,*logTauY,*logTauCurvature;
        double tauMin=1.E-10;
        double tauMax=1.E-4;
        int iBest=tunfoldDensity.ScanLcurve
           (NSTEP,tauMin,tauMax,&lCurve,&logTauX,&logTauY,
            &logTauCurvature);
        Double_t t[1],x[1],y[1];
        logTauX->GetKnot(iBest,t[0],x[0]);
        logTauY->GetKnot(iBest,t[0],y[0]);
        TGraph *bestLcurve=new TGraph(1,x,y);
        TGraph *bestLogTauX=new TGraph(1,t,x);

        lCurve->Write("Lcurve");
        bestLcurve->Write("bestLcurve");
        logTauX->Write("logTauX");
        bestLogTauX->Write("bestLogTauX");
        logTauY->Write("LogTauY");
        logTauCurvature->Write("LogTauCurvature");

        TH1 *unfolded=tunfoldDensity.GetOutput("unfolded");
        TH2 *ematrix=tunfoldDensity.GetEmatrixTotal("ematrixInput");
        unfolded->Write();
        ematrix->Write();

        // extract all distributions, normalize and store
        for(TUnfoldBinning const *node=genBinning;node;) {
           bool skipChild=false;
           if(node->HasUnconnectedBins()) {
              //extract node
              TString path;
              for(TUnfoldBinning const *printNode=node;
                  printNode->GetParentNode();
                  printNode=printNode->GetParentNode()) {
                 path=printNode->GetName()+TString("!")+path;
              }
              if(node->GetChildNode()) {
                 path+=node->GetChildNode()->GetName();
              }
              cout<<"extracting: "<<path<<"\n";
              skipChild=true;
              if(node->GetChildNode()) {
                 int eventBin=node->GetStartBin();
                 // pointer to (2D) binning of multiplicity
                 // mult vs eta
                 TUnfoldBinning const *multHistNode=node->GetChildNode();

                 TVectorD const *etaBins=
                    multHistNode->GetDistributionBinning(1);
                 // loop over eta bins
                 for(int ieta=0;ieta<etaBins->GetNrows()-1;ieta++) {
                    // create 1D histograms for each eta bin
                    TString axisSteering=
                       multHistNode->GetDistributionAxisLabel(1);
                    axisSteering+=TString::Format("[CUO%d]",ieta);
                    TString base_name=path+TString::Format
                       ("_%.1f:eta:%.1f",
                        (*etaBins)(ieta),(*etaBins)(ieta+1));

                    TString title1D=";N_{chg};#LT_{}P(N_{chg})#GT";
                    TString title2D=";N_{chg};N_{chg}";

                    // gen histogram (used for unfolding matrix)
                    TH1 *hist_ntrackMCtruth=multHistNode->CreateHistogram
                       ("MCtruth_"+base_name,true,0,title1D,axisSteering);
                    // truth histograms (if available)
                    TH1 *hist_ntrackTruth=0;
                    if(hist_gen) {
                       hist_ntrackTruth=multHistNode->CreateHistogram
                          ("truth_"+base_name,true,0,title1D,axisSteering);
                    }
                    // unfolded histograms
                    TH2 *hist_ntrackUnfEmatrix=
                       multHistNode->CreateErrorMatrixHistogram
                       ("unfEmatrix_"+base_name,true,0,title2D,axisSteering);
                    TH2 *hist_ntrackUnfRhoIJ=
                       multHistNode->CreateErrorMatrixHistogram
                       ("unfRhoIJ_"+base_name,true,0,title2D,axisSteering);
                    Int_t *binMap;
                    TH1 *hist_ntrackUnf=multHistNode->CreateHistogram
                       ("unfResult_"+base_name,
                        true,&binMap,title1D,axisSteering);
                    // extract multiplicities
                    int i0=multHistNode->GetStartBin();
                    int i1=multHistNode->GetEndBin();
                    extractMultiplicity(eventBin,i0,i1,binMap,
                                        unfolded,hist_ntrackUnf,
                                        ematrix,hist_ntrackUnfEmatrix);
                    extractMultiplicity(eventBin,i0,i1,binMap,
                                        hist_genMC,hist_ntrackMCtruth);
                    extractMultiplicity(eventBin,i0,i1,binMap,
                                        hist_gen,hist_ntrackTruth);
                    delete binMap;

                    for(int i=0;i<=hist_ntrackUnfEmatrix->GetNbinsX()+1;i++) {
                       double v_ii=hist_ntrackUnfEmatrix->GetBinContent(i,i);
                       if(!(v_ii>0.)) continue;
                       for(int j=0;j<=hist_ntrackUnfEmatrix->GetNbinsY()+1;
                           j++) {
                          double v_jj=hist_ntrackUnfEmatrix->GetBinContent(j,j);
                          if(!(v_jj>0.)) continue;
                          hist_ntrackUnfRhoIJ->SetBinContent
                             (i,j,hist_ntrackUnfEmatrix->GetBinContent(i,j)/
                              TMath::Sqrt(v_ii*v_jj));
                       }
                    }

                    hist_ntrackUnf->Write();
                    hist_ntrackUnfEmatrix->Write();
                    hist_ntrackUnfRhoIJ->Write();
                    hist_ntrackMCtruth->Write();
                    if(hist_ntrackTruth) {
                       hist_ntrackTruth->Write();
                    }
                 }
              }
           }
           if(node->GetChildNode() && !skipChild) {
              node=node->GetChildNode();
           } else if(node->GetNextNode()) {
              node=node->GetNextNode();
           } else {
              node=node->GetParentNode();
              if(node) node=node->GetNextNode();
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

void extractMultiplicity(int eventBin,int i0,int i1,Int_t const *binMap,
                         TH1 const *src,TH1 *dest,
                         TH2 const *src_e,TH2 *dest_e) {
   if(!(src && dest)) return;
   double nevent=src->GetBinContent(eventBin);
   double nevent_sq=nevent*nevent;
   for(int i=i0;i<i1;i++) {
      if(binMap[i]<0) continue;
      double e=0.0;
      if(src_e) {
         for(int j=i0;j<i1;j++) {
            if(binMap[j]<0) continue;
            double vij=(src_e->GetBinContent(i,j)
                        +
                        (src_e->GetBinContent(i,eventBin)*src->GetBinContent(j)+
                         src_e->GetBinContent(j,eventBin)*src->GetBinContent(i)
                         )/nevent
                        +
                        src_e->GetBinContent(eventBin,eventBin)*
                        src->GetBinContent(i)*src->GetBinContent(j)/nevent_sq
                        )/nevent_sq
               /dest->GetBinWidth(binMap[i])/dest->GetBinWidth(binMap[j]);
            if(dest_e) {
               dest_e->SetBinContent(binMap[i],binMap[j],vij);
            }
            if(i==j) e=TMath::Sqrt(vij);
         }
      }
      double c=src->GetBinContent(i)/nevent/dest->GetBinWidth(binMap[i]);
      dest->SetBinContent(binMap[i],c);
      if(src_e) dest->SetBinError(binMap[i],e);
   }
}
