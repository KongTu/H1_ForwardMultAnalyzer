#include "VarList.h"
#include "ClassifierBinning.h"

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

using namespace std;


void extractMultiplicity(int i0,int i1,Int_t const *binMap,
                         TH1 const *src,TH1 *dest,
                         TH2 const *src_e=0,
                         TH2 *dest_e=0,TH2 *dest_rhoij=0);

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

  // read steering in XML format                                                
  TDOMParser binningParser;
  Int_t XMLerror=binningParser.ParseFile(xmlBinning);
  if(XMLerror) {
    cout<<"error="<<XMLerror<<" from TDOMParser (while parsing \""
        <<xmlBinning<<"\")\n";
    return 2;
  }
  TXMLDocument const *document=binningParser.GetXMLDocument();

  TString outputFileName;
  TXMLNode *root=document->GetRootNode();

  if(root && (!TString(root->GetNodeName()).CompareTo("TUnfoldBinning")) &&
     (root->GetNodeType()==TXMLNode::kXMLElementNode)) {
    for(TXMLNode *node=root->GetChildren();node;node=node->GetNextNode()) {
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
  TString resultFileName="results_"+outputFileName;

  TFile *outputFile=new TFile(resultFileName,"recreate");
     
  //==================================================
  // read binning schemes and input histograms
  TFile *inputFile=new TFile(outputFileName);

  if(!inputFile->IsOpen()) {
    cout<<"failed to open: "<<outputFileName<<"\n";
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

  // read directories in output file
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
  //genBinning->PrintStream(cout);
  VarList dummyVars;
  ClassifierBinning covClassifier(covBinning,dummyVars);
  ClassifierBinning genClassifier(genBinning,dummyVars);
  ClassifierBinning recClassifier(recBinning,dummyVars);

  // covClassifier:
  //    describes bins in (etaRec,ntrackRec)
  //
  //    the input file contains mant covariance matrices in this binning
  //    each of these covariance matrixes describes the correlations
  //    between the (etaRec,ntrackRec) bins for a fixed (Q2rec,yrec) bin
  //
  //
  // recClassifier
  //    describes binning in (Q2rec,yrec,ntrackRec)
  //
  //    this is the "rec" axis of the unfolding matrix
  //
  //
  // genClassifier:
  //    describes binning in (Q2gen,ygen,ntrackGen)
  //
  //    this is the "gen" axis of the unfolding matrix
  //

  vector<ClassifierBinning const *> covClasses=covClassifier.Enumerate();
  vector<ClassifierBinning const *> recClasses=recClassifier.Enumerate();
  vector<ClassifierBinning const *> genClasses=genClassifier.Enumerate();

  //   covClasses:
  //      these are bins in eta
  //   recClasses
  //      these are bins in (Q2rec,yrec)
  //   genClasses
  //      these are bins in (Q2gen,ygen)
  
  for(size_t src=0;src<source.size();src++) {
     TString name=source[src]->GetName();
     cout<<"plotting "<<name<<"\n";

     // create new output directory
     outputFile->cd();
     TDirectoryFile *outputDir=new TDirectoryFile(name,name);
     outputFile->Add(outputDir);
     outputDir->cd();


     // read unfolded data and pack it to histograms
     //  -> one histogram per (eta,Q2,y) bin
     //  -> one covariance histogram per (eta, (Q2,y)_i, (Q2,y)_j ) bin
     for(size_t iEta=0;iEta<covClasses.size();iEta++) {
        TString etaName=covClasses[iEta]->GetName();
        TH1 *hist_unfold_iEta;
        source[src]->GetObject("hist_unfolded_"+etaName,hist_unfold_iEta);
        TH2 *hist_ematrixTotal_iEta;
        source[src]->GetObject("hist_ematrixTotal_"+etaName,
                               hist_ematrixTotal_iEta);
        TH1 *hist_bias_iEta;
        source[src]->GetObject("hist_bias_"+etaName,hist_bias_iEta);
        for(size_t iQ2y=0;iQ2y<genClasses.size();iQ2y++) {
           TString q2yName=genClasses[iQ2y]->GetName();
           TUnfoldBinning const *binning=
              genClasses[iQ2y]->GetBinningNodeDistribution();
           int i0=binning->GetStartBin();
           int i1=binning->GetEndBin();
           Int_t *binMap=0;
           TH1 *hist_unfold_iEtaQ2y=binning->CreateHistogram
              ("hist_unfolded_"+etaName+q2yName,kTRUE,&binMap);
           TH2 *hist_ematrix_iEtaQ2y=binning->CreateErrorMatrixHistogram
              ("hist_ematrix_"+etaName+q2yName,kTRUE);
           TH2 *hist_rhoij_iEtaQ2y=binning->CreateErrorMatrixHistogram
              ("hist_rhoij_"+etaName+q2yName,kTRUE);
           extractMultiplicity(i0,i1,binMap,
                               hist_unfold_iEta,hist_unfold_iEtaQ2y,
                               hist_ematrixTotal_iEta,
                               hist_ematrix_iEtaQ2y,
                               hist_rhoij_iEtaQ2y);
           TH1 *hist_bias_iEtaQ2y=binning->CreateHistogram
              ("hist_bias_"+etaName+"_"+q2yName,kTRUE);
           extractMultiplicity(i0,i1,binMap,hist_bias_iEta,hist_bias_iEtaQ2y);
           hist_unfold_iEtaQ2y->Write();
           hist_ematrix_iEtaQ2y->Write();
           hist_rhoij_iEtaQ2y->Write();
           hist_bias_iEtaQ2y->Write();
           delete [] binMap;
        }
     }
  }
  delete inputFile;
  delete outputFile;

  return 0;
}

void extractMultiplicity(int i0,int i1,Int_t const *binMap,
                         TH1 const *src,TH1 *dest,
                         TH2 const *src_e,TH2 *dest_e,TH2 *dest_rhoij) {
   if(!(src && dest)) return;
   //cout<<"extractMultiplicity "<<src->GetName()<<" "<<i0<<","<<i1<<" -> "<<dest->GetName()<<"\n";
   double n=0.;
   for(int i=i0;i<i1;i++) {
      if(binMap[i]<0) continue;
      n += src->GetBinContent(i);
   }
   //cout<<"n="<<n<<"\n";
   vector<double> Vin(i1-i0+1);
   double Vnn=0.;
   if(src_e) {
      for(int i=i0;i<i1;i++) {
         if(binMap[i]<0) continue;
         for(int j=i0;j<i1;j++) {
            if(binMap[j]<0) continue;
            Vin[i-i0] += src_e->GetBinContent(i,j);
         }
         Vnn += Vin[i-i0];
      }
   }
   //cout<<"Vnn="<<Vnn<<"\n";
   for(int i=i0;i<i1;i++) {
      if(binMap[i]<0) continue;
      double wi=dest->GetBinWidth(binMap[i]);
      double xi=src->GetBinContent(i)/n;
      dest->SetBinContent(binMap[i],xi/wi);
      if(src_e) {
         double e=0.;
         for(int j=i0;j<i1;j++) {
            if(binMap[j]<0) continue;
            double wj=dest->GetBinWidth(binMap[j]);
            double xj=src->GetBinContent(j)/n;
            double e_ij = 
               (src_e->GetBinContent(i,j)-xi*Vin[j-i0]-xj*Vin[i-i0]+xi*xj*Vnn)/
               (n*n*wi*wj);
            if(dest_e) {
               dest_e->SetBinContent(binMap[i],binMap[j],e_ij); 
            }
            if(i==j) e=TMath::Sqrt(e_ij);
         }
         dest->SetBinError(binMap[i],e);
      }
   }
   if(dest_e) {
      for(int i=0;i<=dest_e->GetNbinsX()+1;i++) {
         double Vii=dest_e->GetBinContent(i,i);
         if(Vii<=0.0) continue;
         for(int j=0;j<=dest_e->GetNbinsY()+1;j++) {
            double Vjj=dest_e->GetBinContent(j,j);
            if(Vjj<=0.) continue;
            dest_rhoij->SetBinContent(i,j,dest_e->GetBinContent(i,j)/
                                      TMath::Sqrt(Vii*Vjj));
         }
      }
   }
}
