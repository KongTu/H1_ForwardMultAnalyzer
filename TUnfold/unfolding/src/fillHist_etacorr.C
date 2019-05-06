
#include <TMath.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include <TXMLAttr.h>
#include <cmath>
#include <set>
#include <sstream>
#include <map>
#include <vector>

#include "VarList.h"
#include "ClassifierBinning.h"

// #define USE_FST
#define LABFRAME

using namespace std;

/*
   program for reading events from a minitree and filling them into
   histograms.

   The histograms are meant to be read by some unfolding code.


   Run this program with one or two arguments

    first argument: XML file which describes the binning and other auxillary
       information (e.g. name of the file which is written)
       if no first argument is given, the program stops

    second argument: XML file which describes the minitrees
       if no second argument is given, the file
          "minitree.xml"
       is tried

    output to a root file (=input to unfolding)
    The output root file is set by the tag "UnfoldingInput" in the XML file
    (name of XML = first argument to this program)

    Content of the root file which is written:
    a new TDirectory for each set of input files from minitree.xml
    with the same name

    in each TDirectory there are histograms

 */

struct MinitreeDescriptor {
   string fFileName;
   string fTreeName;
   double fLuminosity;
};


int main(int argc, char const *argv[]) {
  // switch on histogram errors
  TH1::SetDefaultSumw2();

  TString xmlBinning;
  TString xmlMinitree="minitree.xml";

  if(argc>=2) {
    xmlBinning=argv[1];
    if(argc>=3) {
      xmlMinitree=argv[2];
    }
  } else {
    cout<<"usage: "<<argv[0]<<" binning.xml [minitree.xml]\n";
    return 1;
  }

  cout<<"main: binning is described in: \""<<xmlBinning<<"\"\n";
  cout<<"main: minitrees are described in: \""<<xmlMinitree<<"\"\n";

  // read steering in XML format
  TDOMParser binningParser;
  Int_t error=binningParser.ParseFile(xmlBinning);
  if(error) {
    cout<<"error="<<error<<" from TDOMParser (while parsing \""
	<<xmlBinning<<"\")\n";
    return 2;
  }

  // read minitree steering
  TDOMParser minitreeParser;
  error=minitreeParser.ParseFile(xmlMinitree);
  if(error) {
    cout<<"error="<<error<<" from TDOMParser (while parsing \""
	<<xmlMinitree<<"\")\n";
    return 3;
  }

  // parse steering to get output file name and data luminosity
  TString outputFileName;
  double dataLuminosity=-1.;
  TXMLDocument const *document=binningParser.GetXMLDocument();
  TXMLNode *root=document->GetRootNode();
  if(root && (!TString(root->GetNodeName()).CompareTo("TUnfoldBinning")) &&
     (root->GetNodeType()==TXMLNode::kXMLElementNode)) {
    // locate <UnfoldingInput>
    for(TXMLNode *node=root->GetChildren();node;node=node->GetNextNode()) {
      if(node->GetNodeType()==TXMLNode::kXMLElementNode &&
	 !TString(node->GetNodeName()).CompareTo("UnfoldingInput") &&
	 node->GetAttributes()) {
	TIterator *i=node->GetAttributes()->MakeIterator();
	TXMLAttr *attr;
	while((attr=(TXMLAttr *)i->Next())) {
	  if(!TString(attr->GetName()).CompareTo("dataluminosity")) {
	    dataLuminosity=TString(attr->GetValue()).Atof();
	  }
	}
	outputFileName=node->GetChildren()->GetContent();
	int l;
	do {
	  l=outputFileName.Length();
	  outputFileName=outputFileName.Strip(TString::kLeading,'\n');
	  outputFileName=outputFileName.Strip(TString::kLeading,' ');
	  outputFileName=outputFileName.Strip(TString::kTrailing,'\n');
	  outputFileName=outputFileName.Strip(TString::kTrailing,' ');
	} while(l!=outputFileName.Length());
	cout<<"output file: "<<outputFileName
	    <<" data luminosity: "<<dataLuminosity<<"\n";
	break;
      }
    }
  }

  // open output root file
  TFile *outputFile=new TFile(outputFileName,"recreate");

  // read binning schemes

  // binning schemes from XML
  TUnfoldBinning *recEventBinning,*genEventBinning,
     *recTrackBinning,*genTrackBinning;

  recEventBinning=TUnfoldBinningXML::ImportXML(document,"recEvent");
  genEventBinning=TUnfoldBinningXML::ImportXML(document,"genEvent");
  recTrackBinning=TUnfoldBinningXML::ImportXML(document,"recTrack");
  genTrackBinning=TUnfoldBinningXML::ImportXML(document,"genTrack");

  error=0;
  if(!recEventBinning) {
     cout<<"main: could not read recEvent binning\n";
     error++;
  }
  if(!genEventBinning) {
     cout<<"main: could not read genEvent binning\n";
     error++;
  }
  if(!recTrackBinning) {
     cout<<"main: could not read recTrack binning\n";
     error++;
  }
  if(!genTrackBinning) {
     cout<<"main: could not read genTrack binning\n";
     error++;
  }

  if(error) return 4;

  TUnfoldBinning *recBinning,*genBinning;
  // construct final "big" binning schemes
  recBinning=TUnfoldBinningXML::ImportXML(document,"recEvent");
  recBinning->SetName("recBinning");
  genBinning=TUnfoldBinningXML::ImportXML(document,"genEvent");
  genBinning->SetName("genBinning");

  // add Track binning for each rec node which has one unconnected bin
  for(TUnfoldBinning *node=recBinning;node;) {
     if(node->HasUnconnectedBins()) {
        node->AddBinning(TUnfoldBinningXML::ImportXML(document,"recTrack"));
     }
     // process children
     if(node->GetChildNode()) node=node->GetChildNode();
     // if no children, process next
     else if(node->GetNextNode()) node=node->GetNextNode();
     else {
        // go back to parent
        while((node=node->GetParentNode())) {
           // process next 
           if(node->GetNextNode()) {
              node=node->GetNextNode();
              break;
           }
           // otherwise to grand-parent
        }
     }
  }
  // add Track binning for each gen node which has one unconnected bin
  for(TUnfoldBinning *node=genBinning;node;) {
     if(node->HasUnconnectedBins()) {
        node->AddBinning(TUnfoldBinningXML::ImportXML(document,"genTrack"));
     }
     // process children
     if(node->GetChildNode()) node=node->GetChildNode();
     // if no children, process next
     else if(node->GetNextNode()) node=node->GetNextNode();
     else {
        // go back to parent
        while((node=node->GetParentNode())) {
           // process next 
           if(node->GetNextNode()) {
              node=node->GetNextNode();
              break;
           }
           // otherwise to grand-parent
        }
     }
  }

  // write binning schemes as XML to file
  ostringstream binningXML;
  TUnfoldBinningXML::ExportXML(*recEventBinning,binningXML,true,false);
  TUnfoldBinningXML::ExportXML(*genEventBinning,binningXML,false,false);
  TUnfoldBinningXML::ExportXML(*recTrackBinning,binningXML,false,false);
  TUnfoldBinningXML::ExportXML(*genTrackBinning,binningXML,false,false);
  TUnfoldBinningXML::ExportXML(*recBinning,binningXML,false,false);
  TUnfoldBinningXML::ExportXML(*genBinning,binningXML,false,true);
  TObjString binningXMLroot(binningXML.str().c_str());
  binningXMLroot.Write("binningXML");
  recEventBinning->Write();
  genEventBinning->Write();
  recTrackBinning->Write();
  genTrackBinning->Write();
  recBinning->Write();
  genBinning->Write();
  

  // print binning schemes
  cout<<"\nDetector Binning as seen by TUnfoldBinning::PrintStream\n";
  cout<<"========================================================\n";
  recBinning->PrintStream(cout);

  cout<<"\nGenerator Binning as seen by TUnfoldBinning::PrintStream\n";
  cout<<"========================================================\n";
  genBinning->PrintStream(cout);

  // parse binning schemes and determine variables
  VarList recVariables,genVariables;
  cout<<"\nDetector bininng as seen by ClassifierBinning\n";
  cout<<"========================================================\n";
  ClassifierBinning recClassifier(recBinning,recVariables);

  cout<<"\nGenerator bininng as seen by ClassifierBinning\n";
  cout<<"========================================================\n";
  ClassifierBinning genClassifier(genBinning,genVariables);

  cout<<"\nDefine variables\n";
  recVariables.AddVar("w_mini");
  recVariables.AddVar("eventpass_mini");
  recVariables.AddVar("nRECtrack_mini");
  recVariables.AddVar("nucliaREC_mini");
  recVariables.AddVar("passREC_mini");
#ifdef USE_FST
  recVariables.AddVar("typeChgREC_mini");
#endif
  genVariables.AddVar("nMCtrack_mini");
#ifdef LABFRAME
  recVariables.AddVar("etaREC_mini");
  genVariables.AddVar("pxMC_mini");
  genVariables.AddVar("pyMC_mini");
  genVariables.AddVar("etaMC_mini");
#else
  recVariables.AddVar("etaStarREC_mini");
  genVariables.AddVar("etaStarMC_mini");
#endif

  // binnig information for etaStar
  TVectorD const *etaGenBinning=0,*etaRecBinning=0;
  TVectorD const *ntrackGenBinning=0,*ntrackRecBinning=0;
  for(int i=0;i<genTrackBinning->GetDistributionDimension();i++) {
#ifdef LABFRAME
     if(genTrackBinning->GetDistributionAxisLabel(i)=="etaMC_mini") {
        etaGenBinning=genTrackBinning->GetDistributionBinning(i);
     }
#else
     if(genTrackBinning->GetDistributionAxisLabel(i)=="etaStarMC_mini") {
        etaGenBinning=genTrackBinning->GetDistributionBinning(i);
     }
#endif
     if(genTrackBinning->GetDistributionAxisLabel(i)=="ntrackMC") {
        ntrackGenBinning=genTrackBinning->GetDistributionBinning(i);
     }
  }
  for(int i=0;i<recTrackBinning->GetDistributionDimension();i++) {
#ifdef LABFRAME
     if(recTrackBinning->GetDistributionAxisLabel(i)=="etaREC_mini") {
        etaRecBinning=recTrackBinning->GetDistributionBinning(i);
     }
#else
     if(recTrackBinning->GetDistributionAxisLabel(i)=="etaStarREC_mini") {
        etaRecBinning=recTrackBinning->GetDistributionBinning(i);
     }
#endif
     if(recTrackBinning->GetDistributionAxisLabel(i)=="ntrackREC") {
        ntrackRecBinning=recTrackBinning->GetDistributionBinning(i);
     }
  }
  cout<<"\neta(star) binning\n";
  cout<<"========================================================\n";
  if(etaGenBinning) {
     cout<<"nbin="<<etaGenBinning->GetNrows()-1;
     for(int i=0;i<etaGenBinning->GetNrows();i++) {
        cout<<" "<<(*etaGenBinning)(i);
     }
     cout<<"\n";
  }
  if(etaRecBinning) {
     cout<<"nbin="<<etaRecBinning->GetNrows()-1;
     for(int i=0;i<etaRecBinning->GetNrows();i++) {
        cout<<" "<<(*etaRecBinning)(i);
     }
     cout<<"\n";
  }

  // read information about minitrees
  cout<<"\nParsing file "<<xmlMinitree<<"\n";
  cout<<"========================================================\n";

  map<string,map<string,vector<MinitreeDescriptor > > > minitreeFile;

  document=minitreeParser.GetXMLDocument();
  root=document->GetRootNode();
  if(root && (!TString(root->GetNodeName()).CompareTo("Minitree")) &&
     (root->GetNodeType()==TXMLNode::kXMLElementNode)) {
    // add all <InputFile> nodes
    for(TXMLNode *node=root->GetChildren();node;node=node->GetNextNode()) {
      if(node->GetNodeType()==TXMLNode::kXMLElementNode &&
	 !TString(node->GetNodeName()).CompareTo("InputFile") &&
	 node->GetAttributes()) {
	TIterator *i=node->GetAttributes()->MakeIterator();
	TXMLAttr *attr;
	string type,name;
        MinitreeDescriptor minitreeDescriptor;
        minitreeDescriptor.fLuminosity=-1.;
	while((attr=(TXMLAttr *)i->Next())) {
           TString attrName(attr->GetName());
           if(!attrName.CompareTo("name")) {
              name=attr->GetValue();
           }
           if(!attrName.CompareTo("type")) {
              type=attr->GetValue();
           }
           if(!attrName.CompareTo("luminosity")) {
              minitreeDescriptor.fLuminosity=TString(attr->GetValue()).Atof();
           }
           if(!attrName.CompareTo("treename")) {
              minitreeDescriptor.fTreeName=attr->GetValue();
           }
	}
        TString path=node->GetChildren()->GetContent();
	int l;
	do {
	  l=path.Length();
	  path=path.Strip(TString::kLeading,'\n');
	  path=path.Strip(TString::kLeading,' ');
	  path=path.Strip(TString::kTrailing,'\n');
	  path=path.Strip(TString::kTrailing,' ');
	} while(l!=path.Length());

	minitreeDescriptor.fFileName=path;
	cout<<"adding minitree: \""<<minitreeDescriptor.fFileName
	    <<"\" name: \""<<name<<"\""
	    <<"\" type: \""<<type<<"\""
            <<" treename=\""<<minitreeDescriptor.fTreeName<<"\""
            <<" luminosity="<<minitreeDescriptor.fLuminosity<<"\n";
	minitreeFile[type][name].push_back(minitreeDescriptor);
     }
    }
  }

  cout<<"\nLoop over files\n";
  cout<<"========================================================\n";

  for(map<string,map<string,vector<MinitreeDescriptor > > >::const_iterator
         iType=minitreeFile.begin();iType!=minitreeFile.end();iType++) {

     bool fillRec=true;
     bool fillGen=(*iType).first=="signal";

     for(map<string,vector<MinitreeDescriptor > >::const_iterator
            iFileSet=(*iType).second.begin();iFileSet!=(*iType).second.end();
         iFileSet++) {
        // create output directory
        outputFile->cd();
        TDirectoryFile *outputDir=new TDirectoryFile
           ((*iFileSet).first.c_str(),(*iFileSet).first.c_str());
        outputFile->Add(outputDir);
        for(size_t ifile=0;ifile<(*iFileSet).second.size();ifile++) {
           cout<<"\nOpen type="<<(*iType).first
               <<" name="<<(*iFileSet).first<<" file="<<ifile<<"/"
               <<(*iFileSet).second.size()<<"\n";
           // 
           // access TTree
           TFile *inputFile=TFile::Open
              ((*iFileSet).second[ifile].fFileName.c_str());
           if(!inputFile) {
              cout<<"could not open file "
                  <<(*iFileSet).second[ifile].fFileName<<"\n";
              continue;
           }
           TTree *tree;
           inputFile->GetObject((*iFileSet).second[ifile].fTreeName.c_str(),
                                tree);
           if(!tree) {
              cout<<"failed to read TTree named "
                  <<(*iFileSet).second[ifile].fTreeName.c_str()<<"\n";
              delete inputFile;
              continue;
           }
           outputDir->cd();
           cout<<"reading "<<tree->GetEntries()<<"\n";

           // SetBranchAddress for all variables to be read
           if(fillRec) {
              recVariables.SetBranchAddress(tree);
           }
           if(fillGen) {
              genVariables.SetBranchAddress(tree);
           }

           // book histograms
           TH1 *hist_gen=0;
           TH1 *hist_rec=0;
           TH2 *hist_recCovar=0;
           TH2 *hist_genRec=0;
           TH1 *hist_fake=0;
           if(fillGen) {
              hist_gen=genBinning->CreateHistogram("hist_gen",false);
           }
           if(fillRec) {
              hist_rec=recBinning->CreateHistogram("hist_rec",false);
              hist_recCovar=
                 recBinning->CreateErrorMatrixHistogram("hist_recCovar",false);
           }
           if(fillGen && fillRec) {
              hist_genRec=TUnfoldBinning::CreateHistogramOfMigrations
                 (genBinning,recBinning,"hist_genRec");
              hist_fake=recBinning->CreateHistogram("hist_fake",false);
           }

           // fast access inside event loop
           VarData const *weight=recVariables.FindVar("w_mini");
           VarData const *eventpass_mini=recVariables.FindVar("eventpass_mini");
           VarData const *nRECtrack_mini=recVariables.FindVar("nRECtrack_mini");
           VarData const *nucliaREC_mini=
              recVariables.FindVar("nucliaREC_mini");
           VarData const *passREC_mini=recVariables.FindVar("passREC_mini");

           VarData const *nMCtrack_mini=genVariables.FindVar("nMCtrack_mini");
#ifdef USE_FST
           VarData const *typeChgREC_mini=recVariables.FindVar("typeChgREC_mini");
#endif
#ifdef LABFRAME
           VarData const *pxMC_mini=genVariables.FindVar("pxMC_mini");
           VarData const *pyMC_mini=genVariables.FindVar("pyMC_mini");
           VarData const *etaREC_mini=
              recVariables.FindVar("etaREC_mini");
           VarData const *etaMC_mini=
              genVariables.FindVar("etaMC_mini");
#else
           VarData const *etaStarREC_mini=
              recVariables.FindVar("etaStarREC_mini");
           VarData const *etaStarMC_mini=
              genVariables.FindVar("etaStarMC_mini");
#endif
           double lumiWeight=1.0;
           // loop over events and fill histograms
           if((dataLuminosity>0.)&&
              ((*iFileSet).second[ifile].fLuminosity>0.)) {
              lumiWeight=dataLuminosity/(*iFileSet).second[ifile].fLuminosity;
           }
           int print=0;
           Long_t nEventMax=tree->GetEntries();
           for(Long_t ievt=0;ievt<nEventMax;ievt++) {
              tree->GetEntry(ievt);
              if((ievt%200000)==0) cout<<"ievt="<<ievt<<"\n";
              // skip non-reconstructed events
              // if there are no gen-level histograms
              bool isReconstructed=(eventpass_mini->Int()>0)
#ifdef USE_FST
                 || (typeChgREC_mini->Int()==4)
                 || (typeChgREC_mini->Int()==-4)
#endif
                 ;
              if((!fillGen)&&(!isReconstructed)) continue;
              // event weight
              double w=weight->Double()*lumiWeight;
              //
              if(print) {
                 print--;
                 if(!print) break;
              }
              // find generator-level bins 
              //  this is the bin number to count events
              int genEventBin=0;
              //  these are the bins in track multiplicity
              //    key: eta bin
              //    value.first: global multiplicity bin number
              //    value.second: multiplicity
              map<int,pair<int,int> > genMultBins;
              
              if(fillGen) {
                 ClassifierBinning const *nodeGen=
                    genClassifier.FindNode(genVariables,print?1:0);
                 if(nodeGen) {
                    // count tracks in eta bins
                    //   index: eta bin index
                    //   value: track multiplicity
                    vector<int>
                       genTrackMultiplicity(etaGenBinning->GetNrows()-1);
                    int nGenTrack=nMCtrack_mini->Int();
                    for(int t=0;t<nGenTrack;t++) {
#ifdef LABFRAME
                       double etaLabGen=etaMC_mini->Double(t);
                       double ptLabGen=
                          hypot(pxMC_mini->Double(t),pyMC_mini->Double(t));
#else
                       double etaStarGen=etaStarMC_mini->Double(t);
#endif
                       for(int k=0;k<etaGenBinning->GetNrows()-1;k++) {
#ifdef LABFRAME
                          if(ptLabGen<0.15) continue;
                          if((etaLabGen>=(*etaGenBinning)[k])&&
                             (etaLabGen<(*etaGenBinning)[k+1])) {
                             genTrackMultiplicity[k]++;
                             break;
                          }
#else
                          if((etaStarGen>=(*etaGenBinning)[k])&&
                             (etaStarGen<(*etaGenBinning)[k+1])) {
                             genTrackMultiplicity[k]++;
                             break;
                          }
#endif
                       }
                    }
                    // find bins
                    // event bin
                    if(nodeGen->GetBinningNodeUnconnected()) {
                       genEventBin=nodeGen->GetBinningNodeUnconnected()
                          ->GetStartBin();
                    }
                    // track multiplicity bins
                    TUnfoldBinning const *genDist=
                       nodeGen->GetBinningNodeDistribution();
                    if(genDist) {
                       // for each eta bin with nonzero multiplicity
                       //  locate the multiplicity bin in the distribution
                       //  and store it
                       for(size_t i=0;i<genTrackMultiplicity.size();i++) {
                          int multiplicity=genTrackMultiplicity[i];
                          if(!multiplicity) continue;
                          double eta=0.5*((*etaGenBinning)[i]+
                                          (*etaGenBinning)[i+1]);
                          double multiplicity_max=(*ntrackGenBinning)
                             (ntrackGenBinning->GetNrows()-1);
                          if(multiplicity>multiplicity_max) {
                             multiplicity=(int)multiplicity_max;
                          }
                          int ibin=genDist->GetGlobalBinNumber
                             (multiplicity,eta); // nTrack,eta
                          if(print) {
                             //cout<<"   gen eta="<<eta<<" mult="<<(*i).second<<"\n";
                          }
                          if(ibin) {
                             // remember bin number and multiplicity
                             genMultBins[i]=
                                make_pair(ibin,genTrackMultiplicity[i]);
                          }
                       }
                    }
                 }
              }
              // find reco-level bins
              //   this is the bin number to cound events
              int recEventBin=0;
              // these are the track multiplicities
              //   first key : eta bin number
              //   second key : multiplicity bin
              //        value : probability to have this track multiplicity
              //                given the track-by-track efficiency corrections
              map<int,map<int,double> > recMultBins;
              if(fillRec && isReconstructed) {
                 ClassifierBinning const *nodeRec=
                    recClassifier.FindNode(recVariables,print?1:0);
                 if(nodeRec) {
                    // count tracks in eta bins
                    int nRecTrack=nRECtrack_mini->Int();
                    // vector of eta bins
                    // in each eta bin there is a vector of probabilies
                    //  (probability to have the given multiplicity)
                    vector<vector<double> >
                       recTrackMultiplicity(etaRecBinning->GetNrows()-1);
                    for(size_t k=0;k<recTrackMultiplicity.size();k++) {
                       recTrackMultiplicity[k].reserve(100);
                       // start with 100% probability of zero
                       // multiplicity
                       //   recTrackMultiplicity.size()==1
                       //  and
                       //   recTrackMultiplicity[0]=1.0
                       recTrackMultiplicity[k].push_back(1.0);
                    }
                    for(int t=0;t<nRecTrack;t++) {
                       // reject bad tracks
                       if(!passREC_mini->Int(t)) continue;
#ifdef LABFRAME
                       double etaLabRec=etaREC_mini->Double(t);
#else
                       double etaStarRec=etaStarREC_mini->Double(t);
#endif
                       double trackEff=nucliaREC_mini->Double(t);
                       // locate eta bin
                       for(int k=0;k<etaRecBinning->GetNrows()-1;k++) {
                          if(
#ifdef LABFRAME
                             (etaLabRec>=(*etaRecBinning)[k])&&
                             (etaLabRec<(*etaRecBinning)[k+1])
#else
                             (etaStarRec>=(*etaRecBinning)[k])&&
                             (etaStarRec<(*etaRecBinning)[k+1])
#endif
                             ) {
                             // increase track multiplicity by one with
                             //  probability =trackEff
                             // and keep original multiplicity with
                             //  probability = 1-trackEff
                             vector<double> &mult=recTrackMultiplicity[k];

                             // maximum multiplicity +1
                             mult.resize(mult.size()+1);
                             for(size_t tm=mult.size()-1;tm;tm--) {
                                // probability before adding the new track
                                double &orig=mult[tm-1];
                                // multiplicity+1 has probability
                                // increased by trackEff*orig
                                mult[tm]+=trackEff*orig;
                                // this multiplicity has probabilty
                                // reduced by (1-trackEff)
                                orig *= (1.-trackEff);
                             }
                             break;
                          }
                       }
                    }
                    // find bins
                    // event count
                    if(nodeRec->GetBinningNodeUnconnected()) {
                       recEventBin=nodeRec->GetBinningNodeUnconnected()
                          ->GetStartBin();
                    }
                    // track multiplicity bins
                    TUnfoldBinning const *recDist=
                       nodeRec->GetBinningNodeDistribution();
                    if(recDist) {
                       for(size_t i=0;i<recTrackMultiplicity.size();i++) {
                          double eta=0.5*((*etaRecBinning)[i]+
                                          (*etaRecBinning)[i+1]);
                          // assign probabilities to have
                          // multiplicities >0 to histogram bins
                          for(size_t m=1;m<recTrackMultiplicity[i].size();m++) {
                             int ibin=recDist->GetGlobalBinNumber(m,eta);
                             recMultBins[i][ibin]+=recTrackMultiplicity[i][m];
                          }
                       }
                    }
                 }
              }

              // here we have:
              // (1) pointers to reco and gen event bin
              // (2) pointers to reco and gen track multiplicity bins


              if(fillGen) {
                 // fill gen level distribution
                 if(genEventBin>0) {
                    hist_gen->Fill(genEventBin,w);
                    for(map<int,pair<int,int> >::const_iterator
                           ibin=genMultBins.begin();
                        ibin!=genMultBins.end();ibin++) {
                       hist_gen->Fill((*ibin).second.first,w);
                    }
                 }
              }

              if(fillRec) {
                 // fill reco bins and covariance matrix
                 // determine vector of reco bins
                 map<int,double> binsWithWeight;
                 // event count
                 binsWithWeight[recEventBin]=w;
                 // track multiplicity counts
                 for(map<int,map<int,double> >::const_iterator
                        i=recMultBins.begin();i!=recMultBins.end();i++) {
                    for(map<int,double>::const_iterator im=(*i).second.begin();
                        im!=(*i).second.end();im++) {
                       if((*im).second<=0.) continue;
                       binsWithWeight[(*im).first]=(*im).second * w;
                    }
                 }
                 // fill distribution and covariance matrix
                 // using vector of reco bins
                 for( map<int,double>::const_iterator ibin=
                         binsWithWeight.begin();ibin!=binsWithWeight.end();
                      ibin++) {
                    hist_rec->Fill((*ibin).first,(*ibin).second);
                    for( map<int,double>::const_iterator jbin=
                            binsWithWeight.begin();jbin!=binsWithWeight.end();
                      jbin++) {
                       hist_recCovar->Fill((*ibin).first,(*jbin).first,
                                           (*ibin).second*(*jbin).second);
                    }
                 }
              }

              
              
              if(fillRec && fillGen) {
                 if(!genEventBin) {
                    // fakes (event not in gen binning scheme)
                    if(recEventBin) {
                       hist_fake->Fill(recEventBin,w);
                       // loop iver tmultiplicity eta bins
                       for(map<int,map<int,double> >::const_iterator
                              i=recMultBins.begin();i!=recMultBins.end();i++) {
                          // loop over muultiplicity bins
                          for(map<int,double>::const_iterator
                                 im=(*i).second.begin();im!=(*i).second.end();
                              im++) {
                             // add up multiplicity probabilities
                             hist_fake->Fill((*im).first,w*(*im).second);
                          }
                       }
                    }
                 } else {
                    // fill matrix of migrations
                    // the "reco" contribution
                    // in the genRec matrix 
                    // for each "gen" bin which is set in this event
                    // have to sum up to the total event weight "w"
                    // the difference to "w" is calculated
                    // and is finally added to the reco underflow bin

                    double genEvent_sum=w;
                    // key : bin number
                    // value : sum of weights
                    map<int,double> genMultBins_sum;
                    for(map<int,pair<int,int> >::const_iterator
                           genPtr=genMultBins.begin();
                        genPtr!=genMultBins.end();genPtr++) {
                       genMultBins_sum[(*genPtr).second.first]=w;
                    }
                    if(recEventBin>0) {
                       // count events
                       hist_genRec->Fill(genEventBin,recEventBin,w);
                       genEvent_sum -=w;
                       
                       // count multiplicities
                       // for each non-zero rec multiplicity bin
                       //   make sure the multiplicity weight is filled "somewhere"
                       // loop over recc eta bins
                       for(map<int,map<int,double> >::const_iterator
                              i=recMultBins.begin();i!=recMultBins.end();i++) {
                          // find corresponging gen eta bin
                          map<int,pair<int,int> >::const_iterator genPtr=
                             genMultBins.find((*i).first);
                          if(genPtr!=genMultBins.end()) {
                             // for these (reco) tracks there are (gen) tracks
                             // in the sane eta bin
                             // fill migration matrix using track probabilities
                             for(map<int,double>::const_iterator
                                    im=(*i).second.begin();im!=(*i).second.end();
                                 im++) {
                                if((*im).second<=0.) continue;
                                double wtrack=w*(*im).second;
                                hist_genRec->Fill((*genPtr).second.first, // gen bin
                                                  (*im).first, // rec bin
                                                  wtrack); // event weight
                                genMultBins_sum[(*genPtr).second.first] -= wtrack;
                             }
                          } else {
                             // for these (reco) tracks there is no (gen) track
                             // in the same eta bin
                             // (a) if there are tracks in other eta bins
                             //    -> reco tracks could be from eta migrations
                             //    -> fill off-diagonals in eta
                             //    -> if there are several eta bins,
                             //       share by track multiplicity weight
                             // (b) if there are no other tracks
                             
                             // determine total track multiplicity
                             // across all (gen) eta bins
                             double nTrackOther=0.;
                             for(genPtr=genMultBins.begin();
                                 genPtr!=genMultBins.end();genPtr++) {
                                nTrackOther += (*genPtr).second.second;
                             }
                             if(nTrackOther>0.) {
                                // (a)
                                // loop over reco probabilities
                                for(map<int,double>::const_iterator
                                       im=(*i).second.begin();
                                    im!=(*i).second.end();im++) {
                                   if((*im).second<=0.) continue;
                                   double wTrack=w*(*im).second;
                                   for(genPtr=genMultBins.begin();
                                       genPtr!=genMultBins.end();genPtr++) {
                                      double wTrackBin=
                                         wTrack*(*genPtr).second.second/
                                         nTrackOther;
                                      hist_genRec->Fill
                                         ((*genPtr).second.first, // gen bin
                                          (*im).first, // rec bin
                                          wTrackBin); // event weight
                                      genMultBins_sum[(*genPtr).second.first]-=
                                         wTrackBin;
                                   }
                                }
                             } else {
                                // (b)
                                for(map<int,double>::const_iterator
                                       im=(*i).second.begin();
                                    im!=(*i).second.end();im++) {
                                   if((*im).second<=0.) continue;
                                   double wtrack=w*(*im).second;
                                   hist_genRec->Fill(genEventBin, // gen bin
                                                     (*im).first, // rec bin
                                                     wtrack); // event weight
                                   genEvent_sum -=wtrack;
                                }
                             }
                          }
                       }
                    }
                    // finally, fill underflow bins
                    //   to account for proper nornalisation of gen bins
                    hist_genRec->Fill((double)genEventBin,(double)0,genEvent_sum);
                    for(map<int,double>::const_iterator
                           genPtr=genMultBins_sum.begin();
                        genPtr!=genMultBins_sum.end();genPtr++) {
                       hist_genRec->Fill((double)(*genPtr).first,(double)0,
                                         (*genPtr).second);
                    }
                 }
              }

              if(print) {
                 cout<<"genEventBin "<<genEventBin<<"\n";
                 for(map<int,pair<int,int> >::const_iterator
                        i=genMultBins.begin();
                     i!=genMultBins.end();i++) {
                    cout<<" genMult["<<(*i).first<<"] "<<(*i).second.first
                        <<"("<<(*i).second.second<<")"<<"\n";;
                 }

                 cout<<"recEventBin "<<recEventBin<<"\n";
                 for(map<int,map<int,double> >::const_iterator
                        i=recMultBins.begin();i!=recMultBins.end();i++) {
                    cout<<" recMult["<<(*i).first<<"]";
                    for(map<int,double>::const_iterator im=(*i).second.begin();
                        im!=(*i).second.end();im++) {
                       if((*im).second>0.) { 
                          cout<<" "<<(*im).first<<"="
                              <<(*im).second;
                       }
                    }
                    cout<<"\n";
                 }
              }
           }

           // save histograms to output directory
           outputDir->cd();
           if(hist_gen) hist_gen->Write();
           if(hist_rec) hist_rec->Write();
           if(hist_recCovar) hist_recCovar->Write();
           if(hist_genRec) hist_genRec->Write();
           if(hist_fake) hist_fake->Write();

           delete tree;
        }
     }
  }


  // close output root file
  delete outputFile;

  return 0;
}
