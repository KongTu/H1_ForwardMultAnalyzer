
#include "VarList.h"
#include "ClassifierBinning.h"


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
#include <unistd.h>
#include <cmath>
#include <set>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

using namespace std;

extern char *optarg;
extern int optind, opterr, optopt;

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

TUnfoldBinning *ReadBinning(TXMLDocument const *document,
          char const *nodeNames[],
          char const *distribution);

int main(int argc, char * const argv[]) {
  // switch on histogram errors
  TH1::SetDefaultSumw2();

  TString progname(argv[0]);

  TString selectedDataset;
  TString xmlBinning;
  TString xmlMinitree="minitree.xml";

  bool help=false;
  bool writeDTD=false;
  int opt;
  while( (opt=getopt(argc,argv,"hb:m:s:d"))!=-1) {
     switch(opt) {
     case 'h':
        help=true;
        break;
     case 'b':
        xmlBinning=optarg;
        break;
     case 'm':
        xmlMinitree=optarg;
        break;
     case 's':
        selectedDataset=optarg;
        break;
     case 'd':
        writeDTD=true;
        break;
     }
  }

  if(xmlBinning.Length()==0) help=true;

  if(help) {
     cout<<"usage: "<<progname<<" -b binning.xml "
         <<"[-m minitree.xml] [-d dataset] [-h] [-d]\n";
     cout<<"   -h : print this text\n";
     cout<<"   -b : binning scheme ["<<xmlBinning<<"]\n";
     cout<<"   -m : minitree descriptor ["<<xmlMinitree<<"]\n";
     cout<<"   -s : dataset name ["<<selectedDataset<<"]\n";
     cout<<"   -d : write DTDs binning.dtd,minitree.dtd\n";
     exit(1);
  }

  if(writeDTD) {
     static char const *dtdData=
"<!-- TUnfold Version V17.7 -->\n"
"<!ELEMENT TUnfoldBinning (BinningNode+,Histograms*,UnfoldingInput,UnfoldingOutput) >\n"
"<!ELEMENT BinningNode (BinningNode*|(Binfactorlist?,Axis)|Bins) >\n"
"<!ATTLIST BinningNode name CDATA #REQUIRED firstbin CDATA \"-1\"\n"
"    factor CDATA \"1.\">\n"
"<!ELEMENT Axis ((Bin+,Axis?)|(Axis)) >\n"
"<!ATTLIST Axis name CDATA #REQUIRED lowEdge CDATA #REQUIRED>\n"
"<!ELEMENT Binfactorlist (#PCDATA)>\n"
"<!ATTLIST Binfactorlist length CDATA #REQUIRED>\n"
"<!ELEMENT Bin EMPTY>\n"
"<!ATTLIST Bin width CDATA #REQUIRED location CDATA #IMPLIED\n"
"    center CDATA #IMPLIED repeat CDATA #IMPLIED>\n"
"<!ELEMENT Bins (BinLabel)* >\n"
"<!ATTLIST Bins nbin CDATA #REQUIRED>\n"
"<!ELEMENT BinLabel EMPTY>\n"
"<!ATTLIST BinLabel index CDATA #REQUIRED name CDATA #REQUIRED>\n"
"<!ELEMENT Histograms (Histogram+)>\n"
"<!ELEMENT Histogram (mergeBins)* >\n"
"<!ATTLIST Histogram name CDATA #REQUIRED distribution CDATA #REQUIRED option CDATA \"\">\n"
"<!ELEMENT mergeBins EMPTY >\n"
"<!ATTLIST mergeBins firstCenter CDATA #REQUIRED lastCenter CDATA #REQUIRED>\n"
"<!ELEMENT UnfoldingInput (#PCDATA) >\n"
"<!ATTLIST UnfoldingInput dataluminosity CDATA #REQUIRED>\n"
"<!ELEMENT UnfoldingOutput (#PCDATA) >\n"
"<!ATTLIST UnfoldingOutput>\n";
     ofstream dtd("binning.dtd");
     dtd<<dtdData;
  }
  if(writeDTD) {
     static char const *dtdData=
"<!ELEMENT Minitree (InputFile)+ >\n"
"<!ELEMENT InputFile (#PCDATA) >\n"
"<!ATTLIST InputFile name CDATA #REQUIRED type CDATA #REQUIRED luminosity CDATA #REQUIRED treename CDATA #REQUIRED>\n";
     ofstream dtd("minitree.dtd");
     dtd<<dtdData;
  }

  cout<<"main: binning is described in: \""<<xmlBinning<<"\"\n";
  cout<<"main: minitrees are described in: \""<<xmlMinitree<<"\"\n";
  cout<<"main: selected dataset: \""<<selectedDataset<<"\"\n";

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
  TFile *outputFile;
  if(selectedDataset.Length()) {
     outputFile=new TFile(selectedDataset+"_"+outputFileName,"recreate");
  } else {
     outputFile=new TFile(outputFileName,"recreate");
  }
  
  // read binning schemes

  // construct binning scheme for matrix of migrations
  char const *Q2yREC[]={"recQ2","recY",0};
  char const *Q2yGEN[]={"genQ2","genY",0};
  TUnfoldBinning *recBinning=ReadBinning(document,Q2yREC,"recNtrack");
  TUnfoldBinning *genBinning=ReadBinning(document,Q2yGEN,"genNtrack");

  // pointer to all possible grid points (in Q2,y)

  VarList dummyVars;
  TUnfoldBinning *covBinning=TUnfoldBinningXML::ImportXML(document,"etaBins");
  covBinning->SetName("covBinning");
  for(TUnfoldBinning *node=covBinning;node;) {
     if(!node->GetChildNode()) {
        node->AddBinning(TUnfoldBinningXML::ImportXML(document,"recNtrack"));
        if(node->GetNextNode()) {
           node=node->GetNextNode();
        } else {
           while((node=node->GetParentNode())) {
              if(node->GetNextNode()) {
                 node=node->GetNextNode();
                 break;
              }
           }
        }
     } else {
        node=node->GetChildNode();
     }
  }

  ClassifierBinning covClassifier(covBinning,dummyVars);
  // bins along eta
  vector<ClassifierBinning const *> covClasses=covClassifier.Enumerate();

  covBinning->PrintStream(cout);

  // write binning schemes as XML to file
  ostringstream binningXML;
  TUnfoldBinningXML::ExportXML(*recBinning,binningXML,true,false);
  TUnfoldBinningXML::ExportXML(*genBinning,binningXML,false,false);
  TUnfoldBinningXML::ExportXML(*covBinning,binningXML,false,true);

  TObjString binningXMLroot(binningXML.str().c_str());
  binningXMLroot.Write("binningXML");
  recBinning->Write();
  genBinning->Write();
  covBinning->Write();

  // print binning schemes
  cout<<"\nDetector Binning as seen by TUnfoldBinning::PrintStream\n";
  cout<<"========================================================\n";
  recBinning->PrintStream(cout);

  cout<<"\nGenerator Binning as seen by TUnfoldBinning::PrintStream\n";
  cout<<"========================================================\n";
  genBinning->PrintStream(cout);

  cout<<"\nCovariance Binning as seen by TUnfoldBinning::PrintStream\n";
  cout<<"========================================================\n";
  covBinning->PrintStream(cout);

  // parse binning schemes and determine variables
  VarList recVariables,genVariables;
  cout<<"\nDetector bininng as seen by ClassifierBinning\n";
  cout<<"========================================================\n";
  ClassifierBinning recClassifier(recBinning,recVariables);

  vector<ClassifierBinning const *> recClasses=recClassifier.Enumerate();

  cout<<"\nGenerator bininng as seen by ClassifierBinning\n";
  cout<<"========================================================\n";
  ClassifierBinning genClassifier(genBinning,genVariables);
  
  cout<<"\nDefine variables\n";
  recVariables.AddVar("w_mini");
  recVariables.AddVar("eventpass_mini");
  recVariables.AddVar("nRECtrack_mini");
  recVariables.AddVar("etaREC_mini");
  recVariables.AddVar("nucliaREC_mini");
  recVariables.AddVar("passREC_mini");
  genVariables.AddVar("nMCtrack_mini");
  genVariables.AddVar("etaMC_mini");
  genVariables.AddVar("pxMC_mini");
  genVariables.AddVar("pyMC_mini");
  genVariables.AddVar("isQEDcMC_mini");
  genVariables.AddVar("isQEDbkg_mini");
  genVariables.AddVar("isDaughtersMC_mini");

  cout<<"\neta binning\n";
  cout<<"========================================================\n";
  for(size_t k=0;k<covClasses.size();k++) {
     cout<<" eta bin "<<covClasses[k]->GetName()<<"\n";
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

           if(selectedDataset.Length() &&
              (selectedDataset !=(*iFileSet).first)) continue;

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
           vector<TH1 *> hist_gen,hist_rec,hist_fake,hist_QEDc; // eta bins
           map<ClassifierBinning const *,TH2 *> hist_recCovar;  // Q2,y bins
           vector<TH2 *> hist_genRec; // eta bins
           for(size_t k=0;k<covClasses.size();k++) {
              TString name=covClasses[k]->GetName();
              if(fillGen) {
                 hist_gen.push_back
                    (genBinning->CreateHistogram("hist_gen_"+name,false));
                 if(fillRec) {
                    hist_genRec.push_back
                       (TUnfoldBinning::CreateHistogramOfMigrations
                        (genBinning,recBinning,"hist_genRec_"+name));
                    hist_fake.push_back
                       (recBinning->CreateHistogram("hist_fake_"+name,false));
                    hist_QEDc.push_back
                       (recBinning->CreateHistogram("hist_QEDc_"+name,false));  
                 }
              }
           }
           if(fillRec) {
              for(size_t k=0;k<covClasses.size();k++) {
                 hist_rec.push_back
                    (recBinning->CreateHistogram
                     ("hist_rec_"+covClasses[k]->GetName(),false));
              }
              for(size_t igrid=0;igrid<recClasses.size();igrid++) {
                 hist_recCovar[recClasses[igrid]]=
                    (covBinning->CreateErrorMatrixHistogram
                     (TString("hist_recCovar_")+
                      recClasses[igrid]->GetName(),false));
              }
           }

           // fast access inside event loop
           VarData const *weight=recVariables.FindVar("w_mini");
           VarData const *eventpass_mini=recVariables.FindVar("eventpass_mini");
           VarData const *nRECtrack_mini=recVariables.FindVar("nRECtrack_mini");
           VarData const *etaREC_mini=recVariables.FindVar("etaREC_mini");
           VarData const *nucliaREC_mini=recVariables.FindVar("nucliaREC_mini");
           VarData const *passREC_mini=recVariables.FindVar("passREC_mini");

           VarData const *nMCtrack_mini=genVariables.FindVar("nMCtrack_mini");
           VarData const *etaMC_mini=genVariables.FindVar("etaMC_mini");
           VarData const *pxMC_mini=genVariables.FindVar("pxMC_mini");
           VarData const *pyMC_mini=genVariables.FindVar("pyMC_mini");
           VarData const *isQEDcMC_mini=genVariables.FindVar("isQEDcMC_mini");
           VarData const *isQEDbkg_mini=genVariables.FindVar("isQEDbkg_mini");
           VarData const *isDaughtersMC_mini=genVariables.FindVar("isDaughtersMC_mini");

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
              bool isReconstructed=(eventpass_mini->Int()>0);
              if((!fillGen)&&(!isReconstructed)) continue;
              // event weight
              double w=weight->Double()*lumiWeight;
              //
              if(print) {
                 print--;
                 if(!print) break;
              }

              //gen level signal and background;
              int QEDc = isQEDcMC_mini->Int();
              int QEDbkg = isQEDbkg_mini->Int();
              bool isSignal=true;
              if(QEDc==1||QEDbkg==1) isSignal=false;
              //  these are the bin numbers in (Q2,y,track multiplicity)
              //    index: eta bin
              //    value: global bin number
              vector<int> genMultBins(covClasses.size());
              if(fillGen&&isSignal) {
                 // nodeGen will point to the (Q2,y) bin
                 // or it is zero
                 ClassifierBinning const *nodeGen=
                    genClassifier.FindNode(genVariables,print?1:0);
                 if(nodeGen) {
                    //
                    // loop over eta bins
                    vector<int> genTrackMultiplicity(covClasses.size());
                    int nGenTrack=nMCtrack_mini->Int();
                    for(int t=0;t<nGenTrack;t++) {
                       double etaGen=etaMC_mini->Double(t);
                       double ptGen=TMath::Hypot(pxMC_mini->Double(t),pyMC_mini->Double(t));
                       if(fabs(etaGen)>1.6||ptGen<0.15) continue;
                       if(isDaughtersMC_mini->Int(t)!=0) continue;//add selections on nonV0s on gen
                       for(size_t k=0;k<covClasses.size();k++) {
                          if(covClassifier.IsInside(etaGen,k)) {
                             genTrackMultiplicity[k]++;
                          }
                       }
                    }
                    // track multiplicity bins
                    TUnfoldBinning const *genDist=
                        nodeGen->GetBinningNodeDistribution();
                    if(genDist) {
                       // for each eta bin, locate the multiplicity bin
                       for(size_t i=0;i<genTrackMultiplicity.size();i++) {
                          int multiplicity=genTrackMultiplicity[i];
                          int ibin=genDist->GetGlobalBinNumber(multiplicity);
                          // remember bin number and multiplicity
                          genMultBins[i]=ibin;
                        }
                    }
                  }
              }
              // these are the track multiplicities
              //   index : eta bin number
              //   second key : multiplicity bin
              //        value : probability to have this track multiplicity
              //                given the track-by-track efficiency corrections
              vector<map<int,double> > recMultBins(covClasses.size());
              // index: bin number in covariance histogram
              //         (encodes eta and multiplicity)
              // value : probability
              map<int,double> recMultCovarBins;
              //
              // pointer to reco (Q2,x)
              ClassifierBinning const *nodeRec=0;
              if(fillRec && isReconstructed) {
                 // locate Q2,y bin
                 nodeRec=recClassifier.FindNode(recVariables,print?1:0);
                 if(nodeRec) {
                    // track multiplicity bins
                    TUnfoldBinning const *recDist=
                       nodeRec->GetBinningNodeDistribution();
                    // count tracks in eta bins
                    int nRecTrack=nRECtrack_mini->Int();
                    // vector of eta bins
                    // in each eta bin there is a vector of probabilies
                    //  (probability to have the given multiplicity)
                    vector<vector<double> >
          recTrackMultiplicity(covClasses.size());
                    for(size_t k=0;k<recTrackMultiplicity.size();k++) {
                       recTrackMultiplicity[k].reserve(100);
                       // start with 100% probability of zero
                       recTrackMultiplicity[k].push_back(1.0);
                    }
                    for(int t=0;t<nRecTrack;t++) {
                       // reject bad tracks
                       if(!passREC_mini->Int(t)) continue;
                       double etaRec=etaREC_mini->Double(t);
                       if(fabs(etaRec)>1.6) continue;
                       double trackEff=nucliaREC_mini->Double(t);
                       // locate eta bin
                       for(size_t k=0;k<covClasses.size();k++) {
                          if(covClassifier.IsInside(etaRec,k)) {
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
                          }
                       }
                    }
                    if(recDist) {
                       for(size_t i=0;i<recTrackMultiplicity.size();i++) {
                          // assign probabilities to have
                          // given multiplicity to histogram bins
                          for(size_t m=0;m<recTrackMultiplicity[i].size();m++) {
                             int ibin=recDist->GetGlobalBinNumber(m);
                             recMultBins[i][ibin]+=recTrackMultiplicity[i][m];
                             int covarBin=covClasses[i]->
                                GetBinningNodeDistribution()->
                                GetGlobalBinNumber(m);
                             recMultCovarBins[covarBin]
                                += recTrackMultiplicity[i][m];
                          }
                       }
                    }
                 }
              }

              // here we have pointers to reco and gen track multiplicity bins

              if(fillGen) {
                 // fill gen level distribution
                 for(size_t ieta=0;ieta<genMultBins.size();ieta++) {
                    int imultGenBin=genMultBins[ieta];
                    if(imultGenBin) {
                       hist_gen[ieta]->Fill(imultGenBin,w);
                    } else {
                       // fake -> not counted as generator truth
                    }
                 }
              }

              if(fillRec) {
                 // fill reco bins and covariance matrices
                 // track multiplicity counts
                 for(size_t ieta=0;ieta<recMultBins.size();ieta++) {
                    for(map<int,double>::const_iterator iMultPtr=
                           recMultBins[ieta].begin();
                        iMultPtr!=recMultBins[ieta].end();iMultPtr++) {
                       int iMultRecBin=(*iMultPtr).first;
                       double iMultWeight=(*iMultPtr).second;
                       hist_rec[ieta]->Fill(iMultRecBin,w*iMultWeight);
                       if(!isSignal) hist_QEDc[ieta]->Fill(iMultRecBin, w*iMultWeight);
                    }
                 }
                 // covariance in bins of Q2,y
                 if(nodeRec) {
                    for(map<int,double>::const_iterator iCov=
                           recMultCovarBins.begin();
                        iCov!=recMultCovarBins.end();iCov++) {
                       for(map<int,double>::const_iterator jCov=
                              recMultCovarBins.begin();
                           jCov!=recMultCovarBins.end();jCov++) {
                          //cout<<grid<<"/"<<hist_recCovar.size()
                          //    <<" "<<(*iCov).first<<" "<<(*iCov).second<<"\n";
                          hist_recCovar[nodeRec]->Fill
                             ((*iCov).first,(*jCov).first,
                              w*(*iCov).second*(*jCov).second);
                       }
                    }
                 }
              }

              if(fillRec && fillGen) {
                 for(size_t ieta=0;ieta<genMultBins.size();ieta++) {
                    int iMultGenBin=genMultBins[ieta];
                    if(!iMultGenBin) {
                       // fakes for this eta bin
                       // (event not in gen binning scheme)
                       for(map<int,double>::const_iterator iMultPtr=
                              recMultBins[ieta].begin();
                           iMultPtr!=recMultBins[ieta].end();iMultPtr++) {
                          int iMultRecBin=(*iMultPtr).first;
                          double iMultWeight=(*iMultPtr).second;
                          hist_fake[ieta]->Fill(iMultRecBin,w*iMultWeight);
                       }
                    } else {
                       // fill matrix of migrations (for each ieta)
                       for(map<int,double>::const_iterator iMultPtr=
                              recMultBins[ieta].begin();
                           iMultPtr!=recMultBins[ieta].end();iMultPtr++) {
                          int iMultRecBin=(*iMultPtr).first;
                          double iMultWeight=(*iMultPtr).second;
                          hist_genRec[ieta]->Fill(iMultGenBin,iMultRecBin,
                                                  w*iMultWeight);

                       }
                    }
                 }
              }

              if(print) {
                 for(size_t i=0;i<genMultBins.size();i++) {
                    cout<<" genMultBins["<<i<<"]="<<genMultBins[i]<<"\n";
                 }

                 for(size_t i=0;i<recMultBins.size();i++) {
                    cout<<" recMultBins["<<i<<"]";
                    for(map<int,double>::const_iterator
                           im=recMultBins[i].begin();
                        im!=recMultBins[i].end();im++) {
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
           for(size_t i=0;i<hist_gen.size();i++) hist_gen[i]->Write();
           for(size_t i=0;i<hist_rec.size();i++) hist_rec[i]->Write();
           for( map<ClassifierBinning const *,TH2 *>::const_iterator
                   i=hist_recCovar.begin();
                i!=hist_recCovar.end();i++) (*i).second->Write();
           for(size_t i=0;i<hist_genRec.size();i++) hist_genRec[i]->Write();
           for(size_t i=0;i<hist_fake.size();i++) hist_fake[i]->Write();
           for(size_t i=0;i<hist_QEDc.size();i++) hist_QEDc[i]->Write();
           delete tree;
        }
     }
  }
  
  
  // close output root file
  delete outputFile;

  return 0;
}

TUnfoldBinning *ReadBinning(TXMLDocument const *document,
          char const *nodeNames[],
          char const *distribution) {
   // construct binning scheme from document
   TUnfoldBinning *r=TUnfoldBinningXML::ImportXML(document,nodeNames[0]);
   for(TUnfoldBinning *node=r->GetChildNode();node;
       node=node->GetNextNode()) {
      if(nodeNames[1]) {
         node->AddBinning(ReadBinning(document,nodeNames+1,distribution));
      } else {
         node->AddBinning(TUnfoldBinningXML::ImportXML(document,distribution));
      }
   }
   return r;
}