#include <TProfile2D.h>
#include <TTree.h>
//  #include <iostream>

#include "Correlations.h"
#include "UnfoldingResult.h"
#include "UnfoldingAverage.h"

using namespace std;

UnfoldingAverage::UnfoldingAverage(void) {
   fProfMean=0;
   fProfMeanSQ=0;
   fProfAbserr=0;
   fProfRelerr=0;
   fProfDelta=0;
   fProfDeltarel=0;
   fProfPull=0;
   fHistPull=0;
   fProfPullRMS=0;

   fProfDataEMAT=0;
   fProfMCstatEMAT=0;
   fProfFakesEMAT=0;

   //fToyEMAT=0;
   //fCorrEmat=0;
   //fCorrToy=0;
}

UnfoldingAverage::~UnfoldingAverage(void) {
   // basic histograms and profiles filled from unfolded histogram
   // and associated errors
   DeleteIf(fProfMean);
   DeleteIf(fProfMeanSQ);
   DeleteIf(fProfAbserr);
   DeleteIf(fProfRelerr);
   DeleteIf(fProfDelta);
   DeleteIf(fProfDeltarel);
   DeleteIf(fProfPull);
   DeleteIf(fHistPull);
   DeleteIf(fProfPullRMS);

   // average error matrices
   DeleteIf(fProfDataEMAT);
   DeleteIf(fProfMCstatEMAT);
   DeleteIf(fProfFakesEMAT);

   // analysis of toy data
   //DeleteIf(fToyEMAT);
   //if(fCorrEmat) delete fCorrEmat;
   //if(fCorrToy) delete fCorrToy;
}

void UnfoldingAverage::SaveRoot(void) const {
   WriteIf(fProfMean,"profUnfolded");
   WriteIf(fProfMeanSQ,"profUnfoldedSQ");
   WriteIf(fProfAbserr,"profAbserror");
   WriteIf(fProfRelerr,"profRelerror");
   WriteIf(fProfDelta,"profDelta");
   WriteIf(fProfDeltarel,"profDeltarel");
   WriteIf(fProfPull,"profPull");
   WriteIf(fHistPull,"histPull");
   WriteIf(fProfPullRMS,"profPullRMS");
   WriteIf(fProfDataEMAT,"profDataEMAT");
   WriteIf(fProfMCstatEMAT,"profMCstatEMAT");
   WriteIf(fProfFakesEMAT,"profFakesEMAT");
   //WriteIf(fToyEMAT,"toyEMAT");
   //if(fCorrEmat) fCorrEmat->SaveRoot("corrEmat");
   //if(fCorrToy) fCorrToy->SaveRoot("corrToy");
   //cout<<"Write TTree\n";

   if(fVarsData.size()) {
      //std::cout<<"create TTree size="<<fVarsData.size()<<"\n";
      TTree *tree=new TTree("extra","extra");
      vector<Float_t> vars(fVars.size());
      for(map<TString,int>::const_iterator iv=fVars.begin();iv!=fVars.end();
          iv++) {
         //cout<<"create branch "<<(*iv).first<<"\n";
         tree->Branch((*iv).first,&vars[(*iv).second],(*iv).first);
      }
      for(size_t k=0;k<fVarsData.size();k++) {
         for(size_t l=0;l<fVarsData[k].size();l++) {
            vars[l]=fVarsData[k][l];
         }
         tree->Fill();
      }
      //cout<<"write tree\n";
      tree->Write();
      //cout<<"delete tree\n";
      delete tree;
      //cout<<"done delete\n";
   }
   //cout<<"done write\n";
}

void UnfoldingAverage::Add(UnfoldingResult const *result,
                           const TH1 *hist_truth) {
   if(!fProfMean) {
      static int uniqueID=0;
      TString base=TString::Format("n%d_",uniqueID);
      uniqueID++;
      int n=hist_truth->GetNbinsX();
      double x0=0.5;
      double x1=n+0.5;
      fProfMean=new TProfile(base+"prof_mean",";bin",n,x0,x1);
      fProfMeanSQ=new TProfile2D(base+"prof_meanSQ",";bin;bin",n,x0,x1,n,x0,x1);
      fProfAbserr=new TProfile(base+"prof_abserr",";bin",n,x0,x1);
      fProfRelerr=new TProfile(base+"prof_relerr",";bin",n,x0,x1);
      fProfDelta=new TProfile(base+"prof_delta",";bin",n,x0,x1);
      fProfDeltarel=new TProfile(base+"prof_deltarel",";bin",n,x0,x1);
      fProfPull=new TProfile(base+"prof_pull",";bin",n,x0,x1);
      fHistPull=new TH2D(base+"hist_pull",";bin",n,x0,x1,100,-10.,10.);
      fProfPullRMS=new TProfile(base+"prof_pullRMS",";bin",n,x0,x1,"s");
      //fToyEMAT=new TH2D(base+"hist_toyEMAT",";bin",n,x0,x1,n,x0,x1);
   }
   for(int i=0;i<=hist_truth->GetNbinsX()+1;i++) {
      double x=result->GetUnfoldedData()->GetBinContent(i);
      double e=result->GetUnfoldedData()->GetBinError(i);;
      fProfMean->Fill(i,x);
      for(int j=0;j<=hist_truth->GetNbinsX()+1;j++) {
         double y=result->GetUnfoldedData()->GetBinContent(j);
         fProfMeanSQ->Fill(i,j,x*y);
      }
      fProfAbserr->Fill(i,e);
      double t=hist_truth->GetBinContent(i);
      double delta=x-t;
      fProfDelta->Fill(i,delta);
      if(t>0.) fProfDeltarel->Fill(i,delta/t);
      if(x>0.0) fProfRelerr->Fill(i,e/x);
      double pull=delta/e;
      fProfPull->Fill(i,pull);
      fHistPull->Fill(i,pull);
      fProfPullRMS->Fill(i,pull);
   }
   // fill error matrix from toys
   /*for(int i=0;i<=hist_truth->GetNbinsX()+1;i++) {
      for(int j=0;j<=hist_truth->GetNbinsX()+1;j++) {
         fToyEMAT->SetBinContent
            (i,j,fProfMeanSQ->GetBinContent(i,j)-
             fProfMean->GetBinContent(i)*fProfMean->GetBinContent(j));
      }
      } */

   FillProf(result->GetDataErrorMatrix(),fProfDataEMAT);
   FillProf(result->GetMCstatErrorMatrix(),fProfMCstatEMAT);
   FillProf(result->GetFakesErrorMatrix(),fProfFakesEMAT);

   // handle extra variables, to be saved to a TTree later
   map<TString,double> const &extra=result->GetExtraVars();
   vector<double> rowData(fVars.size());
   for(map<TString,double>::const_iterator iv=extra.begin();iv!=extra.end();
       iv++) {
      map<TString,int >::iterator varPtr=fVars.find((*iv).first);
      if(varPtr==fVars.end()) {
         fVars[(*iv).first]=rowData.size();
         rowData.push_back((*iv).second);
      } else {
         rowData[(*varPtr).second]=(*iv).second;
      }
   }
   fVarsData.push_back(rowData);

}

/*void UnfoldingAverage::UpdateCorrelations(TH1 const *hist_truth) {
   if(fProfDataEMAT) {
      if(!fCorrEmat) {
         fCorrEmat=new Correlations();
      }
      fCorrEmat->Update(fProfDataEMAT,fProfMean,hist_truth);
   }
   if(!fCorrToy) {
      fCorrToy=new Correlations();
   }
   fCorrToy->Update(fToyEMAT,fProfMean,hist_truth);
   } */

void UnfoldingAverage::FillProf(TH2 const *src,TProfile2D *&dest) const {
   if(src) {
      if(!dest) {
         dest=new TProfile2D
            (TString("prof_")+src->GetName(),";bin;bin",
             src->GetNbinsX(),0.5,src->GetNbinsX()+0.5,
             src->GetNbinsY(),0.5,src->GetNbinsY()+0.5);
      }
      for(int i=0;i<=src->GetNbinsX()+1;i++) {
         for(int j=0;j<=src->GetNbinsX()+1;j++) {
            dest->Fill(i,j,src->GetBinContent(i,j));
         }
      }
   }
}

