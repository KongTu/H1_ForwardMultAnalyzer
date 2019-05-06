#include <TH1.h>
#include <TH2.h>

#include "UnfoldingResult.h"
#include "Correlations.h"

using namespace std;

UnfoldingResult::UnfoldingResult(TH1 *unfolded,TH2 *dataEmat,
                                 TH2 *MCstatEmat,TH2 *fakesEmat) {
   fUnfolded=unfolded;
   fDataError=dataEmat;
   fMCstatError=MCstatEmat;
   fFakesError=fakesEmat;
   fCorrelations=0;
   fMCstatCorrelations=0;

   if(fDataError) {
      fCorrelations=new Correlations();
      fCorrelations->Update(fDataError,fUnfolded,0);
      AddVariable("rhoImin",fCorrelations->GetRhoMin());
      AddVariable("rhoImax",fCorrelations->GetRhoMax());
      AddVariable("rhoIavg",fCorrelations->GetRhoAvg());
      AddVariable("rhoIavgSQ",fCorrelations->GetRhoSQavg());
   }
   if(fMCstatError) {
      fMCstatCorrelations=new Correlations();
      fMCstatCorrelations->Update(fMCstatError,fUnfolded,0);
      AddVariable("rhoIminMCstat",fMCstatCorrelations->GetRhoMin());
      AddVariable("rhoImaxMCstat",fMCstatCorrelations->GetRhoMax());
      AddVariable("rhoIavgMCstat",fMCstatCorrelations->GetRhoAvg());
      AddVariable("rhoIavgSQMCstat",fMCstatCorrelations->GetRhoSQavg());
   }
}

UnfoldingResult::~UnfoldingResult() {
   DeleteIf(fUnfolded);
   DeleteIf(fDataError);
   DeleteIf(fMCstatError);
   DeleteIf(fFakesError);
   if(fCorrelations) delete fCorrelations;
   if(fMCstatCorrelations) delete fMCstatCorrelations;
   for(map<TString,TObject *>::iterator io=fExtraObj.begin();
       io!=fExtraObj.end();io++) {
      DeleteIf((*io).second);
   }
}


void UnfoldingResult::AddVariable(char const *name,double value) {
   fExtraVars[name]=value;
}

void UnfoldingResult::AddObject(char const *name,TObject *obj) {
   TObject *&o=fExtraObj[name];
   if(o) delete o;
   o=obj;
}

void UnfoldingResult::SaveRoot(void) const {
   WriteIf(GetUnfoldedData(),"unfolded");
   WriteIf(GetDataErrorMatrix(),"statEMAT");
   WriteIf(GetMCstatErrorMatrix(),"mcstatEMAT");
   WriteIf(GetFakesErrorMatrix(),"fakesEMAT");
   if(fCorrelations) fCorrelations->SaveRoot("stat");
   if(fMCstatCorrelations) fMCstatCorrelations->SaveRoot("MCstat");
   for(map<TString,TObject *>::const_iterator io=fExtraObj.begin();
       io!=fExtraObj.end();io++) {
      if((*io).second) (*io).second->Write((*io).first);
   }
}
