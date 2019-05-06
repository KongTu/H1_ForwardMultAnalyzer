#ifndef __H__UNFOLDINGRESULT__
#define __H__UNFOLDINGRESULT__

#include "UnfoldingBase.h"
#include <map>

class TH1;
class TH2;
class Correlations;

class UnfoldingResult : public UnfoldingBase {
   //====================================================================
   // Results of a single unfolding
   // the following results are set in the constructor
   //    unfolded data (mandatory)
   //    covariance from data stat. uncertainties (recommended)
   //    covariance from MC stat (optional)
   //    covariance from MC fakes (optional)
   // Using AddObject() and AddVariable(), extra data can be saved
   // If the data statistical covariance is given, the corresponding histogram
   // of global correlation coefficiencts is also calculated
   // 
   // the SaveRoot method saves all histograms and extra objects
   // to the current TDirectory
   //
   // the histograms and extra objects are owned by the class UnfoldingResult.

public:
   UnfoldingResult(TH1 *unfolded,TH2 *dataEmat,TH2 *MCstatEmat,TH2 *fakesEmat);
   void AddVariable(char const *name,double value);
   void AddObject(char const *name,TObject *obj);
   ~UnfoldingResult();
   void SaveRoot(void) const;
   inline TH1 const *GetUnfoldedData(void) const { return fUnfolded; }
   inline TH2 const *GetDataErrorMatrix(void) const { return fDataError; }
   inline TH2 const *GetMCstatErrorMatrix(void) const { return fMCstatError; }
   inline TH2 const *GetFakesErrorMatrix(void) const { return fFakesError; }
   inline std::map<TString,double> const &GetExtraVars(void) const { return fExtraVars; }
   inline std::map<TString,TObject*> const &GetExtraObjects(void) const { return fExtraObj; }
   Correlations const *GetDataCorrelations(void) const { return fCorrelations; }
protected:
   std::map<TString,double> fExtraVars;
   std::map<TString,TObject *> fExtraObj;
   // main unfolding result. The unfolding algorithm shall return
   // the central values with (stat) uncertainties from data alone
   TH1 *fUnfolded;
   // other error matices, possibly set by the unfolding algorithm
   TH2 *fDataError,*fMCstatError,*fFakesError;;
   // correlation coefficients calculated from fDataError
   Correlations *fCorrelations,*fMCstatCorrelations;
};

#endif
