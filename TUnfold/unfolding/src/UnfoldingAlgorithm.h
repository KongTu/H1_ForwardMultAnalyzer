#ifndef __H__UNFOLDINGALGORITHM__
#define __H__UNFOLDINGALGORITHM__

#include "UnfoldingAverage.h"

class UnfoldingAlgorithm : public UnfoldingBase {
public:
   UnfoldingAlgorithm(char const *name) : fName(name) { }
   virtual ~UnfoldingAlgorithm();
   std::map<TString,UnfoldingResult *> const &Unfold(const TH1 *data);
   std::map<TString,UnfoldingAverage> const &GetAverage(void);
   void SetTruth(TH1 *hist_truth) { fHistTruth=hist_truth; }
   virtual void InitUnfolding(const TH2 *hist_gen_rec,const TH1 *hist_gen,
                              const TH1 *hist_fakes) =0;
   void ResetAverages(void);
   TString const &GetName(void) const { return fName; }
protected:
   TString fName;
   virtual void RunUnfolding(const TH1 *data)=0;
   const TH1 *fHistTruth;
   std::map<TString,UnfoldingResult *> fUnfoldingResult;
   std::map<TString,UnfoldingAverage> fUnfoldingAverage;
};

#endif
