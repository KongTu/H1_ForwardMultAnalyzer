#ifndef __H__UNFOLDINGAVERAGE__
#define __H__UNFOLDINGAVERAGE__

#include "UnfoldingBase.h"
#include <map>
#include <vector>

class UnfoldingResult;
class Correlations;
class TH1;
class TH2;
class TH2D;
class TProfile;
class TProfile2D;

class UnfoldingAverage : public UnfoldingBase {
   // averages determined from runs over toy experiments
   // provides profile plots of average, pull, uncertainties, etc
   // TODO: include determination of statistical error matrix from toys
public:
   UnfoldingAverage(void);
   ~UnfoldingAverage(void);
   void Add(UnfoldingResult const *result,const TH1 *hist_truth);
   void SaveRoot(void) const;
   //void UpdateCorrelations(TH1 const *hist_truth);
protected:
   void FillProf(TH2 const *src,TProfile2D *&dest) const;
   // basic profiles filled from unfolded histogram and errors
   TH2D *fHistPull;
   TProfile2D *fProfMeanSQ;
   TProfile *fProfDelta,*fProfDeltarel;
   TProfile *fProfMean,*fProfPull,*fProfAbserr,*fProfRelerr,*fProfPullRMS;
   // profiles filled from error matrices
   TProfile2D *fProfDataEMAT,*fProfMCstatEMAT,*fProfFakesEMAT;
   // error matrix from toy experiments
   //TH2D *fToyEMAT;
   // analysis of correlations, either from average error matrix of from toys
   //Correlations *fCorrEmat,*fCorrToy;

   std::map<TString,int > fVars;
   std::vector<std::vector<double> > fVarsData;
};

#endif
