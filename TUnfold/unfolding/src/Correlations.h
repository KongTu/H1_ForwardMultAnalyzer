#ifndef __H__CORRELATIONS__
#define __H__CORRELATIONS__

#include <UnfoldingBase.h>

class TH1;
class TH2;

class Correlations : public UnfoldingBase {
public:
   //====================================================================
   // correlation coefficients, pulls and similar, all derived from
   // a single covariance matrix (plus unfolded and truth information)
   Correlations(void) : fAbsError(0),fRelError(0),fPull(0),fRhoI(0),fRhoIJ(0) {
   }
   virtual ~Correlations();
   void Update(TH2 const *ematrix,TH1 const *data,TH1 const *truth);
   void SaveRoot(TString const &prefix);
   double GetRhoMin(void) const;
   double GetRhoMax(void) const;
   double GetRhoAvg(void) const;
   double GetRhoSQavg(void) const;
protected:
   TH1 *fAbsError,*fRelError,*fPull,*fRhoI;
   TH2 *fRhoIJ;
};

#endif
