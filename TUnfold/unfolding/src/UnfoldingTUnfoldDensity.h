#ifndef __H__UNFOLDINGTUNFOLDDENSITY__
#define  __H__UNFOLDINGTUNFOLDDENSITY__

#include "UnfoldingAlgorithm.h"

#ifdef PRIVATE_TUNFOLD
#include "../TUnfold_V17.7/TUnfoldDensity.h"
#else
#include <TUnfoldDensity.h>
#endif


class Unfolding_TUnfoldDensity : public UnfoldingAlgorithm {
public:
   virtual ~Unfolding_TUnfoldDensity();
   // run unfolding with TUnfold
   // there are many parameters which can be chosen
   // see TUnfold manual for details
   // 
   // by default, this class is doing a template fit without reguularisation,
   //
   Unfolding_TUnfoldDensity
   (char const *name,
    TUnfoldBinning *detectorBinning,TUnfoldBinning *generatorBinning,
    TUnfold::ERegMode regMode=TUnfold::kRegModeCurvature,
    TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone,
    TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth,
    double biasScale=1.0,double oneOverZeroError=0.0)
      : UnfoldingAlgorithm(name),
      fGeneratorBinning(generatorBinning),fDetectorBinning(detectorBinning),
      fRegMode(regMode),
      fConstraintMode(constraintMode),
      fDensityFlags(densityFlags),
      fBiasScale(biasScale),
      fOneOverZeroError(oneOverZeroError),
      fDoSyst(false) {
      fTunfold=0;
      // default: do not unfold with fixed tau
      fFixedTau=0.;
      // default: no L-curve scan
      fLCurveNpoint=0;
      // default: no scan of global correlation
      fScanTauNpoint=0;
      fScanTauDistribution=0;
      fScanTauProjectionMode=0;
   }
   // enable unfolding for fixed tau
   void SetFixedTau(double tau) { fFixedTau=tau; }
   void SetBiasScale(double scale) { fBiasScale=scale; }

   // enable evaulation of MC stat and BGR systematic uncertianties
   void SetDoSyst(bool doSyst) { fDoSyst=doSyst; }

   // enable L-curve scan (set nPoint to zero to disable)
   void SetScanLCurve(int nPoint=30,double logTauMin=0.0,double logTauMax=0.);
   // enable scan of tau (set nPoint to zero to disable)
   void SetScanTau(int nPoint=30,double logTauMin=0.0,double logTauMax=0.,
                   TUnfoldDensity::EScanTauMode mode=
                   TUnfoldDensity::kEScanTauRhoMax,
                   const char *distribution=0,const char *projectionMode=0);
   // specify how to treat bins with zero uncertainty
   void SetOneOverZeroError(double oneOverZeroError=1.) {
      fOneOverZeroError=oneOverZeroError;
   }

   virtual void InitUnfolding(const TH2 *hist_gen_rec,const TH1 *hist_gen,
                              const TH1 *hist_fakes);
   virtual void RunUnfolding(const TH1 *data);   

protected:

   UnfoldingResult *SaveResults(char const *name);
   const TUnfoldBinning *fGeneratorBinning;
   const TUnfoldBinning *fDetectorBinning;
   TUnfold::ERegMode fRegMode;
   TUnfold::EConstraint fConstraintMode;
   TUnfoldDensity::EDensityMode fDensityFlags;
   double fBiasScale;
   double fOneOverZeroError;
   double fFixedTau;
   int fLCurveNpoint;
   double fLCurveLogTauMin,fLCurveLogTauMax;
   int fScanTauNpoint;
   double fScanTauLogTauMin,fScanTauLogTauMax;
   TUnfoldDensity::EScanTauMode fScanTauMode;
   char *fScanTauDistribution,*fScanTauProjectionMode;
   TUnfoldDensity *fTunfold;
   bool fDoSyst;
};

#endif
