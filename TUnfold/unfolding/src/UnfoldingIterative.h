#ifndef __H__UNFOLDINGITERATIVE__
#define __H__UNFOLDINGITERATIVE__

#include "UnfoldingAlgorithm.h"

#include <TVectorD.h>
#include <TMatrixD.h>

class Unfolding_Iterative : public UnfoldingAlgorithm {
public:
   // iterative EM unfolding
   //    L.~A.~Shepp and Y.~Vardi, IEEE Trans.~Med.~Imaging MI-1 (1982) 113
   //    A.~Kondor,  Nucl.\ Instrum.\ Meth.\ 216 (1983) 177
   // (often called D'Agostini, but his paper is much later
   //    G.~D'Agostini, Nucl.\ Instrum.\ Meth.\  A {\bf 362}, 487 (1995))
   Unfolding_Iterative(char const *name="iterative",int maxStep=100,
                       bool poissonErrors=false)
      : UnfoldingAlgorithm(name),fMaxStep(maxStep) {
      fPoissonErrors=poissonErrors;
      fBackgroundSubtractive=false;
      fA=0; fOneOverEps=0; fBgr=0; fX0=0; }
   //virtual ~Unfolding_Iterative();
   void SetBackgroundSubtractive(bool subtractive=true) {
      fBackgroundSubtractive=subtractive;
   }
   void SetPoissonErrors(bool poissonErrors=true) {
      fPoissonErrors=poissonErrors;
   }
   virtual void InitUnfolding(const TH2 *hist_gen_rec,const TH1 *hist_gen,
                              const TH1 *hist_fakes);
   virtual void RunUnfolding(const TH1 *data);
protected:

   void SaveIfRequested
   (int step,
    TVectorD const &x,TMatrixD const &dx_dx0,TMatrixDSym const &Vxx,
    TVectorD const &z,TMatrixD const &dz_dz0,TMatrixDSym const &Vzz);

   void SaveTagged(TString const &tag,int step,TVectorD const &x,
                   TMatrixD const &dx_dx0,TMatrixDSym const &Vxx);

   bool fBackgroundSubtractive,fPoissonErrors;
   TH2 const *fHist_gen_rec;
   TMatrixD fA;
   TVectorD fOneOverEps,fBgr,fX0;
   int fMaxStep;
   bool fSaved_minRhoiSQavg;
   std::vector<double> fRhoI_SQavg;
};

#endif
