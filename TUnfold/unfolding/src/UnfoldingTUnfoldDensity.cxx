#include "UnfoldingTUnfoldDensity.h"
#include "UnfoldingResult.h"


Unfolding_TUnfoldDensity::~Unfolding_TUnfoldDensity() {
   DeleteIf(fTunfold);
   if(fScanTauDistribution) delete fScanTauDistribution;
   if(fScanTauProjectionMode) delete fScanTauProjectionMode;
}

void Unfolding_TUnfoldDensity::SetScanLCurve
(int nPoint,double logTauMin,double logTauMax) {
   fLCurveNpoint=nPoint;
   fLCurveLogTauMin=logTauMin;
   fLCurveLogTauMax=logTauMax;
}

void Unfolding_TUnfoldDensity::SetScanTau
(int nPoint,double logTauMin,double logTauMax,TUnfoldDensity::EScanTauMode mode,
 const char *distribution,const char *projectionMode) {
   fScanTauNpoint=nPoint;
   fScanTauLogTauMin=logTauMin;
   fScanTauLogTauMax=logTauMax;
   fScanTauMode=mode;
   if(fScanTauDistribution) delete fScanTauDistribution;
   if(fScanTauProjectionMode) delete fScanTauProjectionMode;
   fScanTauDistribution=0;
   if(distribution) fScanTauDistribution=strdup(fScanTauDistribution);
   fScanTauProjectionMode=0;
   if(projectionMode) fScanTauProjectionMode=strdup(projectionMode);
}

void Unfolding_TUnfoldDensity::InitUnfolding
(const TH2 *hist_gen_rec,const TH1 *hist_gen,const TH1 *hist_fakes) {
   DeleteIf(fTunfold);
   fTunfold=new TUnfoldDensity
      (hist_gen_rec,TUnfold::kHistMapOutputHoriz,
       fRegMode,fConstraintMode,fDensityFlags,
       fGeneratorBinning,fDetectorBinning);
   if(hist_fakes) {
      fTunfold->SubtractBackground(hist_fakes,"BGR_FAKES");
   }
}

void Unfolding_TUnfoldDensity::RunUnfolding(const TH1 *data) {
   gErrorIgnoreLevel=kWarning+1;
   fTunfold->SetInput(data,fBiasScale,fOneOverZeroError);
   gErrorIgnoreLevel=kInfo+1;

   // get result without regularisation
   fTunfold->DoUnfold(0.);
   SaveResults("tau=0");

   // if fixed tau is requested, repeat unfolding
   if(fFixedTau>0) {
      fTunfold->DoUnfold(fFixedTau);
      SaveResults(TString::Format("tau=%g",fFixedTau));
   }

   // do L-curve scan if requested
   if(fLCurveNpoint>0) {
      TSpline *logTauX,*logTauY;
      TGraph *lCurve=0;

      int iBest=fTunfold->ScanLcurve
	(fLCurveNpoint,fLCurveLogTauMin,fLCurveLogTauMax,
         &lCurve,&logTauX,&logTauY);
      Double_t t[1],x[1],y[1];
      logTauX->GetKnot(iBest,t[0],x[0]);
      logTauY->GetKnot(iBest,t[0],y[0]);
      //TGraph *bestLcurve=new TGraph(1,x,y);
      //TGraph *bestLogTauLogChi2=new TGraph(1,t,x);
      UnfoldingResult *result=SaveResults("LCurveScan");
      result->AddVariable("LCurveScan_iBest",iBest);
      result->AddVariable("LCurveScan_logtau",t[0]);
      result->AddVariable("LCurveScan_logX",x[0]);
      result->AddVariable("LCurveScan_logY",y[0]);
      result->AddObject("LCurveScan_logTauX",logTauX);
      result->AddObject("LCurveScan_logTauY",logTauY);
      result->AddObject("LCurveScan_lCurve",lCurve);
   }

   // do tau scan if requested
   if(fScanTauNpoint>0) {
      TSpline *scanResult;
      TGraph *lCurvePlot;
      TSpline *logTauXPlot,*logTauYPlot;
      int iBest=fTunfold->ScanTau
         (fScanTauNpoint,fScanTauLogTauMin,fScanTauLogTauMax,
          &scanResult,fScanTauMode,fScanTauDistribution,fScanTauProjectionMode,
          &lCurvePlot,&logTauXPlot,&logTauYPlot);
      Double_t t[1],x[1],y[1];
      logTauXPlot->GetKnot(iBest,t[0],x[0]);
      logTauYPlot->GetKnot(iBest,t[0],y[0]);
      UnfoldingResult *result=SaveResults("TauScan");
      result->AddVariable("TauScan_iBest",iBest);
      result->AddVariable("TauScan_logtau",t[0]);
      result->AddVariable("TauScan_logX",x[0]);
      result->AddVariable("TauScan_logY",y[0]);
      result->AddObject("TauScan_scan",scanResult);
      result->AddObject("TauScan_logTauX",logTauXPlot);
      result->AddObject("TauScan_logTauY",logTauYPlot);
      result->AddObject("TauScan_lCurve",lCurvePlot);
   }
}

UnfoldingResult *Unfolding_TUnfoldDensity::SaveResults(const char *name) {
   TString tag=GetName()+"_"+name;
   UnfoldingResult *&result=fUnfoldingResult[tag];
   if(result) delete result;
   TH1 *hist_unfold=fTunfold->GetOutput(tag+"_output",0,0,0,kFALSE);
   TH2 *hist_dataEmatrix=
      fTunfold->GetEmatrixInput(tag+"_dataEMAT",0,0,0,kFALSE);
   // TUnfold includes background errors in its default uncertainties
   // the UnfoldingResult, however, expect the errors to originate 
   // from the data alone
   // -> set errors from diagonals of error matrix
   for(int iGen=0;iGen<=hist_unfold->GetNbinsX()+1;iGen++) {
      hist_unfold->SetBinError
         (iGen,sqrt(hist_dataEmatrix->GetBinContent(iGen,iGen)));
   }
   if(fDoSyst) {
      result=new UnfoldingResult
         (hist_unfold,
          hist_dataEmatrix,
          fTunfold->GetEmatrixSysUncorr(tag+"_MCstatEMAT",0,0,0,kFALSE),
          fTunfold->GetEmatrixSysBackgroundUncorr
          ("BGR_FAKES",tag+"_fakesEMAT",0,0,0,kFALSE)
          );
   } else {
      result=new UnfoldingResult
         (hist_unfold,
          hist_dataEmatrix,
          0,0);
   }
   return result;
}

