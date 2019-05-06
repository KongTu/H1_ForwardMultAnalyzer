#include "UnfoldingIterative.h"
#include "UnfoldingResult.h"

#include <TMatrixDSymEigen.h>
#include <TH2.h>
#include <TH1.h>

using namespace std;

//#define CHECK_VECTOR_MATRIX

#ifdef CHECK_VECTOR_MATRIX

#include <iostream>

#define CHECK_POSITIVE(a) CheckPositive(a,#a)
#define CHECK_FINITE(a) CheckFinite(a,#a)

void CheckPositive(TVectorD const &v,char const *name) {
   int error=0;
   for(int i=0;i<v.GetNrows();i++) {
      if(!(finite(v(i))&&(v(i)>=0))) {
         cout<<"vector "<<name<<"("<<i<<")="<<v(i)<<"\n";
         error++;
      }
   }
   if(error) {
      exit(1);
   }
}
void CheckPositive(TMatrixD const &m,char const *name) {
   int error=0;
   for(int i=0;i<m.GetNrows();i++) {
      for(int j=0;j<m.GetNcols();j++) {
         if(!(finite(m(i,j))&&(m(i,j)>=0))) {
            cout<<"matrix "<<name<<"("<<i<<","<<j<<")="<<m(i,j)<<"\n";
            error++;
         }
      }
   }
   if(error) {
      exit(1);
   }
}
void CheckFinite(TMatrixD const &m,char const *name) {
   int error=0;
   for(int i=0;i<m.GetNrows();i++) {
      for(int j=0;j<m.GetNcols();j++) {
         if(!finite(m(i,j))) {
            cout<<"matrix "<<name<<"("<<i<<","<<j<<")="<<m(i,j)<<"\n";
            error++;
         }
      }
   }
   if(error) {
      exit(1);
   }
}

#else
#define CHECK_POSITIVE(a) { } 
#define CHECK_FINITE(a) { }

#endif

void Unfolding_Iterative::InitUnfolding
(const TH2 *hist_gen_rec,const TH1 *hist_gen,const TH1 *hist_fakes) {
   // initialize iterative unfolding
   // step minus one: set x to truth
   fHist_gen_rec=hist_gen_rec;
   // gen: include underflow and overflow
   int nGen=hist_gen_rec->GetNbinsX()+2;
   // rec excludes underflow and overflow
   int nRec=hist_gen_rec->GetNbinsY();
   fX0.ResizeTo(nGen);
   fA.ResizeTo(nRec,nGen);
   fOneOverEps.ResizeTo(nGen);
   for(int iGen=0;iGen<nGen;iGen++) {
      double nTotal=hist_gen->GetBinContent(iGen);
      fX0(iGen)=nTotal;
      if(nTotal>0.) {
         double nInside=0.;
         for(int iRec=0;iRec<nRec;iRec++) {
            double M_ij= hist_gen_rec->GetBinContent(iGen,iRec+1);
            nInside+=M_ij;
            fA(iRec,iGen)= M_ij / nTotal;
         }
         if(nInside>0.) {
            fOneOverEps(iGen)=nTotal/nInside;
         }
      }
   }
   fBgr.ResizeTo(nRec);
   for(int iRec=0;iRec<nRec;iRec++) {
      fBgr(iRec)=hist_fakes->GetBinContent(iRec+1);
   }
}

void Unfolding_Iterative::RunUnfolding(const TH1 *data) {
   int nRec=fA.GetNrows();
   int nGen=fA.GetNcols();
   /* cout<<"nRec="<<nRec<<" nGen="<<nGen<<"\n";
      exit(0); */

   TVectorD y(nRec);
   for(int iRec=0;iRec<nRec;iRec++) {
      y(iRec)=data->GetBinContent(iRec+1);
   }
   TVectorD bgrAdd(nRec);
   if(fBackgroundSubtractive) {
      y -= fBgr;
   } else {
      bgrAdd=fBgr;
   }
   CHECK_POSITIVE(y);
   CHECK_POSITIVE(bgrAdd);

   fSaved_minRhoiSQavg=false;

   // set start values
   TVectorD x(fX0);
   CHECK_POSITIVE(x);
   // transposed matrix A, normalized by epsilon
   TMatrixD AToverEps(TMatrixD::kTransposed,fA);
   AToverEps.NormByColumn(fOneOverEps,"");
   CHECK_POSITIVE(AToverEps);

   //TMatrixD Vxx(nGen,nGen);
   //TMatrixD Vxy(nGen,nRec);

   // how x depends on the start value
   TMatrixD dx_dx0(nGen,nGen);
   for(int iGen=0;iGen<nGen;iGen++) {
      dx_dx0(iGen,iGen)=1.;
   }
   // how x depends on y
   TMatrixD dx_dy(nGen,nRec);
   // covariance of x (from data)
   TMatrixDSym Vxx(nGen);

   fRhoI_SQavg.resize(0);
   fRhoI_SQavg.reserve(fMaxStep+1);

   for(int step=0;step<=fMaxStep;step++) {

      TVectorD yrec(fA*x+bgrAdd);
      CHECK_POSITIVE(yrec);
      // also called r  below.
      TVectorD yOverYrec(nRec);
      for(int iRec=0;iRec<nRec;iRec++) {
         if(y(iRec)!=0.) {
            yOverYrec(iRec)=y(iRec)/yrec(iRec);
         }
      }
      CHECK_POSITIVE(yOverYrec);
      TVectorD f=AToverEps * yOverYrec;
      CHECK_POSITIVE(f);
      TVectorD z(nGen);
      for(int iGen=0;iGen<nGen;iGen++) {
         if(fOneOverEps(iGen)<=0.0) {
            z(iGen)=x(iGen);
         } else {
            z(iGen)=x(iGen)*f(iGen);
         }
      }
      //   x * df/dr
      TMatrixD xdf_dr(AToverEps);
      for(int i=0;i<nGen;i++) {
         for(int j=0;j<nRec;j++) {
            xdf_dr(i,j) *= x(i);
         }
      }
      CHECK_POSITIVE(xdf_dr);
      TMatrixD dr_dx(nRec,nGen);
      TVectorD dr_dy(nRec);
      for(int j=0;j<nRec;j++) {
         if(yrec(j)!=0.) {
            dr_dy(j)=1.0/yrec(j);
         }
         for(int i=0;i<nGen;i++) {
            double c=yOverYrec(j)*fA(j,i);
            if(c!=0.0) {
               dr_dx(j,i)= -c/yrec(j);
            }
         }
      }
      CHECK_FINITE(dr_dx);
      CHECK_POSITIVE(dr_dy);
      TMatrixD dz_dx(xdf_dr*dr_dx);
      for(int i=0;i<nGen;i++) {
         dz_dx(i,i) +=f(i);
      }
      CHECK_FINITE(dz_dx);
      TMatrixD dz_dy(xdf_dr);
      dz_dy.NormByRow(dr_dy,"");
      dz_dy+=dz_dx*dx_dy;
      CHECK_FINITE(dz_dy);

      // here: 
      //  dz_dx : how z=x[step] depends on x[step-1]
      //  dz_dy : how z=x[step] depends on y

      // now: update x, dx_dx0, dx_dy

      // estimated uncertainties from y errors
      TVectorD ey(nRec);
      if(fPoissonErrors) {
         ey=fA*z+fBgr;
         for(int iRec=0;iRec<nRec;iRec++) {
            ey(iRec)=sqrt(ey(iRec));
         }
      } else {
         for(int iRec=0;iRec<nRec;iRec++) {
            ey(iRec)=data->GetBinError(iRec+1);
         }
      }
      CHECK_POSITIVE(ey);
      TMatrixD dz_dyE_T(TMatrixD::kTransposed,dz_dy);
      dz_dyE_T.NormByColumn(ey,"");
      CHECK_FINITE(dz_dyE_T);
      TMatrixDSym Vzz(TMatrixDSym::kAtA,dz_dyE_T);
      CHECK_FINITE(Vzz);
      // calculate correlation coefficients
      TMatrixDSymEigen evVzz(Vzz);
      TVectorD eVal=evVzz.GetEigenValues();
      TMatrixD eVec=evVzz.GetEigenVectors();
      TVectorD VzzInv_ii(nGen);
      for(int k=0;k<nGen;k++) {
         double ev=eVal(k);
         if(ev>0.0) {
            for(int i=0;i<nGen;i++) {
               VzzInv_ii(i)+= 1./ev*eVec(i,k)*eVec(i,k);
            }
         }
      }
      CHECK_POSITIVE(VzzInv_ii);
      double s0=0.,s1=0.;
      for(int k=0;k<nGen;k++) {
         double Vzz_ii=Vzz(k,k);
         if((Vzz_ii>0.0)&&VzzInv_ii(k)>0.0) {
            s0+=1.;
            s1+=1.-1./(Vzz_ii*VzzInv_ii(k));
         }
      }
      fRhoI_SQavg.push_back(sqrt(s1/s0));

      // decide whether to save this as a result

      TMatrixD dz_dx0(dz_dx,TMatrixD::kMult,dx_dx0);
      CHECK_FINITE(dz_dx0);
      SaveIfRequested(step,x,dx_dx0,Vxx,z,dz_dx0,Vzz);

      // prepare next iteration
      x=z;
      CHECK_POSITIVE(x);
      dx_dx0=dz_dx0;
      CHECK_FINITE(dx_dx0);
      dx_dy =dz_dy;
      CHECK_FINITE(dx_dy);
      Vxx=Vzz;
      CHECK_FINITE(Vxx);
   }
}

void Unfolding_Iterative::SaveIfRequested
(int step,
 TVectorD const &x,TMatrixD const &dx_dx0,TMatrixDSym const &Vxx,
 TVectorD const &z,TMatrixD const &dz_dx0,TMatrixDSym const &Vzz) {
   if(((step>=0)&&(step<=10))
      || ((step<=100)&&(step%10==0))
      || ((step<=1000)&&(step%100==0))
      ) {
      TString tag=GetName()+TString::Format("_iter%d",step);
      SaveTagged(tag,step,z,dz_dx0,Vzz);
   }
   if((step>0)&&(!fSaved_minRhoiSQavg)) {
      if((fRhoI_SQavg[step]>fRhoI_SQavg[step-1])||(step==fMaxStep)) {
         SaveTagged(GetName()+"_minRhoiSQavg",step-1,x,dx_dx0,Vxx);
         fSaved_minRhoiSQavg=true;
      }
   }
   
}

void Unfolding_Iterative::SaveTagged
(TString const &tag,int step,TVectorD const &x,TMatrixD const &dx_dx0,
 TMatrixDSym const &Vxx) {
   UnfoldingResult *&result=fUnfoldingResult[tag];
   int nGen=x.GetNrows();
   TH1D *h_output=new TH1D(tag+"_output",";bin",
                           nGen-2,0.5,nGen-1.5);
   TH2D *h_ematrix=
      new TH2D(tag+"_emat",";bin",
               nGen-2,0.5,nGen-1.5,nGen-2,0.5,nGen-1.5);
   CHECK_POSITIVE(x);
   for(int iGen=0;iGen<nGen;iGen++) {
      h_output->SetBinContent(iGen,x(iGen));
      h_output->SetBinError(iGen,sqrt(Vxx(iGen,iGen)));
      for(int jGen=0;jGen<nGen;jGen++) {
         h_ematrix->SetBinContent(iGen,jGen,Vxx(iGen,jGen));
      }
   }
   result=new UnfoldingResult(h_output,h_ematrix,0,0);
   
   TH1D *hist_D=new TH1D(tag+"_D",";bin",
                         nGen-2,0.5,nGen-1.5);
   TH1D *hist_Dmax=new TH1D(tag+"_Dmax",";bin",
                            nGen-2,0.5,nGen-1.5);
   double dmin=HUGE_VAL,dmax=-HUGE_VAL;
   for(int iGen=0;iGen<nGen;iGen++) {
      double dSQ=0.0;
      double dmax_i=-HUGE_VAL;
      double dmin_i=HUGE_VAL;
      for(int jGen=0;jGen<nGen;jGen++) {
         dSQ+=pow(dx_dx0(iGen,jGen),2.);
         dmin_i=fmin(dmin_i,dx_dx0(iGen,jGen));
         dmax_i=fmax(dmax_i,dx_dx0(iGen,jGen));
      }
      hist_D->SetBinContent(iGen,sqrt(dSQ));
      hist_Dmax->SetBinContent(iGen,dmax_i);
      dmin=fmin(dmin,dmin_i);
      dmax=fmax(dmax,dmax_i);
   }
   result->AddObject("histD",hist_D);
   result->AddObject("histDmax",hist_Dmax);

   result->AddVariable("nstep",(double)step);
   result->AddVariable("dmin",dmin);
   result->AddVariable("dmax",dmax);
   result->AddVariable("rhoSQavg_Nm1",(step>=0) ? fRhoI_SQavg[step-1] : 0.);
   result->AddVariable("rhoSQavg_N",fRhoI_SQavg[step]);
   if(step+1<fRhoI_SQavg.size()) {
      result->AddVariable("rhoSQavg_Np1",fRhoI_SQavg[step+1]);
   }
}
