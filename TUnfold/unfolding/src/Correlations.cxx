#include "Correlations.h"
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixDSymEigen.h>

using namespace std;

Correlations::~Correlations() {
   DeleteIf(fAbsError);
   DeleteIf(fRelError);
   DeleteIf(fPull);
   DeleteIf(fRhoI);
   DeleteIf(fRhoIJ);
}

void Correlations::Update(TH2 const *ematrix,TH1 const *data,TH1 const *truth) {
   static int uniqueID=0;
   if(!fAbsError) {
      if(ematrix) {
         TString base=TString::Format("n%d_",uniqueID)+ematrix->GetName();
         uniqueID++;
         int n=ematrix->GetNbinsX();
         double x0=0.5;
         double x1=n+0.5;
         // with error matrix, can calculate abs errors, rhoI, rhoIJ
         fAbsError=new TH1D(base+"_absError",";bin",n,x0,x1);
         fRhoI=new TH1D(base+"_rhoI",";bin",n,x0,x1);
         fRhoIJ=new TH2D(base+"_rhoIJ",";bin",n,x0,x1,n,x0,x1);
         //
         if(data) {
            // with ematrix & data, can calculate relative errors
            fRelError=new TH1D(base+"_relError",";bin",n,x0,x1);
            if(truth) {
               // with ematrix & data & truth, can calculate pulls
               fPull=new TH1D(base+"_pull",";bin",n,x0,x1);
            }
         }
      }
   }

   // fill histograms
   if(ematrix) {
      int n2=ematrix->GetNbinsX()+2;
      int nNonEmpty=0;
      vector<int> ii;
      ii.reserve(n2);
      for(int i=0;i<n2;i++) {
         double sigma_i=sqrt(ematrix->GetBinContent(i,i));
         // absolute error
         fAbsError->SetBinContent(i,sigma_i);
         if(sigma_i>0.0) {
            ii.push_back(nNonEmpty);
            nNonEmpty++;
         } else {
            ii.push_back(-1);
         }
         if(data && fRelError) {
            // relative error
            double x=data->GetBinContent(i);
            if(x>0.0) {
               fRelError->SetBinContent(i,sigma_i/x);
            } else {
               fRelError->SetBinContent(i,0.0);
            }
            if(truth && fPull) {
               // pull
               double p=x-truth->GetBinContent(i);
               if(sigma_i>0.) {
                  p/=sigma_i;
               } else if (p!=0.0) {
                  p/=sigma_i; // gives +/- infinity
               }
               fPull->SetBinContent(i,p);
            }
         }
      }
      if(nNonEmpty) {
         TMatrixDSym rho(nNonEmpty);
         for(int i=0;i<n2;i++) {
            if(ii[i]<0) continue;
            for(int j=0;j<n2;j++) {
               if(ii[j]<0) continue;
               double rho_ij=ematrix->GetBinContent(i,j)
                  /sqrt(ematrix->GetBinContent(i,i)*
                        ematrix->GetBinContent(j,j));
               rho(ii[i],ii[j])=rho_ij;
               fRhoIJ->SetBinContent(i,j,rho_ij);
            }
         }
         // invert matrix of correlation coefficients
         TMatrixDSymEigen evRho(rho);
         TVectorD eVal=evRho.GetEigenValues();
         TMatrixD eVec=evRho.GetEigenVectors();

         TVectorD rhoInvII(nNonEmpty);
         for(int k=0;k<nNonEmpty;k++) {
            double ev=eVal(k);
            if(ev<=0.0) {
               ev=HUGE_VAL;
            } else {
               ev=1./eVal(k);
            }
            for(int i=0;i<nNonEmpty;i++) {
               rhoInvII(i)+= ev*eVec(i,k)*eVec(i,k);
            }
         }

         for(int i=0;i<n2;i++) {
            double r=-1.0;
            if(ii[i]>=0) {
               r=sqrt(1.-1./rhoInvII(ii[i]));
            }
            fRhoI->SetBinContent(i,r);
         }
      }
   }
}

void Correlations::SaveRoot(TString const &prefix) {
   WriteIf(fAbsError,prefix+"AbsErr");
   WriteIf(fRelError,prefix+"RelErr");
   WriteIf(fPull,prefix+"Pull");
   WriteIf(fRhoI,prefix+"RhoI");
   WriteIf(fRhoIJ,prefix+"RhoIJ");
}

double Correlations::GetRhoMin(void) const {
   double rhoMin=999.;
   static int print=0;
   if(fRhoI) {
      for(int i=0;i<=fRhoI->GetNbinsX()+1;i++) {
         if(fAbsError->GetBinContent(i)<=0.0) continue;
         rhoMin=fmin(rhoMin,fRhoI->GetBinContent(i));
         if(print) {
            cout<<" "<<i<<" "<<fRhoI->GetBinContent(i)<<" "<<rhoMin<<"\n";
         }
      }
   }
   if(print) {
      print--;
      if(!print) exit(0);
   }
   return rhoMin;
}

double Correlations::GetRhoMax(void) const {
   double rhoMax=-1.;
   if(fRhoI) {
      for(int i=0;i<=fRhoI->GetNbinsX()+1;i++) {
         if(fAbsError->GetBinContent(i)<=0.0) continue;
         rhoMax=fmax(rhoMax,fRhoI->GetBinContent(i));
      }
   }
   return rhoMax;
}

double Correlations::GetRhoAvg(void) const {
   double s0=0.,s1=0.;
   if(fRhoI) {
      for(int i=0;i<=fRhoI->GetNbinsX()+1;i++) {
         if(fAbsError->GetBinContent(i)<=0.0) continue;
         s0++;
         s1+=fRhoI->GetBinContent(i);
      }
   }
   return s1/s0;;
}

double Correlations::GetRhoSQavg(void) const {
   double s0=0.,s1=0.;
   if(fRhoI) {
      for(int i=0;i<=fRhoI->GetNbinsX()+1;i++) {
         if(fAbsError->GetBinContent(i)<=0.0) continue;
         s0++;
         s1+=pow(fRhoI->GetBinContent(i),2);
      }
   }
   return sqrt(s1/s0);
}

