#include "ClassifierBinning.h"
#include "VarList.h"
#include <TVectorD.h>

using namespace std;

ClassifierBinning::ClassifierBinning
(TUnfoldBinning const *binning,VarList &vars) {
   fParent=0;
   fBinningClassifier=0;
   fBinningDistribution=0;
   fBinningUnconnected=0;
   AddSubBins(binning,vars,0);
}

ClassifierBinning::ClassifierBinning
(TUnfoldBinning const *binning,VarList &vars,
 ClassifierBinning const *parent,int indentLevel) {
   fParent=parent;
   fBinningClassifier=0;
   fBinningDistribution=0;
   fBinningUnconnected=0;
   AddSubBins(binning,vars,indentLevel);
}

ClassifierBinning::~ClassifierBinning() {
   for(size_t i=0;i<fBins.size();i++) {
      delete fBins[i].fClassifier;
   }
}

ClassifierBinning const *ClassifierBinning::FindNode
(VarList const &vars,int print) const {
   VarList::const_iterator iv=vars.find(fName);
   string indent=string(print,' ');
   if(iv!=vars.end()) {
      double x=(*iv).second.Double();
      for(size_t i=0;i<fBins.size();i++) {
         if(IsInside(x,i)) {
            if(print) {
               cout<<indent<<"ClassifierBinning::VarWithBinning::FindNode "
                   <<fName<<" x is in ["<<fBins[i].fXmin<<","
                   <<fBins[i].fXmax<<"]\n";
            }
            ClassifierBinning const *r=
               fBins[i].fClassifier->FindNode(vars,print?(print+1):0);
            return r;
         }
      }
      if(print) {
         cout<<indent<<"ClassifierBinning::VarWithBinning::FindNode "
             <<fName<<" x="<<x<<" not in binning scheme"
             <<" binBorders["<<fBins.size()<<"]={\n";
         for(size_t i=0;i<fBins.size();i++) {
            cout<<" ["<<fBins[i].fXmin<<","<<fBins[i].fXmax<<"]";
         }
         cout<<"}\n";
      }
   } else {
      if(print && !fName.empty()) {
         cout<<indent<<"not found \""<<fName<<"\"\n";
      }
   }
   ClassifierBinning const *r=this;
   if(print) {
      cout<<indent<<"return node [";
      if(r->GetBinningNodeClassifier()) {
         cout<<" classifier="<<r->GetBinningNodeClassifier()->GetName();
      }
      if(r->GetBinningNodeDistribution()) {
         cout<<" distribution="<<r->GetBinningNodeDistribution()->GetName();
      }
      if(r->GetBinningNodeUnconnected()) {
         cout<<" unconnected="<<r->GetBinningNodeUnconnected()->GetName();
      }
      cout<<" ]\n";
   }
   return r;
}

void ClassifierBinning::AddSubBins(TUnfoldBinning const *binning,
                                   VarList &allVars,int indentLevel) {
   string indent(2*indentLevel,' ');
   //cout<<"AddSubBins "<<indent<<"node="<<binning->GetName()<<"\n";
   // check whether the children are all classifiers
   cout<<indent<<"setting classifier="<<binning->GetName()<<"\n";
   while(binning) {
      int n=0;
      int nClassifier=0;
      for(TUnfoldBinning const *node=binning->GetChildNode();node;
          node=node->GetNextNode()) {
         n++;
      }
      fBins.reserve(n);
      for(TUnfoldBinning const *node=binning->GetChildNode();node;
          node=node->GetNextNode()) {
         TString nodeName=node->GetName();
         Ssiz_t split2=nodeName.Last(':');
         if(split2!=TString::kNPOS) {
            Ssiz_t split1=nodeName.First(':');
            //   VAR:x
            // defines condition  VAR<x
            nClassifier++;
            ClassifierBin bin;
            bin.fXmin=-INFINITY;
            bin.fXmax=TString(nodeName(split2+1,1000)).Atof();
            if(split1==split2) {
               fName=nodeName(0,split1);
               if(fBins.size()) bin.fXmin=(*fBins.rbegin()).fXmax;
            } else {
               bin.fXmin=TString(nodeName(0,split1)).Atof();
               fName=nodeName(split1+1,split2-split1-1);
            }
            allVars.AddVar(fName);
            cout<<indent<<"Classifier "
                <<bin.fXmin<<"<="<<fName<<"<"<<bin.fXmax<<"\n";
            bin.fClassifier=new ClassifierBinning
               (node,allVars,this,indentLevel+1);
            fBins.push_back(bin);
         }
      }
      if(!nClassifier) {
         for(size_t k=0;k<fBins.size();k++) {
            delete fBins[k].fClassifier;
         }
         fBins.resize(0);
      }

      if(binning->HasUnconnectedBins()) {
         fBinningUnconnected=binning;
         cout<<indent<<"node="<<binning->GetName()
             <<" has "<<binning->GetDistributionNumberOfBins()
             <<" unconnected bins"<<"\n";
      }
      if(binning->GetDistributionDimension()) {
         fBinningDistribution=binning;
         cout<<indent<<"node="<<binning->GetName()
             <<" Histogram "<<binning->GetDistributionDimension()
             <<"-dim [";
         for(int i=0;i<binning->GetDistributionDimension();i++) {
            if(binning->HasUnderflow(i)) cout<<"Ufl+";
            cout<<binning->GetDistributionBinning(i)->GetNrows()-1;
            if(binning->HasOverflow(i)) cout<<"+Ofl";
            if(i<binning->GetDistributionDimension()-1)
               cout<<" x ";
            else
               cout<<"]\n";
         }
         for(int i=0;i<binning->GetDistributionDimension();i++) {
            cout<<indent<<"  "<<binning->GetDistributionAxisLabel(i)<<" {";
            if(binning->HasUnderflow(i)) {
               cout<<"-infinity,";
            }
            for(int k=0;k<binning->GetDistributionBinning(i)->GetNrows();
                k++) {
               if(k) {
                  cout<<",";
               }
               cout<<(*binning->GetDistributionBinning(i))[k];
            }
            if(binning->HasOverflow(i)) {
               cout<<",infinity";
            }
            cout<<"}\n";
         }
         cout<<"\n";
      }
      if(fBins.size()) {
         fBinningClassifier=binning;
         break;
      }
      binning=binning->GetChildNode();
   }
}

TString ClassifierBinning::GetName(void) const {
   TString r;
   ClassifierBinning const *node=this;
   for(ClassifierBinning const *parent=node->fParent;parent;
       parent=node->fParent) {
      for(size_t k=0;k<parent->fBins.size();k++) {
         if(parent->fBins[k].fClassifier==node) {
            r=TString::Format("(%g<=",parent->fBins[k].fXmin)+
               parent->fName.c_str()+
               TString::Format("<%g)",parent->fBins[k].fXmax)+r;
         }
      }
      node=parent;
   }
   return r;
}

std::vector<ClassifierBinning const *>
ClassifierBinning::Enumerate(void) const {
   std::vector<ClassifierBinning const *> r;
   Enumerate_r(r);
   return r;
}

void ClassifierBinning::Enumerate_r
(std::vector<ClassifierBinning const *> &clist) const {
   for(size_t k=0;k<fBins.size();k++) {
      if(fBins[k].fClassifier->fBins.size()) {
         fBins[k].fClassifier->Enumerate_r(clist);
      } else {
         clist.push_back(fBins[k].fClassifier);
      }
   }
}

