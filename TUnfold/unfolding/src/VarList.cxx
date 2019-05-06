#include "VarList.h"
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>

using namespace std;


VarData::VarData() : fLeaf(0) { }

void VarData::SetLeaf(TLeaf *l) {
   fLeaf=l;
}


double VarData::Double(int i) const {
   if(fLeaf) return fLeaf->GetValue(i);
   return 0.0;
}

int VarData::Int(int i) const {
   if(fLeaf) return (int)fLeaf->GetValue(i);
   return 0;
}

void VarList::AddVar(std::string const &name) {
   iterator i=find(name);
   if(i==end()) {
      (*this)[name]=VarData();
   }
}

bool VarList::SetBranchAddress(TTree *tree) {
   bool r=true;
   cout<<"VarList::SetBranches\n";
   for(iterator i=begin();i!=end();i++) {
      (*i).second.SetLeaf(tree->FindLeaf((*i).first.c_str()));
      tree->SetBranchStatus((*i).first.c_str());
      if(!(*i).second.GetLeaf()) r=false;
      else cout<<"  "<<(*i).first
               <<" nData="<<(*i).second.GetLeaf()->GetNdata()<<"\n";
   }
   return r;
}


VarData const *VarList::FindVar(string const &name) const {
   const_iterator i=find(name);
   if(i!=end()) {
      return &(*i).second;
   }
   return 0;
}
