#ifndef __H__UNFOLDINGBASE__
#define __H__UNFOLDINGBASE__

#include <TObject.h>
#include <TString.h>


//#include <iostream>
//#include <TNamed.h>

class UnfoldingBase {
 protected:
   void DeleteIf(TObject *o) { if(o) delete o; }
   void WriteIf(TObject const *o,TString const &name) const {
      //std::cout<<"write name: "<<name<<"\n";
      if(o) {
         //TNamed const *tn=dynamic_cast<TNamed const *>(o);
         //if(tn) { std::cout<<"named object\n"<<tn->GetName()<<"\n"; }
         o->Write(name);
      }
   }
};

#endif
