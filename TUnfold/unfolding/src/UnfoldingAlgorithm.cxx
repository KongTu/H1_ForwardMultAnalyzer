#include "UnfoldingAlgorithm.h"
#include "UnfoldingResult.h"

using namespace std;

UnfoldingAlgorithm::~UnfoldingAlgorithm() {
   for( map<TString,UnfoldingResult *>::iterator i=fUnfoldingResult.begin();
        i!=fUnfoldingResult.end();i++) {
      delete (*i).second;
   }
}

void UnfoldingAlgorithm::ResetAverages(void) {
   fUnfoldingAverage.clear();
}

map<TString,UnfoldingAverage> const &UnfoldingAlgorithm::GetAverage(void) {
   /* for(map<TString,UnfoldingAverage>::iterator i=fUnfoldingAverage.begin();
       i!=fUnfoldingAverage.end();i++) {
      (*i).second.UpdateCorrelations(fHistTruth);
      } */
   return fUnfoldingAverage;
}

map<TString,UnfoldingResult*> const &UnfoldingAlgorithm::Unfold
(const TH1 *data) {
   // delete all old unfolding results
   for( map<TString,UnfoldingResult *>::iterator i=fUnfoldingResult.begin();
        i!=fUnfoldingResult.end();i++) {
      delete (*i).second;
   }
   fUnfoldingResult.clear();
   RunUnfolding(data);
   for(map<TString,UnfoldingResult*>::const_iterator
          iResult=fUnfoldingResult.begin();
       iResult!=fUnfoldingResult.end();iResult++) {
      fUnfoldingAverage[(*iResult).first].Add((*iResult).second,fHistTruth);
   }
   return fUnfoldingResult;
}

