#ifdef PRIVATE_TUNFOLD
#include "../TUnfold_V17.7/TUnfoldBinningXML.h"
#else
#include <TUnfoldBinningXML.h>
#endif

#include <string>
#include <vector>

class VarList;

class ClassifierBinning {
  // this helper class is to define "categories" as bins in a given variables
  //
  // in the XML file, the syntax is:
  //   <BinningNode name="[lowEdge:]var:upEdge">
  // where
  //   var : a variable in the minitree
  // upEdge : a number
  // [lowEdge : another number ]
  //
  // this will define a new bin in the range
  //    [lowEdge,upEdge]
  // if lowEdge is not specified, then:
  //   for the first bin, lowEdge defaults to -infinity
  //   for the other bins, previous defaults to upEdge of the preceeding bin
  //
  //  example:
  //    <BinningNode name="cosDPhi_reco:-0.65">
  //      ...
  //    </BinningNode>
  //    <BinningNode name="cosDPhi_reco:1.0">
  //      ...
  //    </BinningNode>
  //
  // this defines two categories
  //   the first category contains events with
  //       cosDPhi_reco<-0.65
  //   the second category contains events with
  //       -0.65<=cosDPhi_reco<1
  //
  // within the categories, define
  //   - further sub-categories (using <BinningNode> )
  //   - N-dimensional histograms (using <Axis> and <Bin> )
  //
public:

  // construct a new "ClassifierBinning"
  // input:  binning schene (from parsing XML file)
  // output: vars (all classifier variables will be added to "vars")
   ClassifierBinning(TUnfoldBinning const *binning,VarList &vars);
   virtual ~ClassifierBinning();
   TString GetName(void) const;
   //std::string const &GetVar(void) const { return fName; }
   //size_t GetNbin(void) const { return fBins.size(); }
   //double GetXmin(size_t k) const { return fBins[k].fXmin; }
   //double GetXmax(size_t k) const { return fBins[k].fXmax; }
   ClassifierBinning const *GetClassifier(size_t k) { return fBins[k].fClassifier; }
   bool IsInside(double x,size_t k) const {
      return (x>=fBins[k].fXmin)&&(x<fBins[k].fXmax); }

   std::vector<ClassifierBinning const *> Enumerate(void) const;
  //
  // locate appropriate BinningNode (unconnected bins,N-dimensional histogram)
  // using classifier information
   ClassifierBinning const *FindNode(VarList const &vars,int print=0) const;

   TUnfoldBinning const *GetBinningNodeClassifier(void) const { return fBinningClassifier; }
   TUnfoldBinning const *GetBinningNodeDistribution(void) const { return fBinningDistribution; }
   TUnfoldBinning const *GetBinningNodeUnconnected(void) const { return fBinningUnconnected; }
protected:
   ClassifierBinning(TUnfoldBinning const *binning,VarList &vars,
                     ClassifierBinning const *parent,int indentLevel);
   void Enumerate_r(std::vector<ClassifierBinning const *> &clist) const;
   

   ClassifierBinning const *fParent;
   std::string fName;
   struct ClassifierBin {
      double fXmin,fXmax;
      ClassifierBinning *fClassifier;
   };
   std::vector<ClassifierBin> fBins;
   TUnfoldBinning const *fBinningClassifier;
   TUnfoldBinning const *fBinningDistribution;
   TUnfoldBinning const *fBinningUnconnected;
   void AddSubBins(TUnfoldBinning const *binning,
                   VarList &allVars,int indentLevel);
};
