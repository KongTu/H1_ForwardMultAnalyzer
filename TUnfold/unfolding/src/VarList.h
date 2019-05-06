#include <map>
#include <string>
#include <map>

#include <Rtypes.h>

class TTree;
class TLeaf;

class VarData {
 public:
   VarData();
   double Double(int i=0) const;
   int Int(int i=0) const;
   void SetLeaf(TLeaf *l);
   TLeaf *GetLeaf(void) const { return fLeaf; }
 protected:
   TLeaf *fLeaf;
};

class VarList : public std::map<std::string,VarData > {
 public:
   // manage a list of variables which are read from a TTree
   // or calculated from other variables
   void AddVar(std::string const &name);
   VarData const *FindVar(std::string const &name) const;
   bool SetBranchAddress(TTree *tree);
};

