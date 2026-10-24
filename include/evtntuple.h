#pragma once

#define BRANCH_B(EXPAND)                        \
  EXPAND(isL1ZDCOr);                            \
  EXPAND(isZeroBias);                           \
  EXPAND(selectedVtxFilter);                    \
  EXPAND(cscTightHalo2015Filter);               \
  EXPAND(ClusterCompatibilityFilter);           \
  EXPAND(ZDCgammaN);                            \
  EXPAND(ZDCNgamma);
#define BRANCH_F(EXPAND)                        \
  EXPAND(HFEMaxPlus_eta5);                      \
  EXPAND(HFEMaxMinus_eta5);                     \
  EXPAND(ZDCsumPlus);                           \
  EXPAND(ZDCsumMinus);
#define BRANCH_I(EXPAND)                        \
  EXPAND(Run);                                  \
  EXPAND(nVtx);                                 \
  EXPAND(nTrackInAcceptanceHP);                 \
  EXPAND(nTrkVtx);                              \
  EXPAND(Dsize);

#define BRANCH_SINGLE(EXPAND)                   \
  BRANCH_B(EXPAND);                             \
  BRANCH_F(EXPAND);                             \
  BRANCH_I(EXPAND);

namespace myvar
{
  class evtntuple {
  public:
    evtntuple(TTree* t) : t_(t) { reset(); setbranches(); }
    bool valid() { return valid_; }
    // float var_value(const std::string& varname, int j = -1);

#define DEC_B(q) bool q;
    BRANCH_B(DEC_B);
#define DEC_F(q) float q;
    BRANCH_F(DEC_F);
#define DEC_I(q) int q;
    BRANCH_I(DEC_I);

  private:
    bool valid_;
    std::map<std::string, bool> good_branch_;
    TTree* t_;

    void reset();
    void setbranches();
  };
}

void myvar::evtntuple::reset() {
  valid_ = false;
};

void myvar::evtntuple::setbranches() {
  valid_ = false;
  if (!t_) {
    __XJJLOG << "!! invalid tree, abort." << std::endl;
    return;
  }
  // t_->SetBranchStatus("*", 0);

#define SETBRANCH(q)                            \
  if (t_->GetBranch(#q)) {                      \
    t_->SetBranchStatus(#q, 1);                 \
    t_->SetBranchAddress(#q, &q);               \
    good_branch_[#q] = true;                    \
  } else {                                      \
    std::cout<<"(x) " #q<<std::endl;            \
    good_branch_[#q] = false;                   \
  }
  BRANCH_SINGLE(SETBRANCH);
  valid_ = true;
}
