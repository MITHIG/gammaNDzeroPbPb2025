#pragma once

#include <memory>

#define private public
#include "dfitter.h"
#undef private

inline TF1* xjjroot::dfitter::f_f(const std::string& name) const {
  if (!fitted_) return nullptr;
  const std::string fname = name.empty() ? Form("%s_total", fun_f_->GetName()) : name;
  return clone_fun(fun_f_, fname);
}

namespace xjjroot {
  struct cached_fit_result { bool fitted = false; double yield = -1; double yieldErr = -1; };

  class cached_dfitter {
  public:
    explicit cached_dfitter(const std::string& option = "") : option_(option) {}
    ~cached_dfitter() { delete total_; }
    bool prepare(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped,
                 const TH1* hmassMCKK, const TH1* hmassMCPiPi) {
      delete total_; total_ = nullptr;
      reference_ = std::make_unique<dfitter>(option_.c_str());
      reference_->fit(hmass, hmassMCSignal, hmassMCSwapped, hmassMCKK, hmassMCPiPi);
      if (!reference_->fitted()) return false;
      nominal_yield_ = reference_->yield(); nominal_yield_err_ = reference_->yieldErr();
      total_ = reference_->f_f(Form("cached_total_%s", xjjc::unique_str().c_str()));
      return total_ != nullptr;
    }
    double nominal_yield() const { return nominal_yield_; }
    double nominal_yield_err() const { return nominal_yield_err_; }
    cached_fit_result fit(TH1* hmass) {
      cached_fit_result result; if (!total_ || !hmass) return result;
      reference_->fit_data_only(hmass);
      if (!reference_->fitted()) return result;
      result.yield = reference_->yield();
      result.yieldErr = reference_->yieldErr();
      result.fitted = true;
      return result;
    }
    TF1* total_model(const std::string& name = "") const { if (!total_) return nullptr; auto* r = new TF1(*total_); if (!name.empty()) r->SetName(name.c_str()); return r; }
    TF1* signal_model(const std::string& name = "") const { return make_match(total_, name); }
  private:
    std::string option_; std::unique_ptr<dfitter> reference_; TF1* total_ = nullptr; double nominal_yield_ = -1, nominal_yield_err_ = -1;
    static TF1* make_match(const TF1* total, const std::string& requested_name = "") {
      if (!total) return nullptr;
      const auto name = requested_name.empty() ? Form("cached_match_%s", xjjc::unique_str().c_str()) : requested_name;
      auto* match = new TF1(name.c_str(), "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))", total->GetXmin(), total->GetXmax());
      const std::map<int,int> params={{0,0},{1,1},{2,2},{3,7},{4,9},{5,10},{6,11},{7,12},{8,13}};
      for (const auto& [target, source] : params) { match->SetParameter(target,total->GetParameter(source)); match->SetParError(target,total->GetParError(source)); }
      return match;
    }
  };
}
