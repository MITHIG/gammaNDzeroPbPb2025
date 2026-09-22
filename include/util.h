#pragma once

namespace util {
  xjjc::info read_info(TDirectory* inf, const std::string& treename = "info", bool verbose = false) {
    auto info = xjjana::getval_regexp(static_cast<TTree*>(inf->Get(treename.c_str())));
    if (verbose) {
      __XJJLOG << "++ info" << std::endl;
      xjjc::print_tab(info, -1);
    }
    return info;
  }
  
  class Writeinfo {
  public:
    Writeinfo() : t_(nullptr) { }
    Writeinfo(const std::string& treename) : t_(nullptr) { init(treename); }
    TTree* init(const std::string& treename = "info") { t_ = new TTree(treename.c_str(), ""); return t_; }
    template<typename T> void cast_branch(const std::string& name, const T& x) {
      t_cont_[name] = xjjc::to_string(x);
      t_->Branch(name.c_str(), &(t_cont_[name]));
    }
    void close() { t_->Fill(); t_->Write(); }
  private:
    TTree* t_;
    std::map<std::string, std::string> t_cont_;
  };

  int mirrorswap_hist(TH1 *h) {
    const auto nbins = h->GetNbinsX();
    if (nbins%2 != 0) {
      __XJJLOG << "!! the histogram can't be mirror swapped, abort." << std::endl;
      return 1;
    }
    for (int i = 1; i <= nbins / 2; ++i) {
      int j = nbins + 1 - i;
      const auto content_i = h->GetBinContent(i);
      const auto error_i   = h->GetBinError(i);
      const auto content_j = h->GetBinContent(j);
      const auto error_j   = h->GetBinError(j);
      h->SetBinContent(i, content_j);
      h->SetBinError(i, error_j);
      h->SetBinContent(j, content_i);
      h->SetBinError(j, error_i);
    }
    return 0;
  }
}

/** solution to use lambda
    auto* t = new TTree("info", "");
    std::map<std::string, std::string> t_cont;
    auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
    t_cont[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_cont[name]));
    };
    cast_branch("input", pinput.content);
    t->Fill();
    t->Write();

    auto* t = new TTree("info", "");
    for (auto& [key, content] : info) {
    t->Branch(key.c_str(), &content);
    }
    t->Fill();
    t->Write();
    xjjroot::closefile(outf);
**/

