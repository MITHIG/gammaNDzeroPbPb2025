namespace globals {
  const std::vector<std::pair<std::string, int>> mask_br = {
    { "DpassCut23*", 0 },
    { "DpassCut23PAS", 1 },
    { "Dtrk*P", 0 },
    { "Dtrk*MassHypo", 0 },
    { "V*", 0 },
    { "clusComp*", 0 },
  };
}

namespace save {
  void mask_branch(TChain* tr) {
    for (const auto& m : globals::mask_br)
      tr->SetBranchStatus(m.first.c_str(), m.second);
  }

  void mask_branch(TTree* tr) {
    for (const auto& m : globals::mask_br)
      tr->SetBranchStatus(m.first.c_str(), m.second);
  }

  std::string cut_adjust_to_mc(const std::string& cut_data) {
    auto cut_r = cut_data;
    for (auto& cut : std::vector<std::string>{ "isL1ZDCOr", "isZeroBias", "cscTightHalo2015Filter" , "ZDCsumPlus < 1100", "ZDCsumMinus < 1000" }) {
      cut_r = xjjc::str_removecut(cut_r, cut);
      // !! regex
      // for (auto& cutstr : { cut+" &&", cut+"&&", "&& "+cut, "&&"+cut }) {
      //   if (xjjc::str_contains(cut_r, cutstr)) {
      //     __XJJLOG << ">> " << cutstr << std::endl;
      //     cut_r = xjjc::str_eraseall(cut_r, cutstr);
      //     break;
      //   }
      // }
    }
    return cut_r;
  }

  template <typename T = int> int tree_is_mc(TTree* nt, const char* brname = "Run") {
    auto* br = nt->GetBranch(brname);
    if (!br) {
      __XJJLOG << "!! bad branch for run number: " << brname << ", abort." << std::endl;
      return -1;
    }
    T Run{};
    nt->SetBranchAddress(brname, &Run);
    if (nt->GetEntry(0) <= 0) {
      __XJJLOG << "!! failed to read entry 0, abort." << std::endl;
      nt->ResetBranchAddress(br);
      return -1;
    }
    auto result = static_cast<int>(Run < 10000);
    nt->ResetBranchAddress(br);
    return result;
  }
}
