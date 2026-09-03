#include <TFile.h>
#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "util.h"
#include "draw.h"

#include "variables.h"

int macro(const std::vector<std::string>& inputnames, const std::string& outputname) {
  TH3D* h3_bins = nullptr;
  std::vector<xjjc::info> infos;
  std::map<std::string, xjjc::array2D<std::vector<TH1D*>>> h1ptys;

  for (int k=0; k<inputnames.size(); k++) {
    const auto inputp = util::parse_input(inputnames[k]);
    auto* inf = TFile::Open(inputp.content.c_str());
    if (!h3_bins) h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins");
    const auto ny = h3_bins->GetXaxis()->GetNbins(), npt = h3_bins->GetZaxis()->GetNbins();

    auto h1s = xjjana::getobj_regexp<TH1D>(inf, "h1_.+_pt-.+_refy-.+");
    for (auto* h : h1s) {
      const auto name = ;
      if (h1ptys.find(name) == h1ptys.end()) {
        h1ptys[name] = xjjc::array2d<std::vector<TH1D*>>(npt, ny);
        for (auto& hhh : h1ptys[name])
          for (auto& hh : hhh)
            hh.resize(inputnames.size(), nullptr);
      }
      const auto i = xjjc::str_extract_index(h->GetName(), "_pt-"),
        j = xjjc::str_extract_index(h->GetName(), "_refy-");
      h1ptys[name][i][j][k] = h;
    }
    auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
    __XJJLOG << "++ info" << std::endl;
    xjjc::print_tab(info, -1);
    infos.push_back(info);
    // const bool has_dreq = !xjjc::str_contains(info.at("dcut"), "Dnoreq");
  }
  // auto this_is = [&info](const std::string& item) {
  //   return static_cast<bool>(std::atoi(info.at(item).c_str()));
  // };
  auto all_are = [&infos](const std::string& item) {
    bool r = true;
    for (const auto& info : infos) {
      if (info.find(item) == info.end()) {
        __XJJLOG << "!! info does not have key: " << item << ", skip." << std::endl;
        r = false;
        continue;
      }
      r = r && static_cast<bool>(std::atoi(info.at(item).c_str()));
    }
    return r;
  };
  auto all_are_same = [&infos]<typename T>(const std::string& item) {
    T r = ;
    if (info.find(item) == info.end()) {
      __XJJLOG << "!! info does not have key: " << item << ", skip." << std::endl;
      continue;
      if ()
    }
  };
  
  const draw::bintex tbins(h3_bins, 0, 2);
  const auto ny = h3_bins->GetXaxis()->GetNbins(), npt = h3_bins->GetZaxis()->GetNbins();
  const auto colors = xjjroot::grayscales_color(ny, xjjroot::mycolor_dark["red"], 1, 0.3);
  std::vector<std::string> ts_common = {
    // std::string(this_is("l1_ZDCOr") ? "L1_ZDCOr + " : (this_is("l1_ZeroBias") ? "L1_ZeroBias + " : "")) + "PV filter + cscTightHalo",
    std::string(all_are("l1_zdcor") ? "L1_ZDCOr + " : (all_are("l1_zerobias") ? "L1_ZeroBias + " : "")) + "PV filter + cscTightHalo",
    (all_are("dir_gammaN") ? "#gamma #rightarrow Z+ (#gammaN)" : "#gamma #rightarrow Z- (N#gamma)"),
    // info["evt_tex"],
  };
}
