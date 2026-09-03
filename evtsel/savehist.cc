#include <TFile.h>
#include <TH3D.h>
#include "xjjanauti.h"
#define __BINS_PTY_PLACEHOLDER__
#include "../include/bins.h"
#include "util.h"
#include "dntuple.h"
#include "variables.h"

enum class Dsel {
  noreq,
  pascut, pascut_sig,
  yrefbdt, yrefbdt_sig,
  yinclbdt, yinclbdt_sig,
};
const std::unordered_map<std::string, Dsel> dsel_maps = {
  { "Dnoreq", Dsel::noreq },
  { "Dpascut", Dsel::pascut },
  { "Dpascut_sig", Dsel::pascut_sig },
  { "Dyrefbdt", Dsel::yrefbdt },
  { "Dyrefbdt_sig", Dsel::yrefbdt_sig },
  { "Dyinclbdt", Dsel::yinclbdt },
  { "Dyinclbdt_sign", Dsel::yinclbdt_sig },
};
bool dcut_sigreg(const Dsel& cut) { return cut == Dsel::pascut_sig || cut == Dsel::yinclbdt_sig || cut == Dsel::yrefbdt_sig; }

int macro(const std::string& inputname, const std::string& evtcat, const std::string& dcutcat,
          const std::string& outputdir)
{
  __XJJLOG << ">> current y binning:" << std::endl;
  xjjc::print_vec_h(bins::ybins, 0);
  __XJJLOG << ">> current pt binning:" << std::endl;
  xjjc::print_vec_h(bins::ptbins, 0);
  
  // parse evtcat
  const auto evtp = util::parse_input(evtcat);
  auto evt_is = [&evtp](const std::string& key) {
    return xjjc::str_contains(evtp.content, key);
  };
  std::map<std::string, bool> want = {
    { "dir_gammaN", evt_is("gammaN") },
    { "dir_Ngamma", evt_is("Ngamma") },
    { "l1_ZDCOr", evt_is("ZDCOr") },
    { "l1_ZeroBias", evt_is("ZeroBias") },
  };

  if (!want.at("dir_gammaN") && !want.at("dir_Ngamma")) {
    __XJJLOG << "!! neither gammaN nor Ngamma, abort." << std::endl;
    return 2;
  }
  const float thre_zdc_gamma = want.at("dir_gammaN") ? 1100 : 1000,
    thre_zdc_N = want.at("dir_gammaN") ? 1000 : 1100,
    thre_gap = evt_is("2023") ? (want.at("dir_gammaN") ? 9.2 : 8.7) : 16.;

  // parse dcutcat
  const auto dcutp = util::parse_input(dcutcat);
  const auto it_dcut = dsel_maps.find(dcutp.content);
  if (it_dcut == dsel_maps.end()) {
    __XJJLOG << "!! bad dsel category: " << dcutp.content << ", abort." << std::endl;
    return 2;
  }
  auto dcut_is = [&it_dcut](const Dsel& dcutref) {
    return it_dcut->second == dcutref;
  };
  want.emplace("Dnoreq", dcut_is(Dsel::noreq));
  const float signal_region_min = dcut_sigreg(it_dcut->second) ? 1.83 : -1,
    signal_region_max = dcut_sigreg(it_dcut->second) ? 1.9 : -1;

  xjjc::print_tab(want, 0);

  // parse inputs
  const auto inputp = util::parse_input(inputname);
  auto* inf = TFile::Open(inputp.content.c_str());
  auto* nt = static_cast<TTree*>(inf->Get("Tree"));
  if (!nt) {
    __XJJLOG << "!! bad input file: " << inputp.content << std::endl;
    return 2;
  }
  myvar::dntuple dnt(nt);
  nt->GetEntry(0);
  const bool is_data = dnt.Run > 10;
  __XJJLOG << ">> This sample is \e[1m" << (is_data ? "data" : "MC") << "\e[0m" << std::endl;
  
  std::map<std::string, TH3D*> h3s;
  for (const auto& the_var : variables) {
    bool have_branches = true;
    for (const auto& br : xjjc::str_divide_trim(the_var.var, ",")) {
      if (!nt->GetBranchStatus(br.c_str())) {
        have_branches = false;
        break;
      }
    }
    if (!have_branches) continue;
    const auto varbins = the_var.bins.empty() ? xjjc::fixedbin_to_edges(the_var.nbin, the_var.varmin, the_var.varmax) : the_var.bins;
    h3s[the_var.varname] = new TH3D(Form("h3_y_%s_pt", the_var.varname.c_str()), Form(";y;%s;p_{T} [GeV]", the_var.vartex.c_str()),
                                    bins::ybins.size()-1, bins::ybins.data(),
                                    varbins.size()-1, varbins.data(),
                                    bins::ptbins.size()-1, bins::ptbins.data());
  }
  xjjroot::print_tab(h3s, 0);
  auto* h3_bins = static_cast<TH3D*>(h3s.begin()->second->Clone("h3_bins"));
  
  std::map<std::string, TH1D*> h1s;
  h1s["nd"] = new TH1D("h1_nd", Form(";Number of %s candidates passing cuts / Event;", xjjroot::CMS::DzDzbar.c_str()), 9, 1, 10);
  // const auto nentries = std::min(nt->GetEntries(), static_cast<long long int>(1.e4));
  const auto nentries = nt->GetEntries();
  for (long long int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries);
    nt->GetEntry(i);

    // basic event selections
    if (is_data && !((want.at("l1_ZDCOr") && dnt.isL1ZDCOr) ||
                     (want.at("l1_ZeroBias") && dnt.isZeroBias))) continue; // trigger
    if (!dnt.selectedVtxFilter) continue; // bkgrej 
    if (is_data && !dnt.cscTightHalo2015Filter) continue;

    std::vector<std::pair<int, int>> to_fill;
    if (want.at("Dnoreq")) {
      to_fill.push_back({ 1, 1 });
    } else {
      for (int j=0; j<dnt.Dsize; j++) {
        // D pre selections
        if (!(TMath::Abs(dnt.Dtrk1PtErr->at(j)/dnt.Dtrk1Pt->at(j))<0.1 && TMath::Abs(dnt.Dtrk2PtErr->at(j)/dnt.Dtrk2Pt->at(j))<0.1 &&
              TMath::Abs(dnt.Dtrk1Eta->at(j)) < 2.4 && TMath::Abs(dnt.Dtrk2Eta->at(j)) < 2.4 &&
              dnt.Dtrk1Pt->at(j) > 0.5 && dnt.Dtrk2Pt->at(j) > 0.5 &&
              dnt.Dchi2cl->at(j) > 0.05 && (dnt.DsvpvDistance->at(j)/dnt.DsvpvDisErr->at(j)) > 1.))
          continue;
        // D topo selections
        if (dcut_is(Dsel::yrefbdt) || dcut_is(Dsel::yrefbdt_sig)) { //
          if (want.at("dir_gammaN") && !((dnt.Dy->at(j)<-1 && dnt.Dmva_BDT->at(j)>0.143) || (dnt.Dy->at(j)>=-1 && dnt.Dy->at(j)<0 && dnt.Dmva_BDT->at(j)>0.142) || (dnt.Dy->at(j)>=0 && dnt.Dy->at(j)<1 && dnt.Dmva_BDT->at(j)>0.123) || (dnt.Dy->at(j)>=1 && dnt.Dmva_BDT->at(j)>0.098))) continue;
          if (want.at("dir_Ngamma") && !((dnt.Dy->at(j)>=1 && dnt.Dmva_BDT->at(j)>0.143) || (dnt.Dy->at(j)<1 && dnt.Dy->at(j)>=0 && dnt.Dmva_BDT->at(j)>0.142) || (dnt.Dy->at(j)<0 && dnt.Dy->at(j)>=-1 && dnt.Dmva_BDT->at(j)>0.123) || (dnt.Dy->at(j)<-1 && dnt.Dmva_BDT->at(j)>0.098))) continue;
        } else if (dcut_is(Dsel::pascut) || dcut_is(Dsel::pascut_sig)) { //
          if (!dnt.DpassCut23PAS->at(j)) continue;
        } else if (dcut_is(Dsel::yinclbdt) || dcut_is(Dsel::yinclbdt_sig)) { //
          if (!(dnt.Dmva_BDT->at(j) > 0.12)) continue;
        }
        if (signal_region_min > 0 && !(dnt.Dmass->at(j) > signal_region_min && dnt.Dmass->at(j) < signal_region_max)) {
          continue;
        }
        // D kinematics
        const std::pair<int, int> ibin = { h3_bins->GetXaxis()->FindBin(dnt.Dy->at(j)),
          h3_bins->GetZaxis()->FindBin(dnt.Dpt->at(j)) };
        if (ibin.first < 1 || ibin.first > h3_bins->GetXaxis()->GetNbins() ||
            ibin.second < 1 || ibin.second > h3_bins->GetZaxis()->GetNbins())
          continue;

        if (std::find(to_fill.begin(), to_fill.end(), ibin) == to_fill.end())
          to_fill.push_back(ibin);
      }
      if (to_fill.size() > 0)
        h1s["nd"]->Fill(to_fill.size());
    } // else of if (dcut_is(Dsel::noreq))

    const float zdcsum_gamma = is_data ? (want.at("dir_gammaN") ? dnt.ZDCsumPlus : dnt.ZDCsumMinus) : 0,
      zdcsum_N = is_data ? (want.at("dir_gammaN") ? dnt.ZDCsumMinus : dnt.ZDCsumPlus) : 4.e3;
    const auto hfemax_eta5_gamma = want.at("dir_gammaN") ? dnt.HFEMaxPlus_eta5 : dnt.HFEMaxMinus_eta5,
      hfemax_eta5_N = want.at("dir_gammaN") ? dnt.HFEMaxMinus_eta5 : dnt.HFEMaxPlus_eta5;
    // const bool is_ZDC_0nXn = (want.at("dir_gammaN") && dnt.ZDCgammaN) || (want.at("dir_Ngamma") && dnt.ZDCNgamma);

    for (const auto& ibin : to_fill) {
      auto fill_value = [&h3s, &ibin](const std::string& name, float value) {
        if (h3s.find(name) != h3s.end()) {
          h3s.at(name)->Fill(h3s.at(name)->GetXaxis()->GetBinCenter(ibin.first), value, h3s.at(name)->GetZaxis()->GetBinCenter(ibin.second));
        }
      };

      fill_value("ZDCsumPlus", dnt.ZDCsumPlus);
      fill_value("ZDCsumMinus", dnt.ZDCsumMinus);
      fill_value("HFEMaxPlus_eta5", dnt.HFEMaxPlus_eta5);
      fill_value("HFEMaxMinus_eta5", dnt.HFEMaxMinus_eta5);
      if (hfemax_eta5_gamma < thre_gap) { // gap
        fill_value("ZDCsumGamma_noxn", zdcsum_gamma);
        if (zdcsum_N > thre_zdc_N) { // gap + Xn
          fill_value("ZDCsumGamma", zdcsum_gamma);
          if (zdcsum_gamma < thre_zdc_gamma) { // gap + Xn + 0n
            fill_value("HFEMaxN_eta5_gap", hfemax_eta5_N);
            fill_value("nTrackInAcceptanceHP", dnt.nTrackInAcceptanceHP);
            fill_value("nTrkVtx", dnt.nTrkVtx);
            fill_value("nVtx", dnt.nVtx);
          }
        }
      } // if (hfemax_eta5_gamma < thre_gap) {
      if (zdcsum_gamma < thre_zdc_gamma) { // 0n
        fill_value("ZDCsumN", zdcsum_N);
        fill_value("ZDCsumN_zoom", zdcsum_N);
        if (zdcsum_N > thre_zdc_N) { // 0n + Xn
          fill_value("HFEMaxGamma_eta5", hfemax_eta5_gamma);
          fill_value("HFEMaxN_eta5", hfemax_eta5_N);
          fill_value("nTrackInAcceptanceHP_nogap", dnt.nTrackInAcceptanceHP);
        }
      }
    }
  }
  xjjc::progressbar_summary(nentries);

  auto* outf = xjjroot::newfile("rootfiles/" + outputdir + ".root");
  for (auto& [_, h] : h3s)
    xjjroot::writehist(h);
  for (auto& [_, h] : h1s)
    xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);
  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_container; // for lifetime of the intermediate string
  auto cast_branch = [&t, &t_container]<typename T>(const std::string& name, const T& x) {
    t_container[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_container[name]));
  };
  cast_branch("input", inputp.content);
  cast_branch("input_tex", inputp.tex);
  cast_branch("input_tag", inputp.tag);
  cast_branch("evt", evtp.content);
  cast_branch("evt_tex", evtp.tex);
  cast_branch("evt_tag", evtp.tag);
  cast_branch("dcut", dcutp.content);
  cast_branch("dcut_tex", dcutp.tex);
  cast_branch("dcut_tag", dcutp.tag);
  cast_branch("is_data", int(is_data));
  cast_branch("thre_zdc_gamma", thre_zdc_gamma);
  cast_branch("thre_zdc_N", thre_zdc_N);
  cast_branch("thre_gap", thre_gap);
  for (auto& [key, val] : want) 
    cast_branch(key, int(val));
  cast_branch("signal_region_min", signal_region_min);
  cast_branch("signal_region_max", signal_region_max);
  t->Fill();
  t->Write();
  outf->cd();
  xjjroot::closefile(outf);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc < 5)
    return 1;
  
  bins::ybins = xjjc::str_convert_vector<double>("-2.4, 2.4", ",");
  bins::ptbins = xjjc::str_convert_vector<double>("0, 100", ",");
  if (!xjjc::str_contains(argv[3], "Dnoreq")) {
    if (argc > 5) {
      bins::ybins = xjjc::str_convert_vector<double>(argv[5], ",");
      if (argc > 6)
        bins::ptbins = xjjc::str_convert_vector<double>(argv[6], ",");
    }
  }

  return macro(argv[1], argv[2], argv[3], argv[4]);
}
