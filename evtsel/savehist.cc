#include <TTreeFormula.h>
#include <TFile.h>
#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"
#define __BINS_PTY_PLACEHOLDER__
#include "bins.h"
#include "util.h"
#include "evtntuple.h"
#include "variables.h"

int macro(const std::string& inputname, const std::string& evtcat, const std::string& dcutstr,
          const std::string& outputdir)
{
  bins::print();
  
  // parse evtcat
  const auto evtp = xjjroot::parse_input(evtcat);
  auto want_evt = [&evtp](const std::string& key) -> bool {
    return xjjc::str_contains(evtp.content, key);
  };
  if (!want_evt("gammaN") && !want_evt("Ngamma")) {
    __XJJLOG << "!! neither gammaN nor Ngamma, abort." << std::endl;
    return 2;
  }
  const bool is_gammaN = want_evt("gammaN");
  const float thre_zdc_gamma = is_gammaN ? 1100 : 1000,
    thre_zdc_N = is_gammaN ? 1000 : 1100,
    thre_gap = want_evt("2023") ? (is_gammaN ? 9.2 : 8.7) : 16.;
    
  // parse inputs
  const auto inputp = xjjroot::parse_input(inputname);
  auto* inf = xjjroot::openfile(inputp.content);
  if (!inf) {
    __XJJLOG << "!! bad inputfile, abort." << std::endl;
    return 2;
  }
  auto* tree = static_cast<TTree*>(inf->Get("Tree"));
  if (!tree) {
    __XJJLOG << "!! bad tree, abort." << std::endl;
    xjjroot::closefile(inf);
    return 2;
  }
  const auto is_data = !xjjana::ismc_runnum(tree);
  myvar::evtntuple ent(tree);

  const std::map<std::string, bool> opt_to_doc = {
    { "is_gammaN", is_gammaN },
    { "is_Ngamma", want_evt("Ngamma") },
    { "is_2023", want_evt("2023") },
    { "is_L1ZDCOr", want_evt("L1ZDCOr") },
    { "is_L1ZeroBias", want_evt("L1ZeroBias") },
    { "is_data", is_data },
  };
  xjjc::print_tab(opt_to_doc, 0);

  // parse dcutstr
  const auto dcutp = xjjroot::parse_input(dcutstr);
  xjjc::print_tab(std::map<std::string, std::string>{
      { "cut_d", dcutp.content } }, -1);
  const bool Dnoreq = dcutp.content == "(1)";
  TTreeFormula form_cut_d("form_cut_d", dcutp.content.c_str(), tree);
  if (form_cut_d.GetNdim() <= 0) {
    __XJJLOG << "!! bad formula, abort." << std::endl;
    xjjroot::closefile(inf);
    return 2;
  }
  int Dsize; tree->SetBranchAddress("Dsize", &Dsize);
  std::vector<float>* Dpt = nullptr; tree->SetBranchAddress("Dpt", &Dpt);
  std::vector<float>* Dy = nullptr; tree->SetBranchAddress("Dy", &Dy);
  
  std::map<std::string, TH3D*> h3s;
  for (const auto& the_var : variables) {
    bool have_branches = true;
    for (const auto& br : xjjc::str_divide_trim(the_var.var, ",")) {
      if (!tree->GetBranchStatus(br.c_str())) {
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
  // xjjroot::print_tab(h3s, 0);
  auto* h3_bins = static_cast<TH3D*>(h3s.begin()->second->Clone("h3_bins"));
  
  std::map<std::string, TH1D*> h1s;
  h1s["nd"] = new TH1D("h1_nd", Form(";Number of %s candidates passing cuts / Event;", xjjroot::CMS::DzDzbar.c_str()), 9, 1, 10);
  // const auto nentries = std::min(tree->GetEntries(), static_cast<long long int>(1.e4));
  const auto nentries = tree->GetEntries();
  for (long long int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries);
    tree->GetEntry(i);

    // basic event selections
    if (is_data && !((want_evt("L1ZDCOr") && ent.isL1ZDCOr) ||
                     (want_evt("L1ZeroBias") && ent.isZeroBias))) continue; // trigger
    if (!ent.selectedVtxFilter) continue; // bkgrej 
    if (is_data && !ent.cscTightHalo2015Filter) continue;

    std::vector<std::pair<int, int>> to_fill;
    if (Dnoreq) {
      to_fill.push_back({ 1, 1 });
    } else {
      // .GetNdata() must be called before EvalInstance
      if (Dsize != form_cut_d.GetNdata()) {
        __XJJLOG << "!! bad GetNdata value, abort." << std::endl
                 << "   Dsize: " << Dsize << std::endl
                 << "   form_cut_d.GetNdata(): " << form_cut_d.GetNdata() << std::endl;
        xjjroot::closefile(inf);
        return 2;
      }
      
      for (int j=0; j<Dsize; j++) {
        if (form_cut_d.EvalInstance(j) == 0.) continue;
        // D kinematics
        const std::pair<int, int> ibin = { h3_bins->GetXaxis()->FindBin(Dy->at(j)),
          h3_bins->GetZaxis()->FindBin(Dpt->at(j)) };
        if (ibin.first < 1 || ibin.first > h3_bins->GetXaxis()->GetNbins() ||
            ibin.second < 1 || ibin.second > h3_bins->GetZaxis()->GetNbins())
          continue;

        if (std::find(to_fill.begin(), to_fill.end(), ibin) == to_fill.end())
          to_fill.push_back(ibin);
      }
      if (to_fill.size() > 0)
        h1s["nd"]->Fill(to_fill.size());
    } // else of if (dcut_is(Dsel::noreq))

    const float zdcsum_gamma = is_data ? (is_gammaN ? ent.ZDCsumPlus : ent.ZDCsumMinus) : 0,
      zdcsum_N = is_data ? (is_gammaN ? ent.ZDCsumMinus : ent.ZDCsumPlus) : 4.e3;
    const auto hfemax_eta5_gamma = is_gammaN ? ent.HFEMaxPlus_eta5 : ent.HFEMaxMinus_eta5,
      hfemax_eta5_N = is_gammaN ? ent.HFEMaxMinus_eta5 : ent.HFEMaxPlus_eta5;

    for (const auto& ibin : to_fill) {
      auto fill_value = [&h3s, &ibin](const std::string& name, float value) {
        if (h3s.find(name) != h3s.end()) {
          h3s.at(name)->Fill(h3s.at(name)->GetXaxis()->GetBinCenter(ibin.first), value, h3s.at(name)->GetZaxis()->GetBinCenter(ibin.second));
        }
      };

      fill_value("ZDCsumPlus", ent.ZDCsumPlus);
      fill_value("ZDCsumMinus", ent.ZDCsumMinus);
      fill_value("HFEMaxPlus_eta5", ent.HFEMaxPlus_eta5);
      fill_value("HFEMaxMinus_eta5", ent.HFEMaxMinus_eta5);
      if (zdcsum_gamma < thre_zdc_gamma) { // 0n
        fill_value("ZDCsumN", zdcsum_N);
        fill_value("ZDCsumN-zoom", zdcsum_N);
        if (zdcsum_N > thre_zdc_N) { // 0n + Xn
          fill_value("HFEMaxGamma_eta5", hfemax_eta5_gamma);
          fill_value("HFEMaxN_eta5", hfemax_eta5_N);
          fill_value("nTrackInAcceptanceHP", ent.nTrackInAcceptanceHP);
          if (hfemax_eta5_gamma < thre_gap) { // 0n + Xn + gap
            fill_value("HFEMaxN_eta5-gap", hfemax_eta5_N);
            fill_value("nTrackInAcceptanceHP-gap", ent.nTrackInAcceptanceHP);
            fill_value("nTrackInAcceptanceHP-gap-extend", ent.nTrackInAcceptanceHP);
            fill_value("nTrkVtx", ent.nTrkVtx);
            fill_value("nVtx", ent.nVtx);
          }
          if (hfemax_eta5_gamma >= thre_gap) { // 0n + Xn + no gap
            fill_value("nTrackInAcceptanceHP-nogap", ent.nTrackInAcceptanceHP);
          }
        }        
      }
      // if (hfemax_eta5_gamma < thre_gap) { // gap
      //   if (zdcsum_N > thre_zdc_N) { // gap + Xn
      //     fill_value("ZDCsumGamma", zdcsum_gamma);
      //   }
      // }
    }
  }
  xjjc::progressbar_summary(nentries);

  auto* outf = xjjroot::newfile("rootfiles/" + outputdir + ".root");
  for (auto& [_, h] : h3s)
    xjjroot::writehist(h);
  for (auto& [_, h] : h1s)
    xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);
  util::Writeinfo winfo("info");
  winfo.cast_branch("input", inputp.content);
  winfo.cast_branch("input_tex", inputp.tex);
  winfo.cast_branch("input_tag", inputp.tag);
  winfo.cast_branch("evt", evtp.content);
  winfo.cast_branch("evt_tex", evtp.tex);
  winfo.cast_branch("evt_tag", evtp.tag);
  winfo.cast_branch("dcut", dcutp.content);
  winfo.cast_branch("dcut_tex", dcutp.tex);
  winfo.cast_branch("dcut_tag", dcutp.tag);
  for (auto& [key, val] : opt_to_doc)
    winfo.cast_branch(key, int(val));
  winfo.cast_branch("thre_zdc_gamma", thre_zdc_gamma);
  winfo.cast_branch("thre_zdc_N", thre_zdc_N);
  winfo.cast_branch("thre_gap", thre_gap);
  winfo.close();
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
