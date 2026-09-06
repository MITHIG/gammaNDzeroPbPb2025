#include <TTreeFormula.h>
#include <TH2D.h>
#include "xjjanauti.h"

#include "../include/save.h"
#define __BINS_PTY_ANA__
#define __BINS_MASS__
#include "../include/bins.h"
#include "../include/util.h"

// #include <cstdlib>

int macro(const std::string& inputstr, const std::string& cutstr_evt, const std::string& cutstr_d,
          const std::string& output) {
  __XJJLOG << ">> current y binning:" << std::endl;
  xjjc::print_vec_h(bins::ybins, 0);
  __XJJLOG << ">> current pt binning:" << std::endl;
  xjjc::print_vec_h(bins::ptbins, 0);
  
  const auto pinput = util::parse_input(inputstr);
  auto* inf = TFile::Open(pinput.content.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! failed to open input file, abort." << std::endl;
    return 1;
  }
  auto* tree = dynamic_cast<TTree*>(inf->Get("Tree"));
  if (!tree) {
    __XJJLOG << "!! failed to get TTree named, abort." << std::endl;
    xjjroot::closefile(inf);
    return 1;
  }
  const auto is_mc = save::tree_is_mc(tree);
  if (is_mc < 0) {
    __XJJLOG << "!! cannot judge if this is MC, abort." << std::endl;
    xjjroot::closefile(inf);
    return 1;
  }

  // branches
  int Dsize; tree->SetBranchAddress("Dsize", &Dsize);
  std::vector<float>* Dpt = nullptr; tree->SetBranchAddress("Dpt", &Dpt);
  std::vector<float>* Dy = nullptr; tree->SetBranchAddress("Dy", &Dy);
  // cuts
  auto pcut_evt = util::parse_input(cutstr_evt);
  if (is_mc) pcut_evt.content = save::cut_adjust_to_mc(pcut_evt.content);
  auto pcut_d = util::parse_input(cutstr_d);
  xjjc::print_tab(std::map<std::string, std::string>{
      { "cut_evt", pcut_evt.content },
      { "cut_d", pcut_d.content } }, -1);
  TTreeFormula form_cut_d("form_cut_d", pcut_d.content.c_str(), tree);
  TTreeFormula form_cut_evt("form_cut_evt", pcut_evt.content.c_str(), tree);
  auto valid_formula = [](const TTreeFormula& formula) -> bool {
    return formula.GetNdim() > 0;
  };
  if (!valid_formula(form_cut_d) || !valid_formula(form_cut_evt)) {
    __XJJLOG << "!! failed to build one or more TTreeFormula, abort." << std::endl
             << "   cut_d: " << pcut_d.content << std::endl 
             << "   cut_evt: " << pcut_evt.content << std::endl;
    xjjroot::closefile(inf);
    return 1;
  }

  std::map<std::string, TH2D*> h2s;
  auto create_hist = [&h2s](const std::string& name) {
    auto* h = new TH2D(Form("h2_%s", name.c_str()), ";y;#it{p}_{T}",
                       bins::ybins.size()-1, bins::ybins.data(),
                       bins::ptbins.size()-1, bins::ptbins.data());
    h->Sumw2();
    h2s[name] = h;
  };
  create_hist("evteff_num");
  create_hist("evteff_den");
  create_hist("bins");

  const auto ny = h2s["bins"]->GetNbinsX(),
    npt = h2s["bins"]->GetNbinsY();
  
  std::vector<size_t> bin_filled(ny * npt, 0);
  const auto nentries = tree->GetEntries();
  for (long long i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    xjjc::progressslide(i, nentries);
    // if (tree->LoadTree(i) < 0) break;

    if (Dsize != form_cut_d.GetNdata() ||
        form_cut_evt.GetNdata() != 1) {
      __XJJLOG << "!! bad GetNdata value, abort." << std::endl
               << "   Dsize: " << Dsize << std::endl
               << "   form_cut_d.GetNdata(): " << form_cut_d.GetNdata() << std::endl
               << "   form_cut_evt.GetNdata()" << form_cut_evt.GetNdata() << std::endl;
      xjjroot::closefile(inf);
      return 3;
    }

    std::fill(bin_filled.begin(), bin_filled.end(), 0);
    for (int j = 0; j < Dsize; ++j) {
      if (form_cut_d.EvalInstance(j) == 0.) continue;

      const auto ibin_y = h2s["bins"]->GetXaxis()->FindBin(Dy->at(j)),
        ibin_pt = h2s["bins"]->GetYaxis()->FindBin(Dpt->at(j));
      if (ibin_y < 1 || ibin_y > ny || ibin_pt < 1 || ibin_pt > npt)
        continue;

      const size_t index_flat = (ibin_pt - 1) * ny + (ibin_y - 1);
      if (bin_filled[index_flat]) continue;

      h2s.at("evteff_den")->Fill(Dy->at(j), Dpt->at(j));
      if (form_cut_evt.EvalInstance() > 0)
        h2s.at("evteff_num")->Fill(Dy->at(j), Dpt->at(j));
      bin_filled[index_flat] = 1;
    }
  }
  xjjc::progressbar_summary(nentries);

  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  for (auto& [_, h] : h2s)
    xjjroot::writehist(h);
  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_cont;
  auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
    t_cont[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_cont[name]));
  };
  cast_branch("input", pinput.content);
  cast_branch("input_tex", pinput.tex);
  cast_branch("input_tag", pinput.tag);
  cast_branch("cut_d", pcut_d.content);
  cast_branch("cut_d_tex", pcut_d.tex);
  cast_branch("cut_d_tag", pcut_d.tag);
  cast_branch("cut_evt", pcut_evt.content);
  cast_branch("cut_evt_tex", pcut_evt.tex);
  cast_branch("cut_evt_tag", pcut_evt.tag);
  cast_branch("is_mc", is_mc);
  t->Fill();
  t->Write();
  xjjroot::closefile(outf);

  return 0;
}

int main (int argc, char* argv[]) {
  if (argc == 7) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[5], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[6], ",");
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  return 1;
}
