#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"

#include "../include/save.h"
#define __BINS_PTY_EFF__
#define __BINS_MULT__
#include "../include/bins.h"

const std::vector<std::string> pnp_labels = { "inclusive", "prompt", "nonprompt" };

int macro(const std::string& inputmcstr, const std::string& cutevtstr, const std::string& cutdstr, const std::string& output, int pnp_flag = 0,
          const std::string& weight_bdt = "") {
  std::map<std::string, TChain*> trs;
  // parse inputmc
  const auto pi_inputmc = xjjroot::parse_input(inputmcstr);
  trs["mc"] = xjjana::chain_files(xjjc::str_divide_trim(pi_inputmc.content, ","), "Tree");
  if (!trs.at("mc")) {
    __XJJLOG << "!! bad inputmc file " << pi_inputmc.content << ", abort." << std::endl;
    return 2;
  }
  save::mask_branch(trs.at("mc"));

  // parse cut
  const auto pi_cutevt = xjjroot::parse_input(cutevtstr);
  const auto cutevt = pi_cutevt.content, cutevt_mc = save::cut_adjust_to_mc(cutevt);
  const auto pi_cutd = xjjroot::parse_input(cutdstr);
  const auto cutd = pi_cutd.content;
  const auto apply_weight = !weight_bdt.empty();
  if (apply_weight)
    __XJJLOG << ">> apply BDT weight: " << weight_bdt << std::endl;
  
  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");

  std::map<std::string, TH3D*> h3;
  auto project = [&h3](TChain* tr, std::string key, std::string vars, std::string icut) {
    h3[key] = new TH3D(Form("h3%s", key.c_str()), ";y;p_{T} [GeV];N_{trk} (highPurity, p_{T} > 0.5 GeV, |#eta| < 2.4)",
                       bins::ny, bins::miny, bins::maxy,
                       bins::npt, bins::minpt, bins::maxpt,
                       bins::nmult, bins::minmult, bins::maxmult
                       );
    __XJJLOG << ">> "<<h3[key]->GetName()<<" ("<<vars<<") \e[2m"<<icut<<"\e[0m"<<std::endl;
    xjjc::saywait();
    tr->Project(h3[key]->GetName(), vars.c_str(), icut.c_str());
    xjjroot::writehist(h3[key]);
    return icut;
  };

  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_cont;
  auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
    t_cont[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_cont[name]));
  };

  __XJJLOG << ">> " << pnp_labels.at(pnp_flag) << std::endl;
  std::string cut_signal_reco = "Dgen==23333", cut_signal_gen = "GisSignalCalc";
  if (pnp_flag == 1) {
    cut_signal_reco += " && DisSignalCalcPrompt";
    cut_signal_gen += " && GisSignalCalcPrompt";
  } else if (pnp_flag == 2) {
    cut_signal_reco += " && DisSignalCalcFeeddown";
    cut_signal_gen += " && GisSignalCalcFeeddown";    
  }

  auto cut_eff_num = project(trs.at("mc"), "_eff_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && " + cut_signal_reco + " && " + cutd);
  cast_branch("cut_eff_num", cut_eff_num);
  auto cut_reco_num = project(trs.at("mc"), "_reco_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && " + cut_signal_reco + " && fabs(Dtrk1Eta) < 2.4 && fabs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5");
  cast_branch("cut_reco_num", cut_reco_num);
  auto cut_acc_num = project(trs.at("mc"), "_acc_num", "nTrackInAcceptanceHP:Gpt:Gy", cutevt_mc + " && " + cut_signal_gen + " && fabs(Gtk1eta) < 2.4 && fabs(Gtk2eta) < 2.4 && Gtk1pt > 0.5 && Gtk2pt > 0.5");
  cast_branch("cut_acc_num", cut_acc_num);
  auto cut_eff_den = project(trs.at("mc"), "_eff_den", "nTrackInAcceptanceHP:Gpt:Gy", cutevt_mc + " && " + cut_signal_gen);
  cast_branch("cut_eff_den", cut_eff_den);

  if (apply_weight) {
    auto cut_eff_num_w = project(trs.at("mc"), "_eff_num_w", "nTrackInAcceptanceHP:Dpt:Dy", "(" + weight_bdt + ")*(" + cut_eff_num + ")");
    cast_branch("cut_eff_num_w", cut_eff_num_w);
    // selection without the BDT cut, unweighted and weighted: eff_calc normalizes the weight with
    // presel_num / presel_num_w in each final bin, so that the weight only changes the BDT cut efficiency
    const auto cutd_nobdt = std::regex_replace(cutd, std::regex(R"(Dmva_BDT\s*>\s*[-+0-9.eE]+)"), "1");
    auto cut_presel_num = project(trs.at("mc"), "_presel_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && " + cut_signal_reco + " && " + cutd_nobdt);
    cast_branch("cut_presel_num", cut_presel_num);
    auto cut_presel_num_w = project(trs.at("mc"), "_presel_num_w", "nTrackInAcceptanceHP:Dpt:Dy", "(" + weight_bdt + ")*(" + cut_presel_num + ")");
    cast_branch("cut_presel_num_w", cut_presel_num_w);

    cast_branch("weight_bdt", weight_bdt);
  }

  cast_branch("inputmc", pi_inputmc.content);
  cast_branch("inputmc_tex", pi_inputmc.tex);
  cast_branch("cutevt", cutevt);
  cast_branch("cutevt_mc", cutevt_mc);
  cast_branch("cutevt_tex", pi_cutevt.tex);
  cast_branch("cutd", cutd);
  cast_branch("cutd_tex", pi_cutd.tex);
  cast_branch("pnp_flag", pnp_flag);

  t->Fill();
  t->Write();
  xjjroot::closefile(outf);
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    return macro(argv[1], argv[2], argv[3], argv[4], std::atoi(argv[5]), argv[6]);
  }
  if (argc == 6) {
    return macro(argv[1], argv[2], argv[3], argv[4], std::atoi(argv[5]));
  }
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  return 1;
}

