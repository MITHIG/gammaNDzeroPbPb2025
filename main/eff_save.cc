#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"

#include "../include/save.h"
#define __BINS_PTY_EFF__
#define __BINS_MULT__
#include "../include/bins.h"

int macro(std::string inputmcstr, std::string cutevtstr, std::string cutdstr, std::string output, std::string inputdatastr = "null") {
  std::map<std::string, TChain*> trs;
  // parse inputmc
  const auto pi_inputmc = xjjroot::parse_input(inputmcstr);
  trs["mc"] = xjjana::chain_files(xjjc::str_divide_trim(pi_inputmc.content, ","), "Tree");
  if (!trs.at("mc")) {
    __XJJLOG << "!! bad inputmc file " << pi_inputmc.content << ", abort." << std::endl;
    return 2;
  }
  save::mask_branch(trs.at("mc"));

  const auto pi_inputdata = xjjroot::parse_input(inputdatastr);
  trs["data"] = xjjana::chain_files(xjjc::str_divide_trim(pi_inputdata.content, ","), "Tree");
  if (!trs.at("data")) {
    __XJJLOG << "?? no inputdata file, only MC is used." << std::endl;
  }

  // parse cut
  const auto pi_cutevt = xjjroot::parse_input(cutevtstr);
  const auto cutevt = pi_cutevt.content, cutevt_mc = save::cut_adjust_to_mc(cutevt);
  const auto pi_cutd = xjjroot::parse_input(cutdstr);
  const auto cutd = pi_cutd.content;
  
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

  auto cut_eff_num = project(trs.at("mc"), "_eff_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && Dgen==23333" + " && " + cutd);
  cast_branch("cut_eff_num", cut_eff_num);
  auto cut_reco_num = project(trs.at("mc"), "_reco_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && Dgen==23333 && fabs(Dtrk1Eta) < 2.4 && fabs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5");
  cast_branch("cut_reco_num", cut_reco_num);
  auto cut_acc_num = project(trs.at("mc"), "_acc_num", "nTrackInAcceptanceHP:Gpt:Gy", cutevt_mc + " && GisSignalCalc && fabs(Gtk1eta) < 2.4 && fabs(Gtk2eta) < 2.4 && Gtk1pt > 0.5 && Gtk2pt > 0.5");
  cast_branch("cut_acc_num", cut_acc_num);
  auto cut_eff_den = project(trs.at("mc"), "_eff_den", "nTrackInAcceptanceHP:Gpt:Gy", cutevt_mc + " && GisSignalCalc");
  cast_branch("cut_eff_den", cut_eff_den);
  std::string cut_data_signalwin = cutevt + " && " + cutd + " && fabs(Dmass-1.8648) < 0.03",
    cut_data_sideband = cutevt + " && " + cutd + " && fabs(Dmass-1.8648) > 0.09 && fabs(Dmass-1.8648) < 0.12";
  if (trs["data"]) {
    cast_branch("cut_data_signalwin", cut_data_signalwin);
    cast_branch("cut_data_sideband", cut_data_sideband);
    project(trs.at("data"), "_data_signalwin", "nTrackInAcceptanceHP:Dpt:Dy", cut_data_signalwin);
    project(trs.at("data"), "_data_sideband", "nTrackInAcceptanceHP:Dpt:Dy", cut_data_sideband);
  }
  
  cast_branch("inputmc", pi_inputmc.content);
  cast_branch("inputmc_tex", pi_inputmc.tex);
  if (trs.at("data")) {
    cast_branch("inputdata", pi_inputdata.content);
    cast_branch("inputdata_tex", pi_inputdata.tex);
  }
  cast_branch("cutevt", cutevt);
  cast_branch("cutevt_mc", cutevt_mc);
  cast_branch("cutevt_tex", pi_cutevt.tex);
  cast_branch("cutd", cutd);
  cast_branch("cutd_tex", pi_cutd.tex);

  t->Fill();
  t->Write();
  xjjroot::closefile(outf);
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 6) {
    return macro(argv[1], argv[2], argv[3], argv[4], argv[5]);
  }
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  return 1;
}

