#include <TH3D.h>
#include "xjjanauti.h"

#include "../include/save.h"
#define __BINS_PTY_EFF__
#define __BINS_MULT__
#include "../include/bins.h"

int macro(std::string inputmcstr, std::string cutevtstr, std::string cutdstr, std::string output, std::string inputdatastr = "null") {
  std::map<std::string, TChain*> trs;
  // parse inputmc
  auto inputmcs = xjjc::str_divide_trim(inputmcstr, ";");
  auto inputmc_tex = inputmcs.size() > 1 ? inputmcs[1] : "";
  trs["mc"] = xjjana::chain_files(xjjc::str_divide_trim(inputmcs[0], ","), "Tree");
  if (!trs.at("mc")) { __XJJLOG<<"!! bad inputmc file "<<inputmcs[0]<<std::endl; return 2; }
  save::mask_branch(trs.at("mc"));

  auto inputdatas = xjjc::str_divide_trim(inputdatastr, ";");
  auto inputdata_tex = inputdatas.size() > 1 ? inputdatas[1] : "";
  trs["data"] = xjjana::chain_files(xjjc::str_divide_trim(inputdatas[0], ","), "Tree");
  if (trs.at("data")) {
    __XJJLOG << "?? no valid inputdata file: " << inputdatas[0] << ", only MC is used." << std::endl;
  }

  // parse cut
  auto cutevts = xjjc::str_divide_trim(cutevtstr, ";");
  auto cutevt = cutevts[0], cutevt_mc = save::cut_adjust_to_mc(cutevt);
  auto cutds = xjjc::str_divide_trim(cutdstr, ";");
  auto cutd = cutds[0];
  
  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  auto* t = new TTree("info", "");

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

  auto cut_eff_num = project(trs.at("mc"), "_eff_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && Dgen==23333" + " && " + cutd);
  t->Branch("cut_eff_num", &cut_eff_num);
  auto cut_reco_num = project(trs.at("mc"), "_reco_num", "nTrackInAcceptanceHP:Dpt:Dy", cutevt_mc + " && Dgen==23333 && fabs(Dtrk1Eta) < 2.4 && fabs(Dtrk2Eta) < 2.4 && Dtrk1Pt > 0.5 && Dtrk2Pt > 0.5");
  t->Branch("cut_reco_num", &cut_reco_num);
  auto cut_eff_den = project(trs.at("mc"), "_eff_den", "nTrackInAcceptanceHP:Gpt:Gy", cutevt_mc + " && GisSignalCalc");
  t->Branch("cut_eff_den", &cut_eff_den);
  std::string cut_data_signalwin = cutevt + " && " + cutd + " && fabs(Dmass-1.8648) < 0.03",
    cut_data_sideband = cutevt + " && " + cutd + " && fabs(Dmass-1.8648) > 0.09 && fabs(Dmass-1.8648) < 0.12";
  if (trs["data"]) {
    t->Branch("cut_data_signalwin", &cut_data_signalwin);
    t->Branch("cut_data_sideband", &cut_data_sideband);
    project(trs.at("data"), "_data_signalwin", "nTrackInAcceptanceHP:Dpt:Dy", cut_data_signalwin);
    project(trs.at("data"), "_data_sideband", "nTrackInAcceptanceHP:Dpt:Dy", cut_data_sideband);
  }
  
  t->Branch("inputmc", &(inputmcs[0]));
  t->Branch("inputmc_tex", &inputmc_tex);
  if (trs.at("data")) {
    t->Branch("inputdata", &(inputdatas[0]));
    t->Branch("inputdata_tex", &inputdata_tex);
  }
  t->Branch("cutevt", &cutevt);
  t->Branch("cutevt_mc", &cutevt_mc);
  t->Branch("cutevt_tex", &(cutevts[1]));
  t->Branch("cutd", &cutd);
  t->Branch("cutd_tex", &(cutds[1]));
  t->Fill();
  t->Write();
  outf->cd();

  outf->Close();
  
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

