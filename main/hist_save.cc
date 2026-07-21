#include <TH3D.h>
#include "xjjanauti.h"

#include "../include/save.h"
#define __BINS_PTY_ANA__
#define __BINS_MASS__
#include "../include/bins.h"

int macro(std::string inputstr, std::string cutstr, std::string output, int isdata = 1) {
  // parse input
  auto inputs = xjjc::str_divide_trim(inputstr, ";");
  auto* trs = xjjana::chain_files(xjjc::str_divide_trim(inputs[0], ","), "Tree");
  if (!trs) { __XJJLOG<<"!! bad input file "<<inputs[0]<<std::endl; return 2; }
  save::mask_branch(trs);

  // parse cut
  auto cuts = xjjc::str_divide_trim(cutstr, ";");
  auto cut = cuts[0];
  if (!isdata) cut = save::cut_adjust_to_mc(cut);

  // parse binning
  __XJJLOG << ">> current y binning:" << std::endl;
  xjjc::print_vec_h(bins::ybins, 0);
  __XJJLOG << ">> current pt binning:" << std::endl;
  xjjc::print_vec_h(bins::ptbins, 0);
  auto massbins = xjjc::fixedbin_to_edges(bins::nmass, bins::minmass, bins::maxmass);

  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  
  std::map<std::string, TH3D*> h3;
  auto project = [&trs, &h3, &massbins](std::string key, std::string icut) {
    h3[key] = new TH3D(Form("h3%s", key.c_str()), ";y;m_{K#pi} [GeV];p_{T} [GeV]",
                       bins::ybins.size()-1, bins::ybins.data(),
                       massbins.size()-1, massbins.data(),
                       bins::ptbins.size()-1, bins::ptbins.data());
    __XJJLOG << ">> "<<h3[key]->GetName()<<" \e[2m"<<icut<<"\e[0m"<<std::endl;
    xjjc::saywait();
    trs->Project(h3[key]->GetName(), "Dpt:Dmass:Dy", icut.c_str());
    xjjroot::writehist(h3[key]);
  };

  project((isdata ? "_data" : ""), cut);
  if (!isdata) {
    project("_match", cut + " && Dgen==23333");
    project("_swap", cut + " && Dgen==23344");
  }
  
  auto* t = new TTree("info", "");
  t->Branch("input", &(inputs[0]));
  auto input_tex = inputs.size() > 1 ? inputs[1] : "";
  t->Branch("input_tex", &input_tex);
  auto input_tag = inputs.size() > 2 ? inputs[2] : "";
  t->Branch("input_tag", &input_tag);
  t->Branch("cut", &cut);
  auto cut_tex = cuts.size() > 1 ? cuts[1] : "";
  t->Branch("cut_tex", &cut_tex);
  auto cut_tag = cuts.size() > 2 ? cuts[2] : "";
  t->Branch("cut_tag", &cut_tag);
  t->Fill();
  t->Write();
  outf->cd();

  outf->Close();
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[5], ",");
    return macro(argv[1], argv[2], argv[3], atoi(argv[6]));
  }
  if (argc == 6) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    return macro(argv[1], argv[2], argv[3], atoi(argv[5]));
  }
  if (argc == 5) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}

