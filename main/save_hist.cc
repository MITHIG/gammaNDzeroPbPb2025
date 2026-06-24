#include <TH3D.h>
#include "xjjanauti.h"

#include "../include/save.h"
#define __BINS_PTY_EQ__
#define __BINS_MASS__
#include "../include/bins.h"

int macro(std::string inputstr, std::string cutstr, std::string output, int isdata = 1) {
  std::cout<<std::endl;

  // parse input
  auto inputs = xjjc::str_divide_trim(inputstr, ";");
  auto input_tex = inputs.size() > 1 ? inputs[1] : "";
  auto* trs = xjjana::chain_files(xjjc::str_divide_trim(inputs[0], ","), "Tree");
  if (!trs) { __XJJLOG<<"!! bad input file "<<inputs[0]<<std::endl; return 2; }
  save::mask_branch(trs);

  // parse cut
  auto cuts = xjjc::str_divide_trim(cutstr, ";");
  auto cut = cuts[0];
  if (isdata) cut += " && isL1ZDCOr && cscTightHalo2015Filter";
  
  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  
  std::map<std::string, TH3D*> h3;
  auto project = [&trs, &h3](std::string key, std::string icut) {
    h3[key] = new TH3D(Form("h3%s", key.c_str()), ";y;m_{K#pi} [GeV];p_{T} [GeV]",
                       bins::ny, bins::miny, bins::maxy,
                       bins::nmass, bins::minmass, bins::maxmass,
                       bins::npt, bins::minpt, bins::maxpt);
    __XJJLOG << ">> "<<h3[key]->GetName()<<" \e[2m"<<icut<<"\e[0m"<<std::endl;
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
  t->Branch("input_tex", &input_tex);
  t->Branch("cut", &cut);
  t->Branch("cut_tex", &(cuts[1]));
  t->Fill();
  t->Write();
  outf->cd();

  outf->Close();
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], atoi(argv[4]));
  }
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}

