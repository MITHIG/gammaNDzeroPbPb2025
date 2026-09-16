#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"

#include "../include/save.h"
#define __BINS_PTY_PLACEHOLDER__
#define __BINS_MASS__
#include "../include/bins.h"

int macro(const std::string& inputstr, const std::string& cutstr, const std::string& output, int is_template = 0) {
  // parse input
  const auto pinput = xjjroot::parse_input(inputstr);
  // __XJJLOG << ">> " << pinput.content << std::endl;
  auto* trs = xjjana::chain_files(xjjc::str_divide_trim(pinput.content, ","), "Tree");
  if (!trs) { __XJJLOG << "!! bad input file, abort." << std::endl; return 2; }
  save::mask_branch(trs);

  const auto is_mc = xjjana::ismc_runnum(trs);
  __XJJLOG << ">> is_mc : " << is_mc << std::endl;

  // parse cut
  const auto pcut = xjjroot::parse_input(cutstr);
  auto cut = pcut.content;
  if (is_mc) cut = save::cut_adjust_to_mc(cut);

  // parse binning
  bins::print();
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

  if (is_template) {
    project("_mc-match", cut + " && Dgen==23333");
    project("_mc-swap", cut + " && Dgen==23344");
    project("_mc-kk", cut + " && DisSignalKK");
    project("_mc-pipi", cut + " && DisSignalpipi");
  } else {
    project("_data", cut);
  }

  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_cont;
  auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
    t_cont[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_cont[name]));
  };
  cast_branch("input", pinput.content);
  cast_branch("input_tex", pinput.tex);
  cast_branch("input_tag", pinput.tag);
  cast_branch("cut", pcut.content);
  cast_branch("cut_tex", pcut.tex);
  cast_branch("cut_tag", pcut.tag);
  cast_branch("is_mc", is_mc);
  cast_branch("is_template", is_template);
  t->Fill();
  t->Write();

  xjjroot::closefile(outf);
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[5], ",");
    return macro(argv[1], argv[2], argv[3], atoi(argv[6]));
  }
  return 1;
}
