#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"

#include "../include/save.h"
#define __BINS_PTY_PLACEHOLDER__
#define __BINS_MASS__
#include "../include/bins.h"

int macro(const std::string& inputstr, const std::string& cutstr, const std::string& output,
          int is_template, const std::string& templatestr = "1") {
  // parse input
  const auto pinput = xjjroot::parse_input(inputstr);
  auto* trs = xjjana::chain_files(xjjc::str_divide_trim(pinput.content, ","), "Tree");
  if (!trs) { __XJJLOG << "!! bad input file, abort." << std::endl; return 2; }
  save::mask_branch(trs);

  const auto is_mc = xjjana::ismc_runnum(trs);
  __XJJLOG << ">> is_mc : " << is_mc << std::endl;

  // template
  const auto ptemplate = xjjroot::parse_input(templatestr);
  const auto alpha = ptemplate.content;
  if (is_template)
    __XJJLOG << ">> alpha : " << alpha << std::endl;
  
  // parse cut
  const auto pcut = xjjroot::parse_input(cutstr);
  auto cut = pcut.content;
  if (is_mc) cut = save::cut_adjust_to_mc(cut);

  // parse binning
  bins::print();
  auto massbins = xjjc::fixedbin_to_edges(bins::nmass, bins::minmass, bins::maxmass);

  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  
  std::map<std::string, TH3D*> h3;
  auto project = [&trs, &h3, &massbins](const std::string& key, const std::string& icut, const std::string& fillstr = "Dpt:Dmass:Dy") {
    h3[key] = new TH3D(Form("h3%s", key.c_str()), ";y;m_{K#pi} [GeV];p_{T} [GeV]",
                       bins::ybins.size()-1, bins::ybins.data(),
                       massbins.size()-1, massbins.data(),
                       bins::ptbins.size()-1, bins::ptbins.data());
    __XJJLOG << ">> " << h3[key]->GetName() << " \e[2m" << icut << "\e[0m [" << fillstr << "]" << std::endl;
    xjjc::saywait();
    trs->Project(h3[key]->GetName(), fillstr.c_str(), icut.c_str());
    xjjroot::writehist(h3[key]);
  };

  if (is_template) {
    std::string scale_mass = "1.865+(Dmass-1.865)*" + alpha;
    project("_mc-match", cut + " && Dgen==23333", "Dpt:" + scale_mass + ":Dy");
    project("_mc-swap", cut + " && Dgen==23344", "Dpt:" + scale_mass + ":Dy");
    project("_mc-kk", cut + " && DisSignalKK");
    project("_mc-pipi", cut + " && DisSignalpipi");
  } else {
    project("_data", cut);
  }

  xjjroot::Writeinfo rinfo("info");
  rinfo.cast_branch("input", pinput.content);
  rinfo.cast_branch("input_tex", pinput.tex);
  rinfo.cast_branch("input_tag", pinput.tag);
  rinfo.cast_branch("cut", pcut.content);
  rinfo.cast_branch("cut_tex", pcut.tex);
  rinfo.cast_branch("cut_tag", pcut.tag);
  rinfo.cast_branch("is_mc", is_mc);
  rinfo.cast_branch("is_template", is_template);
  rinfo.cast_branch("alpha", alpha);
  rinfo.cast_branch("template_tex", ptemplate.tex);
  rinfo.cast_branch("template_tag", ptemplate.tag);
  rinfo.close();

  xjjroot::closefile(outf);
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 8) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[5], ",");
    return macro(argv[1], argv[2], argv[3], atoi(argv[6]), argv[7]);
  }
  return 1;
}
