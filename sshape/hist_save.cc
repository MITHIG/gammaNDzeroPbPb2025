#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"

namespace bins {
  int nmass = 54; double minmass = 1.73, maxmass = 2;
  std::vector<double> ybins = { -1.5, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.5 };
}

int macro(const std::string& inputstr, const std::string& treename,
          const std::string& cutstr, const std::string& output, int is_template = 0) {
  // parse input
  const auto pinput = xjjroot::parse_input(inputstr);
  auto* inf = xjjroot::openfile(pinput.content);
  if (!inf) return 2;
  auto* trs = (TTree*)inf->Get(treename.c_str());
  if (!trs) return 2;

  // parse cut
  const auto pcut = xjjroot::parse_input(cutstr);
  auto cut = pcut.content;
  
  auto massbins = xjjc::fixedbin_to_edges(bins::nmass, bins::minmass, bins::maxmass);

  auto* outf = xjjroot::newfile("rootfiles/" + output + ".root");
  std::map<std::string, TH2D*> h2;
  auto project = [&trs, &h2, &massbins](const std::string& key, const std::string& icut) {
    h2[key] = new TH2D(Form("h2%s", key.c_str()), ";y;m_{K#pi} [GeV];",
                       bins::ybins.size()-1, bins::ybins.data(),
                       massbins.size()-1, massbins.data());
    __XJJLOG << ">> "<<h2[key]->GetName()<<" \e[2m"<<icut<<"\e[0m"<<std::endl;
    xjjc::saywait();
    trs->Project(h2[key]->GetName(), "Dmass:Dy", icut.c_str());
    xjjroot::writehist(h2[key]);
  };

  if (is_template) {
    project("_mc-match", cut + " && Dgen==23333");
    project("_mc-swap", cut + " && Dgen==23344");
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
  rinfo.cast_branch("is_template", is_template);
  rinfo.close();

  xjjroot::closefile(outf);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 6)
    return macro(argv[1], argv[2], argv[3], argv[4], std::atoi(argv[5]));
  return 1;
}
