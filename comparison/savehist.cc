#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include "xjjrootuti.h"
#include "variables.h"

#define __BINS_PTY_EQ__
#include "../include/bins.h"

int macro(std::string inputname, std::string cutstr,
          std::string varname, std::string output)
{
  // parse var
  auto it_var = std::find_if(vars.begin(), vars.end(), [&varname](const xjjana::variable& v) {
    return v.varname == varname;
  });
  if (it_var == vars.end()) {
    __XJJLOG << "!! bad varname: " << varname << std::endl;
    return 2;
  }
  const auto& the_var = *it_var;
  auto isDvar = std::regex_match(the_var.varname, std::regex("-?D.+"));
  __XJJLOG << ">> isDvar: " << isDvar << std::endl;

  // parse cut
  auto cuts = xjjc::str_divide_trim(cutstr, ";");
  auto cut = cuts[0], cut_tex = cuts[1];
  // cut = cut + " && isL1ZDCOr";

  // parse input
  auto inputs = xjjc::str_divide_trim(inputname, ";");
  auto input = inputs[0], input_tex = inputs[1];
  auto* nt = xjjana::chain_files(xjjc::str_divide_trim(input, ","), "Tree");
  if (!nt) {
    __XJJLOG << "!! bad input file: " << input << std::endl; 
    return 2;
  }

  // output
  auto* outf = xjjroot::newfile(output + ".root");
  // /eos/user/c/cmsdqm/www/CAF/certification/Collisions23HI/Cert_Collisions2023HI_374288_375823_Good_ZDC_Golden.json
  TH3D* h3;
  if (isDvar) {
    h3 = new TH3D("h3_y_var_pt", Form(";y;%s;p_{T} [GeV]", the_var.vartex.c_str()),
                  bins::ny, bins::miny, bins::maxy,
                  the_var.nbin, the_var.varmin, the_var.varmax,
                  bins::npt, bins::minpt, bins::maxpt);
  } else {
    h3 = new TH3D("h3_l1_var_vz", Form(";isL1ZDCOr_Min400_Max10000;%s;v_{z} [cm]", the_var.vartex.c_str()),
                  2, 0, 2,
                  the_var.nbin, the_var.varmin, the_var.varmax,
                  60, -15, 15);
  }
  __XJJLOG << ">> "<<h3->GetName()<<" [ "<<the_var.varname<<" ] \e[2m"<<cut<<"\e[0m"<<std::endl;
  std::string str_proj = isDvar ? ("Dpt:" + the_var.var + ":Dy") : ("VZ:" + the_var.var + ":isL1ZDCOr_Min400_Max10000");
  nt->Project(h3->GetName(), str_proj.c_str(), cut.c_str());
  xjjroot::writehist(h3);

  auto* tinfo = new TTree("info", "");
  tinfo->Branch("varname", &varname);
  tinfo->Branch("cut", &cut);
  tinfo->Branch("cut_tex", &cut_tex);
  tinfo->Branch("input", &input);
  tinfo->Branch("input_tex", &input_tex);
  tinfo->Fill();
  tinfo->Write();

  outf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  return 1;
}
