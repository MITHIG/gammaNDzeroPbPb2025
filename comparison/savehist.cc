#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include "xjjrootuti.h"
#include "variables.h"

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

  // parse cut
  auto cuts = xjjc::str_divide_trim(cutstr, ";");
  auto cut = cuts[0], cut_tex = cuts[1];
  cut = cut + " && (isL1ZDCOr_Max400_Pixel || isL1ZDCOr_Min400_Max10000)";

  // parse input
  auto inputs = xjjc::str_divide_trim(inputname, ";");
  auto input = inputs[0], input_tex = inputs[1];
  auto* nt = xjjana::chain_files(xjjc::str_divide_trim(input, ","), "Tree");
  if (!nt) {
    __XJJLOG << "!! bad input file: " << input << std::endl; 
    return 2;
  }

  // output
  auto* outf = xjjroot::newfile("rootfiles/" + varname + "/save_" + output + ".root");
  // /eos/user/c/cmsdqm/www/CAF/certification/Collisions23HI/Cert_Collisions2023HI_374288_375823_Good_ZDC_Golden.json
  auto* h3 = new TH3D("h3_run_var_l1", Form(";Run;%s;isL1ZDCOr_Min400_Max10000", the_var.vartex.c_str()), 943, 374804, 375747, the_var.nbin, the_var.varmin, the_var.varmax, 2, 0, 2);
  __XJJLOG << ">> "<<h3->GetName()<<" [ "<<the_var.varname<<" ] \e[2m"<<cut<<"\e[0m"<<std::endl;
  nt->Project(h3->GetName(), Form("isL1ZDCOr_Min400_Max10000:%s:Run", the_var.var.c_str()), cut.c_str());
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
