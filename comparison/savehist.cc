#include <TFile.h>
#include <TTree.h>
#include <TH3D.h>
#include "xjjrootuti.h"
#include "xjjstruct.h"
#include "variables.h"

#define __BINS_PTY_INCL__
#include "bins.h"

int macro(const std::string& inputstr, const std::string& cutstr,
          const std::string& varname, const std::string& output)
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
  auto isGvar = std::regex_match(the_var.varname, std::regex("-?G.+"));
  __XJJLOG << ">> isDvar: " << isDvar << std::endl;
  __XJJLOG << ">> isGvar: " << isGvar << std::endl;
  const auto varbins = xjjc::fixedbin_to_edges(the_var.nbin, the_var.varmin, the_var.varmax);

  // parse cut
  auto const pcut = xjjroot::parse_input(cutstr);
  auto cut = pcut.content;

  // parse input
  const auto pinput = xjjroot::parse_input(inputstr);
  auto* nt = xjjana::chain_files(xjjc::str_divide_trim(pinput.content, ","), "Tree");
  if (!nt) {
    __XJJLOG << "!! bad input file, abort." << std::endl; 
    return 2;
  }
  const auto nentries = nt->GetEntries();

  // output
  auto* outf = xjjroot::newfile(output + ".root");
  TH3D* h3;
  if (isDvar || isGvar) {
    h3 = new TH3D("h3_y_var_pt", Form(";y;%s;p_{T} [GeV]", the_var.vartex.c_str()),
                  bins::ybins.size()-1, bins::ybins.data(),
                  varbins.size()-1, varbins.data(),
                  bins::ptbins.size()-1, bins::ptbins.data());
  } else {
    h3 = new TH3D("h3_vz_var_l1", Form(";v_{z} [cm];%s;isL1ZDCOr_Min400_Max10000", the_var.vartex.c_str()),
                  8, -20, 20,
                  the_var.nbin, the_var.varmin, the_var.varmax,
                  2, 0, 2);
  }
  __XJJLOG << ">> "<<h3->GetName()<<" [ "<<the_var.varname<<" ] \e[2m"<<cut<<"\e[0m"<<std::endl;
  std::string str_proj = isDvar ? ("Dpt:" + the_var.var + ":Dy") : ( isGvar ? ("Gpt:" + the_var.var + ":Gy") : ("isL1ZDCOr_Min400_Max10000:" + the_var.var + ":VZ"));
  nt->Project(h3->GetName(), str_proj.c_str(), cut.c_str());
  xjjroot::writehist(h3);

  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_cont;
  auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
    t_cont[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_cont[name]));
  };
  cast_branch("varname", varname);
  cast_branch("cut", pcut.content);
  cast_branch("cut_tex", pcut.tex);
  cast_branch("input", pinput.content);
  cast_branch("input_tex", pinput.tex);
  cast_branch("nentries", nentries);
  t->Fill();
  t->Write();
  
  xjjroot::closefile(outf);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  return 1;
}
