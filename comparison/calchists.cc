#include "xjjanauti.h"

int macro(std::string input)
{
  auto* inf = TFile::Open(input.c_str());
  if (!inf) {
    __XJJLOG << "warning: bad input " << input << ", abort."<< std::endl;
    return 2;
  }

  // 
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  info["input_tag"] = xjjc::str_eraseall(xjjc::str_tag_from_file(input), "save_");
  xjjc::print_tab(info, -1);

  //
  auto* h3 = (TH3D*)inf->Get("h3_run_var_l1");

  //
  std::map<std::string, TH1D*> h1s;
  h1s["var"] = (TH1D*)h3->ProjectionY(Form("%s", xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_var").c_str()), // Y : var
                                      0, -1, // X : run
                                      0, -1, "e"); // Z : l1
  h1s["var__Max400"] = (TH1D*)h3->ProjectionY(Form("%s__Max400", xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_var").c_str()), // Y : var
                                              0, -1, // X : run
                                              1, 1, "e"); // Z : l1
  h1s["var__Min400"] = (TH1D*)h3->ProjectionY(Form("%s__Min400", xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_var").c_str()), // Y : var
                                              0, -1, // X : run
                                              2, 2, "e"); // Z : l1
  TH1D* h1 = nullptr;
  h1 = (TH1D*)h3->ProjectionX("h1", // X : run
                              0, -1, // Y : var
                              0, -1, "e"); // Z : l1
  h1s["run"] = xjjana::rmthemptybins(h1, xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_run").c_str());
  delete h1;
  h1 = (TH1D*)h3->ProjectionX("h1", // X : run
                              0, -1, // Y : var
                              1, 1, "e"); // Z : l1
  h1s["run__Max400"] = xjjana::rmthemptybins(h1, Form("%s__Max400", xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_run").c_str()));
  delete h1;
  h1 = (TH1D*)h3->ProjectionX("h1", // X : run
                              0, -1, // Y : var
                              2, 2, "e"); // Z : l1
  h1s["run__Min400"] = xjjana::rmthemptybins(h1, Form("%s__Min400", xjjc::str_replaceall(h3->GetName(), "h3_run_var_l1", "h1_run").c_str()));
  delete h1;

  auto* outf = xjjroot::newfile(xjjc::str_replaceall(input, "save_", "calc_"));
  xjjroot::writehist(h3);
  for (auto& h1 : h1s) {
    xjjroot::sethempty(h1.second, 0, 0);
    xjjroot::writehist(h1.second);
  }
  auto* tinfo = xjjana::write_info(info);
  tinfo->Fill();
  tinfo->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc==2) {
    return macro(argv[1]);
  }
  return 1;
}
