#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"

#include "xjjanauti.h"

int macro(const std::string& inputname_data,
          const std::string& inputname_prompt, const std::string& inputname_nonprompt,
          const std::string& outputname,
          const std::string& var) {

  std::map<std::string, xjjc::info> infos;
  std::map<std::string, std::vector<TH1D*>> h1ys;

  // data inputs
  auto* inf_data = TFile::Open(inputname_data.c_str());
  if (!inf_data) {
    __XJJLOG << "!! bad data input: " << inputname_data << ", abort."<< std::endl;
    return 2;
  }
  infos["data"] = xjjana::get_info(inf_data, "data/info");
  auto* h3_bins = xjjana::getobj<TH3D>(inf_data, "h3_bins_y-mass-pt");
  float sfstep = 0.02, sfmin = 0.9, sfmax = 1.6; // ip3D
  if (xjjc::str_contains(var, "sig")) { // ip3Dsig
    sfstep = 0.02, sfmin = 0.7; sfmax = 1.2;
  }
  int nsf = (sfmax-sfmin)/sfstep + 1;
  auto* h1_bins_sf = new TH1D("h1_bins_sf", ";Resolution scale factor#scale[0.5]{ }#alpha_{reso};", nsf, sfmin - sfstep*0.5, sfmax + sfstep*0.5);

  auto* dir = (TDirectory*)inf_data->GetDirectory(Form("dir_%s", var.c_str()));
  if (!dir) {
    __XJJLOG << "?? no dir for variable: dir_" << var << ", skip."<< std::endl;
    return 2;
  }
  TH1D* h1_dump = nullptr;
  for (const std::string type : { "data-sub", "data-sigswap" }) {
    auto h1s = xjjana::getobj_regexp<TH1D>(dir, ".+" + type + ".+");
    if (h1s.empty()) continue;
    if (!h1_dump) {
      h1_dump = (TH1D*)h1s.front()->Clone(Form("h1_dump_%s", var.c_str()));
      h1_dump->Reset("ICESM");
      // !! need to sort !!
    }
    h1ys[type] = h1s;
  }

  auto ny = h3_bins->GetXaxis()->GetNbins();
  
  // mc inputs
  auto read_file = [&infos, &ny, &h1ys, &h1_dump, &var, &nsf, &sfmin, &sfstep](std::string inputname, std::string name) {
    auto* inf = TFile::Open(inputname.c_str());
    if (!inf) {
      __XJJLOG << "!! bad mc input: " << inputname << ", abort."<< std::endl;
      return 2;
    }
    auto datasets = xjjana::getobj_regexp<RooDataSet>(inf, ".*mc_match__y-.+");
    // !! need to sort
    if (datasets.size() != ny) {
      __XJJLOG << "!! inconsistent number of y bins : " << ny << " vs " << datasets.size() << ", abort." << std::endl;
      return 3;
    }
    for (auto* ds : datasets) { // loop y
      const auto index_y = std::atoi(xjjc::str_erasestar(ds->GetName(), "*__y-").c_str());

      auto* h = (TH1D*)h1_dump->Clone(Form("h1_%s_mc-%s__y-%d", var.c_str(), name.c_str(), index_y));
      // h->Sumw2();
      std::vector<TH1D*> hs_sf(nsf, nullptr);
      for (int k=0; k<nsf; k++)
        hs_sf[k] = (TH1D*)h->Clone(Form("h1_%s_mc-%s_sf-%d__y-%d", var.c_str(), name.c_str(), k, index_y));

      __XJJLOG << "++ " << h->GetName() << std::endl;
      xjjc::saywait();
      for (int i = 0; i < ds->numEntries(); ++i) {
        const auto* row = ds->get(i);
        auto value = row->getRealValue(var.c_str());
        h->Fill(value);
        for (int k=0; k<nsf; k++)
          hs_sf[k]->Fill(value * (sfmin + k*sfstep));
      }
      h1ys[name].push_back(h);
      for (int k=0; k<nsf; k++) {
        h1ys[Form("%s_sf-%d", name.c_str(), k)].push_back(hs_sf[k]);
      }
    }
    infos[name] = xjjana::get_info(inf, "info");
    
    return 0;
  };

  if (read_file(inputname_prompt, "prompt")) return 2;
  if (read_file(inputname_nonprompt, "nonprompt")) return 2;

  xjjroot::print_tab(h1ys, 0);

  auto* outf = xjjroot::newfile(outputname + ".root");
  auto* dir_var = outf->mkdir(Form("dir_%s", var.c_str()));
  dir_var->cd();
  for (int i=0; i<ny; i++) {
    dir_var->mkdir(Form("dir__y-%d", i))->cd();
    for (auto& [_, hh] : h1ys) {
      auto* h = hh[i];
      auto* h_norm = (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "h1_", "h1_norm_").c_str());
      h_norm->Scale(1./h->Integral(), "width");
      xjjroot::writehist(h_norm);
    }
    outf->cd();
  }
  outf->cd();
  xjjroot::writehist(h3_bins);
  xjjroot::writehist(h1_bins_sf);
  for (auto& [iname, info] : infos) {
    outf->mkdir(iname.c_str())->cd();
    auto* t_data = new TTree("info", "");
    for (auto& [key, content] : info) {
      t_data->Branch(key.c_str(), &content);
    }
    t_data->Fill();
    t_data->Write();
    outf->cd();
  }
  xjjroot::closefile(outf);
  return 0;
}

int main(int argc, char* argv[]) {
  std::vector<std::string> inputs;
  if (argc == 6) {
    for (int i=1; i<argc; i++)
      inputs.push_back(argv[i]);
  }
  else if (argc == 4) {
    inputs = xjjc::str_divide_trim(argv[1], ",");
    inputs.push_back(argv[2]);
    inputs.push_back(argv[3]);
  }
  if (inputs.size() == 5) {
    return macro(inputs[0], inputs[1], inputs[2], inputs[3], inputs[4]);
  }
  return 1;
}
