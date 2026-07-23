#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"

#include "xjjanauti.h"

#define __COOK_NAME__
#include "style.h"

int macro(const std::string& inputname_data,
          const std::string& inputname_prompt, const std::string& inputname_nonprompt,
          const std::string& outputname,
          const std::string& var) {

  std::map<std::string, xjjc::info> infos;
  std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, std::vector<std::vector<TH1D*>>> h1ysfs;

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
    sfmin = 0.7; sfmax = 1.2;
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
    auto h1v = xjjana::getobj_regexp<TH1D>(dir, ".+" + type + ".+");
    std::sort(h1v.begin(), h1v.end(),
              [](const TH1D* ha, const TH1D* hb) {
                return index_y(ha) < index_y(hb);
              });
    if (h1v.empty()) continue;
    if (!h1_dump) {
      h1_dump = (TH1D*)h1v.front()->Clone(Form("h1_dump_%s", var.c_str()));
      h1_dump->Reset("ICESM");
    }
    h1ys[type] = h1v;
  }

  auto ny = h3_bins->GetXaxis()->GetNbins();
  
  // mc inputs
  auto read_file = [&infos, &ny, &h1ys, &h1ysfs, &h1_dump, &var, &nsf, &sfmin, &sfstep](std::string inputname, std::string name) {
    auto* inf = TFile::Open(inputname.c_str());
    if (!inf) {
      __XJJLOG << "!! bad mc input: " << inputname << ", abort."<< std::endl;
      return 2;
    }
    auto datasets = xjjana::getobj_regexp<RooDataSet>(inf, ".*mc_match__y-.+");
    std::sort(datasets.begin(), datasets.end(),
              [](const RooDataSet* ha, const RooDataSet* hb) {
                return index_y(ha) < index_y(hb);
              });
    if (datasets.size() != ny) {
      __XJJLOG << "!! inconsistent number of y bins : " << ny << " vs " << datasets.size() << ", abort." << std::endl;
      return 3;
    }
    for (auto* ds : datasets) { // loop y
      auto* h = (TH1D*)h1_dump->Clone(Form("h1_%s_mc-%s__y-%d", var.c_str(), name.c_str(), index_y(ds)));
      // h->Sumw2();
      std::vector<TH1D*> hs_sf(nsf, nullptr);
      for (int k=0; k<nsf; k++)
        hs_sf[k] = (TH1D*)h->Clone(add_suffix__y(h->GetName(), Form("_sf-%d", k)).c_str());

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
      h1ysfs[name].push_back(hs_sf);
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
    for (auto& [_, hhh] : h1ysfs) {
      auto hh = hhh[i];
      for (auto& h_sf : hh) {
        auto* h_sf_norm = (TH1D*)h_sf->Clone(xjjc::str_replaceall(h_sf->GetName(), "h1_", "h1_norm_").c_str());
        h_sf_norm->Scale(1./h_sf->Integral(), "width");
        xjjroot::writehist(h_sf_norm);
      }
    }
    outf->cd();
  }
  outf->cd();
  xjjroot::writehist(h3_bins);
  xjjroot::writehist(h1_bins_sf);
  outf->mkdir("info")->cd();
  for (auto& [iname, info] : infos) {
    auto* t_data = new TTree(iname.c_str(), "");
    for (auto& [key, content] : info) {
      t_data->Branch(key.c_str(), &content);
    }
    t_data->Fill();
    t_data->Write();
  }
  outf->cd();
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
