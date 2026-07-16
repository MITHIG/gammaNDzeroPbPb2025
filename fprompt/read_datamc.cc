#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"

#include "xjjanauti.h"

int macro(std::string inputname_data,
          std::string inputname_prompt, std::string inputname_nonprompt,
          std::string outputname) {
  std::map<std::string, xjjc::info> infos;
  std::map<std::string, std::map<std::string, std::vector<TH1D*>>> h1ys;

  // data inputs
  auto* inf_data = TFile::Open(inputname_data.c_str());
  if (!inf_data) {
    __XJJLOG << "!! bad data input: " << inputname_data << ", abort."<< std::endl;
    return 2;
  }
  infos["data"] = xjjana::get_info(inf_data, "data/info");
  auto* h3_bins = xjjana::getobj<TH3D>(inf_data, "h3_bins_y-mass-pt");
  std::map<std::string, TH1D*> h1s_dump; // vs var
  for (const std::string var : { "Dip3D", "Dip3Dsig" }) {
    auto* dir = (TDirectory*)inf_data->GetDirectory(Form("dir_%s", var.c_str()));
    if (!dir) {
      __XJJLOG << "?? no dir for variable: dir_" << var << ", skip."<< std::endl;
      continue;
    }
    for (const std::string type : { "data-sub", "data-sigswap" }) {
      auto h1s = xjjana::getobj_regexp<TH1D>(dir, ".+" + type + ".+");
      if (h1s.empty()) continue;
      if (!h1s_dump[var]) {
        h1s_dump[var] = (TH1D*)h1s.front()->Clone(Form("h1_dump_%s", var.c_str()));
        h1s_dump[var]->Reset("ICESM");
        // !! need to sort !!
      }
      h1ys[var][type] = h1s;
    }
  }

  // mc inputs
  auto read_file = [&infos, &h1ys, &h1s_dump](std::string inputname, std::string name) {
    auto* inf = TFile::Open(inputname.c_str());
    if (!inf) {
      __XJJLOG << "!! bad mc input: " << inputname << ", abort."<< std::endl;
      return 2;
    }
    auto datasets = xjjana::getobj_regexp<RooDataSet>(inf, ".*mc_match__y-.+");
    // !! need to sort
    for (auto* ds : datasets) {
      const auto index_y = std::atoi(xjjc::str_erasestar(ds->GetName(), "*__y-").c_str());
      for (auto& [var, _] : h1ys) {
        auto* h = (TH1D*)h1s_dump[var]->Clone(Form("h1_%s_mc-%s__y-%d", var.c_str(), name.c_str(), index_y));
        // h->Sumw2();
        __XJJLOG << "++ " << h->GetName() << std::endl;
        xjjc::saywait();
        for (int i = 0; i < ds->numEntries(); ++i) {
          const auto* row = ds->get(i);
          h->Fill(row->getRealValue(var.c_str()));
        }
        h1ys[var][name].push_back(h);
      }
    }
    infos[name] = xjjana::get_info(inf, "info");
    
    return 0;
  };

  if (read_file(inputname_prompt, "prompt")) return 2;
  if (read_file(inputname_nonprompt, "nonprompt")) return 2;

  for (auto& [var, h1s] : h1ys) {
    __XJJLOG << ">> " << var << std::endl;
    xjjroot::print_tab(h1s, 0);
  }

  auto* outf = xjjroot::newfile(outputname + ".root");
  for (auto& [var, h1s] : h1ys) {
    outf->mkdir(Form("dir_%s", var.c_str()))->cd();
    for (auto& [_, hh] : h1s)
      for (auto& h : hh) {
        auto* h_norm = (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "h1_", "h1_norm_").c_str());
        h_norm->Scale(1./h->Integral(), "width");
        xjjroot::writehist(h_norm);
      }
  }
  outf->cd();
  xjjroot::writehist(h3_bins);
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
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  if (argc == 3) {
    auto inputs = xjjc::str_divide_trim(argv[1], ",");
    if (inputs.size() < 3) return 1;
    return macro(inputs[0], inputs[1], inputs[2], argv[2]);
  }
  return 1;
}
