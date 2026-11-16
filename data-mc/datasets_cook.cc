#include "RooDataSet.h"

#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"

#define __VARIABLES_ROOSPLOT__
#include "variables.h"
#define __BINS_PTY_PLACEHOLDER__
#include "../include/bins.h"

int macro(const std::string& inputname, const std::string& outputname) {
  // parse binning
  bins::print();
  auto* h2_bins = new TH2D("h2_bins_y-pt", ";y;#it{p}_{T}",
                           bins::ybins.size()-1, bins::ybins.data(),
                           bins::ptbins.size()-1, bins::ptbins.data());

  const auto inputfile = xjjroot::parse_input(inputname).content;
  auto* inf = TFile::Open(inputfile.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
    return 2;
  }
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
 
  auto datasets = xjjana::getobj_regexp<RooDataSet>(inf);
  if (datasets.empty()) {
    __XJJLOG << "!! no RooDataSet, abort." << std::endl;
    xjjroot::closefile(inf);
    // throw std::runtime_error(Form("Could not find RooDataSet %s", name));
    return 2;
  }
  std::map<std::string, RooDataSet*> datas;
  for (auto& ds : datasets) {
    __XJJLOG << ">> " << std::left << std::setw(10) << ds->GetName() << ": "<< ds->numEntries() << std::endl;
    if (ds->numEntries() == 0) return 2;
    datas[ds->GetName()] = ds;
  }
  xjjroot::print_tab(datas, 0);
  
  __XJJLOG << "++ split y bins" << std::endl;
  std::map<std::string, std::vector<RooDataSet*>> datays;
  for (int i=0; i<bins::ybins.size()-1; i++) {
    for (auto& [t, ds] : datas) {
      // ds->Print("v");
      // std::cout << Form("Dy >= %f && Dy < %f", bins::ybins[i], bins::ybins[i+1]) << std::endl;
      auto* ds_reduced = static_cast<RooDataSet*>(ds->reduce(Form("Dy >= %f && Dy < %f", bins::ybins[i], bins::ybins[i+1])));
      ds_reduced->SetName(Form("%s__y-%d", ds->GetName(), i));
      datays[t].push_back(ds_reduced);
      __XJJLOG << "   >> " << ds_reduced->GetName() << " (" << ds_reduced->numEntries() << ")" << std::endl;
    }
  }
  
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h2_bins);
  for (auto& [_, vds] : datays)
    for (auto& ds : vds)
      xjjroot::writehist(ds);

  auto* t_data = new TTree("info", "");
  for (auto& [key, content] : info) {
    t_data->Branch(key.c_str(), &content);
  }
  t_data->Fill();
  t_data->Write();
  
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[3], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[4], ",");
    return macro(argv[1], argv[2]);
  }
  return 1;
}
