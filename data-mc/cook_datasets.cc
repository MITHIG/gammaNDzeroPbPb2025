#include "RooDataSet.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/util.h"

#define __VARIABLES_ROOINCL__
#include "variables.h"
#define __BINS_PTY_FPROMPT__
#include "../include/bins.h"

int macro(std::string inputname_data, std::string inputname_template, std::string outputname) {

  std::map<std::string, RooDataSet*> datas;
  std::map<std::string, xjjc::info> infos;
  auto read_file = [&datas, &infos](std::string inputname, std::string dname, std::string iname) {
    const auto inputfile = util::parse_input(inputname).content;
    auto* inf = TFile::Open(inputfile.c_str());
    if (!inf || inf->IsZombie()) {
      __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
      return 2;
    }
    auto* dataset = dynamic_cast<RooDataSet*>(inf->Get(dname.c_str()));
    if (!dataset) {
      __XJJLOG << "!! no RooDataSet: " << dname << ", abort." << std::endl;
      // throw std::runtime_error(Form("Could not find RooDataSet %s", name));
      return 2;
    }
    if (datas.find(dname) != datas.end()) {
      __XJJLOG << "!! repleated keyname: " << dname << ", abort." << std::endl;
      return 2;
    }
    __XJJLOG << ">> " << std::left << std::setw(10) << dname << ": "<< dataset->numEntries() << std::endl;
    if (dataset->numEntries() == 0) return 2;
    
    datas[dname] = dataset;
    if (infos.find(iname) == infos.end())
      infos[iname] = xjjana::getval_regexp((TTree*)inf->Get("info"));
    return 0;
  };
  
  if (read_file(inputname_data, "data_main", "data")) return 2;
  if (read_file(inputname_template, "mc_match", "template")) return 2;
  if (read_file(inputname_template, "mc_swap", "template")) return 2;

  xjjroot::print_tab(datas, 0);
  
  // parse binning
  auto* h2_bins = new TH2D("h2_bins_y-pt", ";y;#it{p}_{T}",
                           bins::ybins.size()-1, bins::ybins.data(),
                           bins::npt, xjjc::fixedbin_to_edges(bins::npt, bins::minpt, bins::maxpt).data());
  xjjc::print_vec_h(bins::ybins, 0);

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
  
  auto* outf = xjjroot::newfile(outputname + ".root");
  xjjroot::writehist(h2_bins);
  for (auto& [_, vds] : datays)
    for (auto& ds : vds)
      xjjroot::writehist(ds);

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
  
  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}
