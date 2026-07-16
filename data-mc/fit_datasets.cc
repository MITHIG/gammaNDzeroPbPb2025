#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
// #include "RooStats/SPlot.h"
#include "RooHist.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

#define __BINS_MASS__
#include "../include/bins.h"

#include "../include/util.h"
#include "../include/droofitter.h"
#include "../include/draw.h"

int macro(std::string inputname_data, std::string inputname_template, std::string outputname) {

  std::map<std::string, std::vector<RooDataSet*>> datays;
  TH2D* h2_bins = nullptr;
  std::map<std::string, xjjc::info> infos;
  auto read_file = [&datays, &infos, &h2_bins](std::string inputname, std::string dname, std::string infoname = "") {
    const auto inputfile = util::parse_input(inputname).content;
    auto* inf = TFile::Open(inputfile.c_str());
    if (!inf || inf->IsZombie()) {
      __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
      return 2;
    }
    auto info = xjjana::getval_regexp(dynamic_cast<TTree*>(inf->Get("info")));
    if (!infoname.empty() && infos.find(infoname) == infos.end()) {
      xjjc::print_tab(info, -1);
      infos[infoname] = info;
    }
    auto dsy = xjjana::getobj_regexp<RooDataSet>(inf, dname + "__y-.+");
    if (dsy.empty()) {
      __XJJLOG << "!! no RooDataSet: " << dname << "__y*, abort." << std::endl;
      return 2;
    }
    datays[dname] = dsy;

    if (!h2_bins) h2_bins = xjjana::getobj<TH2D>(inf, "h2_bins_y-pt");
    
    return 0;
  };
  
  if (read_file(inputname_data, "data_main", "data")) return 2;
  if (read_file(inputname_template, "mc_match", "template")) return 2;
  if (read_file(inputname_template, "mc_swap")) return 2;
    
  // RooRealVar Dmass("Dmass", v_by_name("Dmass").vartex.c_str(), bins::minmass, bins::maxmass);
  auto Dmass = dynamic_cast<RooRealVar*>(datays.at("data_main").front()->get()->find("Dmass"));

  std::vector<droofitter*> fitterys;
  for (int i=0; i<datays.at("data_main").size(); i++) {
    auto* ft = new droofitter(datays.at("data_main")[i], datays.at("mc_match")[i], datays.at("mc_swap")[i], *Dmass); 
    ft->fit(bins::minmass, bins::maxmass);
    fitterys.push_back(ft);
  }

  draw::bintex btex(h2_bins, 0, 1);
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  for (int i=0; i<fitterys.size(); i++) {
    pdf->prepare();
    auto* frame_data_y = fitterys[i]->draw_data(bins::nmass);
    frame_data_y->Draw();
    fitterys[i]->leg()->Draw();
    xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt(), infos.at("data").at("cut_tex") }, 0.035, 13);
    xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
    pdf->write();

    pdf->prepare();
    auto* frame_mc_y = fitterys[i]->draw_mc_swap(bins::nmass);
    fitterys[i]->draw_mc_sig(bins::nmass, frame_mc_y);
    frame_mc_y->Draw();
    xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt() }, 0.035, 13);
    xjjroot::drawCMS(xjjroot::CMS::simulation, infos.at("template").at("input_tex") + " (5.36 TeV)");
    pdf->write();
  }

  pdf->close();

  auto* h3_bins = new TH3D("h3_bins_y-mass-pt", ";y;m_{#piK} [GeV/c^{2}];#it{p}_{T}",
                           h2_bins->GetXaxis()->GetNbins(), h2_bins->GetXaxis()->GetXbins()->GetArray(), // y 
                           bins::nmass, xjjc::fixedbin_to_edges(bins::nmass, bins::minmass, bins::maxmass).data(), // mass 
                           h2_bins->GetYaxis()->GetNbins(), h2_bins->GetYaxis()->GetXbins()->GetArray()); // pt
  
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h3_bins);
  for (int i=0; i<fitterys.size(); i++) {
    auto* ws = fitterys[i]->make_ws(Form("ws__y-%d", i));
    ws->Write();
  }
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
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}
