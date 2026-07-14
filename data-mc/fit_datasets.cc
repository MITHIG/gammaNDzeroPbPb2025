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

int macro(std::string inputname, std::string outputname) {
  const auto inputfile = util::parse_input(inputname).content;
  auto* inf = TFile::Open(inputfile.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
    return 2;
  }
  std::map<std::string, xjjc::info> infos;
  for (std::string tr : { "data/info", "template/info" }) {
    auto info = xjjana::getval_regexp((TTree*)inf->Get(tr.c_str()));
    infos[xjjc::str_eraseall(tr, "/info")] = info;
  }
  for (auto& [key, info] : infos) {
    __XJJLOG << "++ infos [" << key << "]" << std::endl;
    xjjc::print_tab<std::string, std::string>(info, -1);
  }

  std::map<std::string, std::vector<RooDataSet*>> datays;
  TH2D* h2_bins = nullptr;
  auto read_file = [&inf, &datays, &h2_bins](std::string dname) {
    auto dsy = xjjana::getobj_regexp<RooDataSet>(inf, dname + "__y-.+");
    if (dsy.empty()) {
      __XJJLOG << "!! no RooDataSet: " << dname << "__y*, abort." << std::endl;
      return 2;
    }
    datays[dname] = dsy;

    if (!h2_bins) h2_bins = xjjana::getobj<TH2D>(inf, "h2_bins_y-pt");
    
    return 0;
  };
  
  if (read_file("data_main")) return 2;
  if (read_file("mc_match")) return 2;
  if (read_file("mc_swap")) return 2;
    
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
  outf->Close();
  
  // RooStats::SPlot sData("sData", "sData", *datas.at("data_main"), &pdf_total,
  //                       RooArgList(n_sigswap, n_bkg));

  // TH1D hSigDpt("hSigDpt", ";D^{0} p_{T} (GeV);Signal sWeighted candidates",
  //              ptBins, ptMin, ptMax);
  // TH1D hSigSwapDpt("hSigSwapDpt",
  //                  ";D^{0} p_{T} (GeV);Signal+swap sWeighted candidates",
  //                  ptBins, ptMin, ptMax);
  // hSigDpt.Sumw2();
  // hSigSwapDpt.Sumw2();
  // RooArgSet massObs(Dmass);
  // for (int i = 0; i < datas.at("data_main")->numEntries(); ++i) {
  //   const RooArgSet *row = datas.at("data_main")->get(i);
  //   const auto *ptValue = dynamic_cast<const RooRealVar *>(row->find("Dpt"));
  //   const auto *massValue =
  //     dynamic_cast<const RooRealVar *>(row->find("Dmass"));
  //   const auto *weight =
  //     dynamic_cast<const RooRealVar *>(row->find("n_sigswap_sw"));
  //   if (!ptValue || !massValue || !weight) {
  //     throw std::runtime_error("Could not read Dpt/Dmass/n_sigswap_sw.");
  //   }

  //   Dmass.setVal(massValue->getVal());
  //   const double sigDensity = pdf_sig.getVal(&massObs);
  //   const double swapDensity = pdf_swap.getVal(&massObs);
  //   const double signalPart =
  //     signalFractionValue * sigDensity /
  //     (signalFractionValue * sigDensity + swapFractionValue * swapDensity);
  //   const double sigSwapWeight = weight->getVal();

  //   hSigSwapDpt.Fill(ptValue->getVal(), sigSwapWeight);
  //   hSigDpt.Fill(ptValue->getVal(), sigSwapWeight * signalPart);
  // }

  // TCanvas sPlotCanvas(Form("%s_sPlotDpt_canvas", outPrefix),
  //                     "Signal sPlot Dpt", 850, 750);
  // hSigDpt.SetMarkerStyle(kFullCircle);
  // hSigDpt.SetMarkerSize(0.9);
  // hSigDpt.SetLineColor(kBlue + 1);
  // hSigDpt.SetMarkerColor(kBlue + 1);
  // hSigDpt.Draw("E1");
  // sPlotCanvas.SaveAs(Form("%s_sPlotDpt.pdf", outPrefix));
  // sPlotCanvas.SaveAs(Form("%s_sPlotDpt.png", outPrefix));

  // TFile output(Form("%s_workspaceInputs.root", outPrefix), "RECREATE");
  // sigFit->Write("signalMcFitResult");
  // swapFit->Write("swapMcFitResult");
  // dataFit->Write("dataFitResult");
  // hSigDpt.Write("hSigDpt");
  // hSigSwapDpt.Write("hSigSwapDpt");
  // output.Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
