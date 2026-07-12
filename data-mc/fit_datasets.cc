#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
// #include "RooStats/SPlot.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

// #define __VARIABLES_ROOINCL__
// #include "variables.h"
#define __BINS_MASS__
#include "../include/bins.h"

#include "../include/util.h"
#include "../include/droofitter.h"

int macro(std::string inputname, std::string outputname) {

  std::map<std::string, RooDataSet*> datas;
  std::map<std::string, std::vector<RooDataSet*>> datays;
  TH2D* h2_bins = nullptr;
  auto read_file = [&datas, &datays, &h2_bins](std::string inputname, std::string dname) {
    if (datas.find(dname) != datas.end()) {
      __XJJLOG << "!! repleated keyname: " << dname << ", abort." << std::endl;
      return 2;
    }
    const auto inputfile = util::parse_input(inputname).content;
    auto* inf = TFile::Open(inputfile.c_str());
    if (!inf || inf->IsZombie()) {
      __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
      return 2;
    }

    auto* ds = dynamic_cast<RooDataSet*>(inf->Get(dname.c_str()));
    if (!ds) {
      __XJJLOG << "!! no RooDataSet: " << dname << ", abort." << std::endl;
      // throw std::runtime_error(Form("Could not find RooDataSet %s", name));
      return 2;
    }
    xjjroot::print_obj(ds);
    if (ds->numEntries() == 0) return 2;
    datas[dname] = ds;

    auto dsy = xjjana::getobj_regexp<RooDataSet>(inf, dname + "__y-.+");
    if (dsy.empty()) {
      __XJJLOG << "!! no RooDataSet: " << dname << "__y*, abort." << std::endl;
      return 2;
    }
    datays[dname] = dsy;

    if (!h2_bins) h2_bins = xjjana::getobj<TH2D>(inf, "h1_bins");
    
    return 0;
  };
  
  if (read_file(inputname, "data_main")) return 2;
  if (read_file(inputname, "mc_match")) return 2;
  if (read_file(inputname, "mc_swap")) return 2;
    
  // auto v_by_name = [](std::string name) {
  //   for (const auto& v : variables) {
  //     if (v.varname == name) return v;
  //   }
  //   throw std::runtime_error(Form("bad varname %s", name.c_str()));
  // };
 
  // RooRealVar Dmass("Dmass", v_by_name("Dmass").vartex.c_str(), bins::minmass, bins::maxmass);
  // Dmass.setBins(bins::nmass); // did anything?
  // Dmass.setRange("range_fit", bins::minmass, bins::maxmass);

  auto* fitter = new droofitter(datas.at("data_main"), datas.at("mc_match"), datas.at("mc_swap"));
  fitter->fit(bins::minmass, bins::maxmass);
  
  // xjjroot::setgstyle(1);
  // auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  // pdf->prepare();
  // auto *frame_sig = Dmass.frame();
  // xjjroot::sethempty(frame_sig);
  // datas.at("mc_match")->plotOn(frame_sig, RooFit::Binning(bins::nmass));
  // pdf_sig.plotOn(frame_sig, RooFit::LineColor(kRed + 1));
  // pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g1), RooFit::LineStyle(kDashed),
  //                RooFit::LineColor(kBlue + 1));
  // pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g2), RooFit::LineStyle(kDashed),
  //                RooFit::LineColor(kGreen + 2));
  // pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g3), RooFit::LineStyle(kDashed),
  //                RooFit::LineColor(kMagenta + 1));
  // // frame_save(frame_sig, Form("%s_signal", outPrefix));
  // frame_sig->Draw();
  // pdf->write();

  // pdf->prepare();
  // auto *frame_swap = Dmass.frame();
  // xjjroot::sethempty(frame_swap);
  // datas.at("mc_swap")->plotOn(frame_swap, RooFit::Binning(bins::nmass));
  // pdf_swap.plotOn(frame_swap, RooFit::LineColor(kRed + 1));
  // // frame_save(frame_swap, Form("%s_swap", outPrefix));
  // frame_swap->Draw();
  // pdf->write();
  
  // pdf->prepare();
  // auto *frame_data = Dmass.frame();
  // xjjroot::sethempty(frame_sig);
  // datas.at("data_main")->plotOn(frame_data, RooFit::Binning(bins::nmass));
  // pdf_total.plotOn(frame_data, RooFit::LineColor(kRed + 1));
  // pdf_total.plotOn(frame_data, RooFit::Components(pdf_sig), RooFit::LineColor(kBlue + 1),
  //                  RooFit::LineStyle(kDashed));
  // pdf_total.plotOn(frame_data, RooFit::Components(pdf_swap), RooFit::LineColor(kGreen + 2),
  //                  RooFit::LineStyle(kDashed));
  // pdf_total.plotOn(frame_data, RooFit::Components(pdf_bkg), RooFit::LineColor(kMagenta + 1),
  //                  RooFit::LineStyle(kDashed));
  // pdf_total.paramOn(frame_data, RooFit::Layout(0.58, 0.88, 0.88),
  //                   RooFit::Parameters(RooArgSet(n_sigswap, n_bkg, frac_sig, par_mean,
  //                                                par_sig_sigma1, par_sig_sigma2, par_sig_sigma3, par_swap_sigma,
  //                                                par_bkg_c1, par_bkg_c2, par_bkg_c3)));
  // // frame_save(frame_data, Form("%s_data", outPrefix));
  // frame_data->Draw();
  // pdf->write();

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

  // pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
