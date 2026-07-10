#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
// #include "RooStats/SPlot.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

#define __VARIABLES_ROOINCL__
#include "variables.h"

#define __BINS_MASS__
#include "../include/bins.h"
#include "../include/util.h"

void set_list_constant(RooArgList vars, bool constant = true) {
  for (int i = 0; i < vars.getSize(); ++i) {
    auto *var = dynamic_cast<RooRealVar *>(vars.at(i));
    if (var) var->setConstant(constant);
  }
}

int macro(std::string inputname_data, std::string inputname_template, std::string outputname) {

  std::map<std::string, RooDataSet*> datas;
  auto read_file = [&datas](std::string inputname, std::string dname) {
    const auto inputfile = util::parse_input(inputname).file;
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
    return 0;
  };
  
  if (read_file(inputname_data, "data_main")) return 2;
  if (read_file(inputname_template, "mc_match")) return 2;
  if (read_file(inputname_template, "mc_swap")) return 2;
  
  // datas["data_reduced"] = static_cast<RooDataSet*>(datas.at("data_main")->reduce(Form("Dmass >= %f && Dmass < %f", bins::minmass, bins::maxmass)));

  auto v_by_name = [](std::string name) {
    for (const auto& v : variables) {
      if (v.varname == name) return v;
    }
    throw std::runtime_error(Form("bad varname %s", name.c_str()));
  };
  
  RooRealVar Dmass("Dmass", v_by_name("Dmass").vartex.c_str(), bins::minmass, bins::maxmass);
  Dmass.setBins(bins::nmass);

  RooRealVar mean("mean", "mean", 1.865, 1.82, 1.91);
  RooRealVar sig_sigma1("sig_sigma1", "signal core width", 0.006, 0.001, 0.060);
  RooRealVar sig_sigma2("sig_sigma2", "signal middle width", 0.015, 0.002, 0.090);
  RooRealVar sig_sigma3("sig_sigma3", "signal tail width", 0.035, 0.004, 0.150);
  RooGaussian pdf_sig_g1("pdf_sig_g1", "signal Gaussian 1", Dmass, mean, sig_sigma1);
  RooGaussian pdf_sig_g2("pdf_sig_g2", "signal Gaussian 2", Dmass, mean, sig_sigma2);
  RooGaussian pdf_sig_g3("pdf_sig_g3", "signal Gaussian 3", Dmass, mean, sig_sigma3);
  RooRealVar sig_frac1("sig_frac1", "signal Gaussian 1 fraction", 0.50, 0.0, 1.0);
  RooRealVar sig_frac2("sig_frac2", "signal Gaussian 2 fraction", 0.30, 0.0, 1.0);
  RooAddPdf pdf_sig("pdf_sig", "triple-Gaussian signal",
                    RooArgList(pdf_sig_g1, pdf_sig_g2, pdf_sig_g3),
                    RooArgList(sig_frac1, sig_frac2), true);
  // fir signal
  auto sigFit = std::unique_ptr<RooFitResult>(pdf_sig.fitTo(*datas.at("mc_match"), RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1)));

  RooRealVar swap_sigma("swap_sigma", "swap width", 0.035, 0.002, 0.200);
  RooGaussian pdf_swap("pdf_swap", "swap Gaussian", Dmass, mean, swap_sigma);
  mean.setConstant(true);
  // fit swap
  auto swapFit = std::unique_ptr<RooFitResult>(pdf_swap.fitTo(*datas.at("mc_swap"), RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1)));

  // prepare to fit data
  set_list_constant(RooArgList(mean, sig_sigma1, sig_sigma2, sig_sigma3, sig_frac1, sig_frac2,
                               swap_sigma));
  const double signalFractionValue = 0.5;
  const double swapFractionValue = 1.0 - signalFractionValue;
  RooRealVar frac_sig("frac_sig", "signal fraction in signal+swap", signalFractionValue);
  frac_sig.setConstant(true);
  RooAddPdf pdf_sigswap("pdf_sigswap", "fixed signal + swap mixture",
                        RooArgList(pdf_sig, pdf_swap),
                        RooArgList(frac_sig));

  RooRealVar bkg_c1("bkg_c1", "background c1", 0.0, -2.0, 2.0);
  RooRealVar bkg_c2("bkg_c2", "background c2", 0.0, -2.0, 2.0);
  RooRealVar bkg_c3("bkg_c3", "background c3", 0.0, -2.0, 2.0);
  RooChebychev pdf_bkg("pdf_bkg", "third-order combinatorial background",
                       Dmass, RooArgList(bkg_c1, bkg_c2, bkg_c3));

  const double ncount_data = datas.at("data_main")->numEntries();
  RooRealVar n_sigswap("n_sigswap", "signal+swap yield", 0.30 * ncount_data, 0.0,
                       2.0 * ncount_data);
  RooRealVar n_bkg("n_bkg", "background yield", 0.70 * ncount_data, 0.0, 3.0 * ncount_data);
  RooAddPdf pdf_total("pdf_total", "fixed(signal + swap) + background",
                      RooArgList(pdf_sigswap, pdf_bkg),
                      RooArgList(n_sigswap, n_bkg));

  auto dataFit = std::unique_ptr<RooFitResult>(pdf_total.fitTo(*datas.at("data_main"), RooFit::Extended(true), RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1)));

  __XJJLOG << "++ sig fit:" << std::endl;
  sigFit->Print("v");
  __XJJLOG << "++ swap fit:" << std::endl;
  swapFit->Print("v");
  __XJJLOG << "++ data fit:" << std::endl;
  dataFit->Print("v");
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  pdf->prepare();
  auto *frame_sig = Dmass.frame();
  xjjroot::sethempty(frame_sig);
  datas.at("mc_match")->plotOn(frame_sig, RooFit::Binning(bins::nmass));
  pdf_sig.plotOn(frame_sig, RooFit::LineColor(kRed + 1));
  pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g1), RooFit::LineStyle(kDashed),
                 RooFit::LineColor(kBlue + 1));
  pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g2), RooFit::LineStyle(kDashed),
                 RooFit::LineColor(kGreen + 2));
  pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g3), RooFit::LineStyle(kDashed),
                 RooFit::LineColor(kMagenta + 1));
  // frame_save(frame_sig, Form("%s_signal", outPrefix));
  frame_sig->Draw();
  pdf->write();

  pdf->prepare();
  auto *frame_swap = Dmass.frame();
  xjjroot::sethempty(frame_swap);
  datas.at("mc_swap")->plotOn(frame_swap, RooFit::Binning(bins::nmass));
  pdf_swap.plotOn(frame_swap, RooFit::LineColor(kRed + 1));
  // frame_save(frame_swap, Form("%s_swap", outPrefix));
  frame_swap->Draw();
  pdf->write();
  
  pdf->prepare();
  auto *frame_data = Dmass.frame();
  xjjroot::sethempty(frame_sig);
  datas.at("data_main")->plotOn(frame_data, RooFit::Binning(bins::nmass));
  pdf_total.plotOn(frame_data, RooFit::LineColor(kRed + 1));
  pdf_total.plotOn(frame_data, RooFit::Components(pdf_sig), RooFit::LineColor(kBlue + 1),
                   RooFit::LineStyle(kDashed));
  pdf_total.plotOn(frame_data, RooFit::Components(pdf_swap), RooFit::LineColor(kGreen + 2),
                   RooFit::LineStyle(kDashed));
  pdf_total.plotOn(frame_data, RooFit::Components(pdf_bkg), RooFit::LineColor(kMagenta + 1),
                   RooFit::LineStyle(kDashed));
  pdf_total.paramOn(frame_data, RooFit::Layout(0.58, 0.88, 0.88),
                    RooFit::Parameters(RooArgSet(n_sigswap, n_bkg, frac_sig, mean,
                                                 sig_sigma1, sig_sigma2, sig_sigma3, swap_sigma,
                                                 bkg_c1, bkg_c2, bkg_c3)));
  // frame_save(frame_data, Form("%s_data", outPrefix));
  frame_data->Draw();
  pdf->write();

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

  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}
