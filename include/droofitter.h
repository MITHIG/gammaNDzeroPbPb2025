#pragma once

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"

class droofitter
{
public:
  droofitter(RooDataSet* data_main, RooDataSet* mc_match, RooDataSet* mc_swap);
  droofitter(RooDataSet* data_main, RooDataSet* mc_match, RooDataSet* mc_swap, RooRealVar aDmass);
  void fit(float minmass, float maxmass);
  // void draw_mc_sig(int nmass);
  // void draw_mc_swap(int nmass);
  // void draw_data(int nmass);
private:
  RooRealVar Dmass;
  std::unique_ptr<RooFitResult> sigFit, swapFit, dataFit;
  RooDataSet *ds_data_main, *ds_mc_match, *ds_mc_swap;
  void set_list_constant(RooArgList vars, bool constant = true) {
    for (int i = 0; i < vars.getSize(); ++i) {
      auto *var = dynamic_cast<RooRealVar *>(vars.at(i));
      if (var) var->setConstant(constant);
    }
  }
};

droofitter::droofitter(RooDataSet* data_main, RooDataSet* mc_match, RooDataSet* mc_swap)
  : Dmass(*dynamic_cast<RooRealVar*>(data_main->get()->find("Dmass"))),
    ds_data_main(data_main), ds_mc_match(mc_match), ds_mc_swap(mc_swap) {
  // Dmass = dynamic_cast<RooRealVar*>(ds_data_main->get()->find("Dmass"));
}

droofitter::droofitter(RooDataSet* data_main, RooDataSet* mc_match, RooDataSet* mc_swap, RooRealVar aDmass)
  : Dmass(aDmass),
    ds_data_main(data_main), ds_mc_match(mc_match), ds_mc_swap(mc_swap) {}

void droofitter::fit(float minmass, float maxmass) {
  Dmass.setRange("range_fit", minmass, maxmass);
  // Dmass.setBins(nmass);

  RooRealVar par_mean("par_mean", "mean", 1.865, 1.82, 1.91);
  RooRealVar par_sig_sigma1("par_sig_sigma1", "signal core width", 0.006, 0.001, 0.060);
  RooRealVar par_sig_sigma2("par_sig_sigma2", "signal middle width", 0.015, 0.002, 0.090);
  RooRealVar par_sig_sigma3("par_sig_sigma3", "signal tail width", 0.035, 0.004, 0.150);
  RooGaussian pdf_sig_g1("pdf_sig_g1", "signal Gaussian 1", Dmass, par_mean, par_sig_sigma1);
  RooGaussian pdf_sig_g2("pdf_sig_g2", "signal Gaussian 2", Dmass, par_mean, par_sig_sigma2);
  RooGaussian pdf_sig_g3("pdf_sig_g3", "signal Gaussian 3", Dmass, par_mean, par_sig_sigma3);
  RooRealVar par_sig_frac1("par_sig_frac1", "signal Gaussian 1 fraction", 0.50, 0.0, 1.0);
  RooRealVar par_sig_frac2("par_sig_frac2", "signal Gaussian 2 fraction", 0.30, 0.0, 1.0);
  RooAddPdf pdf_sig("pdf_sig", "triple-Gaussian signal",
                    RooArgList(pdf_sig_g1, pdf_sig_g2, pdf_sig_g3),
                    RooArgList(par_sig_frac1, par_sig_frac2), true);
  // fit signal
  sigFit = std::unique_ptr<RooFitResult>(pdf_sig.fitTo(*ds_mc_match, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  RooRealVar par_swap_sigma("par_swap_sigma", "swap width", 0.035, 0.002, 0.200);
  RooGaussian pdf_swap("pdf_swap", "swap Gaussian", Dmass, par_mean, par_swap_sigma);
  par_mean.setConstant(true);
  // fit swap
  swapFit = std::unique_ptr<RooFitResult>(pdf_swap.fitTo(*ds_mc_swap, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  // prepare to fit data
  set_list_constant(RooArgList(par_mean, par_sig_sigma1, par_sig_sigma2, par_sig_sigma3, par_sig_frac1, par_sig_frac2,
                               par_swap_sigma));

  const double signalFractionValue = ds_mc_match->numEntries()*1. / (ds_mc_match->numEntries() + ds_mc_swap->numEntries()), // in fitting range
    swapFractionValue = 1.0 - signalFractionValue;
  RooRealVar frac_sig("frac_sig", "signal fraction in signal+swap", signalFractionValue);
  frac_sig.setConstant(true);
  RooAddPdf pdf_sigswap("pdf_sigswap", "fixed signal + swap mixture",
                        RooArgList(pdf_sig, pdf_swap),
                        RooArgList(frac_sig));

  RooRealVar par_bkg_c1("par_bkg_c1", "background c1", 0.0, -2.0, 2.0);
  RooRealVar par_bkg_c2("par_bkg_c2", "background c2", 0.0, -2.0, 2.0);
  RooRealVar par_bkg_c3("par_bkg_c3", "background c3", 0.0, -2.0, 2.0);
  RooChebychev pdf_bkg("pdf_bkg", "third-order combinatorial background",
                       Dmass, RooArgList(par_bkg_c1, par_bkg_c2, par_bkg_c3));

  const double ncount_data = ds_data_main->numEntries();
  RooRealVar n_sigswap("n_sigswap", "signal+swap yield", 0.30 * ncount_data, 0.0,
                       2.0 * ncount_data);
  RooRealVar n_bkg("n_bkg", "background yield", 0.70 * ncount_data, 0.0, 3.0 * ncount_data);
  RooAddPdf pdf_total("pdf_total", "fixed(signal + swap) + background",
                      RooArgList(pdf_sigswap, pdf_bkg),
                      RooArgList(n_sigswap, n_bkg));

  dataFit = std::unique_ptr<RooFitResult>(pdf_total.fitTo(*ds_data_main, RooFit::Extended(true), // Extended(true) tells RooFit to fit both the shape and the total number of events, instead of only the shape
                                                          RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  __XJJLOG << "++ sig fit:" << std::endl;
  sigFit->Print("v");
  __XJJLOG << "++ swap fit:" << std::endl;
  swapFit->Print("v");
  __XJJLOG << "++ data fit:" << std::endl;
  dataFit->Print("v");
}

// void droofitter::draw_mc_sig(int nmass) {
//   auto *frame_sig = Dmass.frame();
//   xjjroot::sethempty(frame_sig);
//   ds_mc_match->plotOn(frame_sig, RooFit::Binning(nmass));
//   pdf_sig.plotOn(frame_sig, RooFit::LineColor(kRed + 1));
//   pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g1), RooFit::LineStyle(kDashed),
//                  RooFit::LineColor(kBlue + 1));
//   pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g2), RooFit::LineStyle(kDashed),
//                  RooFit::LineColor(kGreen + 2));
//   pdf_sig.plotOn(frame_sig, RooFit::Components(pdf_sig_g3), RooFit::LineStyle(kDashed),
//                  RooFit::LineColor(kMagenta + 1));
//   frame_sig->Draw();
// }

// void droofitter::draw_mc_swap(int nmass) {
//   auto *frame_swap = Dmass.frame();
//   xjjroot::sethempty(frame_swap);
//   ds_mc_match->plotOn(frame_swap, RooFit::Binning(nmass));
//   pdf_swap.plotOn(frame_swap, RooFit::LineColor(kRed + 1));
//   frame_swap->Draw();
// }

// void droofitter::draw_data(int nmass) {
//   auto *frame_data = Dmass.frame();
//   xjjroot::sethempty(frame_sig);
//   ds_data_main->plotOn(frame_data, RooFit::Binning(nmass));
//   pdf_total.plotOn(frame_data, RooFit::LineColor(kRed + 1));
//   pdf_total.plotOn(frame_data, RooFit::Components(pdf_sig), RooFit::LineColor(kBlue + 1),
//                    RooFit::LineStyle(kDashed));
//   pdf_total.plotOn(frame_data, RooFit::Components(pdf_swap), RooFit::LineColor(kGreen + 2),
//                    RooFit::LineStyle(kDashed));
//   pdf_total.plotOn(frame_data, RooFit::Components(pdf_bkg), RooFit::LineColor(kMagenta + 1),
//                    RooFit::LineStyle(kDashed));
//   pdf_total.paramOn(frame_data, RooFit::Layout(0.58, 0.88, 0.88),
//                     RooFit::Parameters(RooArgSet(n_sigswap, n_bkg, frac_sig, par_mean,
//                                                  par_sig_sigma1, par_sig_sigma2, par_sig_sigma3, par_swap_sigma,
//                                                  par_bkg_c1, par_bkg_c2, par_bkg_c3)));
//   frame_data->Draw();
// }
