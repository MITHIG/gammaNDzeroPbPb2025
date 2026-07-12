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
  RooPlot* draw_mc_sig(int nmass, RooPlot *frame = nullptr);
  RooPlot* draw_mc_swap(int nmass, RooPlot *frame = nullptr);
  RooPlot* draw_data(int nmass, RooPlot *frame = nullptr);
  static void draw_ds(RooDataSet*, RooPlot*);
private:
  RooRealVar Dmass;
  RooDataSet *ds_data_main, *ds_mc_match, *ds_mc_swap;
  std::unique_ptr<RooFitResult> fitr_sig, fitr_swap, fitr_data;
  std::unique_ptr<RooAddPdf> pdf_total, pdf_sig, pdf_sigswap;
  std::unique_ptr<RooGaussian> pdf_sig_g1, pdf_sig_g2, pdf_sig_g3, pdf_swap;
  std::unique_ptr<RooChebychev> pdf_bkg;
  std::unique_ptr<RooRealVar> par_mean,
    par_sig_sigma1, par_sig_sigma2, par_sig_sigma3, par_sig_frac1, par_sig_frac2,
    par_swap_sigma1, par_sigswap_frac,
    par_bkg_c1, par_bkg_c2, par_bkg_c3,
    n_sigswap, n_bkg;
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

  par_mean = std::make_unique<RooRealVar>("par_mean", "mean", 1.865, 1.82, 1.91);
  par_sig_sigma1 = std::make_unique<RooRealVar>("par_sig_sigma1", "signal core width", 0.006, 0.001, 0.060);
  par_sig_sigma2 = std::make_unique<RooRealVar>("par_sig_sigma2", "signal middle width", 0.015, 0.002, 0.090);
  par_sig_sigma3 = std::make_unique<RooRealVar>("par_sig_sigma3", "signal tail width", 0.035, 0.004, 0.150);
  pdf_sig_g1 = std::make_unique<RooGaussian>("pdf_sig_g1", "signal Gaussian 1", Dmass, *par_mean, *par_sig_sigma1);
  pdf_sig_g2 = std::make_unique<RooGaussian>("pdf_sig_g2", "signal Gaussian 2", Dmass, *par_mean, *par_sig_sigma2);
  pdf_sig_g3 = std::make_unique<RooGaussian>("pdf_sig_g3", "signal Gaussian 3", Dmass, *par_mean, *par_sig_sigma3);
  par_sig_frac1 = std::make_unique<RooRealVar>("par_sig_frac1", "signal Gaussian 1 fraction", 0.50, 0.0, 1.0);
  par_sig_frac2 = std::make_unique<RooRealVar>("par_sig_frac2", "signal Gaussian 2 fraction", 0.30, 0.0, 1.0);
  pdf_sig = std::make_unique<RooAddPdf>("pdf_sig", "triple-Gaussian signal",
                                        RooArgList(*pdf_sig_g1, *pdf_sig_g2, *pdf_sig_g3),
                                        RooArgList(*par_sig_frac1, *par_sig_frac2), true);
  // fit signal
  fitr_sig = std::unique_ptr<RooFitResult>(pdf_sig->fitTo(*ds_mc_match, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  par_swap_sigma1 = std::make_unique<RooRealVar>("par_swap_sigma1", "swap width", 0.035, 0.002, 0.200);
  // pdf_swap_g1 = std::make_unique<RooGaussian>("pdf_swap_g1", "swap Gaussian 1", Dmass, *par_mean, *par_swap_sigma1);
  pdf_swap = std::make_unique<RooGaussian>("pdf_swap", "swap Gaussian", Dmass, *par_mean, *par_swap_sigma1);
  par_mean->setConstant(true);
  // fit swap
  fitr_swap = std::unique_ptr<RooFitResult>(pdf_swap->fitTo(*ds_mc_swap, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  // prepare to fit data
  set_list_constant(RooArgList(*par_mean, *par_sig_sigma1, *par_sig_sigma2, *par_sig_sigma3, *par_sig_frac1, *par_sig_frac2,
                               *par_swap_sigma1));

  const double signalFractionValue = ds_mc_match->numEntries()*1. / (ds_mc_match->numEntries() + ds_mc_swap->numEntries()), // in fitting range
    swapFractionValue = 1.0 - signalFractionValue;
  par_sigswap_frac = std::make_unique<RooRealVar>("par_sigswap_frac", "signal fraction in signal+swap", signalFractionValue);
  par_sigswap_frac->setConstant(true);
  pdf_sigswap = std::make_unique<RooAddPdf>("pdf_sigswap", "fixed signal + swap mixture",
                                            RooArgList(*pdf_sig, *pdf_swap),
                                            RooArgList(*par_sigswap_frac));

  par_bkg_c1 = std::make_unique<RooRealVar>("par_bkg_c1", "background c1", 0.0, -2.0, 2.0);
  par_bkg_c2 = std::make_unique<RooRealVar>("par_bkg_c2", "background c2", 0.0, -2.0, 2.0);
  par_bkg_c3 = std::make_unique<RooRealVar>("par_bkg_c3", "background c3", 0.0, -2.0, 2.0);
  pdf_bkg = std::make_unique<RooChebychev>("pdf_bkg", "third-order combinatorial background",
                                           Dmass, RooArgList(*par_bkg_c1, *par_bkg_c2, *par_bkg_c3));

  const double ncount_data = ds_data_main->numEntries();
  n_sigswap = std::make_unique<RooRealVar>("n_sigswap", "signal+swap yield", 0.30 * ncount_data, 0.0,
                                           2.0 * ncount_data);
  n_bkg = std::make_unique<RooRealVar>("n_bkg", "background yield", 0.70 * ncount_data, 0.0, 3.0 * ncount_data);
  pdf_total = std::make_unique<RooAddPdf>("pdf_total", "fixed(signal + swap) + background",
                                          RooArgList(*pdf_sigswap, *pdf_bkg),
                                          RooArgList(*n_sigswap, *n_bkg));

  fitr_data = std::unique_ptr<RooFitResult>(pdf_total->fitTo(*ds_data_main, RooFit::Extended(true), // Extended(true) tells RooFit to fit both the shape and the total number of events, instead of only the shape
                                                             RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1), RooFit::Range("range_fit")));

  __XJJLOG << "++ sig fit:" << std::endl;
  fitr_sig->Print("v");
  __XJJLOG << "++ swap fit:" << std::endl;
  fitr_swap->Print("v");
  __XJJLOG << "++ data fit:" << std::endl;
  fitr_data->Print("v");
}

void droofitter::draw_ds(RooDataSet* ds, RooPlot* frame) {
  xjjroot::sethempty(frame, 0, 0.1);
  // ds->plotOn(frame, RooFit::Binning(nmass), RooFit::MarkerStyle(20), RooFit::MarkerSize(1.3), RooFit::XErrorSize(0));
  ds->plotOn(frame, RooFit::MarkerStyle(20), RooFit::MarkerSize(1.3), RooFit::XErrorSize(0));
  frame->SetMinimum(0);
  frame->SetMaximum(xjjana::gethmaximum(frame->getHist())*1.4*1.2);
  auto* x = frame->getPlotVar();
  int nbins = x->getBins();
  auto xmin = x->getMin(), xmax = x->getMax();
  double binWidth = (xmax - xmin) / nbins;
  frame->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", binWidth*1.e+3));
}

RooPlot* droofitter::draw_mc_sig(int nmass, RooPlot *frame) {
  // pdf_total->Print("t");
  if (!frame) {
    Dmass.setBins(nmass);
    frame = Dmass.frame();
  }
  draw_ds(ds_mc_match, frame);
  pdf_sig->plotOn(frame, RooFit::LineColor(kRed + 1));
  pdf_sig->plotOn(frame, RooFit::Components(*pdf_sig_g1), RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kBlue + 1));
  pdf_sig->plotOn(frame, RooFit::Components(*pdf_sig_g2), RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kGreen + 2));
  pdf_sig->plotOn(frame, RooFit::Components(*pdf_sig_g3), RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kMagenta + 1));
  return frame;
}

RooPlot* droofitter::draw_mc_swap(int nmass, RooPlot *frame) {
  if (!frame) {
    Dmass.setBins(nmass);
    frame = Dmass.frame();
  }
  draw_ds(ds_mc_swap, frame);
  pdf_swap->plotOn(frame, RooFit::LineColor(kRed + 1));
  return frame;
}

RooPlot* droofitter::draw_data(int nmass, RooPlot *frame) {
  if (!frame) {
    Dmass.setBins(nmass);
    frame = Dmass.frame();
  }
  draw_ds(ds_data_main, frame);
  pdf_total->plotOn(frame, RooFit::LineColor(kRed + 1));
  pdf_total->plotOn(frame, RooFit::Components("pdf_sig"), RooFit::LineColor(kBlue + 1),
                    RooFit::LineStyle(kDashed));
  pdf_total->plotOn(frame, RooFit::Components("pdf_swap"), RooFit::LineColor(kGreen + 2),
                    RooFit::LineStyle(kDashed));
  pdf_total->plotOn(frame, RooFit::Components("pdf_bkg"), RooFit::LineColor(kMagenta + 1),
                    RooFit::LineStyle(kDashed));
  // pdf_total->paramOn(frame, RooFit::Layout(0.58, 0.88, 0.88));
  // pdf_total->paramOn(frame, RooFit::Layout(0.58, 0.88, 0.88),
  //                   RooFit::Parameters(RooArgSet(n_sigswap, n_bkg, par_sigswap_frac, par_mean,
  //                                                par_sig_sigma1, par_sig_sigma2, par_sig_sigma3, par_swap_sigma,
  //                                                par_bkg_c1, par_bkg_c2, par_bkg_c3)));
  // frame->Draw();
  return frame;
}
