#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"
#include "dfitter_cached.h"
#include "draw.h"

namespace {
  TH1D* get_mass_hist(TFile* file, const std::string& component, int ipt, int iy) {
    return static_cast<TH1D*>(file->Get(Form("h1_mass_%s__pt-%d__y-%d",
                                             component.c_str(), ipt, iy)));
  }

  TH1D* make_toy(const TH1D* source, TF1* truth_total, TRandom3& rng,
                 int toy, int ipt, int iy) {
    auto* result = static_cast<TH1D*>(source->Clone(Form("closure_mass__toy-%d__pt-%d__y-%d",
                                                         toy, ipt, iy)));
    result->SetDirectory(nullptr);
    for (int ibin = 1; ibin <= source->GetNbinsX(); ++ibin) {
      const auto low = source->GetXaxis()->GetBinLowEdge(ibin);
      const auto high = source->GetXaxis()->GetBinUpEdge(ibin);
      const auto expected = std::max(0., truth_total->Integral(low, high) /
                                     source->GetBinWidth(ibin));
      const auto count = rng.Poisson(expected);
      result->SetBinContent(ibin, count);
      result->SetBinError(ibin, std::sqrt(static_cast<double>(count)));
    }
    return result;
  }
}

int macro(const std::string& inputname, const std::string& outputname,
          const std::string& fit_opt = "3G-Peaky", int ntoys = 100,
          unsigned int seed = 12345) {
  if (ntoys <= 0) return 1;
  gROOT->SetBatch(true);
  xjjroot::setgstyle(1);
  TCanvas canvas("closure_canvas", "closure_canvas");

  auto* input = TFile::Open(inputname.c_str(), "READ");
  if (!input || input->IsZombie()) return 2;
  auto* h3_bins = static_cast<TH3D*>(input->Get("h3_bins"));
  if (!h3_bins) return 3;

  const draw::bintex bins(h3_bins, 0, 2);
  const auto fitopt = xjjroot::parse_input(fit_opt);
  auto* fout = xjjroot::newfile("rootfiles/" + outputname + ".root");
  auto* tree = new TTree("toys", "Fitted-model closure-test results");

  int toy = -1, pt = -1, y = -1;
  double yield = -1, yield_err = -1, truth_yield = -1;
  int fit_ok = 0;
  tree->Branch("toy", &toy);
  tree->Branch("pt", &pt);
  tree->Branch("y", &y);
  tree->Branch("yield", &yield);
  tree->Branch("yield_err", &yield_err);
  tree->Branch("truth_yield", &truth_yield);
  tree->Branch("fit_ok", &fit_ok);

  TRandom3 rng(seed);
  for (int ipt = 0; ipt < bins.npt(); ++ipt) {
    for (int iy = 0; iy < bins.ny(); ++iy) {
      auto* hdata = get_mass_hist(input, "data", ipt, iy);
      auto* hmatch = get_mass_hist(input, "mc-match", ipt, iy);
      auto* hswap = get_mass_hist(input, "mc-swap", ipt, iy);
      auto* hkk = get_mass_hist(input, "mc-kk", ipt, iy);
      auto* hpipi = get_mass_hist(input, "mc-pipi", ipt, iy);
      if (!hdata || !hmatch || !hswap || !hkk || !hpipi) return 4;

      xjjroot::cached_dfitter cached_fit(fitopt.content);
      if (!cached_fit.prepare(hdata, hmatch, hswap, hkk, hpipi)) return 5;
      auto* truth_total = cached_fit.total_model(Form("closure_total__pt-%d__y-%d", ipt, iy));
      auto* truth_signal = cached_fit.signal_model(Form("closure_signal__pt-%d__y-%d", ipt, iy));
      if (!truth_total || !truth_signal) return 6;
      truth_yield = truth_signal->Integral(hdata->GetXaxis()->GetXmin(),
                                           hdata->GetXaxis()->GetXmax()) /
                    hdata->GetBinWidth(1);
      pt = ipt;
      y = iy;

      for (toy = 0; toy < ntoys; ++toy) {
        auto* htoy = make_toy(hdata, truth_total, rng, toy, ipt, iy);
        const auto fit = cached_fit.fit(htoy);
        fit_ok = fit.fitted ? 1 : 0;
        yield = fit_ok ? fit.yield : -1;
        yield_err = fit_ok ? fit.yieldErr : -1;
        tree->Fill();
        delete htoy;
      }
      delete truth_total;
      delete truth_signal;
      __XJJLOG << "++ closure complete pt=" << ipt << ", y=" << iy
               << " (" << ntoys << " toys)" << std::endl;
    }
  }

  tree->Write();
  auto* h3_bins_out = static_cast<TH3D*>(h3_bins->Clone("h3_bins"));
  h3_bins_out->Write();
  auto fitopt_content = fitopt.content;
  auto input_content = inputname;
  auto* info = new TTree("info", "closure-test configuration");
  info->Branch("input", &input_content);
  info->Branch("fitopt", &fitopt_content);
  info->Branch("ntoys", &ntoys);
  info->Branch("seed", &seed);
  int npt = bins.npt(), ny = bins.ny();
  info->Branch("npt", &npt);
  info->Branch("ny", &ny);
  info->Fill();
  info->Write();
  xjjroot::closefile(fout);
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc < 3 || argc > 6) return 1;
  const auto fitopt = argc > 3 ? argv[3] : "3G-Peaky";
  const auto ntoys = argc > 4 ? std::atoi(argv[4]) : 100;
  const auto seed = argc > 5 ? static_cast<unsigned int>(std::atoi(argv[5])) : 12345u;
  return macro(argv[1], argv[2], fitopt, ntoys, seed);
}
