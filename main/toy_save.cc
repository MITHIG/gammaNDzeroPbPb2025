#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include "xjjanauti.h"
#include "xjjstruct.h"
#include "dfitter_cached.h"
#include "draw.h"
#include "util.h"

namespace {
  const std::vector<std::string> type_names = { "Bootstrap", "Closure" };
  TH1D* make_toy_boostrap(const TH1D* source, TRandom3& rng, const std::string& tunique) {
    auto* result = static_cast<TH1D*>(source->Clone(Form("h1_mass__%s", tunique.c_str())));
    result->SetDirectory(nullptr);
    for (int ibin = 1; ibin <= source->GetNbinsX(); ++ibin) {
      const auto count = rng.Poisson(std::max(0., source->GetBinContent(ibin)));
      result->SetBinContent(ibin, count);
      result->SetBinError(ibin, std::sqrt(static_cast<double>(count)));
    }
    return result;
  }
  TH1D* make_toy_closure(const TH1D* source, TF1* truth_total, TRandom3& rng, const std::string& tunique) {
    auto* result = static_cast<TH1D*>(source->Clone(Form("h1_mass__%s", tunique.c_str())));
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
          int toytype = 0,
          int ntoys = 500, unsigned int seed = 12345,
          const std::string& fit_opt_force = "") {

  if (ntoys <= 0) return 2;
  gROOT->SetBatch(true);

  auto* inf = xjjroot::readfile(inputname);
  if (!inf) return 2;
  const auto info = util::read_info(inf, "fit/info");
  const auto fitopt = fit_opt_force.empty() ? info.at("fitopt") : fit_opt_force;
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins");
  if (!h3_bins) return 3;
  const draw::bintex bins(h3_bins, 0, 2);

  auto* fout = xjjroot::newfile("rootfiles/" + outputname + ".root");
  auto* tree = new TTree("toys", "");
  int itoy = -1; tree->Branch("itoy", &itoy);
  int pt = -1; tree->Branch("ipt", &pt);
  int y = -1; tree->Branch("iy", &y);
  double yield = -1; tree->Branch("yield", &yield);
  double yield_err = -1; tree->Branch("yield_err", &yield_err);
  double nominal_yield = -1; tree->Branch("nominal_yield", &nominal_yield); 
  double nominal_yield_err = -1; tree->Branch("nominal_yield_err", &nominal_yield_err);
  // double truth_yield = -1; if (toytype == 0) { tree->Branch("truth_yield", &truth_yield); }

  TRandom3 rng(seed);
  for (int ipt = 0; ipt < bins.npt(); ++ipt) {
    for (int iy = 0; iy < bins.ny(); ++iy) {
      auto get_mass_hist = [&inf, &ipt, &iy](const std::string& comp) -> TH1D* {
        return xjjana::getobj<TH1D>(inf, Form("h1_mass_%s__pt-%d__y-%d", comp.c_str(), ipt, iy));
      };

      auto* hdata = get_mass_hist("data");
      auto* hmatch = get_mass_hist("mc-match");
      auto* hswap = get_mass_hist("mc-swap");
      auto* hkk = get_mass_hist("mc-kk");
      auto* hpipi = get_mass_hist("mc-pipi");
      if (!hdata || !hmatch || !hswap || !hkk || !hpipi) {
        __XJJLOG << "!! missing mass histogram for pt=" << ipt << ", y=" << iy << std::endl;
        return 4;
      }

      xjjroot::cached_dfitter cached_fit(fitopt);
      if (!cached_fit.prepare(hdata, hmatch, hswap, hkk, hpipi)) return 5;
      auto* truth_total = cached_fit.total_model(Form("closure_total__pt-%d__y-%d", ipt, iy));
      auto* truth_signal = cached_fit.signal_model(Form("closure_signal__pt-%d__y-%d", ipt, iy));
      nominal_yield = cached_fit.nominal_yield();
      nominal_yield_err = cached_fit.nominal_yield_err();
      // truth_yield = truth_signal->Integral(hdata->GetXaxis()->GetXmin(), hdata->GetXaxis()->GetXmax()) /
      //   hdata->GetBinWidth(1);
      pt = ipt;
      y = iy;

      for (itoy = 0; itoy < ntoys; ++itoy) {
        TH1D* htoy = nullptr;
        if (toytype == 0)
          htoy = make_toy_boostrap(hdata, rng, Form("%d_%d_%d", ipt, iy, itoy));
        else if (toytype == 1)
          htoy = make_toy_closure(hdata, truth_total, rng, Form("%d_%d_%d", ipt, iy, itoy));
        if (!htoy) continue;
        
        const auto fit = cached_fit.fit(htoy);
        yield = fit.fitted ? fit.yield : -1;
        yield_err = fit.fitted ? fit.yieldErr : -1;
        tree->Fill();
        delete htoy;
      }
    }
  }
  xjjc::progressbar_summary(bins.npt()*bins.ny());

  tree->Write();
  xjjroot::writehist(h3_bins);
  
  util::Writeinfo tinfo("info");
  tinfo.cast_branch("fitopt", fitopt);
  tinfo.cast_branch("inputname", inputname);
  tinfo.cast_branch("toytype", toytype);
  tinfo.cast_branch("toytype_tex", type_names[toytype]);
  tinfo.cast_branch("ntoys", ntoys);
  tinfo.cast_branch("seed", seed);
  tinfo.close();

  xjjroot::closefile(fout);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 7) return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]), argv[6]);
  if (argc == 6) return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]));
  if (argc == 5) return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]));

  return 1;
}
