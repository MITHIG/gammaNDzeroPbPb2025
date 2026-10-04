#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"
#include "dfitter.h"
#include "draw.h"

namespace {
  TH1D* get_mass_hist(TFile* file, const std::string& component, int ipt, int iy) {
    return static_cast<TH1D*>(file->Get(Form("h1_mass_%s__pt-%d__y-%d",
                                             component.c_str(), ipt, iy)));
  }

  TH1D* make_toy(const TH1D* source, TRandom3& rng, int toy, int ipt, int iy) {
    auto* result = static_cast<TH1D*>(source->Clone(Form("toy_mass__toy-%d__pt-%d__y-%d",
                                                         toy, ipt, iy)));
    result->SetDirectory(nullptr);
    for (int ibin = 1; ibin <= source->GetNbinsX(); ++ibin) {
      const auto count = rng.Poisson(std::max(0., source->GetBinContent(ibin)));
      result->SetBinContent(ibin, count);
      result->SetBinError(ibin, std::sqrt(static_cast<double>(count)));
    }
    return result;
  }
}

int macro(const std::string& inputname, const std::string& outputname,
          const std::string& fit_opt = "3G-Peaky", int ntoys = 1000,
          unsigned int seed = 12345) {
  if (ntoys <= 0) return 1;

  gROOT->SetBatch(true);
  xjjroot::setgstyle(1);
  TCanvas canvas("toy_canvas", "toy_canvas");

  auto* input = TFile::Open(inputname.c_str(), "READ");
  if (!input || input->IsZombie()) {
    __XJJLOG << "!! unable to open input ROOT file: " << inputname << std::endl;
    return 2;
  }

  auto* h3_bins = static_cast<TH3D*>(input->Get("h3_bins"));
  if (!h3_bins) {
    __XJJLOG << "!! missing h3_bins in " << inputname << std::endl;
    return 3;
  }
  const draw::bintex bins(h3_bins, 0, 2);
  const auto fitopt = xjjroot::parse_input(fit_opt);

  auto* fout = xjjroot::newfile("rootfiles/" + outputname + ".root");
  TRandom3 rng(seed);

  int toy = -1;
  int pt = -1, y = -1;
  double yield = -1, yield_err = -1;
  double nominal_yield = -1, nominal_yield_err = -1;
  auto* tree = new TTree("toys", "Poisson pseudo-data fit results");
  tree->Branch("toy", &toy);
  tree->Branch("pt", &pt);
  tree->Branch("y", &y);
  tree->Branch("yield", &yield);
  tree->Branch("yield_err", &yield_err);
  tree->Branch("nominal_yield", &nominal_yield);
  tree->Branch("nominal_yield_err", &nominal_yield_err);

  for (int ipt = 0; ipt < bins.npt(); ++ipt) {
    for (int iy = 0; iy < bins.ny(); ++iy) {
      auto* hdata = get_mass_hist(input, "data", ipt, iy);
      auto* hmatch = get_mass_hist(input, "mc-match", ipt, iy);
      auto* hswap = get_mass_hist(input, "mc-swap", ipt, iy);
      auto* hkk = get_mass_hist(input, "mc-kk", ipt, iy);
      auto* hpipi = get_mass_hist(input, "mc-pipi", ipt, iy);
      if (!hdata || !hmatch || !hswap || !hkk || !hpipi) {
        __XJJLOG << "!! missing mass histogram for pt=" << ipt << ", y=" << iy << std::endl;
        return 4;
      }

      xjjroot::dfitter nominal_fit(fitopt.content.c_str());
      nominal_fit.fit(hdata, hmatch, hswap, hkk, hpipi);
      nominal_yield = nominal_fit.fitted() ? nominal_fit.yield() : -1;
      nominal_yield_err = nominal_fit.fitted() ? nominal_fit.yieldErr() : -1;
      pt = ipt;
      y = iy;

      for (toy = 0; toy < ntoys; ++toy) {
        auto* htoy = make_toy(hdata, rng, toy, ipt, iy);
        xjjroot::dfitter fit(fitopt.content.c_str());
        fit.fit(htoy, hmatch, hswap, hkk, hpipi);
        yield = fit.fitted() ? fit.yield() : -1;
        yield_err = fit.fitted() ? fit.yieldErr() : -1;
        tree->Fill();
        delete htoy;
      }
      __XJJLOG << "++ completed pt=" << ipt << ", y=" << iy
               << " (" << ntoys << " toys)" << std::endl;
    }
  }

  tree->Write();
  auto* h3_bins_out = static_cast<TH3D*>(h3_bins->Clone("h3_bins"));
  h3_bins_out->Write();
  auto fitopt_content = fitopt.content;
  auto input_content = inputname;
  auto* info = new TTree("info", "toy study configuration");
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
  const auto ntoys = argc > 4 ? std::atoi(argv[4]) : 1000;
  const auto seed = argc > 5 ? static_cast<unsigned int>(std::atoi(argv[5])) : 12345u;
  return macro(argv[1], argv[2], fitopt, ntoys, seed);
}
