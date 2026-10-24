#pragma once

#include "../include/variable.h"

std::vector<global::variable> variables = {
  { .varname = "ZDCsumPlus", .var = "ZDCsumPlus", .vartex = "Plus ZDC E_{sum} [GeV]", .varmin = 0, .varmax = 9000, .nbin = 450, .logy = 1 },
  { .varname = "ZDCsumMinus", .var = "ZDCsumMinus", .vartex = "Minus ZDC E_{sum} [GeV]", .varmin = 0, .varmax = 9000, .nbin = 450, .logy = 1 },
  // { .varname = "ZDCsumGamma", .var = "ZDCsumPlus,ZDCsumMinus", .vartex = "#gamma-going ZDC E_{sum} [GeV]", .varmin = 0, .varmax = 3000, .nbin = 500, .logy = 1, .note = "Gap + ZDC^{oppo} Xn" },
  { .varname = "ZDCsumN", .var = "ZDCsumPlus,ZDCsumMinus", .vartex = "N-going ZDC E_{sum} [GeV]", .varmin = 0, .varmax = 9000, .nbin = 450, .logy = 1, .note = "ZDC^{oppo} 0n" },
  { .varname = "ZDCsumN-zoom", .var = "ZDCsumPlus,ZDCsumMinus", .vartex = "N-going ZDC E_{sum} [GeV]", .varmin = 900, .varmax = 9000, .nbin = 405, .logy = 1, .note = "ZDC^{oppo} 0n" },
  { .varname = "HFEMaxPlus_eta5", .var = "HFEMaxPlus_eta5", .vartex = "Plus HF+ PF E_{max}^{3 < |#eta| < 5} [GeV]", .varmin = 0, .varmax = 40, .nbin = 100, .logy = 1 },
  { .varname = "HFEMaxMinus_eta5", .var = "HFEMaxMinus_eta5", .vartex = "Minus HF- PF E_{max}^{3 < |#eta| < 5} [GeV]", .varmin = 0, .varmax = 40, .nbin = 100, .logy = 1 },
  { .varname = "HFEMaxGamma_eta5", .var = "HFEMaxPlus_eta5,HFEMaxMinus_eta5", .vartex = "0n (#gamma-going) side HF PF E_{max}^{3 < |#eta| < 5} [GeV]", .varmin = 0, .varmax = 40, .nbin = 100, .logy = 1, .note = "ZDC 1nXOR" },
  // { .varname = "HFEMaxN_eta5", .var = "HFEMaxPlus_eta5,HFEMaxMinus_eta5", .vartex = "Xn (N-going) side HF PF E_{max}^{3 < |#eta| < 5} [GeV]", .varmin = 0, .varmax = 40, .nbin = 100, .logy = 1, .note = "ZDC 1nXOR" },
  { .varname = "HFEMaxN_eta5-gap", .var = "HFEMaxPlus_eta5,HFEMaxMinus_eta5", .vartex = "Xn (N-going) side HF PF E_{max}^{3 < |#eta| < 5} [GeV]", .varmin = 0, .varmax = 40, .nbin = 100, .logy = 1, .note = "ZDC 1nXOR + Gap" },
  { .varname = "nTrackInAcceptanceHP", .var = "nTrackInAcceptanceHP", .vartex = "N_{trk} (highPurity, |#eta| < 2.4,#scale[0.4]{ }#it{p}_{T} > 0.5 GeV)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR" },
  { .varname = "nTrackInAcceptanceHP-gap", .var = "nTrackInAcceptanceHP", .vartex = "N_{trk} (highPurity, |#eta| < 2.4,#scale[0.4]{ }#it{p}_{T} > 0.5 GeV)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR + Gap" },
  { .varname = "nTrackInAcceptanceHP-nogap", .var = "nTrackInAcceptanceHP", .vartex = "N_{trk} (highPurity, |#eta| < 2.4,#scale[0.4]{ }#it{p}_{T} > 0.5 GeV)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR + No gap" },
  { .varname = "nTrackInAcceptanceHP-gap-low", .var = "nTrackInAcceptanceHP", .vartex = "N_{trk} (highPurity, |#eta| < 2.4,#scale[0.4]{ }#it{p}_{T} > 0.5 GeV)", .varmin = 0, .varmax = 80, .nbin = 80, .logy = 0, .note = "ZDC 1nXOR + Gap" },
  { .varname = "nTrkVtx", .var = "nTrkVtx", .vartex = "N_{trk} asso to PV (weight > 0.5)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR" },
  { .varname = "nTrkVtx-gap", .var = "nTrkVtx", .vartex = "N_{trk} asso to PV (weight > 0.5)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR + Gap" },
  { .varname = "nTrkVtx-nogap", .var = "nTrkVtx", .vartex = "N_{trk} asso to PV (weight > 0.5)", .varmin = 0, .varmax = 250, .nbin = 125, .logy = 1, .note = "ZDC 1nXOR + No gap" },
  { .varname = "nVtx", .var = "nVtx", .vartex = "Number of vertices", .varmin = 0, .varmax = 10, .nbin = 10, .logy = 1, .note = "ZDC 1nXOR + Gap" },
};

global::variable var_by_name(const std::string& name) {
  global::variable r;
  for (const auto& v : variables) {
    if (v.varname == name)
      return v;
  }
  return r;
}
