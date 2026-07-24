#pragma once

#include "../include/variable.h"

#ifdef __VARIABLES_ROOSPLOT__

std::vector<global::variable> variables = {
  //
  { .varname = "Dmass", .var = "Dmass", .vartex = "m_{#piK} [GeV/c^{2}]", .varmin = 1.66, .varmax = 2.06, .nbin = 80, .logy = 0, .isbranch = 1 },
  { .varname = "Dpt", .var = "Dpt", .vartex = "#it{p}_{T} [GeV/c]", .varmin = 2., .varmax = 5., .nbin = 10, .logy = 0, .isbranch = 1 },
  { .varname = "Dy", .var = "Dy", .vartex = "y", .varmin = -2., .varmax = 2., .nbin = 20, .logy = 0, .isbranch = 1 },
  { .varname = "Dmva_BDT", .var = "Dmva_BDT", .vartex = "BDT", .varmin = 0, .varmax = 0.4, .nbin = 20, .logy = 0, .isbranch = 1 },
  { .varname = "Dalpha", .var = "Dalpha", .vartex = "#alpha (3D pointing angle)", .varmin = 0, .varmax = 0.6, .nbin = 15, .logy = 0, .isbranch = 1 },
  // { .varname = "Dalpha", .var = "Dalpha", .vartex = "#alpha (3D pointing angle)", .varmin = 0, .varmax = 3.2, .nbin = 64, .logy = 0, .isbranch = 1 },
  // { .varname = "Dalpha_zoom", .var = "Dalpha", .vartex = "#alpha (3D pointing angle)", .varmin = 0, .varmax = 0.6, .nbin = 60, .logy = 0, .isbranch = -1 },
  { .varname = "DsvpvDistance", .var = "DsvpvDistance", .vartex = "3D decay length [cm]", .varmin = 0, .varmax = 1., .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1 } },    
  { .varname = "DsvpvDisErr", .var = "DsvpvDisErr", .vartex = "3D decay length error [cm]", .varmin = 0, .varmax = 0.2, .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1 } },
  { .varname = "Ddls", .var = "DsvpvDistance/DsvpvDisErr", .vartex = "3D decay length significance", .varmin = 0, .varmax = 30, .nbin = 30, .logy = 1, .isbranch = 0, //
    .bins = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 13, 17, 22, 30 } },
  { .varname = "DsvpvDistance_2D", .var = "DsvpvDistance_2D", .vartex = "2D decay length [cm]", .varmin = 0, .varmax = 1., .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1 } },
  { .varname = "DsvpvDisErr_2D", .var = "DsvpvDisErr_2D", .vartex = "2D decay length error [cm]", .varmin = 0, .varmax = 0.2, .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1 } },
  { .varname = "Ddls_2D", .var = "DsvpvDistance_2D/DsvpvDisErr_2D", .vartex = "2D decay length significance", .varmin = 0, .varmax = 30, .nbin = 30, .logy = 1, .isbranch = 0, //
    .bins = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 13, 17, 22, 30 } },
  { .varname = "Dchi2cl", .var = "Dchi2cl", .vartex = "Secondary vertex prob", .varmin = 0, .varmax = 1, .nbin = 10, .logy = 0, .isbranch = 1 },
  { .varname = "Dip3D", .var = "Dip3D", .vartex = "DCA [cm]", .varmin = 0, .varmax = 0.1, .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.001, 0.00227, 0.0038829, 0.00593128, 0.00853273, 0.0118366, 0.0160324, 0.0213612, 0.0281287, 0.0367235, 0.0476388, 0.0615013, 0.0791066, 0.1 } },
  { .varname = "Dip3derr", .var = "Dip3derr", .vartex = "#sigma(DCA) [cm]", .varmin = 0, .varmax = 0.1, .nbin = 20, .logy = 1, .isbranch = 1 },
  { .varname = "Dip3Dsig", .var = "Dip3D/Dip3derr", .vartex = "DCA /#scale[0.4]{ }#sigma_{DCA}", .varmin = 0, .varmax = 25, .nbin = 25*5, .logy = 1, .isbranch = 0, // 
    .bins = { 0, 0.2, 0.5, 1., 2, 3, 4, 5, 7, 9, 11, 14, 18, 25 } },
  { .varname = "Dtrk1Pt", .var = "Dtrk1Pt", .vartex = "Track 1#scale[0.5]{ }#it{p}_{T} [GeV/c]", .varmin = 0, .varmax = 5, .nbin = 25, .logy = 0, .isbranch = 1 },
  { .varname = "Dtrk2Pt", .var = "Dtrk2Pt", .vartex = "Track 2#scale[0.5]{ }#it{p}_{T} [GeV/c]", .varmin = 0, .varmax = 5, .nbin = 25, .logy = 0, .isbranch = 1 },
  { .varname = "Dtrk1PtErr", .var = "Dtrk1PtErr", .vartex = "Track 1#scale[0.5]{ }#sigma(p_{T}) [GeV/c]", .varmin = 0, .varmax = 0.4, .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.15, 0.25 } },
  { .varname = "Dtrk2PtErr", .var = "Dtrk2PtErr", .vartex = "Track 2#scale[0.5]{ }#sigma(p_{T}) [GeV/c]", .varmin = 0, .varmax = 0.4, .nbin = 20, .logy = 1, .isbranch = 1,
    .bins = { 0, 0.01, 0.02, 0.03, 0.04, 0.06, 0.09, 0.15, 0.25 } },
  { .varname = "Dtrk1Eta", .var = "Dtrk1Eta", .vartex = "Track 1#scale[0.5]{ }#eta", .varmin = -2.4, .varmax = 2.4, .nbin = 24, .logy = 0, .isbranch = 1 },
  { .varname = "Dtrk2Eta", .var = "Dtrk2Eta", .vartex = "Track 2#scale[0.5]{ }#eta", .varmin = -2.4, .varmax = 2.4, .nbin = 24, .logy = 0, .isbranch = 1 },
  // 
};

#endif

#ifdef __SPECTATORS__

std::vector<global::variable> spectators = {
  { .varname = "isL1ZDCOr", .var = "isL1ZDCOr", .varmin = 0, .varmax = 2 },
  { .varname = "cscTightHalo2015Filter", .var = "cscTightHalo2015Filter", .varmin = 0, .varmax = 2 },
  { .varname = "selectedVtxFilter", .var = "selectedVtxFilter", .varmin = 0, .varmax = 2 },
  { .varname = "ZDCgammaN", .var = "ZDCgammaN", .varmin = 0, .varmax = 2 },
  { .varname = "ZDCNgamma", .var = "ZDCNgamma", .varmin = 0, .varmax = 2 },
  { .varname = "HFEMaxPlus_eta5", .var = "HFEMaxPlus_eta5", .varmin = 0, .varmax = 100 },
  { .varname = "HFEMaxMinus_eta5", .var = "HFEMaxMinus_eta5", .varmin = 0, .varmax = 100 },
  { .varname = "Dgen", .var = "Dgen", .varmin = -1, .varmax = 100000 },
  //
  // { .varname = "", .var = "", .varmin = 0, .varmax = 2 },
};

#endif

global::variable var_by_name(const std::string& name, const std::vector<global::variable>& vars = variables) {
  global::variable result;
  for (auto v : vars) {
    if (v.varname == name) {
      result = v;
      break;
    }
  }
  if (result.varname.empty()) {
    __XJJLOG << "!! bad variable: " << name << std::endl;
  }
  return result;
}

global::variable var_by_exp(const std::string& name, const std::vector<global::variable>& vars = variables) {
  global::variable result;
  for (auto v : vars) {
    if (v.var == name) {
      result = v;
      break;
    }
  }
  if (result.varname.empty()) {
    __XJJLOG << "!! bad variable: " << name << std::endl;
  }
  return result;
}

