#include "xjjanauti.h"

std::vector<xjjana::variable> vars = {
  { .varname = "ZDCsumPlus",           .var = "ZDCsumPlus",                  .vartex = "Offline ZDC Plus Energy [GeV]",                   .varmin = 0,   .varmax = 9000, .nbin = 1000, .logy = 1 },
  { .varname = "ZDCsumMinus",          .var = "ZDCsumMinus",                 .vartex = "Offline ZDC Minus Energy [GeV]",                  .varmin = 0,   .varmax = 9000, .nbin = 1000, .logy = 1 },
  { .varname = "nTrackInAcceptanceHP", .var = "nTrackInAcceptanceHP",        .vartex = "N_{trk} (highPurity, |#eta|<2.4, p_{T}>0.5 GeV)", .varmin = 0,   .varmax = 50,   .nbin = 50,   .logy = 1 },
  { .varname = "HFEMaxPlusforest",     .var = "HFEMaxPlus_forest",           .vartex = "Leading HF+ PF (3<|#eta|<6) [GeV]",               .varmin = 0,   .varmax = 50,   .nbin = 100,  .logy = 0 },
  { .varname = "HFEMaxPlusforest-Low", .var = "HFEMaxPlus_forest",           .vartex = "Leading HF+ PF (3<|#eta|<6) [GeV]",               .varmin = 0,   .varmax = 10,   .nbin = 100,  .logy = 0 },
  { .varname = "HFEMaxMinusforest",    .var = "HFEMaxMinus_forest",          .vartex = "Leading HF- PF (3<|#eta|<6) [GeV]",               .varmin = 0,   .varmax = 50,   .nbin = 100,  .logy = 0 },
  //
  { .varname = "Dmass",                .var = "Dmass",                       .vartex = "m_{#piK} [GeV] ",                                 .varmin = 1.7, .varmax = 2.0,  .nbin = 60,   .logy = 0 },
  { .varname = "Dalpha",               .var = "Dalpha",                      .vartex = "#alpha (3D pointing angle)",                      .varmin = 0,   .varmax = 3.2,  .nbin = 64,   .logy = 0 },
  { .varname = "Dalpha-Low",           .var = "Dalpha",                      .vartex = "#alpha (3D pointing angle)",                      .varmin = 0,   .varmax = 0.5,  .nbin = 50,   .logy = 0 },
  { .varname = "Ddls",                 .var = "DsvpvDistance/DsvpvDisErr",   .vartex = "3D decay length significance",                    .varmin = 0,   .varmax = 20,   .nbin = 50,   .logy = 1 },
  { .varname = "Dchi2cl",              .var = "Dchi2cl",                     .vartex = "Secondary vertex prob",                           .varmin = 0,   .varmax = 1,    .nbin = 50,   .logy = 0 },
  { .varname = "Dtrk1Pt",              .var = "Dtrk1Pt",                     .vartex = "Track 1 p_{T} [GeV]",                             .varmin = 0,   .varmax = 5,    .nbin = 50,   .logy = 1 },
  { .varname = "Dtrk2Pt",              .var = "Dtrk2Pt",                     .vartex = "Track 2 p_{T} [GeV]",                             .varmin = 0,   .varmax = 5,    .nbin = 50,   .logy = 1 },
  { .varname = "Dtrk1ptrel",           .var = "Dtrk1PtErr/Dtrk1Pt",          .vartex = "Track 1 #sigma(p_{T})/p_{T}",                     .varmin = 0,   .varmax = 0.4,  .nbin = 40,   .logy = 1 },
  { .varname = "Dtrk2ptrel",           .var = "Dtrk2PtErr/Dtrk2Pt",          .vartex = "Track 2 #sigma(p_{T})/p_{T}",                     .varmin = 0,   .varmax = 0.4,  .nbin = 40,   .logy = 1 },
  { .varname = "Dtrk1nhit",            .var = "Dtrk1PixelHit+Dtrk1StripHit", .vartex = "Track 1 number of hits",                          .varmin = 2,   .varmax = 32,   .nbin = 30,   .logy = 0 },
  { .varname = "Dtrk2nhit",            .var = "Dtrk2PixelHit+Dtrk2StripHit", .vartex = "Track 2 number of hits",                          .varmin = 2,   .varmax = 32,   .nbin = 30,   .logy = 0 },
};
