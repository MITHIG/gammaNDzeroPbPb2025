#include "xjjanauti.h"

#ifdef __VARIABLES_ROOINCL__

std::vector<xjjana::variable> variables = {
  //
  { .varname = "Dmass", .var = "Dmass", .vartex = "m_{#piK} [GeV]", .varmin = 1.66, .varmax = 2.06, .nbin = 60, .logy = 0 },
  { .varname = "Dpt", .var = "Dpt", .vartex = "#it{p}_{T} [GeV]", .varmin = 2., .varmax = 5., .nbin = 80, .logy = 0 },
  { .varname = "Dy", .var = "Dy", .vartex = "y", .varmin = 2., .varmax = 5., .nbin = 80, .logy = 0 },
  { .varname = "Dmva_BDT", .var = "Dmva_BDT", .vartex = "BDT", .varmin = -0.3, .varmax = 0.4, .nbin = 70, .logy = 0 },
  { .varname = "Dalpha", .var = "Dalpha", .vartex = "#alpha (3D pointing angle)", .varmin = 0, .varmax = 3.2, .nbin = 64, .logy = 0 },
  { .varname = "DsvpvDistance", .var = "DsvpvDistance", .vartex = "3D decay length", .varmin = 0, .varmax = 20, .nbin = 50, .logy = 1 },
  { .varname = "DsvpvDisErr", .var = "DsvpvDisErr", .vartex = "3D decay length error", .varmin = 0, .varmax = 20, .nbin = 50, .logy = 1 },
  { .varname = "Dchi2cl", .var = "Dchi2cl", .vartex = "Secondary vertex prob", .varmin = 0, .varmax = 1, .nbin = 50, .logy = 0 },
  { .varname = "Dip3D", .var = "Dip3D", .vartex = "DCA [cm]", .varmin = 0, .varmax = 0.1, .nbin = 20, .logy = 1 },
  { .varname = "Dip3derr", .var = "Dip3derr", .vartex = "#sigma(DCA) [cm]", .varmin = 0, .varmax = 0.1, .nbin = 20, .logy = 1 },
  { .varname = "Dtrk1Pt", .var = "Dtrk1Pt", .vartex = "Track 1#scale[0.5]{ }#it{p}_{T} [GeV]", .varmin = 0, .varmax = 5, .nbin = 50, .logy = 1 },
  { .varname = "Dtrk2Pt", .var = "Dtrk2Pt", .vartex = "Track 2#scale[0.5]{ }#it{p}_{T} [GeV]", .varmin = 0, .varmax = 5, .nbin = 50, .logy = 1 },
  { .varname = "Dtrk1PtErr", .var = "Dtrk1PtErr", .vartex = "Track 1 #sigma(p_{T}) [GeV]", .varmin = 0, .varmax = 1, .nbin = 50, .logy = 1 },
  { .varname = "Dtrk2PtErr", .var = "Dtrk2PtErr", .vartex = "Track 2 #sigma(p_{T}) [GeV]", .varmin = 0, .varmax = 1, .nbin = 50, .logy = 1 },
  { .varname = "Dtrk1Eta", .var = "Dtrk1Eta", .vartex = "Track 1#scale[0.5]{ }#eta", .varmin = -2.4, .varmax = 2.4, .nbin = 48, .logy = 0 },
  { .varname = "Dtrk2Eta", .var = "Dtrk2Eta", .vartex = "Track 2#scale[0.5]{ }#eta", .varmin = -2.4, .varmax = 2.4, .nbin = 48, .logy = 0 },
};

#endif
