#ifndef __BINS__
#define __BINS__

namespace bins {

#ifdef __BINS_PTY_ANA__
  std::vector<double> ybins; // to be defined in .cc
  const int npt = 1; const double minpt = 2, maxpt = 5;
#endif

#ifdef __BINS_PTY_EFF__
  const int ny = 48; const float miny = -2.4, maxy = 2.4;
  const int npt = 30; const float minpt = 2, maxpt = 5;
#endif

#ifdef __BINS_MULT__
  const int nmult = 50; const float minmult = 0, maxmult = 50;
#endif

#ifdef __BINS_DCA__
  const std::vector<float> ip3dbins = { 0, 0.001, 0.00227, 0.0038829, 0.00593128, 0.00853273, 0.0118366, 0.0160324, 0.0213612, 0.0281287, 0.0367235, 0.0476388, 0.0615013, 0.0791066, 0.101465 };
  // const std::vector<float> ip3dsigbins = {  };
#endif

#ifdef __BINS_MASS__
  const int nmass = 80; const float minmass = 1.66, maxmass = 2.06;
#endif

#ifdef __BINS_PTY_EQ__
  const int ny = 4; const float miny = -2, maxy = 2;
  const int npt = 1; const float minpt = 2, maxpt = 5;
#endif
  
}

#endif
