#pragma once

namespace global {
  struct variable
  {
    std::string varname;
    std::string var;
    std::string vartex;
    float varmin;
    float varmax;
    int nbin = 50;
    int logy = 0;
    int isbranch = 1;
  };
}

