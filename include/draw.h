#include "xjjcuti.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace draw {
  class bintex {
  public:
    explicit bintex(TH1* h, int xyz_y = 0,  int xyz_pt = 0)
      : h_bins(h), axis_pt(decide_axis(xyz_pt)), axis_y(decide_axis(xyz_y)) { h_dump.resize(3, nullptr); }

    std::string label_y(int i = -1) const {
      const auto ymin = (i >= 0 ? axis_y->GetBinLowEdge(i+1) : axis_y->GetBinLowEdge(1));
      const auto ymax = (i >= 0 ? axis_y->GetBinUpEdge(i+1) : axis_y->GetBinUpEdge(axis_y->GetNbins()));
      return xjjc::number_range_string(ymin, ymax, "y", -1.e1);
    }
    
    std::string label_pt(int i = -1) const {
      const auto ptmin = (i >= 0 ? axis_pt->GetBinLowEdge(i+1) : axis_pt->GetBinLowEdge(1));
      const auto ptmax = (i >= 0 ? axis_pt->GetBinUpEdge(i+1) : axis_pt->GetBinUpEdge(axis_pt->GetNbins()));
      return xjjc::number_range_string(ptmin, ptmax, "#it{p}_{T}") + " GeV";
    }

    template<class T> T* make_h1_y(std::string name) {
      auto* h = new T(name.c_str(), Form(";%s;", axis_y->GetTitle()), axis_y->GetNbins(), axis_y->GetXbins()->GetArray());
      return h;
    }
    
  private:
    TH1* h_bins;
    std::vector<TH1D*> h_dump;
    const TAxis *axis_pt, *axis_y;
    const TAxis* decide_axis(int xyz) {
      const TAxis* axis = nullptr;
      if (xyz == 0) {
        axis = h_bins->GetXaxis();
      }
      else if (xyz == 1) {
        if (!dynamic_cast<TH2*>(h_bins))
          throw std::runtime_error("Y axis requested but histogram is not TH2/TH3");
        axis = h_bins->GetYaxis();
      }
      else if (xyz == 2) {
        if (!dynamic_cast<TH3*>(h_bins))
          throw std::runtime_error("Z axis requested but histogram is not TH3");
        axis = h_bins->GetZaxis();
      }
      else {
        throw std::runtime_error("xyz must be 0, 1, or 2");
      }
      return axis;
    }
  };

}
