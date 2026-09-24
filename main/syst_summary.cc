#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "xjjstruct.h"

#include "draw.h"
#include "util.h"

#include "theory.h"

struct Property {
  std::string tex;
};

TGraphAsymmErrors* CombineErrors(
                                 const std::vector<TGraphAsymmErrors*>& gs)
{
  if (gs.empty())
    return nullptr;

  const int n = gs[0]->GetN();
  auto *gtotal = new TGraphAsymmErrors(n);
  for (int i = 0; i < n; ++i) {
    double x, y;
    gs[0]->GetPoint(i, x, y);
    double errLow2  = 0.0;
    double errHigh2 = 0.0;
    for (auto *g : gs) {
      double xi, yi;
      g->GetPoint(i, xi, yi);
      // Assuming same x and same central value
      errLow2  += std::pow(g->GetErrorYlow(i),  2);
      errHigh2 += std::pow(g->GetErrorYhigh(i), 2);
    }
    const double errLow  = std::sqrt(errLow2);
    const double errHigh = std::sqrt(errHigh2);
    gtotal->SetPoint(i, x, y);
    gtotal->SetPointError(
                          i,
                          gs[0]->GetErrorXlow(i),
                          gs[0]->GetErrorXhigh(i),
                          errLow,
                          errHigh
                          );
  }
  return gtotal;
}

const float ytop = 0.84, lspace = 1.2, tsize = 0.038;
float ypos(float i = 0, float margin = 0) { return ytop+margin-i*(tsize*lspace); }

int macro(const std::vector<std::string>& inputnames, const std::string& outputname, const std::string& tags) {
  __XJJLOG << ">> " << inputnames.size() << std::endl;
  const std::vector<Color_t> colors = { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["magenta"], xjjroot::mycolor_middle["cyan"], xjjroot::mycolor_middle["violet"] };
  // const auto alphas = xjjroot::grayscales_alpha(inputnames.size());

  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", xjjroot::parse_input(inputnames.front()).content.c_str()));
  draw::bintex tbins(h3_bins, 0, 2);
  const auto colors_y = xjjroot::grayscales_color(tbins.ny());

  const std::map<std::string, Property> texs = {
    { "gammaN", { .tex = "Xn0n (#gammaN)" } },
    { "NgammaRef", { .tex = "Xn0n (N#gamma) (#it{y} #rightarrow -#it{y})" } },
    { "sum", { .tex = "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" } }
  };
  TLegend *leg = nullptr;

  std::map<std::string, std::vector<TH1D*>> hs;
  std::map<std::string, std::vector<TGraphAsymmErrors*>> gs_syst;
  for (const auto& [key, prop] : texs) {
    hs[key].resize(tbins.npt(), nullptr); //
    for (int j=0; j<tbins.npt(); j++) {

      std::vector<TGraphAsymmErrors*> gs_err;
      for (int i=0; i<inputnames.size(); i++) {
        const auto inputp = xjjroot::parse_input(inputnames[i]);
        auto* inf = TFile::Open(inputp.content.c_str());
        const auto cc = colors[i%colors.size()];
        auto* dir = (TDirectory*)inf->Get(Form("dir__pt-%d", j));

        auto* h = xjjana::getobj<TH1D>(dir, "h1_y_xsec_" + key);
        xjjroot::setthgrstyle(h, kBlack, xjjroot::mstylelist_solid(i), 1.4, kBlack, 1, 1);
        if (!hs.at(key)[j]) hs.at(key)[j] = h; //
        auto* g_err = xjjana::getobj<TGraphAsymmErrors>(dir, "gr-err_y_xsec_" + key);
        xjjroot::setthgrstyle(g_err, cc, xjjroot::mstylelist_solid(i), 1.4, cc, 1, 1);
        gs_err.push_back(g_err); //
      }
      gs_syst[key].push_back(CombineErrors(gs_err));
      xjjroot::setthgrstyle(gs_syst.at(key)[j], kBlack, 20, 1.4, kBlack, 1, 1, 0, 0, 0);
    }
  } //

  auto drawh1_X0 = [](TH1D* h) {
    auto* g = xjjana::shifthistcenter(h, Form("gX0_%s", h->GetName()), 0);
    xjjroot::setthgrstyle(g, h->GetMarkerColor(), h->GetMarkerStyle(), h->GetMarkerSize(),
                          h->GetLineColor(), h->GetLineStyle(), h->GetLineWidth());
    g->Draw("pe1 same");
  };
  
  xjjroot::setgstyle(1);
  TLegend* leg_fonll = nullptr;
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (int j=0; j<tbins.npt(); j++) {
    fonll::DrawSets ds("../data/FONLL_all_predictions_vs_y_21Sept2026_CTEQ18_with_pt_2_5.root", tbins.edgelow_pt(j), tbins.edgeup_pt(j));
    const auto ymax = ds.ymaximum();
    for (const auto& [key, prop] : texs) {
      if (key == "sum")
        hs.at(key)[j]->SetMaximum(ymax*1.3);
      pdf->prepare();
      hs.at(key)[j]->Draw("axis");
      if (ds.valid() && key == "sum") ds.draw()->Draw();
      drawh1_X0(hs.at(key)[j]);
      gs_syst.at(key)[j]->Draw("5 same");
      // gs_syst.at(key)[j]->Draw("[] same");
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex, tbins.label_pt(j) }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { xjjroot::CMS::DzDzbar2 }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      gPad->RedrawAxis();
      pdf->write();
    }
  }

  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3]);
  }
  return 1;
}

// // clean which hists to draw
// void check_clean_consist(std::map<std::string, std::vector<TH1D*>> hs) {
//   for (const auto& [key, vh] : hs) {
//     int bad = 0;
//     for (auto& h : vh)
//       if (!h) {
//         bad++;
//         break;
//       }
//     if (bad) {
//       __XJJLOG << "!! not all input files have: " << key << ", erase it." << std::endl;
//       hs.erase(key);
//     }
//   }
// }
// check_clean_consist(hs);
