#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "xjjstruct.h"

#include "draw.h"
#include "util.h"

#include "theory.h"

struct Property {
  std::string tex;
};

TGraphAsymmErrors* combine_errors(const std::vector<TGraphAsymmErrors*>& gs) {
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
TGraphAsymmErrors* get_rel_unc(const TGraphAsymmErrors* g);

const float ytop = 0.84, lspace = 1.2, tsize = 0.038;
float ypos(float i = 0, float margin = 0) { return ytop+margin-i*(tsize*lspace); }

int macro(const std::vector<std::string>& inputnames, const std::string& inputprompt,
          const std::string& outputname, int save_png) {
  __XJJLOG << ">> " << inputnames.size() << std::endl;
  const std::vector<Color_t> colors = { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["magenta"], xjjroot::mycolor_middle["cyan"], xjjroot::mycolor_middle["violet"] };
  // const auto alphas = xjjroot::grayscales_alpha(inputnames.size());

  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", xjjroot::parse_input(inputnames.front()).content.c_str()));
  draw::bintex tbins(h3_bins, 0, 2);
  const auto colors_y = xjjroot::grayscales_color(tbins.ny());

  const auto info = util::read_info(TFile::Open(xjjroot::parse_input(inputnames.front()).content.c_str()), "info", true);
  const auto lumi = std::atof(info.at("lumi").c_str());

  const std::map<std::string, Property> texs = {
    { "gammaN", { .tex = "Xn0n (#gammaN)" } },
    { "NgammaRef", { .tex = "Xn0n (N#gamma) (#it{y} #rightarrow -#it{y})" } },
    { "sum", { .tex = "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" } },
    { "sum-prompt", { .tex = "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" } }
  };
  TLegend *leg = nullptr;

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });

  std::vector<std::map<std::string, TGraphAsymmErrors*>> gs_syst(tbins.npt()); // [npt][type]
  std::vector<std::map<std::string, TH1D*>> hs(tbins.npt());
  for (int j=0; j<tbins.npt(); j++) {
    for (const std::string key : { "sum", "gammaN", "NgammaRef" }) {
      // hs[key].resize(tbins.npt(), nullptr);
      std::vector<TGraphAsymmErrors*> gs_err, gs_rel;
      auto* leg = new TLegend(0.45, ypos(2. + inputnames.size()), 0.45+0.3, ypos(2.));
      xjjroot::setleg(leg, tsize);
      for (int i=0; i<inputnames.size(); i++) {
        const auto inputp = xjjroot::parse_input(inputnames[i]);
        auto* inf = TFile::Open(inputp.content.c_str());
        const auto cc = colors[i%colors.size()];
        auto* dir = (TDirectory*)inf->Get(Form("dir__pt-%d", j));

        auto* g_err = xjjana::getobj<TGraphAsymmErrors>(dir, "gr-err_y_xsec_" + key);
        xjjroot::setthgrstyle(g_err, cc, xjjroot::mstylelist_solid(i), 1.4, cc, 1, 1);
        gs_err.push_back(g_err); //
        gs_rel.push_back(get_rel_unc(g_err)); // only for drawing
        xjjroot::addentrybystyle(leg, inputp.tex.c_str(), "f", { .lcolor = 0, .lwidth = 0, .fcolor = colors[i], .falpha = 0.5, .fstyle = 1001 });

        if (i == 0) {
          auto* h = xjjana::getobj<TH1D>(dir, "h1_y_xsec_" + key);
          hs[j][key] = h;
        }
      }
      gs_syst[j][key] = combine_errors(gs_err);
      xjjroot::setthgrstyle(gs_syst[j][key], kBlack, 20, 1.4, kBlack, 1, 1, 0, 0, 0);
      auto* gr_rel_total = get_rel_unc(gs_syst[j][key]); // only for drawing
      xjjroot::setthgrstyle(gr_rel_total, 0, 0, 0, kBlack, 1, 1, 0, 0, 0);

      // draw uncertainty
      if (gs_rel.empty()) continue;
      auto* hempty = new TH2D("hempty", Form(";%s;Relative uncertainty", hs[j][key]->GetXaxis()->GetTitle()),
                              10, hs[j][key]->GetXaxis()->GetXmin(), hs[j][key]->GetXaxis()->GetXmax(),
                              10, xjjana::gethwerrminimum(gr_rel_total)*1.1, xjjana::gethwerrmaximum(gr_rel_total)*2);
      xjjroot::sethempty(hempty, 0, 0.2);
      pdf->prepare();
      hempty->Draw("axis");
      xjjroot::drawline_horizon(0, hempty);
      for (int k=0; k<hs[j][key]->GetNbinsX(); k++) {
        const auto xleft = hs[j][key]->GetBinCenter(k+1) - hs[j][key]->GetBinWidth(k+1) * 0.5,
          step = hs[j][key]->GetBinWidth(k+1) / gs_rel.size();
        for (int i=0; i<gs_rel.size(); i++) {
          xjjroot::drawbox(xleft + i*step, 0 - gs_rel[i]->GetErrorYlow(k),
                           xleft + (i+1)*step, gs_rel[i]->GetErrorYhigh(k),
                           colors[i], 0.5);
        }
      }
      gr_rel_total->Draw("5 same");
      leg->Draw();
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { texs.at(key).tex, tbins.label_pt(j) }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { xjjroot::CMS::DzDzbar2 }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, Form("%.1f #mub^{-1} (PbPb 5.36 TeV)", lumi));
      pdf->write(Form("%s_relerr_%s_pt-%d.pdf", name_png.c_str(), key.c_str(), j), save_png && key == "sum" ? "" : "X");
      delete hempty;
    }

    auto* inf0 = TFile::Open(inputprompt.c_str());
    auto* dir0 = (TDirectory*)inf0->Get(Form("dir__pt-%d", j));
    auto* h1_xsec_prompt = xjjana::getobj<TH1D>(dir0, "h1_y_xsec_sum-prompt");
    if (h1_xsec_prompt) {
      hs[j]["sum-prompt"] = h1_xsec_prompt;
      auto* g_syst_prompt = (TGraphAsymmErrors*)gs_syst[j]["sum"]->Clone(xjjc::str_replaceall(gs_syst[j]["sum"]->GetName(), "sum", "sum-prompt").c_str());
      for (int i=0; i<g_syst_prompt->GetN(); i++) {
        const auto scale = hs[j]["sum-prompt"]->GetBinContent(i+1) / hs[j]["sum"]->GetBinContent(i+1);
        g_syst_prompt->GetY()[i] *= scale;
        g_syst_prompt->GetEYhigh()[i] *= scale;
        g_syst_prompt->GetEYlow()[i] *= scale;
      }
      gs_syst[j]["sum-prompt"] = g_syst_prompt;
      auto* gr_prompt_syst = xjjana::getobj<TGraphErrors>(dir0, "gr_y_xsec_sum-prompt");
      auto* ga = new TGraphAsymmErrors(gr_prompt_syst->GetN());
      for (int i = 0; i < gr_prompt_syst->GetN(); ++i) {
        ga->SetPoint(i, gr_prompt_syst->GetPointX(i), gr_prompt_syst->GetPointY(i));
        ga->SetPointError(i, gr_prompt_syst->GetErrorX(i),
                          gr_prompt_syst->GetErrorX(i),
                          gr_prompt_syst->GetErrorY(i),
                          gr_prompt_syst->GetErrorY(i));
      }
      gs_syst[j]["sum-prompt-syst"] = ga;
      // gs_syst[j]["sum-prompt-syst"] = new TGraphAsymmErrors(*gr_prompt_syst);
      xjjroot::setthgrstyle(gs_syst[j]["sum-prompt-syst"], -1, -1, -1, kBlack, 1, 1);
    }
    for (auto& [_, h] : hs[j])
      xjjroot::setthgrstyle(h, kBlack, 20, 1.5, kBlack, 1, 1);
  } //

  auto drawh1_X0 = [](TH1D* h) {
    auto* g = xjjana::shifthistcenter(h, Form("gX0_%s", h->GetName()), 0);
    xjjroot::setthgrstyle(g, h->GetMarkerColor(), h->GetMarkerStyle(), h->GetMarkerSize(),
                          h->GetLineColor(), h->GetLineStyle(), h->GetLineWidth());
    g->Draw("pe1 same");
  };
  
  TLegend* leg_fonll = nullptr;
  TLegend* leg_data = new TLegend(0.22, 0.76-0.038, 0.40, 0.76);
  xjjroot::setleg(leg_data, 0.035);
  leg_data->AddEntry(hs.front()["sum"], "Data", "p");

  for (int j=0; j<tbins.npt(); j++) {
    fonll::DrawSets ds("../theory/FONLL_all_predictions_vs_y_21Sept2026_CTEQ18_with_pt_2_5.root", tbins.edgelow_pt(j), tbins.edgeup_pt(j));
    const auto ymax = ds.ymaximum();

    for (const auto key : { "sum", "sum-prompt", "gammaN", "NgammaRef" } ) {
      if (hs[j].find(key) == hs[j].end()) continue;
      
      // auto prop = texs[key];
      const auto draw_theory = key == "sum" || key == "sum-prompt";
      if (draw_theory)
        hs[j][key]->SetMaximum(ymax*1.3);
      pdf->prepare();
      hs[j][key]->Draw("axis");
      if (ds.valid() && draw_theory) ds.draw();
      drawh1_X0(hs[j][key]);
      gs_syst[j][key]->Draw("5 same");
      if (key == "sum-prompt")
        gs_syst[j]["sum-prompt-syst"]->Draw("[] same");
      // gs_syst[j][key]->Draw("[] same");
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { texs.at(key).tex, tbins.label_pt(j) }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { xjjroot::CMS::DzDzbar2 }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, Form("%.1f #mub^{-1} (PbPb 5.36 TeV)", lumi));
      leg_data->Draw();
      gPad->RedrawAxis();
      pdf->write(Form("%s_%s_pt-%d.pdf", name_png.c_str(), key, j), save_png && draw_theory ? "" : "X");
    }
  }

  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], std::atoi(argv[4]));
  }
  return 1;
}

TGraphAsymmErrors* get_rel_unc(const TGraphAsymmErrors* g) {
  auto* gr = new TGraphAsymmErrors(g->GetN());
  gr->SetName(Form("%s_relerr", gr->GetName()));

  for (int i = 0; i < g->GetN(); ++i) {
    const double x = g->GetPointX(i);
    const double y = g->GetPointY(i);

    gr->SetPoint(i, x, 0.);

    if (y != 0.) {
      gr->SetPointError(i,
                        g->GetErrorXlow(i),
                        g->GetErrorXhigh(i),
                        g->GetErrorYlow(i)  / std::abs(y),
                        g->GetErrorYhigh(i) / std::abs(y)
                        );
    } else {
      gr->SetPointError(i,
                        g->GetErrorXlow(i),
                        g->GetErrorXhigh(i),
                        0., 0.
                        );
    }
  }

  return gr;
}
