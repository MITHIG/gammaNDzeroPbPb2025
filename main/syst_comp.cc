#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "xjjstruct.h"

#include "draw.h"
#include "util.h"

struct Property {
  std::string tex;
  // int isNgamma;
};

const float ytop = 0.84, lspace = 1.2, tsize = 0.036;
float ypos(float i = 0, float margin = 0) { return ytop+margin-i*(tsize*lspace); }

int macro(const std::vector<std::string>& inputnames, const std::string& outputname, const std::string& tags,
          float legdx = 0, float legdy = 0,
          int numscan = 0, const std::string& varscan = "") { //
  __XJJLOG << ">> " << inputnames.size() << std::endl;
  std::vector<Color_t> colors = { kBlack };
  // if (numscan) xjjc::vec_append(colors, xjjroot::grayscales_color(inputnames.size(), xjjroot::mycolor_middle["red"]));
  if (numscan) xjjc::vec_append(colors, std::vector<Color_t>(inputnames.size()-1, xjjroot::mycolor_middle["blue"]));
  else  xjjc::vec_append(colors, { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["magenta"], xjjroot::mycolor_middle["cyan"], xjjroot::mycolor_middle["violet"] });
  const auto alphas = xjjroot::grayscales_alpha(inputnames.size());

  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", xjjroot::parse_input(inputnames.front()).content.c_str()));
  draw::bintex tbins(h3_bins, 0, 2);
  const auto colors_y = xjjroot::grayscales_color(tbins.ny());

  const std::map<std::string, Property> texs = {
    { "gammaN", { .tex = "Xn0n (#gammaN)" } },
    { "Ngamma", { .tex = "Xn0n (N#gamma) (#it{y} #rightarrow -#it{y})" } },
    { "sum", { .tex = "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" } }
  };
  std::map<std::string, xjjc::array2D<TH1D*>> hys;
  std::map<std::string, xjjc::array2D<TGraphErrors*>> gs_scan;
  for (const auto& [key, prop] : texs) {
    hys[key].resize(tbins.npt()); //
    gs_scan[key] = xjjc::array2d<TGraphErrors*>(tbins.npt(), tbins.ny(), nullptr); //
    for (int j=0; j<tbins.npt(); j++) {
      for (int k=0; k<tbins.ny(); k++) {
        const auto k_style = k;
        const auto ccy = colors_y[k_style];
        gs_scan[key][j][k] = new TGraphErrors();
        gs_scan[key][j][k]->SetName(Form("gr-scan_y_xsec_%s__pt-%d__y-%d", key.c_str(), j, k));
        xjjroot::setthgrstyle(gs_scan[key][j][k], ccy, xjjroot::mstylelist_solid(k_style), 1.4, ccy, 1, 1);
      }
    }
  }
  auto *leg_y = new TLegend(0.21, ypos(0.7 + std::ceil(tbins.ny()/3.)), 0.91, ypos(0.7));
  xjjroot::setleg(leg_y, tsize);
  leg_y->SetNColumns(3);
  for (int k=0; k<tbins.ny(); k++)
    leg_y->AddEntry(gs_scan.at("sum")[0][k], tbins.label_y(k).c_str(), "p");
  
  std::vector<double> xvals; 
  TLegend *leg = nullptr;
  for (int i=0; i<inputnames.size(); i++) {
    const auto inputp = xjjroot::parse_input(inputnames[i]);
    if (numscan && inputp.pars.size() < 4) {
      __XJJLOG << "!! numscan is true but no scan numbers are provided, switch off numscan." << std::endl;
      numscan = 0;
    }
    auto* inf = TFile::Open(inputp.content.c_str());
    const auto cc = colors[i%colors.size()];
    for (int j=0; j<tbins.npt(); j++) {
      auto* dir = (TDirectory*)inf->Get(Form("dir__pt-%d", j));
      for (const auto& [key, _] : texs) {
        auto* h = xjjana::getobj<TH1D>(dir, "h1_y_xsec_" + key);
        if (numscan)
          xjjroot::setthgrstyle(h, cc, xjjroot::mstylelist_solid(i), 1.4, cc, 1, 1, 0, 0, 0, alphas[i], alphas[i]);
        else
          xjjroot::setthgrstyle(h, cc, xjjroot::mstylelist_solid(i), 1.4, cc, 1, 1);
        hys.at(key)[j].push_back(h); //
      }
    }
    if (hys.at("sum").front().empty()) {
      __XJJLOG << "!! no xsec histogram, abort." << std::endl;
      return 2;
    }
    auto* hdump = hys.at("sum").front().back();
    // legend
    if (!leg) {
      leg = new TLegend(0.60+legdx, ypos(1+inputnames.size(), legdy), 0.85+legdx, ypos(1, legdy));
      xjjroot::setleg(leg, tsize);
    }
    leg->AddEntry(hdump, inputp.tex.c_str(), "p");

    if (!numscan) continue;
    const auto xval = static_cast<float>(std::atof(inputp.pars[3].c_str()));
    xvals.push_back(xval); //
  }

  std::map<std::string, std::vector<TGraphAsymmErrors*>> gs_err, gs_rel;
  for (auto& [key, hs] : hys) {
    gs_err[key].resize(tbins.npt(), nullptr);
    gs_rel[key].resize(tbins.npt(), nullptr);
    for (int j=0; j<tbins.npt(); j++) {
      auto* h0 = hs[j][0];
      std::vector<double> xs(tbins.ny(), 0), xserr(tbins.ny(), 0),
        ys(tbins.ny(), 0), maxdev_up(tbins.ny(), 0), maxdev_low(tbins.ny(), 0),
        y0(tbins.ny(), 0), maxrel_up(tbins.ny(), 0), maxrel_low(tbins.ny(), 0);
      for (int i=0; i<hs[j].size(); i++) {
        auto* h = hs[j][i];
        for (int k=0; k<tbins.ny(); k++) {
          auto* gr = gs_scan[key][j][k];
          const int n = gr->GetN();
          gr->SetPoint(n, xvals[i], h->GetBinContent(k+1));
          gr->SetPointError(n, 0, h->GetBinError(k+1));

          const double dev = h->GetBinContent(k+1) - h0->GetBinContent(k+1);
          maxdev_up[k] = std::max(maxdev_up[k], dev);
          maxdev_low[k] = std::min(maxdev_low[k], dev);
          if (!i) {
            xs[k] = h->GetBinCenter(k+1);
            ys[k] = h0->GetBinContent(k+1);
            xserr[k] = h->GetBinWidth(k+1) / 2.;
          }
        }
      }
      xjjana::sethsabsmin(hs[j], 0);
      xjjana::sethsmax(hs[j], 1.8);
      for (int k=0; k<tbins.ny(); k++) {
        maxdev_low[k] = std::abs(maxdev_low[k]);
        maxrel_up[k] = maxdev_up[k] / ys[k];
        maxrel_low[k] = maxdev_low[k] / ys[k];        
      }
      auto* gr_err = new TGraphAsymmErrors(tbins.ny(), xs.data(), ys.data(), xserr.data(), xserr.data(), maxdev_low.data(), maxdev_up.data());
      gr_err->SetName(Form("gr-err_y_xsec_%s", key.c_str()));
      xjjroot::setthgrstyle(gr_err, kBlack, 20, 1.4, kBlack, 1, 2, 0, 0, 0);
      gs_err[key][j] = gr_err;
      auto* gr_rel = new TGraphAsymmErrors(tbins.ny(), xs.data(), y0.data(), xserr.data(), xserr.data(), maxrel_low.data(), maxrel_up.data());
      gr_rel->SetName(Form("gr-rel_y_xsec_%s", key.c_str()));
      xjjroot::setthgrstyle(gr_rel, kBlack, 20, 1.4, kBlack, 1, 1, xjjroot::color_alpha(kBlack, 0.2), 1, 1001);
      gs_rel[key][j] = gr_rel;
    }
  }
  xjjroot::print_tab(hys.at("sum"), 0);
  std::sort(xvals.begin(), xvals.end());
  if (numscan && xvals.size() < 2) numscan = 0;

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (int j=0; j<tbins.npt(); j++) {
    for (const auto& [key, prop] : texs) {
      auto& vh = hys.at(key)[j];
      pdf->prepare();
      vh.front()->Draw("axis");
      for (auto& h : vh)
        h->Draw("pe1 same");
      vh.front()->Draw("pe1 same");      
      gs_err.at(key)[j]->Draw("5 same");
      if (inputnames.size() < 10) leg->Draw();
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex,  }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      pdf->write();
    }
    for (const auto& [key, prop] : texs) {
      auto& g = gs_rel.at(key)[j];
      const auto ymax = xjjana::gethmaximum(g), ymin = xjjana::gethminimum(g);
      // std::cout<<ymin<<" -> "<<ymax<<std::endl;
      auto* hempty = new TH2F(Form("hempty-rel_xsec_%s__pt-%d", key.c_str(), j), Form(";%s;Relative Uncertainty", hys.at(key).front().front()->GetXaxis()->GetTitle()),
                              10, hys.at(key).front().front()->GetXaxis()->GetXmin(), hys.at(key).front().front()->GetXaxis()->GetXmax(),
                              10, -0.2, 0.2);
                              // 10, ymin*1.5, ymax*1.5);
      xjjroot::sethempty(hempty, 0, 0.2);
      pdf->prepare();
      hempty->Draw("axis");
      g->Draw("2 same");
      xjjroot::drawline_horizon(0, hempty, kBlack, 2, 1);
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex  }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      gPad->RedrawAxis();
      pdf->write();
    }
  }

  if (numscan) {
    for (int j=0; j<tbins.npt(); j++) {
      for (const auto& [key, prop] : texs) {
        auto& vg = gs_scan.at(key)[j];
        const auto ymax = xjjana::gethsmaximum(vg);
        auto* hempty = new TH2F(Form("hempty-scan_xsec_%s__pt-%d", key.c_str(), j), Form(";%s;%s", varscan.c_str(), hys.at(key).front().front()->GetYaxis()->GetTitle()),
                                10, xvals.front() - std::abs(xvals[1]-xvals[0]), xvals.back() + std::abs(xvals[1]-xvals[0]),
                                10, 0, ymax * 1.8);
        xjjroot::sethempty(hempty, 0, 0.2);
        pdf->prepare();
        hempty->Draw("axis");
        for (auto& g : vg) {
          auto* g2 = (TGraphErrors*)g->Clone(Form("%s_2", g->GetName()));
          xjjana::sortgrx(g2);
          xjjroot::setthgrstyle(g2, -1, -1, -1, -1, 2, 1);
          g2->Draw("l same");
          g->Draw("pe same");
        }
        leg_y->Draw();
        // xjjroot::drawtexgroup(0.23, ypos(-0.2), xjjc::str_divide_trim(tags, ","), tsize, 13);
        xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex }, tsize, 13);
        xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
        xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
        pdf->write();
      }
    }
  }
  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h3_bins);
  for (int j=0; j<tbins.npt(); j++) {
    auto* dir = outf->mkdir(Form("dir__pt-%d", j));
    dir->cd();
    for (const auto& [key, _] : texs) {
      xjjroot::writehist(hys.at(key)[j].front());
      xjjroot::writehist(gs_err.at(key)[j]);
      xjjroot::writehist(gs_rel.at(key)[j]);
    }
    outf->cd();
  }
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 8) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], std::atof(argv[4]), std::atof(argv[5]), std::atoi(argv[6]), argv[7]);
  }
  if (argc == 6) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], std::atof(argv[4]), std::atof(argv[5]));
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
