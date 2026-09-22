#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "xjjstruct.h"

#include "draw.h"
#include "util.h"

// const float cr = 2./3, ytop = 0.84, lspace = 1.2, tsize = 0.038, tsize_up = 0.038/cr;
const float ytop = 0.84, lspace = 1.2, tsize = 0.036;
float ypos(float i = 0, float margin = 0) { return ytop+margin-i*(tsize*lspace); }

int macro(const std::vector<std::string>& inputnames, const std::string& outputname, const std::string& tags,
          float legdx = 0, float legdy = 0,
          int numscan = 0, const std::string& varscan = "") { //
  __XJJLOG << ">> " << inputnames.size() << std::endl;
  std::vector<Color_t> colors = { kBlack };
  if (numscan) xjjc::vec_append(colors, xjjroot::grayscales_color(inputnames.size(), xjjroot::mycolor_middle["red"]));
  else  xjjc::vec_append(colors, { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["magenta"], xjjroot::mycolor_middle["cyan"], xjjroot::mycolor_middle["violet"] });

  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", xjjroot::parse_input(inputnames.front()).content.c_str()));
  TH2D* h2_bins = nullptr;
  if (!h3_bins) xjjana::getobj<TH2D>(Form("%s::h2_bins", xjjroot::parse_input(inputnames.front()).content.c_str()));
  draw::bintex tbins;
  if (h3_bins) {
    tbins.seth(h3_bins, 0, 2);
  } else if (h2_bins) {
    tbins.seth(h2_bins, 0, 1);
  }
  const auto colors_y = xjjroot::grayscales_color(tbins.ny());

  std::map<std::string, xjjc::array2D<TH1D*>> hys;
  std::map<std::string, xjjc::array2D<TGraphErrors*>> gs_scan;
  auto *leg_y = new TLegend(0.21, ypos(0.7 + std::ceil(tbins.ny()/3.)), 0.91, ypos(0.7));
  xjjroot::setleg(leg_y, tsize);
  leg_y->SetNColumns(3);
  auto init_keys = [&hys, &gs_scan, &tbins, &leg_y, &colors_y](TFile* inf) {
    const auto hsdump = xjjana::getobj_regexp<TH1D>(inf, "h1_y_.+__pt-0");
    const auto info = util::read_info(inf, "info", false);
    const auto isNgamma = std::atoi(info.at("event_is").c_str());
    std::vector<std::string> names;
    for (const auto h : hsdump)
      names.push_back(xjjc::str_eraseall(h->GetName(), { "h1_y_", "__pt-0" }));
    xjjc::print_vec_h(names, 1);
    for (int i=0; i < names.size(); i++) {
      const auto key = names[i];
      hys[key].resize(tbins.npt()); //
      gs_scan[key] = xjjc::array2d<TGraphErrors*>(tbins.npt(), tbins.ny(), nullptr); //
      for (int j=0; j<tbins.npt(); j++) {
        for (int k=0; k<tbins.ny(); k++) {
          const auto k_style = isNgamma>0 ? (tbins.ny()-k-1) : k;
          const auto ccy = colors_y[k_style];
          gs_scan[key][j][k] = new TGraphErrors();
          gs_scan[key][j][k]->SetName(Form("gr_%s__pt-%d__y-%d", key.c_str(), j, k));
          xjjroot::setthgrstyle(gs_scan[key][j][k], ccy, xjjroot::mstylelist_solid(k_style), 1.4, ccy, 1, 1);
          if (!i)
            leg_y->AddEntry(gs_scan[key][j][k], tbins.label_y(k).c_str(), "p");
        }
      }
    }
  };
  
  std::vector<double> xs; 
  TLegend *leg = nullptr;
  for (int i=0; i<inputnames.size(); i++) {
    const auto inputp = xjjroot::parse_input(inputnames[i]);
    if (numscan && inputp.pars.size() < 4) {
      __XJJLOG << "!! numscan is true but no scan numbers are provided, switch off numscan." << std::endl;
      numscan = 0;
    }
    auto* inf = TFile::Open(inputp.content.c_str());
    if (hys.empty())
      init_keys(inf);

    // histograms
    const auto cc = colors[i%colors.size()];
    for (auto& [key, vhs] : hys) {
      auto hs = xjjana::getobj_regexp<TH1D>(inf, "h1_y_" + key + "__pt-[0-9]+");
      if (hs.size() != tbins.npt()) {
        __XJJLOG << "!! number of histograms and number of pT bins are inconsistent, abort." << std::endl;
        return 2;
      }
      for (auto& h : hs) {
        xjjroot::setthgrstyle(h, cc, xjjroot::mstylelist_solid(i), 1.4, cc, 1, 1);
        const auto j = xjjc::str_extract_index(h->GetName(), "__pt-");
        vhs[j].push_back(h);
      }
    }
    if (hys.empty() || hys.find("xsec") == hys.end() || hys.at("xsec").empty() || hys.at("xsec").front().empty()) {
      __XJJLOG << "!! no xsec histogram, abort." << std::endl;
      return 2;
    }
    auto* hdump = hys.at("xsec").front().back();

    // legend
    if (!leg) {
      leg = new TLegend(0.60+legdx, ytop+legdy-tsize*lspace*inputnames.size(), 0.85+legdx, ytop+legdy);
      xjjroot::setleg(leg, tsize);
    }
    leg->AddEntry(hdump, inputp.tex.c_str(), "p");

    // numscan
    if (!numscan) continue;
    const auto xval = static_cast<float>(std::atof(inputp.pars[3].c_str()));
    xs.push_back(xval);
    for (auto& [key, vhs] : hys) {
      for (int j=0; j<tbins.npt(); j++) {
        auto* h = vhs[j].back();
        for (int k=0; k<tbins.ny(); k++) {
          auto* gr = gs_scan[key][j][k];
          const int n = gr->GetN();
          gr->SetPoint(n, xval, h->GetBinContent(k+1));
          gr->SetPointError(n, 0, h->GetBinError(k+1));
        }
      }
    } // for (auto& [key, vhs] : hys) {
  }
  xjjroot::print_tab(hys.at("xsec"), 0);
  std::sort(xs.begin(), xs.end());
  if (numscan && xs.size() < 2) numscan = 0;

  // std::vector<TH1D*> hsrelerr;
  for (auto& [_, vhs] : hys) {
    for (auto& vh : vhs) {
      // for (auto& h : vh) {
      //   auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      //   hratio->Divide(vh.front());
      //   hratio->GetYaxis()->SetTitle("Ratio");
      //   // xjjana::sethabsminmax(hratio, 0.71, 1.29);
      //   hratios[key].push_back(hratio);
      // }
      xjjana::sethsabsmin(vh, 0);
      xjjana::sethsmax(vh, 1.8);
      // xjjana::sethsmin(hratios[key], 0.71);
      // xjjana::sethsmax(hratios[key], 1.29);
    }
  }
  // xjjroot::print_tab(hratios, 0);
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  if (numscan) {
    for (int j=0; j<gs_scan.at("xsec").size(); j++) {
      auto& vg = gs_scan.at("xsec")[j];
      const auto ymax = xjjana::gethsmaximum(vg);
      auto* hempty = new TH2F(Form("hempty_xsec__pt-%d", j), Form(";%s;%s", varscan.c_str(), hys.at("xsec").front().front()->GetYaxis()->GetTitle()),
                              10, xs.front() - std::abs(xs[1]-xs[0]), xs.back() + std::abs(xs[1]-xs[0]),
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
      xjjroot::drawtexgroup(0.23, ypos(-0.2), xjjc::str_divide_trim(tags, ","), tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      pdf->write();
    }
  }
  
  // for (auto& [key, vh] : hs) {
  //   pdf->prepare();
  //   auto pads = xjjroot::twopads(pdf->getc(), vh.front(), hratios[key].front(), cr, 0.3);
  //   pads.front()->cd();
  //   vh.front()->Draw("axis");
  //   for (auto& h : vh)
  //     h->Draw("pe1 same");
  //   leg->Draw();
  //   auto labels = xjjc::str_divide_trim(tags, ",");
  //   if (tbins.valid()) {
  //     const auto index_pt = xjjc::str_extract_index(key, "__pt-");
  //     labels.push_back(tbins.label_pt(index_pt));
  //     const auto index_y = xjjc::str_extract_index(key, "__y-");
  //     labels.push_back(tbins.label_y(index_y));
  //   }
  //   xjjroot::drawtexgroup(0.24, ytop-(lspace-1)*tsize/2, labels, tsize, 13, 42, lspace);
  //   xjjroot::drawCMS(xjjroot::CMS::internal, title, 1./cr);

  //   pads.back()->cd();
  //   hratios[key].front()->Draw("axis");
  //   xjjroot::drawline(hratios[key].front()->GetXaxis()->GetXmin(), 1, hratios[key].front()->GetXaxis()->GetXmax(), 1,
  //                     hratios[key].front()->GetLineColor(), 2, 1);
  //   xjjroot::drawbox(hratios[key].front()->GetXaxis()->GetXmin(), 1.05, hratios[key].front()->GetXaxis()->GetXmax(), 0.95,
  //                    hratios[key].front()->GetLineColor(), 0.05, 1001);
  //   for (auto& h : hratios[key])
  //     if (h != hratios[key].front())
  //       h->Draw("pe1 same");
  //   pdf->getc()->cd();
  //   pdf->write();
  // }

  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  __XJJLOG << ">> argc " << argc << std::endl;
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
