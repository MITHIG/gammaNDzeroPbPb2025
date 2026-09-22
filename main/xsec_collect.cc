#include "xjjanauti.h"
#include "xjjmypdf.h"

// #include "measurements.h"
#include "draw.h"
#include "util.h"

int macro(const std::vector<std::string>& inputnames, const std::string& outputname,
          int incl_syst = 0) {
  if (inputnames.size() != 2) {
    __XJJLOG << "!! inputnames should have two files, abort." << std::endl;
    return 2;
  }

  const std::vector<Color_t> colors = { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], kBlack };
  const std::vector<std::string> labels = { "gammaN", "Ngamma", "sum" };
  const std::vector<std::string> texs = { "Xn0n (#gammaN)", "0nXn (N#gamma)", "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" };
  
  xjjc::array2D<TH1D*> hsxsec;
  std::vector<std::string> tagspt;
  float lumi = -1.;
  TH3D* h3_bins = nullptr;
  for (const auto& input : inputnames) {
    __XJJLOG << "++ " << input << std::endl;
    auto* inf = TFile::Open(input.c_str());
    const auto info = util::read_info(inf);
    lumi = std::atof(info.at("lumi").c_str());
    const auto isNgamma = std::atoi(info.at("event_is").c_str());
    if (isNgamma < 0 || isNgamma > 1) continue;
    
    __XJJLOG << ">> [" << labels[isNgamma] << "]" << std::endl;
    const auto cc = colors[isNgamma];
    if (!h3_bins) h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins");
    auto vh = xjjana::getobj_regexp<TH1D>(inf, "h1_y_xsec__pt-[0-9]+");
    if (hsxsec.empty()) {
      hsxsec = xjjc::array2d<TH1D*>(vh.size(), 3, nullptr);
      tagspt.resize(vh.size());
    }
    for (auto& h : vh) {
      const auto j = xjjc::str_extract_index(h->GetName(), "__pt-");
      tagspt[j] = xjjc::str_extract_regex(h->GetName(), "(__pt-[0-9]+)").front();
      auto name_base = xjjc::str_replaceall_regex(h->GetName(), "__pt-[0-9]+", "");
      h->SetName(xjjc::str_replaceall(name_base, "_xsec", "_xsec_" + labels[isNgamma]).c_str());
      if (isNgamma) {
        if (util::mirrorswap_hist(h)) return 2;
      }
      xjjroot::setthgrstyle(h, cc, 20, 1.5, cc, 1, 1);
      hsxsec[j][isNgamma] = h;
    }
  }
  if (!h3_bins) {
    __XJJLOG << "!! cannot read binning histogram, abort." << std::endl;
    return 2;
  }
  draw::bintex tbins(h3_bins, 0, 2);

  for (auto& hs : hsxsec) {
    auto *h_gammaN = hs[0], *h_Ngamma = hs[1];
    auto* h_sum = (TH1D*)h_gammaN->Clone(xjjc::str_replaceall(h_gammaN->GetName(), "gammaN", "sum").c_str());
    h_sum->Add(h_Ngamma);
    hs[2] = h_sum;
    for (int i=0; i<hs.size(); i++)
      xjjroot::setthgrstyle(hs[i], colors[i], 20, 1.5, colors[i], 1, 1);
  }

  auto draw_global = [&lumi]() {
    xjjroot::drawtexgroup(0.88, 0.86, { xjjroot::CMS::DzDzbar2 }, 0.043, 33, 42, 1.25);
    xjjroot::drawCMS(xjjroot::CMS::internal, Form("%.1f #mub^{-1} (PbPb 5.36 TeV)", lumi*1.e3));
  };

  auto draw_hs = [](std::vector<TH1D*> vh) {
    if (vh.empty()) return;
    xjjana::sethsmax(vh, 1.6);
    xjjana::sethsmin(vh, 0.);
    vh.front()->Draw("axis");
    for (auto h : vh) h->Draw("pe1 same");
  };

  auto* leg = new TLegend(0.65, 0.70-0.038*1.25*2, 0.85, 0.70);
  xjjroot::setleg(leg, 0.038);
  for (const int i : {0, 1})
    leg->AddEntry(hsxsec.front()[i], texs[i].c_str(), "p");
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (int j=0; j<hsxsec.size(); j++) {
    pdf->prepare();
    draw_hs({ hsxsec[j][0], hsxsec[j][1] });
    draw_global();
    leg->Draw();
    xjjroot::drawtexgroup(0.23, 0.85, { tbins.label_pt(j) }, 0.04, 13, 42, 1.25);
    pdf->write();

    pdf->prepare();
    draw_hs({ hsxsec[j][2] });
    draw_global();
    xjjroot::drawtexgroup(0.23, 0.85, { texs[2], tbins.label_pt(j) }, 0.04, 13, 42, 1.25);
    pdf->write();
  }

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h3_bins);
  for (int j=0; j<hsxsec.size(); j++) {
    auto* dir = outf->mkdir(Form("dir%s", tagspt[j].c_str()));
    dir->cd();
    for (auto& h : hsxsec[j])
      xjjroot::writehist(h);

    outf->cd();
  }
  xjjroot::closefile(outf);
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], std::atoi(argv[3]));
  }
}
