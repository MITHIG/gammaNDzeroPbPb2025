#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "measurements.h"
#include "draw.h"
#include "util.h"

int macro(const std::vector<std::string>& inputnames, const std::string& outputname, int save_png = 0) {
  if (inputnames.size() != 2) {
    __XJJLOG << "!! inputnames should have two files, abort." << std::endl;
    return 2;
  }

  enum class Cat { gammaN, Ngamma, NgammaRef, sum };
  const std::vector<std::string> labels = { "gammaN", "Ngamma", "Ngamma-ref", "sum" };
  const std::vector<std::string> texs = { "Xn0n (#gammaN)", "0nXn (N#gamma)", "0nXn (N#gamma)#scale[0.4]{ }#it{y}#scale[0.4]{ }#rightarrow -#it{y}", "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" };
  const std::vector<Color_t> colors = { xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["red"], kBlack };
  const auto ncat = labels.size();
  
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
      hsxsec = xjjc::array2d<TH1D*>(vh.size(), ncat, nullptr);
      tagspt.resize(vh.size());
    }
    for (auto& h : vh) {
      const auto j = xjjc::str_extract_index(h->GetName(), "__pt-");
      tagspt[j] = xjjc::str_extract_regex(h->GetName(), "(__pt-[0-9]+)").front();
      auto name_base = xjjc::str_replaceall_regex(h->GetName(), "__pt-[0-9]+", "");
      h->SetName(xjjc::str_replaceall(name_base, "_xsec", "_xsec_" + labels[isNgamma]).c_str());
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
    auto* h_NgammaRef = (TH1D*)h_Ngamma->Clone(xjjc::str_replaceall(h_Ngamma->GetName(), "Ngamma", "NgammaRef").c_str());
    if (util::mirrorswap_hist(h_NgammaRef)) return 2;
    hs[2] = h_NgammaRef;
    auto* h_sum = (TH1D*)h_gammaN->Clone(xjjc::str_replaceall(h_gammaN->GetName(), "gammaN", "sum").c_str());
    h_sum->Add(h_NgammaRef);
    hs[3] = h_sum;
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

  const float x1 = 0.54, y1 = 0.755, tsize = 0.038, lspace = 1.25, tlsize = tsize*lspace;
  auto y_pos = [&y1, &tlsize](float i = 0, float margin = 0) { return (y1 + margin - i*tlsize); };

  auto* leg = new TLegend(x1, y_pos(2), x1 + 0.3, y_pos());
  xjjroot::setleg(leg, tsize);
  for (const int i : { int(Cat::gammaN), int(Cat::NgammaRef) })
    leg->AddEntry(hsxsec.front()[i], texs[i].c_str(), "p");

  auto draw_HIN_25_002 = [&tbins, &tsize, &y1, &tlsize](int isNgamma, TH1D* hdata, const std::string& tdata) {
    measurement::draw_HIN_25_002(static_cast<Event>(isNgamma));
    auto xx = isNgamma ? 0.222 : 0.54;
    auto* leg_vs23 = new TLegend(xx, y1-tlsize*0.9 - tlsize*(tbins.npt()+2), xx+0.3, y1-tlsize*0.9);
    xjjroot::setleg(leg_vs23, tsize);
    leg_vs23->AddEntry(hdata, tdata.c_str(), "p");
    leg_vs23->AddEntry((TObject*)0, "", NULL);
    leg_vs23->AddEntry(measurement::get_style(), xjjroot::str_fixspace(xjjc::number_range_string(float(2), float(5), "#it{p}_{T}") + " GeV").c_str(), "pf");
    xjjroot::drawtex(leg_vs23->GetX1() + 0.008, leg_vs23->GetY2()+tlsize*0.5 - 0.005, "This analysis", tsize, 12);
    leg_vs23->Draw();
    xjjroot::drawtex(leg_vs23->GetX1() + 0.008, leg_vs23->GetY1()+tlsize*1.5 - 0.005, "2023 PbPb (HIN-25-002)", tsize, 12);;
  };
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  const auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });
  for (int j=0; j<hsxsec.size(); j++) {
    pdf->prepare();
    draw_hs({ hsxsec[j][int(Cat::gammaN)], hsxsec[j][int(Cat::NgammaRef)] });
    draw_global();
    leg->Draw();
    xjjroot::drawtexgroup(0.23, 0.85, { tbins.label_pt(j) }, 0.04, 13, 42, 1.25);
    pdf->write(Form("%s_symm_pt-%d.pdf", name_png.c_str(), j), save_png ? "" : "X");

    xjjroot::setthgrstyle(hsxsec[j][int(Cat::gammaN)], xjjroot::mycolor_middle["red"], -1, -1, xjjroot::mycolor_middle["red"]);

    for (const int i : { int(Cat::gammaN), int(Cat::Ngamma) }) {
      pdf->prepare();
      hsxsec[j][i]->Draw("axis");
      draw_HIN_25_002(i, hsxsec[j][i], tbins.label_pt(j));
      hsxsec[j][i]->Draw("pe1 same");
      draw_global();
      xjjroot::drawtexgroup(0.23, 0.85, { texs[i] }, 0.04, 13, 42, 1.25);
      pdf->write(Form("%s_%s_pt-%d.pdf", name_png.c_str(), labels[i].c_str(), j), save_png ? "" : "X");
    }
    
    pdf->prepare();
    draw_hs({ hsxsec[j][int(Cat::sum)] });
    draw_global();
    xjjroot::drawtexgroup(0.23, 0.85, { texs[int(Cat::sum)], tbins.label_pt(j) }, 0.04, 13, 42, 1.25);
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
  util::Writeinfo tinfo;
  tinfo.init("info");
  tinfo.cast_branch("lumi", lumi);
  tinfo.close();
  xjjroot::closefile(outf);
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], std::atoi(argv[3]));
  }
}
