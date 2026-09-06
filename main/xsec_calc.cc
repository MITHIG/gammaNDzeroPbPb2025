#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/measurements.h"
#include "../include/draw.h"

namespace global {
  float BR_DtoKpi = 0.03936, err_BR_DtoKpi = 0.030*1.e-2;
}

int macro(const std::string& inputname_raw, const std::string& inputname_effd,
          const std::string& inputname_effevent, const std::string inputname_fprompt = "null",
          float lumi = 1., const std::string& outputdir = "") {

  __XJJLOG << ">> lumi: " << lumi << " nb-1" << std::endl;
  std::map<std::string, std::vector<TH1D*>> h1pts;
  // std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, std::map<std::string, std::string>> infos;
  auto* h3_bins = xjjana::getobj<TH3D>(inputname_raw + "::h3_bins");
  if (!h3_bins)
    return 2;
  draw::bintex tbins(h3_bins, 0, 2);
  //
  std::string tag = "xsec";
  auto get_h1pts = [&h1pts, &infos, &tag](const std::string &inputname, const std::string &category,
                                                    const std::vector<std::string>& h1names,
                                                    const std::vector<std::string>& infots) {
    __XJJLOG << "[" << category << "] " << inputname << std::endl;
    auto* inf = TFile::Open(inputname.c_str());
    auto itag = xjjc::str_tag_from_file(inputname);
    tag += ("_" + itag);
    if (!inf) {
      __XJJLOG << "?? no " << category << " input file: " << inputname << ", skip." << std::endl;
      return 2;
    }
    // TH1D
    std::vector<Color_t> colors = { xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["red"] };
    for (auto& name : h1names) {
      // auto vh = xjjana::getobj_regexp<TH1D>(inf, name + ".*__pt-.+");
      auto vh = xjjana::getobj_regexp<TH1D>(inf, name + "__pt-.+");
      if (vh.empty()) {
        __XJJLOG << "?? TH1D " << name << " not found, skip." << std::endl;
        continue;
      }
      for (int i=0; i<vh.size(); i++) {
        xjjroot::sethempty(vh[i], 0, 0.4);
        auto cc = vh.size() > 1 ? colors[i] : kBlack;
        xjjroot::setthgrstyle(vh[i], cc, 21, 1.6, cc, 1, 1);
      }
      h1pts[xjjc::str_eraseall(name, "h1_")] = vh;
    }
    // info
    for (auto& tr : infots) {
      auto info = xjjana::getval_regexp((TTree*)inf->Get(tr.c_str()));
      auto name = category;
      if (name != "info")
        name += ("_" + xjjc::str_eraseall(tr, "/info"));
      infos[name] = info;
    }
    return 0;
  };

  get_h1pts(inputname_raw, "raw", { "h1_y_yield" }, { "data/info", "template/info" });
  get_h1pts(inputname_effd, "effd", { "h1_y_eff__rebin" }, { "info" });
  get_h1pts(inputname_effevent, "effevent", { "h1_y_evteff" }, { "info" });
  get_h1pts(inputname_fprompt, "fprompt", {  }, {  });

  Event event_is = xjjc::str_contains(infos.at("raw_data").at("cut_tex"), "#gammaN") ? Event::gammaN :
    (xjjc::str_contains(infos.at("raw_data").at("cut_tex"), "N#gamma") ? Event::Ngamma : Event::Other);
  // auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  // __XJJLOG << "++ info" << std::endl;
  // xjjc::print_tab(info, -1);

  xjjroot::print_tab(h1pts, 0);
  for (int j=0; j<tbins.npt(); j++) {
    auto* h_corr = (TH1D*)h1pts.at("y_yield")[j]->Clone(xjjc::str_replaceall(h1pts.at("y_yield")[j]->GetName(), "yield", "corr").c_str());
    h_corr->Divide(h1pts.at("y_eff__rebin")[j]);
    h_corr->GetYaxis()->SetTitle("Corrected Yield");
    h1pts["y_corr"].push_back(h_corr);
    auto* h_xsec = (TH1D*)h_corr->Clone(xjjc::str_replaceall(h_corr->GetName(), "corr", "xsec").c_str());
    h_xsec->Divide(h1pts.at("y_evteff")[j]);
    h_xsec->Scale(0.5/global::BR_DtoKpi/(lumi*1.e6)/tbins.binwidth_pt(j)/*dpt*/, "width");
    h_xsec->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{y}d#it{p}_{T}} [mb/GeV]");
    xjjroot::sethempty(h_xsec, 0, 0.2);
    if (tbins.npt() == 1)
      xjjroot::setthgrstyle(h_xsec, xjjroot::mycolor_middle["red"], -1, -1, xjjroot::mycolor_middle["red"]);
    h1pts["y_xsec"].push_back(h_xsec);
  }

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputdir + "/" + tag + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" }});

  auto draw_global = [&infos, &lumi](bool divide2) {
    xjjroot::drawCMS(xjjroot::CMS::internal, xjjc::str_replaceall(infos.at("raw_data").at("input_tex"), ")", Form(", %.1f#scale[0.3]{ }#mub^{-1})", lumi*1.e3)));
    xjjroot::drawtexgroup(0.23, 0.85, { infos.at("raw_data").at("cut_tex") }, 0.04, 13, 42, 1.25);
    xjjroot::drawtexgroup(0.88, 0.86, { divide2 ? xjjroot::CMS::DzDzbar2 : xjjroot::CMS::DznDzbar }, 0.043, 33, 42, 1.25);
  };
  
  xjjana::sethsmin(h1pts.at("y_xsec"), 0.);
  auto ymax = xjjana::sethsmax(h1pts.at("y_xsec"), 1.8);
  if (ymax < 2.5) xjjana::sethsabsmax(h1pts.at("y_xsec"), 3.2);
  else if (ymax >= 2.5 && ymax < 5.) xjjana::sethsabsmax(h1pts.at("y_xsec"), 6.6);

  const auto x1 = (event_is==Event::Ngamma ? 0.222 : 0.54), y1 = 0.755,
    tsize = 0.038, lspace = 1.25, tlsize = tsize*lspace;
  auto* legvs23 = new TLegend(x1, y1-tlsize*0.9 - tlsize*(tbins.npt()+2), x1+0.3, y1-tlsize*0.9);
  xjjroot::setleg(legvs23, tsize);
  for (int i=0; i<tbins.npt(); i++)
    legvs23->AddEntry(h1pts.at("y_xsec")[i], tbins.label_pt(i).c_str(), "p");
  legvs23->AddEntry((TObject*)0, "", NULL);
  legvs23->AddEntry(measurement::get_style(), xjjroot::str_fixspace(xjjc::number_range_string(float(2), float(5), "#it{p}_{T}") + " GeV").c_str(), "pf");

  pdf->prepare();
  h1pts.at("y_xsec").front()->Draw("axis");
  measurement::draw_HIN_25_002(event_is);
  for (auto& h : h1pts.at("y_xsec"))
    h->Draw("pe1 same");
  draw_global(true);
  legvs23->Draw();
  xjjroot::drawtex(legvs23->GetX1() + 0.008, legvs23->GetY2()+tlsize*0.5 - 0.005, "This analysis", tsize, 12);
  xjjroot::drawtex(legvs23->GetX1() + 0.008, legvs23->GetY1()+tlsize*1.5 - 0.005, "2023 PbPb (HIN-25-002)", tsize, 12);
  pdf->write(name_png + "_xsec_vs23.pdf");

  for (auto& h : h1pts.at("y_xsec"))
    xjjroot::print_th(h);
  
  auto* leg = new TLegend(x1, y1-tlsize*tbins.npt(), x1+0.3, y1);
  xjjroot::setleg(leg, tsize);
  for (int i=0; i<tbins.npt(); i++)
    leg->AddEntry(h1pts.at("y_yield")[i], tbins.label_pt(i).c_str(), "p");
  leg->Draw();

  for (auto& t : { "y_yield", "y_eff__rebin", "y_evteff", "y_corr", "y_xsec" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts.at(t), 0.);
    xjjana::sethsmax(h1pts.at(t), 1.8);
    h1pts.at(t).front()->Draw("axis");
    for (auto& h : h1pts.at(t)) {
      if (h1pts.size() == 1)
        xjjroot::setthgrstyle(h, kBlack, -1, -1, kBlack, -1, -1);
      h->Draw("pe1 same");
    }
    leg->Draw();
    draw_global(xjjc::str_contains(t, "xsec"));
    pdf->write(name_png + xjjc::str_replaceall(t, "y_", "_") + ".pdf");
  }

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputdir + "/" + tag + ".root");
  // for (auto& [_, h] : h2s) xjjroot::writehist(h);
  for (auto& [_, hh] : h1pts)
    for (auto& h : hh)
      xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);
  // auto* t = new TTree("info", "");
  // for (auto& [key, content] : info) {
  //   t->Branch(key.c_str(), &content);
  // }
  // t->Fill();
  // t->Write();
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    return macro(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), argv[6]);
  }
  return 1;
}

