#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/measurements.h"
#include "../include/draw.h"

namespace global {
  float BR_DtoKpi = 0.03936, err_BR_DtoKpi = 0.030*1.e-2;
}

int macro(const std::string& inputname_raw, const std::string& inputname_effd,
          const std::string& inputname_effevent = "null", const std::string inputname_fprompt = "null",
          float lumi = 1., const std::string& outputdir = "") {

  __XJJLOG << ">> lumi: " << lumi << " nb-1" << std::endl;
  std::map<std::string, std::vector<TH1D*>> h1pts;
  std::map<std::string, std::map<std::string, std::string>> infos;
  std::string tag = "xsec";
  auto* h3_bins = xjjana::getobj<TH3D>(inputname_raw + "::h3_bins");
  if (!h3_bins)
    return 2;
  draw::bintex tbins(h3_bins, 0, 2);
  auto get_h1pts = [&h1pts, &infos, &tag](const std::string &inputname, const std::string &category,
                                          const std::vector<std::string>& h1names,
                                          const std::vector<std::string>& infots) {
    __XJJLOG << "[" << category << "] " << inputname << std::endl;
    auto* inf = TFile::Open(inputname.c_str());
    tag += ("_" + xjjc::str_tag_from_file(inputname));
    if (!inf) {
      __XJJLOG << "?? no " << category << " input file: " << inputname << ", skip." << std::endl;
      return 2;
    }
    // TH1D
    std::vector<Color_t> colors = { xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["red"] };
    for (auto& name : h1names) {
      auto vh = xjjana::getobj_regexp<TH1D>(inf, name+".*__pt-.+");
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
  get_h1pts(inputname_effd, "effd", { "h1_y_eff__rebin", "h1_y_effdata__rebin" }, { "info" });
  get_h1pts(inputname_effevent, "effevent", {  }, {  });
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
    xjjroot::drawtexgroup(0.24, 0.86, { infos.at("raw_data").at("cut_tex") }, 0.04, 13, 42, 1.25);
    xjjroot::drawtexgroup(0.88, 0.87, { divide2 ? xjjroot::CMS::DzDzbar2 : xjjroot::CMS::DznDzbar }, 0.043, 33, 42, 1.25);
  };
  
  xjjana::sethsmin(h1pts.at("y_xsec"), 0.);
  xjjana::sethsmax(h1pts.at("y_xsec"), 2.);
  // for (auto& h : h1pts.at("y_xsec")) h->SetMaximum(5.);
  pdf->prepare();
  h1pts.at("y_xsec").front()->Draw("axis");
  auto* g_HIN_25_002 = measurement::draw_HIN_25_002(event_is);
  for (auto& h : h1pts.at("y_xsec"))
    h->Draw("pe1 same");
  draw_global(true);
  auto* legvs23 = new TLegend((event_is==Event::Ngamma ? 0.24 : 0.55), 0.75-0.042*(tbins.npt()+2), (event_is==Event::Ngamma ? 0.24 : 0.55)+0.3, 0.75);
  xjjroot::setleg(legvs23, 0.038);
  legvs23->SetHeader("This analysis");
  for (int i=0; i<tbins.npt(); i++)
    legvs23->AddEntry(h1pts.at("y_xsec")[i], tbins.label_pt(i).c_str(), "p");
  legvs23->AddEntry(g_HIN_25_002, "PAS-HIN-25-002", "pf");
  legvs23->Draw();
  pdf->write(name_png + "_xsec_vs23.pdf");

  auto* leg = new TLegend((event_is==Event::Ngamma ? 0.24 : 0.55), 0.75-0.042*tbins.npt(), (event_is==Event::Ngamma ? 0.24 : 0.55)+0.3, 0.75);
  xjjroot::setleg(leg, 0.038);
  for (int i=0; i<tbins.npt(); i++)
    leg->AddEntry(h1pts.at("y_yield")[i], tbins.label_pt(i).c_str(), "p");
  leg->Draw();

  for (auto& t : { "y_yield", "y_eff__rebin", "y_corr", "y_xsec" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts.at(t), 0.);
    xjjana::sethsmax(h1pts.at(t), 1.8);
    h1pts.at(t).front()->Draw("axis");
    for (auto& h : h1pts.at(t)) {
      xjjroot::setthgrstyle(h, kBlack, -1, -1, kBlack, -1, -1);
      h->Draw("pe1 same");
    }
    leg->Draw();
    draw_global(xjjc::str_contains(t, "xsec"));
    pdf->write(name_png + xjjc::str_replaceall(t, "y_", "_") + ".pdf");
  }

  pdf->close();

  // auto* outf = xjjroot::newfile(xjjc::str_replaceall(inputname, "saveeff", "calceff"));
  // for (auto& [_, h] : h2s) xjjroot::writehist(h);
  // for (auto& [_, h] : h1s) xjjroot::writehist(h);
  // auto* t = new TTree("info", "");
  // for (auto& [key, content] : info) {
  //   t->Branch(key.c_str(), &content);
  // }
  // t->Fill();
  // t->Write();
  // outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    return macro(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), argv[6]);
  }
  return 1;
}

