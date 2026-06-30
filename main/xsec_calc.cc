#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/measurements.h"

namespace global {
  float BR_DtoKpi = 0.03936, err_BR_DtoKpi = 0.030*1.e-2;
}
// enum Event { gammaN, Ngamma, Other };
int macro(const std::string& inputname_raw, const std::string& inputname_effd,
          const std::string& inputname_effevent = "null", const std::string inputname_fprompt = "null",
          float lumi = 1., const std::string& outputdir = "") { // nb-1

  __XJJLOG << ">> lumi: " << lumi << std::endl;
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, std::map<std::string, std::string>> infos;
  std::string tag;
  auto get_h1s = [&h1s, &infos, &tag](const std::string &inputname, const std::string &category,
                                      const std::vector<std::string>& h1names,
                                      const std::vector<std::string>& infots) {
    auto* inf = TFile::Open(inputname.c_str());
    tag += ((tag.empty() ? "" : "_") + xjjc::str_tag_from_file(inputname));
    if (!inf) {
      __XJJLOG << "!! no " << category << " input file: " << inputname << ", skip." << std::endl;
      return 2;
    }
    // TH1D
    for (auto& name : h1names) {
      auto* h = xjjana::getobj<TH1D>(inf, name);
      if (!h) {
        __XJJLOG << "!! TH1D " << name << " not found, skip." << std::endl;
        continue;
      }
      xjjroot::sethempty(h, 0, 0.4);
      xjjroot::setthgrstyle(h, kBlack, 21, 1.6, kBlack, 1, 1);
      h1s[xjjc::str_eraseall(name, "h1_")] = h;
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

  get_h1s(inputname_raw, "raw", { "h1_y_yield" }, { "data/info", "template/info" });
  get_h1s(inputname_effd, "effd", { "h1_y_eff__rebin", "h1_y_effdata__rebin" }, { "info" });
  get_h1s(inputname_effevent, "effevent", {  }, {  });
  get_h1s(inputname_fprompt, "fprompt", {  }, {  });
    
  Event event_is = xjjc::str_contains(infos.at("raw_data").at("cut_tex"), "#gammaN") ? Event::gammaN :
    (xjjc::str_contains(infos.at("raw_data").at("cut_tex"), "N#gamma") ? Event::Ngamma : Event::Other);
  // auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  // __XJJLOG << "++ info" << std::endl;
  // xjjc::print_tab(info, -1);

  xjjroot::print_tab(h1s, 0);

  h1s["y_corr"] = (TH1D*)h1s.at("y_yield")->Clone(xjjc::str_replaceall(h1s.at("y_yield")->GetName(), "yield", "corr").c_str());
  h1s["y_corr"]->Divide(h1s.at("y_eff__rebin"));
  h1s["y_corr"]->GetYaxis()->SetTitle("Corrected Yield");
  xjjana::sethminmax(h1s.at("y_corr"), 0., 1.8);
  h1s["y_xsec"]	= (TH1D*)h1s.at("y_corr")->Clone(xjjc::str_replaceall(h1s.at("y_corr")->GetName(), "corr", "xsec").c_str());
  h1s["y_xsec"]->Scale(0.5/global::BR_DtoKpi/(lumi*1.e6)/3/*dpt: 5-2*/, "width");
  h1s["y_xsec"]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{y}d#it{p}_{T}} [mb/GeV]");
  xjjana::sethminmax(h1s.at("y_xsec"), 0., 1.7);
  xjjroot::sethempty(h1s.at("y_xsec"), 0, 0.2);

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputdir + "/" + tag + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" }});

  auto draw_global = [&infos, &lumi]() {
    xjjroot::drawCMS(xjjroot::CMS::internal, xjjc::str_replaceall(infos.at("raw_data").at("input_tex"), ")", Form(", %.1f nb^{-1})", lumi*1.e3)));
    xjjroot::drawtexgroup(0.24, 0.86, { "2 < p_{T} < 5 GeV", infos.at("raw_data").at("cut_tex") }, 0.04, 13);
  };
  
  for (auto& t : { "y_yield", "y_corr", "y_xsec" }) {
    pdf->prepare();
    h1s.at(t)->Draw("pe1");
    draw_global();
    pdf->write(name_png + xjjc::str_replaceall(t, "y_", "_") + ".pdf");
  }

  xjjana::sethminmax(h1s.at("y_xsec"), 0, 2.5);
  xjjroot::setthgrstyle(h1s.at("y_xsec"), xjjroot::mycolor_satmiddle.at("red"), -1, -1, xjjroot::mycolor_satmiddle.at("red"));
  pdf->prepare();
  h1s.at("y_xsec")->Draw("axis");
  auto* g_HIN_25_002 = measurement::draw_HIN_25_002(event_is);
  h1s.at("y_xsec")->Draw("pe1 same");
  draw_global();
  auto* legvs23 = new TLegend(0.55, 0.75-2*0.042, 0.85, 0.75);
  xjjroot::setleg(legvs23, 0.04);
  legvs23->AddEntry(h1s.at("y_xsec"), "This analysis", "p");
  legvs23->AddEntry(g_HIN_25_002, "PAS-HIN-25-002", "pf");
  legvs23->Draw();
  pdf->write(name_png + "_xsec_vs23.pdf");

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

