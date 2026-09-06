#include <TH2D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/draw.h"

int macro(const std::string& inputname, const std::string& outputname, int save_png = 0) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! failed to open input file: " << inputname << ", abort." << std::endl;
    return 1;
  }
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);

  std::map<std::string, TH2D*> h2s;
  std::map<std::string, std::vector<TH1D*>> h1pts;

  for (auto* h2 : xjjana::getobj_regexp<TH2D>(inf, "h2_.+")) {
    const auto name = xjjc::str_eraseall(h2->GetName(), "h2_");
    h2->Sumw2();
    h2s[name] = h2;
    if (!xjjc::str_contains(name, "eff")) continue;
    
    for (int j=0; j<h2->GetYaxis()->GetNbins(); j++) {
      auto* h_y = h2->ProjectionX(Form("h1_y_%s__pt-%d", name.c_str(), j),
                                  j+1, j+1, "e");
      h1pts[name + "-y"].push_back(h_y);
    }
  }
  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1pts, 0);

  const auto tbins = draw::bintex(h2s["bins"], 0, 1);
  const auto npt = tbins.npt(), ny = tbins.ny();
  
  auto make_eff = [&h2s, &h1pts](const std::string& name_new, const std::string& name_num, const std::string& name_den,
                                 const std::string& title) {
    __XJJLOG << ">> " << name_new << std::endl;
    auto *h2eff_num = h2s.at(name_num), *h2eff_den = h2s.at(name_den);
    auto* h2eff = (TH2D*)h2eff_num->Clone(xjjc::str_replaceall(h2eff_num->GetName(), name_num, name_new).c_str());
    h2eff->Divide(h2eff_den);
    h2eff->GetZaxis()->SetTitle(title.c_str());
    h2s[name_new] = h2eff;
    
    for (int j=0; j<h1pts.at(name_num + "-y").size(); j++) {
      auto *h1eff_num = h1pts.at(name_num + "-y")[j], *h1eff_den = h1pts.at(name_den + "-y")[j];
      auto* h1eff = (TH1D*)h1eff_num->Clone(xjjc::str_replaceall(h1eff_num->GetName(), name_num, name_new).c_str());
      h1eff->Divide(h1eff_den);
      h1eff->GetYaxis()->SetTitle(title.c_str());
      h1pts[name_new + "-y"].push_back(h1eff);
    }
  };
  make_eff("evteff", "evteff_num", "evteff_den", "#epsilon_{event selections}");

  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1pts, 0);

  for (auto& [_, h] : h2s) {
    xjjroot::sethempty(h, 0, -0.1);
    h->GetZaxis()->SetTitleSize(0);
  }

  // const auto colors = xjjroot::grayscales_color(npt, xjjroot::mycolor_dark["red"]);
  const auto alphas = xjjroot::grayscales_alpha(npt);
  for (auto& [name, hh] : h1pts) {
    for (int i=0; i<hh.size(); i++) {
      const auto cc = kBlack;
      // const auto cc = hh.size()==1 ? kBlack : xjjroot::mycolor_dark["red"];
      xjjroot::setthgrstyle(hh[i], cc, xjjroot::markerlist_solid[i%xjjroot::markerlist_solid.size()], 1.5, cc, 1, 1, -1, -1, -1, alphas[i], alphas[i]);
      xjjroot::sethempty(hh[i], 0, 0.3);
    }
  }

  auto draw_global = [&info]() {
    xjjroot::drawCMS(std::atoi(info.at("is_mc").c_str()) ? xjjroot::CMS::simulation : xjjroot::CMS::internal, info.at("input_tex"));
  };
  
  xjjroot::setgstyle(1, 2, xjjroot::Colz);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto png_name = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" } });
  
  // TH2
  for (std::string name : { "evteff" } ) {
    pdf->prepare();
    h2s[name]->Draw("colz");
    xjjroot::drawtexgroup(0.18, 0.85, {
        h2s.at(name)->GetZaxis()->GetTitle(),
      }, 0.04, 13);
    draw_global();
    pdf->write(png_name + "_" + name + "-pt-y.pdf", save_png ? "" : "X");
  }

  // TH1
  xjjroot::setcstyle(pdf->getc(), 1, xjjroot::Standard);
  xjjroot::setgstyle(1, 2, xjjroot::Standard);
  auto* leg = new TLegend(0.20, 0.85-npt*0.042, 0.6, 0.85);
  xjjroot::setleg(leg, 0.038);
  for (int i=0; i<npt; i++)
    leg->AddEntry(h1pts.at("evteff-y")[i], tbins.label_pt(i).c_str(), "p");
  
  pdf->prepare();
  gPad->Modified();
  gPad->Update();

  for (auto& name : std::vector<std::string>{ "evteff" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts[name + "-y"], 0);
    xjjana::sethsmax(h1pts[name + "-y"], 1.5);
    h1pts[name + "-y"].front()->Draw("axis");
    // xjjroot::drawtexgroup(xjjroot::get_pad_center(), 0.24, {
    //     ts_formula.at(name)
    //   }, 0.038, 21, 42, 1.2, 1, 1, { 16 });
    for (auto& h : h1pts[name + "-y"])
      h->Draw("pe1 same");
    draw_global();
    leg->Draw();
    pdf->write(png_name + "_" + name + "-y.pdf", save_png ? "" : "X");
  }

  // pdf->draw_cover( {
  //     "#bf{eff num} " + info.at("cut_eff_num"),
  //     "#bf{effreco num} " + info.at("cut_reco_num"),
  //     "#bf{acc num} " + info.at("cut_acc_num"),
  //     "#bf{eff den} " + info.at("cut_eff_den"),
  //     "#bf{File} " + info.at("inputmc")
  //   }, 0.03);

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (auto& [_, h] : h2s) xjjroot::writehist(h);
  for (auto& [_, hh] : h1pts)
    for (auto& h : hh)
      xjjroot::writehist(h);
  auto* t = new TTree("info", "");
  for (auto& [key, content] : info) {
    t->Branch(key.c_str(), &content);
  }
  t->Fill();
  t->Write();
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4)
    return macro(argv[1], argv[2], atoi(argv[3]));
  if (argc == 3)
    return macro(argv[1], argv[2]);
  return 1;
}

