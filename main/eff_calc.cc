#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#define __BINS_PTY_ANA__
#include "../include/bins.h"

#include "../include/draw.h"

int macro(const std::string& inputname, const std::string& outputname, int save_png = 1) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);

  __XJJLOG << "++ y binning" << std::endl;
  xjjc::print_vec_h(bins::ybins);
  __XJJLOG << "++ pt binning" << std::endl;
  xjjc::print_vec_h(bins::ptbins);
  
  std::map<std::string, TH3D*> h3s;
  std::map<std::string, TH2D*> h2s;
  std::map<std::string, std::vector<TH1D*>> h1pts, h1ys;

  for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_.+")) {
    auto name = xjjc::str_eraseall(h3->GetName(), "h3_");
    h3->Sumw2();
    h3s[name] = h3;
    h3->GetZaxis()->SetRange(0, h3->GetZaxis()->GetNbins() + 1); // include multiplicity overflow
    h2s[name + "-y-pt"] = (TH2D*)h3->Project3D("yx");
    h2s[name + "-y-pt"]->SetName(Form("h2_y-pt_%s", name.c_str()));
    h2s[name + "-y-pt__rebin"] = xjjana::rebin(h2s[name + "-y-pt"],
                                               bins::ybins.size()-1, bins::ybins.data(),
                                               bins::ptbins.size()-1, bins::ptbins.data(),
                                               Form("%s__rebin", h2s[name + "-y-pt"]->GetName()));
    for (int j=0; j<h2s.at(name+"-y-pt__rebin")->GetYaxis()->GetNbins(); j++) {
      auto* h_y = h2s.at(name + "-y-pt__rebin")->ProjectionX(Form("h1_y_%s__rebin__pt-%d", name.c_str(), j),
                                                             j+1, j+1, "e");
      h1pts[name + "-y__rebin"].push_back(h_y);
    }
    for (int i=0; i<h2s.at(name+"-y-pt__rebin")->GetXaxis()->GetNbins(); i++) {
      auto* h_pt = h2s.at(name + "-y-pt__rebin")->ProjectionY(Form("h1_pt_%s__rebin__y-%d", name.c_str(), i),
                                                              i+1, i+1, "e");
      h1ys[name + "-pt__rebin"].push_back(h_pt);
    }
  }
  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1pts, 0);
  if (h1ys.begin()->second.size() > 1)
    xjjroot::print_tab(h1ys, 0);

  const auto npt = bins::ptbins.size() - 1, ny = bins::ybins.size() - 1;
  auto nbin_consist = [](const std::map<std::string, std::vector<TH1D*>>& h1s, int nbin) {
    for (const auto& [_, hh] : h1s) {
      if (hh.size() != nbin)
        return 3;
    }
    return 0;
  };
  if (nbin_consist(h1pts, npt)) return 3;
  if (nbin_consist(h1ys, ny)) return 3;
  
  std::map<std::string, std::string> ts_formula;
  auto make_eff = [&h2s, &h1pts, &h1ys, &ts_formula](const std::string& name_new, const std::string& name_num, const std::string& name_den,
                                                     const std::string& title, const std::string& formula) {
    __XJJLOG << ">> " << name_new << std::endl;
    for (const std::string& suffix : { "-y-pt", "-y-pt__rebin" }) {
      h2s[name_new + suffix] = (TH2D*)h2s.at(name_num + suffix)->Clone(xjjc::str_replaceall(h2s.at(name_num + suffix)->GetName(), name_num, name_new).c_str());
      h2s[name_new + suffix]->Divide(h2s.at(name_den + suffix));
      h2s[name_new + suffix]->GetZaxis()->SetTitle(title.c_str());
    }
    for (int j=0; j<h1pts.at(name_num + "-y__rebin").size(); j++) {
      auto* heff_j = (TH1D*)h1pts.at(name_num+"-y__rebin")[j]->Clone(xjjc::str_replaceall(h1pts.at(name_num+"-y__rebin")[j]->GetName(), name_num, name_new).c_str());
      heff_j->Divide(h1pts.at(name_den+"-y__rebin")[j]);
      heff_j->GetYaxis()->SetTitle(title.c_str());
      h1pts[name_new + "-y__rebin"].push_back(heff_j);
    }
    for (int j=0; j<h1ys.at(name_num + "-pt__rebin").size(); j++) {
      auto* heff_j = (TH1D*)h1ys.at(name_num+"-pt__rebin")[j]->Clone(xjjc::str_replaceall(h1ys.at(name_num+"-pt__rebin")[j]->GetName(), name_num, name_new).c_str());
      heff_j->Divide(h1ys.at(name_den+"-pt__rebin")[j]);
      heff_j->GetYaxis()->SetTitle(title.c_str());
      h1ys[name_new + "-pt__rebin"].push_back(heff_j);
    }
    ts_formula[name_new] = formula;
  };

  const std::string Nreco = "N_{reco}^{gen-matched}#scale[0.2]{ }", Ngen = "N_{gen}^{signal}#scale[0.2]{ }";
  make_eff("eff", "eff_num", "eff_den", xjjroot::CMS::DzDzbar + "#scale[0.4]{ }#LT#alpha#scale[0.4]{ }#times#scale[0.4]{ }#epsilon_{reco}#scale[0.4]{ }#times#scale[0.4]{ }#epsilon_{sel}#GT",
           "#frac{"+Nreco+"(All selections)}{"+Ngen+"}");
  make_eff("effsel", "eff_num", "reco_num", xjjroot::CMS::DzDzbar + "#scale[0.4]{ }#LT#epsilon_{sel}#GT",
           "#frac{"+Nreco+"(All selections)}{"+Nreco+"(#it{p}_{T}^{trk} > 0.5 GeV, |#eta^{trk}| < 2.4)}");
  make_eff("effreco", "reco_num", "acc_num", xjjroot::CMS::DzDzbar + "#scale[0.4]{ }#LT#epsilon_{reco}#GT",
           "#frac{"+Nreco+"(#it{p}_{T}^{trk} > 0.5 GeV, |#eta^{trk}| < 2.4)}{"+Ngen+"(#it{p}_{T}^{trk} > 0.5 GeV, |#eta^{trk}| < 2.4)}");
  make_eff("acc", "acc_num", "eff_den", xjjroot::CMS::DzDzbar + "#scale[0.4]{ }#LT#alpha#GT",
           "#frac{"+Ngen+"(#it{p}_{T}^{trk} > 0.5 GeV, |#eta^{trk}| < 2.4)}{"+Ngen+"}");

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

  auto* h2_bins = (TH2D*)h2s["eff-y-pt__rebin"]->Clone("h2_bins");
  h2_bins->Reset("ICESM");
  const auto tbins = draw::bintex(h2_bins, 0, 1);
  auto draw_global = [&h2s, &info](bool drawpt = false) {
    xjjroot::drawCMS(xjjroot::CMS::simulation, info.at("inputmc_tex"));
  };
  
  xjjroot::setgstyle(1, 2, xjjroot::Colz);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto png_name = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" } });
  
  // TH2
  for (std::string name : { "eff", "acc", "effreco", "effsel" } ) {
    pdf->prepare();
    h2s[name + "-y-pt"]->Draw("colz");
    xjjroot::drawtexgroup(0.18, 0.85, {
        h2s.at(name + "-y-pt")->GetZaxis()->GetTitle(),
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
    leg->AddEntry(h1pts.at("eff-y__rebin")[i], tbins.label_pt(i).c_str(), "p");
  
  pdf->prepare();
  gPad->Modified();
  gPad->Update();

  for (auto& name : std::vector<std::string>{ "eff", "acc", "effreco", "effsel" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts[name + "-y__rebin"], 0);
    xjjana::sethsmax(h1pts[name + "-y__rebin"], 1.5);
    h1pts[name + "-y__rebin"].front()->Draw("axis");
    xjjroot::drawtexgroup(xjjroot::get_pad_center(), 0.24, {
        ts_formula.at(name)
      }, 0.038, 21, 42, 1.2, 1, 1, { 16 });
    for (auto& h : h1pts[name + "-y__rebin"])
      h->Draw("pe1 same");
    draw_global();
    leg->Draw();
    pdf->write(png_name + "_" + name + "-y.pdf", save_png ? "" : "X");
  }

  pdf->draw_cover( {
      "#bf{eff num} " + info.at("cut_eff_num"),
      "#bf{effreco num} " + info.at("cut_reco_num"),
      "#bf{acc num} " + info.at("cut_acc_num"),
      "#bf{eff den} " + info.at("cut_eff_den"),
      "#bf{File} " + info.at("inputmc")
    }, 0.03);

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (auto& [_, h] : h2s) xjjroot::writehist(h);
  for (auto& [_, hh] : h1pts)
    for (auto& h : hh)
      xjjroot::writehist(h);
  xjjroot::writehist(h2_bins);
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
  if (argc >= 5) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[3], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[4], ",");
    if (argc == 6)
      return macro(argv[1], argv[2], atoi(argv[5]));
    if (argc == 5)
      return macro(argv[1], argv[2]);
  }
  return 1;
}

