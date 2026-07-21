#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#define __BINS_PTY_ANA__
#include "../include/bins.h"

#include "../include/draw.h"

namespace eff {
  std::pair<TH1D*, TH1D*> sum_norm_byweight(TH2D* hreveff, TH2D* hweight, TH1D* hrefbin);
}

int macro(std::string inputname, int save_png = 1) {
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
    
    // h1s[name + "-y"] = h3->ProjectionX(Form("h1_y_%s", name.c_str()),
    //                                    1, h3->GetYaxis()->GetNbins(), // pt only in analysis range
    //                                    0, h3->GetZaxis()->GetNbins()+1, // multiplicity overflow
    //                                    "e");
    // h1s[name + "-y__rebin"] = (TH1D*)h1s[name + "-y"]->Rebin(bins::ybins.size()-1, Form("%s__rebin", h1s[name + "-y"]->GetName()), bins::ybins.data());
    
  }
  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1pts, 0);

  auto make_eff = [&h2s, &h1pts](std::string name_new, std::string name_num, std::string name_den, std::string title) {
    __XJJLOG << ">> " << name_new << std::endl;
    for (const std::string& suffix : { "-y-pt", "-y-pt__rebin" }) {
      std::cout<<suffix<<std::endl;
      h2s[name_new + suffix] = (TH2D*)h2s.at(name_num + suffix)->Clone(xjjc::str_replaceall(h2s.at(name_num + suffix)->GetName(), name_num, name_new).c_str());
      h2s[name_new + suffix]->Divide(h2s.at(name_den + suffix));
      h2s[name_new + suffix]->GetZaxis()->SetTitle(title.c_str());
    }
    // h1s[name_new + "-y"] = (TH1D*)h1s.at(name_num + "-y")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y")->GetName(), name_num, name_new).c_str());
    // h1s[name_new + "-y"]->Divide(h1s.at(name_den + "-y"));
    // h1s[name_new + "-y"]->GetYaxis()->SetTitle(title.c_str());
    // h1s[name_new + "-y__rebin"] = (TH1D*)h1s.at(name_num + "-y__rebin")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y__rebin")->GetName(), name_num, name_new).c_str());
    // h1s[name_new + "-y__rebin"]->Divide(h1s.at(name_den + "-y__rebin"));
    // h1s[name_new + "-y__rebin"]->GetYaxis()->SetTitle(title.c_str());
    for (int j=0; j<h1pts.at(name_num + "-y__rebin").size(); j++) {
      auto* heff_j = (TH1D*)h1pts.at(name_num+"-y__rebin")[j]->Clone(xjjc::str_replaceall(h1pts.at(name_num+"-y__rebin")[j]->GetName(), name_num, name_new).c_str());
      heff_j->Divide(h1pts.at(name_den+"-y__rebin")[j]);
      heff_j->GetYaxis()->SetTitle(title.c_str());
      h1pts[name_new + "-y__rebin"].push_back(heff_j);
    }
    // for (int i=0; i<h1ys.at(name_num+"-pt").size(); i++) {
    //   auto* heff_i = (TH1D*)h1ys.at(name_num+"-pt")[i]->Clone(xjjc::str_replaceall(h1ys.at(name_num+"-pt")[i]->GetName(), name_num, name_new).c_str());
    //   heff_i->Divide(h1ys.at(name_den+"-pt")[i]);
    //   heff_i->GetYaxis()->SetTitle(title.c_str());
    //   h1ys[name_new + "-pt"].push_back(heff_i);
    // }
  };

  make_eff("eff", "eff_num", "eff_den", xjjroot::CMS::DzDzbar + "#scale[0.5]{ }#LT#alpha#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{reco}#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{sel}#GT");
  make_eff("reveff", "eff_den", "eff_num", xjjroot::CMS::DzDzbar + " 1 /#scale[0.5]{ }#LT#alpha#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{reco}#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{sel}#GT");
  make_eff("effreco", "reco_num", "eff_den", xjjroot::CMS::DzDzbar + "#scale[0.5]{ }#LT#alpha#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{reco}#GT");
  make_eff("effsel", "eff_num", "reco_num", xjjroot::CMS::DzDzbar + "#scale[0.5]{ }#LT#epsilon_{sel}#GT");

  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1pts, 0);

  for (auto& [_, hh] : h1pts) {
    if (hh.size() != bins::ptbins.size()-1)
      return 3;
  }
  const auto npt = bins::ptbins.size() - 1;
  for (auto& [_, h] : h2s) {
    xjjroot::sethempty(h, 0, -0.1);
    h->GetZaxis()->SetTitleSize(0);
  }

  for (auto& [name, hh] : h1pts) {
    for (int i=0; i<hh.size(); i++) {
      xjjroot::setthgrstyle(hh[i], xjjroot::color_alpha(kBlack, 1-i*0.15), 21, 1.4, xjjroot::color_alpha(kBlack, 1-i*0.15), 1, 1);
      xjjroot::sethempty(hh[i], 0, 0.3);
    }
  }

  const auto tbins = draw::bintex(h2s["eff-y-pt"], 0, 1);
  auto draw_global = [&h2s, &info](bool drawpt = false) {
    xjjroot::drawCMS(xjjroot::CMS::simulation, info.at("inputmc_tex"));
  };
  
  xjjroot::setgstyle(1, 2, xjjroot::Colz);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(inputname, { { "saveeff", "calceff" }, { "rootfiles/", "figspdf/" }, { ".root", ".pdf" } }));
  auto png_name = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" } });

  // TH2
  for (std::string name : { "eff", "effreco", "effsel" } ) {
    pdf->prepare();
    h2s[name + "-y-pt"]->Draw("colz");
    xjjroot::drawtexgroup(0.18, 0.85, {
        h2s.at(name + "-y-pt")->GetZaxis()->GetTitle()
      }, 0.04, 13);
    draw_global();
    pdf->write(png_name + "_" + name + "-pt-y.pdf", save_png ? "" : "X");
  }

  // TH1
  xjjroot::setcstyle(pdf->getc(), 1, xjjroot::Standard);
  xjjroot::setgstyle(1, 2, xjjroot::Standard);
  auto* leg = new TLegend(0.20, 0.3-npt*0.042, 0.6, 0.3);
  xjjroot::setleg(leg, 0.038);
  for (int i=0; i<npt; i++)
    leg->AddEntry(h1pts.at("eff-y__rebin")[i], tbins.label_pt(i).c_str(), "p");
  
  pdf->prepare();
  gPad->Modified();
  gPad->Update();

  for (auto& name : std::vector<std::string>{ "eff", "effreco", "effsel" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts[name + "-y__rebin"], 0);
    xjjana::sethsmax(h1pts[name + "-y__rebin"], 1.7);
    h1pts[name + "-y__rebin"].front()->Draw("axis");
    for (auto& h : h1pts[name + "-y__rebin"])
      h->Draw("pe1 same");
    draw_global();
    leg->Draw();
    pdf->write(png_name + "_" + name + ".pdf", save_png ? "" : "X");
  }

  pdf->close();

  auto* outf = xjjroot::newfile(xjjc::str_replaceall(inputname, "saveeff", "calceff"));
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
  outf->Close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[2], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[3], ",");
    return macro(argv[1], atoi(argv[4]));
  }
  return 1;
}

std::pair<TH1D*, TH1D*> eff::sum_norm_byweight(TH2D* hreveff, TH2D* hweight, TH1D* hrefbin) {
  auto nbinx = hrefbin->GetXaxis()->GetNbins();
  std::vector<double> sum_y(nbinx, 0), norm_y(nbinx, 0);
  auto* hsum = (TH1D*)hrefbin->Clone(Form("%s__sum", hweight->GetName()));
  auto* hnorm = (TH1D*)hrefbin->Clone(Form("%s__norm", hweight->GetName()));
  hsum->Reset(); hnorm->Reset();
  for (int i=0; i<hreveff->GetXaxis()->GetNbins(); i++) {
    auto ybin = hrefbin->FindBin(hreveff->GetXaxis()->GetBinCenter(i+1));
    if (ybin < 1 || ybin > nbinx) {
      __XJJLOG << "?? bad FindBin (" << ybin << "). skip." << std::endl;
      __XJJLOG << "   >> " << hreveff->GetXaxis()->GetBinCenter(i+1) << " -> " << ybin <<std::endl;
      continue;
    }
    for (int j=0; j<hreveff->GetYaxis()->GetNbins(); j++) {
      auto reveff = hreveff->GetBinContent(i+1, j+1),
        reveff_e = hreveff->GetBinError(i+1, j+1),
        nweight = hweight->GetBinContent(i+1, j+1),
        nweight_e = hweight->GetBinError(i+1, j+1);

      if (reveff == 0) {
        __XJJLOG << "!! reveff is 0 for the bin (" << i+1 <<", " << j+1 << ")" << std::endl;
        continue;
      }

      double sume = hsum->GetBinError(ybin);
      hsum->SetBinContent(ybin, hsum->GetBinContent(ybin) + reveff*nweight);
      hsum->SetBinError(ybin, std::sqrt(sume*sume + std::pow(nweight * reveff_e, 2) + std::pow(reveff * nweight_e, 2)));
      double norme = hnorm->GetBinError(ybin);
      hnorm->SetBinContent(ybin, hnorm->GetBinContent(ybin) + nweight);
      hnorm->SetBinError(ybin, std::sqrt(norme*norme + nweight_e*nweight_e));
    }
  }
  return std::pair<TH1D*, TH1D*>{ hsum, hnorm };
}
