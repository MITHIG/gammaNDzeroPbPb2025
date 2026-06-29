#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"

#define __BINS_PTY_EQ__
#include "../include/bins.h"

int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);

  auto ybin_analysis = xjjc::fixedbin_to_edges(bins::ny, bins::miny, bins::maxy);
  __XJJLOG << "++ analysis binning" << std::endl;
  xjjc::print_vec_h(ybin_analysis);
  
  std::map<std::string, TH3D*> h3s;
  // std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, TH2D*> h2s;

  for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_.+")) {
    auto name = xjjc::str_eraseall(h3->GetName(), "h3_");
    h3->Sumw2();
    h3s[name] = h3;
    h3->GetZaxis()->SetRange(0, h3->GetZaxis()->GetNbins() + 1); // include overflow
    h2s[name + "-y-pt"] = (TH2D*)h3->Project3D("yx");
    h2s[name + "-y-pt"]->SetName(Form("h2_y-pt_%s", name.c_str()));
    h1s[name + "-y"] = h3->ProjectionX(Form("h1_y_%s", name.c_str()),
                                       1, h3->GetYaxis()->GetNbins(),
                                       0, h3->GetZaxis()->GetNbins()+1, // overflow
                                       "e");
    h1s[name + "-y__rebin"] = (TH1D*)h1s[name + "-y"]->Rebin(ybin_analysis.size()-1, Form("%s__rebin", h1s[name + "-y"]->GetName()), ybin_analysis.data());
    // h2s[name + "-y-pt__rebin"] = (TH2D*)h2s[name + "-y-pt"]->Rebin(ybin_analysis.size()-1, Form("%s__rebin", h2s[name + "-y-pt"]->GetName()), ybin_analysis.data());
    // for (int i=0; i<h2s.at(name+"-y-pt__rebin")->GetXaxis()->GetNbins(); i++) {
    //   auto* hpt = h2s.at(name + "-y-pt")->ProjectionY(Form("h1_pt_%s__rebin__y-%d", name.c_str(), i),
    //                                                   i+1, i+1, "e");
    //   h1ys[name + "-pt"].push_back(hpt);
    // }
  }
  
  auto make_eff = [&h2s, &h1s](std::string name_new, std::string name_num, std::string name_den, std::string title) {
    __XJJLOG << ">> " << name_new << std::endl;
    h2s[name_new + "-y-pt"] = (TH2D*)h2s.at(name_num + "-y-pt")->Clone(xjjc::str_replaceall(h2s.at(name_num + "-y-pt")->GetName(), name_num, name_new).c_str());
    h2s[name_new + "-y-pt"]->Divide(h2s.at(name_den + "-y-pt"));
    h2s[name_new + "-y-pt"]->GetZaxis()->SetTitle(title.c_str());
    h1s[name_new + "-y"] = (TH1D*)h1s.at(name_num + "-y")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y")->GetName(), name_num, name_new).c_str());
    h1s[name_new + "-y"]->Divide(h1s.at(name_den + "-y"));
    h1s[name_new + "-y"]->GetYaxis()->SetTitle(title.c_str());
    h1s[name_new + "-y__rebin"] = (TH1D*)h1s.at(name_num + "-y__rebin")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y__rebin")->GetName(), name_num, name_new).c_str());
    h1s[name_new + "-y__rebin"]->Divide(h1s.at(name_den + "-y__rebin"));
    h1s[name_new + "-y__rebin"]->GetYaxis()->SetTitle(title.c_str());
    // for (int i=0; i<h1ys.at(name_num+"-pt").size(); i++) {
    //   auto* heff_i = (TH1D*)h1ys.at(name_num+"-pt")[i]->Clone(xjjc::str_replaceall(h1ys.at(name_num+"-pt")[i]->GetName(), name_num, name_new).c_str());
    //   heff_i->Divide(h1ys.at(name_den+"-pt")[i]);
    //   heff_i->GetYaxis()->SetTitle(title.c_str());
    //   h1ys[name_new + "-pt"].push_back(heff_i);
    // }
  };

  make_eff("eff", "eff_num", "eff_den", xjjroot::CMS::D0 + "#scale[0.5]{ }#LT#alpha#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{reco}#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{sel}#GT");
  make_eff("effreco", "reco_num", "eff_den", xjjroot::CMS::D0 + "#scale[0.5]{ }#LT#alpha#scale[0.5]{ }#times#scale[0.5]{ }#epsilon_{reco}#GT");
  make_eff("effsel", "eff_num", "reco_num", xjjroot::CMS::D0 + "#scale[0.5]{ }#LT#epsilon_{sel}#GT");

  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1s, 0);
  // xjjroot::print_tab(h1ys, 0);

  if (h3s.find("data_signalwin") != h3s.end()) {
    __XJJLOG << "++ data driven eff" << std::endl;
  
    auto *heff = h2s.at("eff-y-pt"), *hdata = h2s.at("data_signalwin-y-pt");

    auto make_effdata = [&h1s, &heff, &hdata](std::string name_refbin) {
      auto name_new = xjjc::str_replaceall(name_refbin, "eff", "effdata");
      h1s[name_new] = (TH1D*)h1s.at(name_refbin)->Clone(xjjc::str_replaceall(h1s.at(name_refbin)->GetName(), "eff", "effdata").c_str());
      h1s[name_new]->Reset();
      std::vector<double> sum_y(h1s.at(name_new)->GetXaxis()->GetNbins(), 0), norm_y(sum_y.size(), 0);
      for (int i=0; i<heff->GetXaxis()->GetNbins(); i++) {
        for (int j=0; j<heff->GetYaxis()->GetNbins(); j++) {
          auto eff = heff->GetBinContent(i+1, j+1),
            eff_e = heff->GetBinError(i+1, j+1),
            ndata = hdata->GetBinContent(i+1, j+1),
            ndata_e = hdata->GetBinError(i+1, j+1);
            
          if (eff == 0) {
            __XJJLOG << "!! eff is 0 for the bin (" << i+1 <<", " << j+1 << ")" << std::endl;
            continue;
          }
          auto ybin = h1s.at(name_new)->FindBin(heff->GetXaxis()->GetBinCenter(i+1));
          if (ybin < 1 || ybin > h1s.at(name_new)->GetXaxis()->GetNbins()) {
            __XJJLOG << "?? bad FindBin (" << ybin << "). skip." << std::endl;
            __XJJLOG << "   >> " << heff->GetXaxis()->GetBinCenter(i+1) << " -> " << ybin <<std::endl;
            continue;
          }
          sum_y[ybin-1] += (1./eff)*ndata;
          norm_y[ybin-1] += ndata;
        }
      }
      
      for (int i=0; i<h1s.at(name_new)->GetXaxis()->GetNbins(); i++) {
        if (sum_y[i] == 0) {
          __XJJLOG << "?? ndata is 0 for the bin (" << i+1 << "). skip." << std::endl;
          continue;
        }
        
        h1s.at(name_new)->SetBinContent(i+1, 1./(sum_y[i]/norm_y[i]));
        h1s.at(name_new)->SetBinError(i+1, 0.001); // ! how to get the uncertainty ?
      }
    };

    make_effdata("eff-y");
    make_effdata("eff-y__rebin");
  }

  for (auto& [_, h] : h2s) {
    xjjroot::sethempty(h, 0, -0.1);
    h->GetZaxis()->SetTitleSize(0);
  }
  for (auto& [name, h] : h1s) {
    auto is_data = xjjc::str_contains(name, "data");
    int color = is_data ? xjjroot::mycolor_satmiddle["red"] : kBlack;
    xjjroot::setthgrstyle(h, color, is_data ? 20 : 21, 1.4, color, 1, 1);
    xjjroot::sethempty(h, 0, 0.3);
  }
  // for (auto& [_, hh] : h1ys) {
  //   for (int i=0; i<hh.size(); i++) {
  //     xjjroot::sethempty(hh[i], 0, 0.3);
  //     xjjroot::setthgrstyle(hh[i], xjjroot::colorlist_middle[i], 21, 1.3, xjjroot::colorlist_middle[i], 1, 1);
  //   }
  // }

  auto draw_global = [&h2s, &info](bool drawpt = true) {
    xjjroot::drawCMS(xjjroot::CMS::internal, info.at("inputmc_tex"));
    if (drawpt)
      xjjroot::drawtexgroup(1-gStyle->GetPadRightMargin()-0.06, 0.85, {
          xjjc::number_range_string(h2s.at("eff-y-pt")->GetYaxis()->GetXmin(), h2s.at("eff-y-pt")->GetYaxis()->GetXmax(), "p_{T}") + " GeV",
        }, 0.04, 33);
  };
  
  xjjroot::setgstyle(1, 2, xjjroot::Colz);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(inputname, { { "saveeff", "calceff" }, { "rootfiles/", "figspdf/" }, { ".root", ".pdf" } }));
  auto png_name = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" } });

  for (std::string name : { "eff", "effreco", "effsel" } ) {
    pdf->prepare();
    h2s[name + "-y-pt"]->Draw("colz");
    xjjroot::drawtexgroup(0.18, 0.85, {
        h2s.at(name + "-y-pt")->GetZaxis()->GetTitle()
      }, 0.04, 13);
    draw_global(false);
    pdf->write(png_name + "_" + name + "-pt-y.pdf");
  }

  xjjroot::setcstyle(pdf->getc(), 1, xjjroot::Standard);
  xjjroot::setgstyle(1, 2, xjjroot::Standard);
  
  auto* leg = new TLegend(0.20, 0.3-2*0.042, 0.6, 0.3);
  xjjroot::setleg(leg, 0.04);
  leg->AddEntry(h1s["eff-y__rebin"], "Directly from MC", "p");
  leg->AddEntry(h1s["effdata-y__rebin"], "MC eff + data signal region kinematics", "p");
  
  pdf->prepare();
  gPad->Modified();
  gPad->Update();
  // xjjana::sethabsminmax(h1s["eff-y__rebin"], 0, .2);
  xjjana::sethminmax(h1s["eff-y__rebin"], 0, 1.4);
  h1s["eff-y__rebin"]->Draw("pe1");
  h1s["effdata-y__rebin"]->Draw("pe1 same");
  leg->Draw();
  draw_global();
  pdf->write(png_name + "_eff.pdf");

  pdf->prepare();
  // xjjana::sethabsminmax(h1s["effreco-y__rebin"], 0, 1.);
  xjjana::sethminmax(h1s["effreco-y__rebin"], 0, 1.4);
  h1s["effreco-y__rebin"]->Draw("pe1");
  draw_global();
  pdf->write(png_name + "_effreco.pdf");

  pdf->prepare();
  // xjjana::sethabsminmax(h1s["effsel-y__rebin"], 0, .2);
  xjjana::sethminmax(h1s["effsel-y__rebin"], 0, 1.4);
  h1s["effsel-y__rebin"]->Draw("pe1");
  draw_global();
  pdf->write(png_name + "_effsel.pdf");

  pdf->close();

  auto* outf = xjjroot::newfile(xjjc::str_replaceall(inputname, "saveeff", "calceff"));
  for (auto& [_, h] : h2s) xjjroot::writehist(h);
  for (auto& [_, h] : h1s) xjjroot::writehist(h);
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
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

