#include <TH3D.h>
#include "xjjanauti.h"

#define __BINS_PTY_EQ__
#include "../include/bins.h"

int macro(std::string inputname) {
  std::cout<<std::endl;

  auto* inf = TFile::Open(inputname.c_str());
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);

  auto ybin_analysis = xjjana::fixedbin_to_edges(bins::ny, bins::miny, bins::maxy);
  
  std::map<std::string, TH3D*> h3s;
  std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, TH1D*> h2s;

  for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_.+")) {
    auto name = xjjc::str_eraseall(h3->GetName(), "h3_");
    h3->Sumw2();
    h3s[name] = h3;
    h3->GetZaxis()->SetRange(0, h3->GetZaxis()->GetNbins() + 1); // include overflow
    h2s[name + "-y-pt"] = (TH2D*)h3->Project3D("xy");
    h2s[name + "-y-pt"]->SetName(Form("h1_y-pt_%s", name.c_str()));
    h1s[name + "-y"] = h3->ProjectionX(Form("h1_y_%s", name.c_str()),
                                       1, h3->GetYaxis()->GetNbins(),
                                       0, h3->GetZaxis()->GetNbins()+1, // overflow
                                       "e");
    h1s[name + "-y__rebin"] = (TH1D*)h1s[name + "-y"]->Rebin(ybin_analysis.size()-1, Form("%s__rebin", h1s[name + "-y"]->GetName()), ybin_analysis.data());
  }
  
  auto make_eff = [&h2s, &h1s](std::string name_new, std::string name_num, std::string name_den, std::string title) {
    h2s[name_new + "-y-pt"] = (TH2D*)h2s.at(name_num + "-y-pt")->Clone(xjjc::str_replaceall(h2s.at(name_num + "-y-pt")->GetName(), name_num, name_new).c_str());
    h2s[name_new + "-y-pt"]->Divide(h2s.at(name_den + "-y-pt"));
    h2s[name_new + "-y-pt"]->GetZaxis()->SetTitle(title.c_str());
    h1s[name_new + "-y"] = (TH2D*)h1s.at(name_num + "-y")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y")->GetName(), name_num, name_new).c_str());
    h1s[name_new + "-y"]->Divide(h1s.at(name_den + "-y"));
    h1s[name_new + "-y"]->GetYaxis()->SetTitle(title.c_str());
    h1s[name_new + "-y__rebin"] = (TH2D*)h1s.at(name_num + "-y__rebin")->Clone(xjjc::str_replaceall(h1s.at(name_num + "-y__rebin")->GetName(), name_num, name_new).c_str());
    h1s[name_new + "-y__rebin"]->Divide(h1s.at(name_den + "-y__rebin"));
    h1s[name_new + "-y__rebin"]->GetYaxis()->SetTitle(title.c_str());
  };

  make_eff("eff", "eff_num", "eff_den", "#alpha #times #epsilon_{reco} #times #epsilon_{sel}");
  make_eff("effreco", "reco_num", "eff_den", "#alpha #times #epsilon_{reco}");
  make_eff("effsel", "eff_num", "reco_num", "#epsilon_{sel}");

  if (h3s.find("data_signalwin") != h3s.end()) {
    auto *heff = h2s.at("eff-y-pt"), *hdata = h2s.at("data_signalwin");

    auto make_effdata = [&h1s, &heff, &hdata](std::string name_refbin) {
      auto name_new = xjjc::str_replaceall(name_refbin, "eff", "effdata");
      h1s[name_new] = (TH1D*)h1s.at(name_refbin)->Clone(xjjc::str_replaceall(h1s.at(name_refbin)->GetName(), "eff", "effdata").c_str());
      h1s[name_new]->Reset();
      std::vector<double> num_y(h1s.at(name_new)->GetXaxis()->GetNbins(), 0), den_y(num_y.size(), 0);
      for (int i=0; i<heff->GetXaxis()->GetNbins(); i++) {
        for (int j=0; j<heff->GetYaxis()->GetNbins(); j++) {
          auto eff = heff->GetBinContent(i+1, j+1), // ! uncertainty
            ndata = hdata->GetBinContent(i+1, j+1);
          if (eff == 0) {
            __XJJLOG << "!! eff is 0 for the bin (" << i+1 <<", " << j+1 << ")" << std::endl;
            continue;
          }
          ybin = h1s[name_new]->FindBin(heff->GetBinCenter(i+1));
          num_y[ybin-1] += (1./eff)*ndata;
          den_y[ybin-1] += ndata;
        }
      }
      for (int i=0; i<h1s.at(name_new)->GetXaxis()->GetNbins(); i++) {
        if (num_y[i] == 0) {
          __XJJLOG << "?? ndata is 0 for the bin (" << i+1 << ")" << std::endl;
          continue;
        }
        h1s.at(name_new)->SetBinContent(i+1, 1./(num_y[i]/den_y[i]));
      }
    };

    make_effdata("eff-y");
    make_effdata("eff-y__rebin");
  }

  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1s, 0);
  for (auto& h : h2s)
    xjjroot::sethempty(h.second, 0, 0);
  for (auto& h : h1s) {
    int color = xjjc::str_contains(h.first, "data") ? xjjroot::mycolor_satmiddle["red"] : kBlack;
    xjjroot::setthgrstyle(h.second, color, 21, 1.2, color, 1, 1);
    xjjroot::sethempty(h.second, 0, 0);
    // xjjana::sethminmax(h.second, 0, 1.5);
  }

  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(inputname, { { "saveeff", "calceff" }, { ".root", ".pdf" } }));

  pdf->prepare();
  h2s["eff-y-pt"]->Draw("colz");
  pdf->write();

  pdf->prepare();
  h1s["eff-y__rebin"]->Draw("pe");
  h1s["effdata-y__rebin"]->Draw("pe same");
  pdf->write();

  pdf->prepare();
  h1s["effreco-y__rebin"]->Draw("pe");
  pdf->write();

  pdf->prepare();
  h1s["effsel-y__rebin"]->Draw("pe");
  pdf->write();

  pdf->close();

  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) {
    return macro(argv[1]);
  }
  return 1;
}

