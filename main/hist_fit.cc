#include <TH3D.h>

#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/dfitter.h"

int macro(std::string input_data, std::string input_template, int fit_type = 0) {
  std::cout<<std::endl;

#define READ_FILE(q)                                                    \
  auto* inf##q = TFile::Open(input##q.c_str());                         \
  auto info##q = xjjana::getval_regexp((TTree*)inf##q->Get("info"));    \
  __XJJLOG << "++ info" << std::endl;                                   \
  xjjc::print_tab(info##q, -1);

  READ_FILE(_data);
  READ_FILE(_template);

  std::map<std::string, TH3D*> h3s;
  std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, TH1D*> h1s;

  auto read_hists = [&h3s, &h1ys, &h1s](TFile* inf) {
    for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_.+")) { 
      auto name = xjjc::str_eraseall(h3->GetName(), "h3_"); 
      if (h3s.find(name) != h3s.end()) { 
        __XJJLOG << "!! error: name " << name << " already in the map." << std::endl; 
        xjjroot::print_tab(h3s, 0); 
        return 3; 
      } 
      auto ny = h3->GetXaxis()->GetNbins(); 
      h1ys[name].reserve(ny); 
      for (int i=0; i<ny; i++) { 
        auto* h1_mass = h3->ProjectionY(Form("h1_mass_%s__y-%d", name.c_str(), i), 
                                        i+1, i+1, 
                                        1, h3->GetZaxis()->GetNbins(), // pt: only in analysis range
                                        "e"); 
        h1ys[name].push_back(h1_mass); 
      } 
      if (h1s.find("dump-y") == h1s.end()) { 
        h1s["dump-y"] = h3->ProjectionX(Form("h1_y_dump"), 0, -1, 0, -1); // get binning of y
        h1s["dump-y"]->Reset();
      } 
      h3s[name] = h3; 
    }
    return 0;
  };

  if (read_hists(inf_data)) { return 3; }
  if (read_hists(inf_template)) { return 3; }

  // prepare TH1
  auto make_h1_y = [&h1s](const std::string& name, const std::string& ytitle) {
    // h1_y_yield
    h1s[name + "-y"] = (TH1D*)h1s.at("dump-y")->Clone(Form("h1_y_%s", name.c_str()));
    h1s[name + "-y"]->Reset();
    h1s[name + "-y"]->GetYaxis()->SetTitle(ytitle.c_str());
  };
  make_h1_y("yield", "Raw Yield");
  make_h1_y("width68mc", "Signal Effective#scale[0.5]{ }#sigma in MC [GeV]");
  make_h1_y("width95mc", "Signal Effective 2#sigma in MC [GeV]");
    
  xjjroot::print_tab(h3s, 0);
  xjjroot::print_tab(h1ys, 0);
  xjjroot::print_tab(h1s, 0);

  auto* h3_dump = static_cast<TH3D*>(h3s.at("data")->Clone("h3_dump")); h3_dump->Reset();
  auto label_y = [&h3_dump](int i) {
    const auto ymin = (i>=0 ? h3_dump->GetXaxis()->GetBinLowEdge(i+1) : h3_dump->GetXaxis()->GetBinLowEdge(1)),
      ymax = (i>=0 ? h3_dump->GetXaxis()->GetBinUpEdge(i+1) : h3_dump->GetXaxis()->GetBinUpEdge(h3_dump->GetXaxis()->GetNbins()));
    return xjjc::number_range_string(ymin, ymax, "y", -1.e1);
  };
  auto label_pt = [&h3_dump](int i = -1) {
    const auto ptmin = (i>=0 ? h3_dump->GetZaxis()->GetBinLowEdge(i+1) : h3_dump->GetZaxis()->GetBinLowEdge(1)),
      ptmax = (i>=0 ? h3_dump->GetZaxis()->GetBinUpEdge(i+1) : h3_dump->GetZaxis()->GetBinUpEdge(h3_dump->GetZaxis()->GetNbins()));
    return xjjc::number_range_string(ptmin, ptmax, "p_{T}") + " GeV";
  };

  auto* df = new xjjroot::dfitter("S3");
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf(xjjc::str_replaceall(input_data, { { "rootfiles/", "figspdf/" }, { ".root", ".pdf" }, { "savehist_", "fithist_" } }));
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });
  for (int i=0; i<h1ys["data"].size(); i++) {
    const auto &h = h1ys.at("data").at(i),
      &hmc = h1ys.at("match").at(i), &hmcswap = h1ys.at("swap").at(i);
    pdf->prepare();
    df->fit(h, hmc, hmcswap);
    xjjroot::drawtexgroup(0.25, 0.86, { label_y(i), label_pt(), "#bf{" + info_data.at("cut_tex") + "}" }, 0.035, 13);
    df->draw_leg();
    df->draw_result(0.25, 0.86-3*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::internal, info_data.at("input_tex"));
    pdf->write(Form("%s_y-%d.pdf", name_png.c_str(), i));

    h1s["yield-y"]->SetBinContent(i+1, df->yield());
    h1s["yield-y"]->SetBinError(i+1, df->yieldErr());
    auto w68mc = df->width_mc_mass(xjjana::frac_1sigma),
      w95mc = df->width_mc_mass(xjjana::frac_2sigma);
    h1s["width68mc-y"]->SetBinContent(i+1, w68mc.first);
    h1s["width68mc-y"]->SetBinError(i+1, w68mc.second);
    h1s["width95mc-y"]->SetBinContent(i+1, w95mc.first);
    h1s["width95mc-y"]->SetBinError(i+1, w95mc.second);
    
    pdf->prepare();
    df->set_hist(hmc);
    hmc->Draw("pe1"); 
    df->set_hist(hmcswap);
    hmcswap->Draw("pe1 same");
    df->draw_fmc();
    xjjroot::drawtexgroup(0.25, 0.86, { label_y(i), label_pt() }, 0.035, 13);
    df->draw_params(0.25, 0.86-2*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::simulation, info_template.at("input_tex"));
    pdf->write(Form("%s_mc_y-%d.pdf", name_png.c_str(), i));
  }

  xjjroot::print_th(h1s.at("yield-y"));

  for (const std::string& p : { "yield-y", "width68mc-y", "width95mc-y" }) {
    xjjroot::sethempty(h1s.at(p), 0, 0.1);
    xjjroot::setthgrstyle(h1s.at(p), kBlack, 21, 1.5, kBlack, 1, 1);
    xjjana::sethminmax(h1s.at(p), 0, 1.5);
    pdf->prepare();
    h1s.at(p)->Draw("pe1");
    xjjroot::drawtexgroup(0.25, 0.86, { label_pt(), info_data.at("cut_tex") }, 0.035, 13, 42, 1.25);
    xjjroot::drawCMS(xjjc::str_contains(p, "mc") ? xjjroot::CMS::simulation : xjjroot::CMS::internal, info_data.at("input_tex"));
    pdf->write();
  }
  
  pdf->draw_cover( {
      "#bf{Data} " + info_data.at("input"),
      "#bf{Cut} " + info_data.at("cut"),
      "#bf{Template} " + info_template.at("input"),
      "#bf{Cut} " + info_template.at("cut"),
    }, 0.03);

  pdf->close();
  
  h1s.erase("dump-y");
  auto* outf = xjjroot::newfile(xjjc::str_replaceall(input_data, { { "savehist_", "fithist_" } }));

  for (const auto& [_, h] : h1s) xjjroot::writehist(h);
  for (const auto& [_, h] : h3s) xjjroot::writehist(h);
  for (const auto& [_, hh] : h1ys)
    for (const auto& h : hh) xjjroot::writehist(h);
  
  outf->mkdir("data")->cd();
  auto* t_data = new TTree("info", "");
  for (auto& [key, content] : info_data) {
    t_data->Branch(key.c_str(), &content);
  }
  t_data->Fill();
  t_data->Write();
  outf->cd();

  outf->mkdir("template")->cd();
  auto* t_template = new TTree("info", "");
  for (auto& [key, content] : info_template) {
    t_template->Branch(key.c_str(), &content);
  }
  t_template->Fill();
  t_template->Write();
  outf->cd();

  outf->Close();
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}

