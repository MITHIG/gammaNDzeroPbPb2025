#include <TH3D.h>

#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/dfitter.h"
#include "../include/draw.h"
#include "../include/util.h"

int macro(const std::string& input_data, const std::string& input_template, const std::string& outputname, const std::string& fit_opt = "3", int save_png = 1) {
  std::cout<<std::endl;

  const auto fitopt = util::parse_input(fit_opt);
  __XJJLOG << ">> fitting option : " << fit_opt << " -> " << fitopt.content << std::endl;
  xjjc::info info_fit = {
    { "fitopt", fitopt.content },
    { "fitopt_tex", fitopt.tex },
    { "fitopt_tag", fitopt.tag },
  };

#define READ_FILE(q)                                                    \
  auto* inf##q = TFile::Open(input##q.c_str());                         \
  auto info##q = xjjana::getval_regexp((TTree*)inf##q->Get("info"));    \
  __XJJLOG << "++ info" << std::endl;                                   \
  xjjc::print_tab(info##q, -1);

  READ_FILE(_data);
  READ_FILE(_template);

  std::map<std::string, TH3D*> h3s;
  std::map<std::string, xjjc::array2D<TH1D*>> h1ptys;
  std::map<std::string, std::vector<TH1D*>> h1pts;

  auto read_hists = [&h3s, &h1ptys](TFile* inf) {
    for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_.+")) { 
      auto name = xjjc::str_eraseall(h3->GetName(), "h3_"); 
      if (h3s.find(name) != h3s.end()) { 
        __XJJLOG << "!! error: name " << name << " already in the map." << std::endl; 
        xjjroot::print_tab(h3s, 0); 
        return 3; 
      } 
      auto ny = h3->GetXaxis()->GetNbins(), npt = h3->GetZaxis()->GetNbins();
      h1ptys[name] = xjjc::array2d<TH1D*>(npt, ny);
      for (int i=0; i<ny; i++) { 
        for (int j=0; j<npt; j++) { 
          auto* h1_mass = h3->ProjectionY(Form("h1_mass_%s__pt-%d__y-%d", name.c_str(), j, i), 
                                          i+1, i+1, // y bin
                                          j+1, j+1, // pt bin
                                          "e"); 
          h1ptys.at(name)[j][i] = h1_mass; 
        }
      }
      h3s[name] = h3; 
    }
    return 0;
  };

  if (read_hists(inf_data)) { return 3; }
  if (read_hists(inf_template)) { return 3; }

  auto* h3_bins = static_cast<TH3D*>(h3s.at("data")->Clone("h3_bins")); h3_bins->Reset();
  const draw::bintex tbins(h3_bins, 0, 2);
  
  // prepare TH1
  auto make_h1_y = [&tbins, &h1pts](const std::string& name, const std::string& ytitle) {
    // h1_y_yield
    for (int j=0; j<tbins.npt(); j++) {
      auto* h = tbins.make_h1_y<TH1D>(Form("h1_y_%s__pt-%d", name.c_str(), j));
      h->GetYaxis()->SetTitle(ytitle.c_str());
      h1pts[name + "-y"].push_back(h);
    }
  };
  make_h1_y("yield", "Raw Yield");
  make_h1_y("width68mc", "Signal Effective#scale[0.5]{ }#sigma in MC [GeV]");
  make_h1_y("width95mc", "Signal Effective 2#sigma in MC [GeV]");
    
  xjjroot::print_tab(h3s, 0);
  xjjroot::print_tab(h1pts, 0);

  xjjroot::setgstyle(1);
  auto dfs = xjjc::array2d<xjjroot::dfitter*>(tbins.npt(), tbins.ny());
  // auto ndiv = std::ceil(std::sqrt(tbins.npt()));
  // gStyle->SetLineScalePS(3./ndiv);
  // auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf", "c", ndiv, ndiv);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });
  for (int i=0; i<tbins.ny(); i++) {
    // pdf->prepare();
    // xjjroot::divide(ndiv, ndiv);
    for (int j=0; j<tbins.npt(); j++) {
      pdf->prepare();
      // auto* p = (TPad*)pdf->getc()->cd(j+1);
      const auto &h = h1ptys.at("data")[j][i],
        &hmc = h1ptys.at("mc-match")[j][i], &hmcswap = h1ptys.at("mc-swap")[j][i],
        &hmckk = h1ptys.at("mc-kk")[j][i], &hmcpipi = h1ptys.at("mc-pipi")[j][i];

      auto* df = new xjjroot::dfitter(fitopt.content.c_str());
      df->fit(h, hmc, hmcswap, hmckk, hmcpipi);
      xjjroot::drawtexgroup(0.24, 0.86, { tbins.label_y(i), tbins.label_pt(j), "#bf{" + info_data.at("cut_tex") + "}" }, 0.035, 13);
      df->draw_leg();
      df->draw_result(0.24, 0.86-3*0.035*1.15, 0.035);
      xjjroot::drawCMS(xjjroot::CMS::internal, info_data.at("input_tex"));

      h1pts.at("yield-y")[j]->SetBinContent(i+1, df->yield());
      h1pts.at("yield-y")[j]->SetBinError(i+1, df->yieldErr());
      auto w68mc = df->width_mc_match(xjjana::frac_1sigma),
        w95mc = df->width_mc_match(xjjana::frac_2sigma);
      h1pts.at("width68mc-y")[j]->SetBinContent(i+1, w68mc.first);
      h1pts.at("width68mc-y")[j]->SetBinError(i+1, w68mc.second);
      h1pts.at("width95mc-y")[j]->SetBinContent(i+1, w95mc.first);
      h1pts.at("width95mc-y")[j]->SetBinError(i+1, w95mc.second);
    
      dfs[j][i] = df;
      gPad->RedrawAxis();
      pdf->write(Form("%s_pt-%d_y-%d.pdf", name_png.c_str(), j, i), save_png ? "" : "X");
    }
  }

  for (int i=0; i<tbins.ny(); i++) {
    // pdf->prepare();
    // xjjroot::divide(ndiv, ndiv);
    for (int j=0; j<tbins.npt(); j++) {
      pdf->prepare();
      // auto* p = (TPad*)pdf->getc()->cd(j+1);
      auto* df = dfs[j][i];

      for (const std::string s : { "match", "swap", "kk", "pipi" }) {
        if (!df->haskkpipi() && (s == "kk" || s == "pipi")) continue; 
        auto* hmc = h1ptys.at("mc-" + s)[j][i];
        df->set_hist(hmc);
        if (xjjroot::fstyle.find(s) != xjjroot::fstyle.end() && s != "match")
          xjjroot::setthgrstyle(hmc, xjjroot::fstyle.at(s).lcolor, -1, -1, xjjroot::fstyle.at(s).lcolor, -1, -1, -1, -1, -1, 0.6, 0.6);
        hmc->Draw(s == "match" ? "pe1" : "pe1 same");
      }
      df->draw_fmc();
      df->draw_legmc();
      xjjroot::drawtexgroup(0.24, 0.86, { tbins.label_y(i), tbins.label_pt(j) }, 0.035, 13);
      df->draw_params(0.24, 0.86-2*0.035*1.15, 0.035);
      // xjjroot::drawtexgroup(0.88, 0.86, { info_fit.at("fitopt_tex") }, 0.035, 33);
      xjjroot::drawCMS(xjjroot::CMS::simulation, info_template.at("input_tex"));
      gPad->RedrawAxis();
      pdf->write(Form("%s_mc_pt-%d_y-%d.pdf", name_png.c_str(), j, i), save_png ? "" : "X");
    }
  }

  for (const std::string& p : { "yield-y", "width68mc-y", "width95mc-y" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts.at(p), 0);
    xjjana::sethsmax(h1pts.at(p), 1.6);
    auto* leg = new TLegend(0.60, 0.86-0.035*1.25*tbins.npt(), 0.80, 0.86);
    xjjroot::setleg(leg, 0.035);
    auto cc = xjjroot::grayscales_alpha(tbins.npt());
    for (int j=0; j<tbins.npt(); j++) {
      xjjroot::sethempty(h1pts.at(p)[j], 0, 0.5);
      xjjroot::setthgrstyle(h1pts.at(p)[j], kBlack, 21, 1.5, kBlack, 1, 1, -1, -1, -1, cc[j], cc[j]);
      leg->AddEntry(h1pts.at(p)[j], tbins.label_pt(j).c_str(), "p");
    }
    h1pts.at(p).front()->Draw("axis");
    for (auto& h : h1pts.at(p)) h->Draw("pe1 same");
    leg->Draw();
    xjjroot::drawtexgroup(0.24, 0.85, {
        info_data.at("cut_tex"),
        info_fit.at("fitopt_tex"),
      }, 0.035, 13, 42, 1.25);
    xjjroot::drawCMS(xjjc::str_contains(p, "mc") ? xjjroot::CMS::simulation : xjjroot::CMS::internal,
                     xjjc::str_contains(p, "mc") ? info_template.at("input_tex") : info_data.at("input_tex"));
    auto save_indi = save_png && (p == "width68mc-y");
    pdf->write(Form("%s_%s.pdf", name_png.c_str(), p.c_str()), save_indi ? "" : "X");
  }
  
  pdf->draw_cover( {
      "#bf{Data} " + info_data.at("input"),
      "#bf{Cut} " + info_data.at("cut"),
      "#bf{Template} " + info_template.at("input"),
      "#bf{Cut} " + info_template.at("cut"),
    }, 0.03);

  pdf->close();
  
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");

  for (const auto& [_, h] : h3s) xjjroot::writehist(h);
  for (const auto& [_, hh] : h1pts)
    for (const auto& h : hh)
      xjjroot::writehist(h);
  for (const auto& [_, hhh] : h1ptys)
    for (const auto& hh : hhh)
      for (const auto& h : hh)
        xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);

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

  outf->mkdir("fit")->cd();
  auto* t_fit = new TTree("info", "");
  for (auto& [key, content] : info_fit) {
    t_fit->Branch(key.c_str(), &content);
  }
  t_fit->Fill();
  t_fit->Write();
  outf->cd();

  xjjroot::closefile(outf);
  
  return 0; 
}

int main(int argc, char* argv[]) {
  if (argc == 6) {
    return macro(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]));
  }
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], argv[4]);
  }
  if (argc == 4) {
    return macro(argv[1], argv[2], argv[3]);
  }
  return 1;
}
