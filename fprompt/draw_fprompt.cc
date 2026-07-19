#include "xjjanauti.h"
#include "xjjmypdf.h"

// #include "fpfitter.h"
#include "../include/draw.h"

int index_sf(const TH1D* h) {
  if (!h) { return -1; }
  auto str_name = xjjc::str_erasestar(h->GetName(), "*_sf-");
  str_name = xjjc::str_eraseall(str_name, "__y-*");
  return std::stoi(str_name.c_str());
};

struct Style {
  std::string title = "";
  Color_t color = kBlack;
  Style_t mstyle = 20;
  Style_t lstyle = 1;
};

const std::map<std::string, Style> m_style_data = {
  { "sub", Style{ .title = "Sideband sub", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "sigswap", Style{ .title = "sPlot", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
  { "fix", Style{ .title = "MC template", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "best", Style{ .title = "Best#scale[0.5]{ }#chi^{2}", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
  { "Dip3D", Style{ .title = "DCA", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "Dip3Dsig", Style{ .title = "DCA significance", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
};

Style style_data(const std::string& name, int exact = 0) {
  Style empty;
  for (const auto& [key, style] : m_style_data) {
    if (exact && name == key)
      return style;
    if (!exact && xjjc::str_contains(name, key))
      return style;
  }
  return empty;
}

struct Fprompt {
  std::map<std::string, TH1D*> h1s;
  int ibin_best = -1;
  int ibin_fix = -1;
};

Fprompt init(TDirectory* dir, std::string name) {
  Fprompt fp;
  fp.h1s["chi2"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_chi2_.+");
  fp.h1s["fprompt"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_fprompt_.+");
  fp.h1s["data"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_norm_.+_data-.+");
  fp.h1s["nonprompt-fitted"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_norm_.+_mc-nonprompt.+");
  fp.h1s["total-fitted"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_norm_.+_total.+");
  fp.h1s["pull"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_norm_.+_pull-.+");
  xjjroot::print_tab(fp.h1s, 0);
  fp.ibin_best = index_sf(fp.h1s["total-fitted"]) + 1;
  if (fp.ibin_best != fp.h1s["chi2"]->GetMinimumBin()) {
    __XJJLOG << "!! bad ibin_best, " << fp.ibin_best << " (index_sf) vs. " << fp.h1s["chi2"]->GetMinimumBin() << " (GetMinimumBin())" << std::endl;
    fp.ibin_best = -1;
  }
  fp.ibin_fix = fp.h1s["chi2"]->FindBin(1);
  if (fp.ibin_fix < 1 || fp.ibin_fix > fp.h1s["chi2"]->GetXaxis()->GetNbins())
    fp.ibin_fix = -1;
  for (auto& [key, h] : fp.h1s) {
    if (!h) {
      fp.ibin_best = -1;
      fp.ibin_fix = -1;
      continue;
    }
    h->SetName(Form("h1_%s_%s", key.c_str(), name.c_str()));
  }
  return fp;
}

int macro(const std::string& inputname, const std::string& outputname) {
  //
  auto* inf = TFile::Open(inputname.c_str());
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
  draw::bintex tbins(h3_bins, 0, 2);
  auto ny = h3_bins->GetXaxis()->GetNbins();
  auto* h1_bins_sf = xjjana::getobj<TH1D>(inf, "h1_bins_sf");
  auto nsf = h1_bins_sf->GetXaxis()->GetNbins();

  std::map<std::string, xjjc::info> infos;
  for (const std::string name : { "data", "prompt", "nonprompt", "fit" } ) {
    infos[name] = xjjana::get_info(inf, Form("%s/info", name.c_str()));
  }
  const auto var = infos.at("fit").at("var");
  const auto var_tex = style_data(var, 1).title;

  auto make_hfprompt = [&tbins](const char* name) {
    auto* h = tbins.make_h1_y<TH1D>(name);
    h->GetYaxis()->SetTitle("#it{f}_{prompt}");
    xjjroot::sethempty(h);
    xjjana::sethabsminmax(h, 0, 1.5);
    return h;
  };
  auto make_halpha = [&tbins](const char* name) {
    auto* h = tbins.make_h1_y<TH1D>(name);
    h->GetYaxis()->SetTitle("Resolution scale factor#scale[0.5]{ }#alpha_{reso}");
    xjjroot::sethempty(h);
    xjjana::sethabsminmax(h, 0.5, 2.2);
    return h;
  };

  auto ss_best = style_data("best"), ss_fix = style_data("fix");
  
  std::map<std::string, std::vector<Fprompt>> fpys;
  std::map<std::string, std::vector<TH1D*>> h1sfs_fprompt;
  std::map<std::string, TH1D*> h1s_fprompt_best, h1s_fprompt_fix, h1s_alpha_best;
  for (const std::string type_data : { "data-sub", "data-sigswap" }) {
    auto ss = style_data(type_data);
    for (int k=0; k<nsf; k++) {
      auto* h = make_hfprompt(Form("h1_fprompt_%s_sf-%d", type_data.c_str(), k));
      xjjroot::setthgrstyle(h, xjjroot::color_alpha(kBlack, 0.2), ss.mstyle, 0, xjjroot::color_alpha(kBlack, 0.2), 1, 1);
      h1sfs_fprompt[type_data].push_back(h);
    }
    h1s_fprompt_best[type_data] = make_hfprompt(Form("h1_fprompt_%s_best", type_data.c_str()));
    h1s_fprompt_fix[type_data] = make_hfprompt(Form("h1_fprompt_%s_fix", type_data.c_str()));
    xjjroot::setthgrstyle(h1s_fprompt_best[type_data], ss_best.color, ss.mstyle, 1.7, ss_best.color, 1, 1);
    xjjroot::setthgrstyle(h1s_fprompt_fix[type_data], ss_fix.color, ss.mstyle, 1.7, ss_fix.color, 1, 1);
    h1s_alpha_best[type_data] = make_halpha(Form("h1_alpha_%s_best", type_data.c_str()));
    h1s_alpha_best[type_data]->GetYaxis()->SetTitle("Best scale factor#scale[0.5]{ }#alpha_{reso}");
    xjjroot::setthgrstyle(h1s_alpha_best[type_data], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);

    // set bin content
    for (int i=0; i<ny; i++) {
      auto* dir = inf->GetDirectory(Form("dir__y-%d/dir_%s", i, type_data.c_str()));
      if (!dir) {
        __XJJLOG << "!! bad dir, abort." << std::endl;
        return 2;
      }
      auto fp = init(dir, Form("%s__y-%d", type_data.c_str(), i));
      if (fp.ibin_best < 0) 
        return 2;
      if (nsf != fp.h1s["chi2"]->GetXaxis()->GetNbins()) {
        __XJJLOG << "!! inconsistent binning, abort." << std::endl;
        return 2;
      }
      xjjroot::setthgrstyle(fp.h1s["chi2"], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);
      xjjroot::setthgrstyle(fp.h1s["fprompt"], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);
      for (int k=0; k<nsf; k++) {
        const auto fprompt_k = fp.h1s["fprompt"]->GetBinContent(k+1), fprompterr_k = fp.h1s["fprompt"]->GetBinError(k+1);
        h1sfs_fprompt[type_data][k]->SetBinContent(i+1, fprompt_k);
        h1sfs_fprompt[type_data][k]->SetBinError(i+1, 0.001);
        // h1sfs_fprompt[type_data][k]->SetBinError(i+1, fprompterr_k);
        if (k == fp.ibin_best-1) {
          h1s_fprompt_best[type_data]->SetBinContent(i+1, fprompt_k);
          h1s_fprompt_best[type_data]->SetBinError(i+1, fprompterr_k);
          h1s_alpha_best[type_data]->SetBinContent(i+1, h1_bins_sf->GetBinCenter(k));
          h1s_alpha_best[type_data]->SetBinError(i+1, 0);
        }
        if (k == fp.ibin_fix-1) {
          h1s_fprompt_fix[type_data]->SetBinContent(i+1, fprompt_k);
          h1s_fprompt_fix[type_data]->SetBinError(i+1, fprompterr_k);
        }
      }
      fpys[type_data].push_back(fp);
    }
  }

  // auto var_tex = xjjc::str_erasestar(h1s_fprompt_best["data-sub"]->GetXaxis()->GetTitle(), " [*");
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  pdf->prepare();
  h1s_fprompt_best["data-sub"]->Draw("axis");
  for (auto& [_, vhs] : h1sfs_fprompt) {
    for (auto& h : vhs)
      h->Draw("pe1 same");
  }
  for (auto& [_, h] : h1s_fprompt_fix)
    h->Draw("pe1 same");
  for (auto& [_, h] : h1s_fprompt_best)
    h->Draw("pe1 same");
  xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
  xjjroot::drawtexgroup(0.24, 0.86, {
      tbins.label_pt(-1),
      infos.at("data").at("cut_tex"),
      "#bf{Fit on " + var_tex + "}",
    }, 0.038, 13, 42, 1.2);
  auto* leg1 = new TLegend(0.45, 0.35-2*0.038*1.2, 0.8, 0.35);
  xjjroot::setleg(leg1, 0.038);
  leg1->SetNColumns(2);
  for (auto& [key, h] : h1s_fprompt_fix)
    leg1->AddEntry(h, style_data(key).title.c_str(), "p");
  for (auto& [key, h] : h1s_fprompt_best)
    leg1->AddEntry(h, style_data(key).title.c_str(), "p");
  leg1->Draw();
  xjjroot::drawtex(0.45-0.01, 0.35-0.005, ss_fix.title.c_str(), 0.038, 33, 62, ss_fix.color);
  xjjroot::drawtex(0.45-0.01, 0.35-0.005-0.038*1.2, ss_best.title.c_str(), 0.038, 33, 62, ss_best.color);
  pdf->getc()->RedrawAxis();
  pdf->write();

  auto* leg2 = new TLegend(0.24, 0.86-0.01-0.038*1.2*(3+2), 0.5, 0.86-0.01-0.038*1.2*3);
  xjjroot::setleg(leg2, 0.038);
  for (auto& [key, h] : h1s_alpha_best)
    leg2->AddEntry(h, style_data(key).title.c_str(), "pl");
  
  pdf->prepare();
  h1s_alpha_best["data-sub"]->Draw("axis");
  xjjroot::drawbox(h1s_alpha_best["data-sub"]->GetXaxis()->GetXmin(), h1_bins_sf->GetBinCenter(1),
                   h1s_alpha_best["data-sub"]->GetXaxis()->GetXmax(), h1_bins_sf->GetBinCenter(h1_bins_sf->GetNbinsX()), kBlack, 0.05);
  for (auto& [_, h] : h1s_alpha_best)
    h->Draw("pl same");
  xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
  xjjroot::drawtexgroup(0.24, 0.86, {
      tbins.label_pt(-1),
      infos.at("data").at("cut_tex"),
      "#bf{Fit on " + var_tex + "}",
    }, 0.038, 13, 42, 1.2);
  leg2->Draw();
  pdf->getc()->RedrawAxis();
  pdf->write();

  for (int i=0; i<ny; i++) {
    std::vector<TH1D*> hlist;
    for (auto& [_, fp] : fpys) {
      hlist.push_back(fp[i].h1s["chi2"]);
    }
    xjjana::sethsmin(hlist, 0);
    xjjana::sethsmax(hlist, 1.5);
    
    pdf->prepare();
    hlist.front()->Draw("axis");
    for (auto& h : hlist) h->Draw("p same");
    xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
    xjjroot::drawtexgroup(0.24, 0.86, {
        tbins.label_pt(-1),
        tbins.label_y(i),
        infos.at("data").at("cut_tex"),
        "#bf{Fit on " + var_tex + "}",
      }, 0.038, 13, 42, 1.2);
    xjjroot::moveleg_n_draw(leg2, -1, 0.86-0.01-0.038*1.2*4);
    leg2->Draw();
    pdf->write();
  }
  
  pdf->close();

  // auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");

  // outf->cd();
  // xjjroot::writehist(h3_bins);
  // xjjroot::closefile(outf);
  // inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
