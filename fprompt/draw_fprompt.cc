#include "TFitResult.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "fpfitter.h"
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
  { "Dip3Dsig", Style{ .title = "DCA /#scale[0.5]{ }#sigma(DCA)", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
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
  std::map<std::string, TGraphAsymmErrors*> grs;
  std::map<std::string, fpfitter*> frs = { { "fix", nullptr }, { "best", nullptr } };
  std::map<std::string, int> ibin_wp = { { "fix", -1 }, { "best", -1 } };
};

Fprompt init(TDirectory* dir, std::string name) {
  Fprompt fp;
  fp.h1s["chi2"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_chi2_.+");
  fp.h1s["fprompt"] = xjjana::getobj_regexp_first<TH1D>(dir, "h1_fprompt_.+");
  fp.grs["fprompt"] = xjjana::getobj_regexp_first<TGraphAsymmErrors>(dir, "gr_fprompt_.+");
  for (const std::string wp : { "best", "fix" }) {
    auto* dir_fit = dir->GetDirectory(Form("fit_%s", wp.c_str()));
    fp.h1s[wp + "-data"] = xjjana::getobj_regexp_first<TH1D>(dir_fit, "h1_norm_.+_data-.+");
    fp.h1s[wp + "-nonprompt-fitted"] = xjjana::getobj_regexp_first<TH1D>(dir_fit, "h1_norm_.+_mc-nonprompt.+");
    fp.h1s[wp + "-total-fitted"] = xjjana::getobj_regexp_first<TH1D>(dir_fit, "h1_norm_.+_total.+");
    fp.h1s[wp + "-pull"] = xjjana::getobj_regexp_first<TH1D>(dir_fit, "h1_norm_.+_pull-.+");
    fp.h1s[wp + "-ratio"] = xjjana::getobj_regexp_first<TH1D>(dir_fit, "h1_norm_.+_ratio-.+");
    auto* fitresult = xjjana::getobj_regexp_first<TFitResult>(dir_fit, "TFitResult-.+");
    fp.frs[wp] = new fpfitter(fp.h1s.at(wp + "-data"), fp.h1s.at(wp + "-total-fitted"), fp.h1s.at(wp + "-nonprompt-fitted"), fitresult, fp.h1s.at(wp + "-pull"), fp.h1s.at(wp + "-ratio"));
  }
  fp.ibin_wp.at("best") = index_sf(fp.h1s.at("best-total-fitted")) + 1;
  if (fp.ibin_wp.at("best") != fp.h1s.at("chi2")->GetMinimumBin()) {
    __XJJLOG << "!! bad ibin_best, " << fp.ibin_wp.at("best") << " (index_sf) vs. " << fp.h1s.at("chi2")->GetMinimumBin() << " (GetMinimumBin())" << std::endl;
    fp.ibin_wp.at("best") = -1;
  } 
  fp.ibin_wp.at("fix") = index_sf(fp.h1s.at("fix-total-fitted")) + 1;
  if (fp.ibin_wp.at("fix") != fp.h1s.at("chi2")->FindBin(1)) {
    __XJJLOG << "!! bad ibin_fix, " << fp.ibin_wp.at("fix") << " (index_sf) vs. " << fp.h1s.at("chi2")->FindBin(1) << " (FindBin(1))" << std::endl;
    fp.ibin_wp.at("fix") = -1;
  }
  for (auto& [key, h] : fp.h1s) {
    if (!h) {
      fp.ibin_wp.at("best") = -1;
      fp.ibin_wp.at("fix") = -1;
      continue;
    }
    h->SetName(Form("h1_%s_%s", key.c_str(), name.c_str()));
  }
  xjjroot::print_tab(fp.h1s, 0);
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
    infos[name] = xjjana::get_info(inf, Form("info/%s", name.c_str()));
    xjjc::print_tab(infos.at(name), -1);
  }
  const auto var = infos.at("fit").at("var");
  const auto var_tex = style_data(var, 1).title;

  auto make_hvsy = [&tbins](const char* name, const char* ytitle,
                            double hmin, double hmax) {
    auto* h = tbins.make_h1_y<TH1D>(name);
    h->GetYaxis()->SetTitle(ytitle);
    xjjroot::sethempty(h);
    xjjana::sethabsminmax(h, hmin, hmax);
    return h;
  };

  auto ss_best = style_data("best"), ss_fix = style_data("fix");
  
  std::map<std::string, std::vector<Fprompt>> fpys;
  std::map<std::string, std::vector<TH1D*>> h1sfs_fprompt;
  std::map<std::string, TGraphAsymmErrors*> grs_fprompt_best, grs_fprompt_fix;
  std::map<std::string, TH1D*> h1s_alpha_best, h1s_chi2_best;
  for (const std::string type_data : { "data-sub", "data-sigswap" }) {
    auto ss = style_data(type_data);
    for (int k=0; k<nsf; k++) {
      auto* h = make_hvsy(Form("h1_fprompt_%s_sf-%d", type_data.c_str(), k), "#it{f}_{prompt}", 0, 1.5);
      xjjroot::setthgrstyle(h, xjjroot::color_alpha(kBlack, 0.2), ss.mstyle, 0, xjjroot::color_alpha(kBlack, 0.2), 1, 1);
      h1sfs_fprompt[type_data].push_back(h);
    }
    grs_fprompt_best[type_data] = new TGraphAsymmErrors(ny);
    grs_fprompt_best[type_data]->SetName(Form("gr_fprompt_%s_best", type_data.c_str()));
    grs_fprompt_fix[type_data] = new TGraphAsymmErrors(ny);
    grs_fprompt_fix[type_data]->SetName(Form("gr_fprompt_%s_fix", type_data.c_str()));
    xjjroot::setthgrstyle(grs_fprompt_best[type_data], ss_best.color, ss.mstyle, 1.7, ss_best.color, 1, 1);
    xjjroot::setthgrstyle(grs_fprompt_fix[type_data], ss_fix.color, ss.mstyle, 1.7, ss_fix.color, 1, 1);
    h1s_alpha_best[type_data] = make_hvsy(Form("h1_alpha_%s_best", type_data.c_str()), "Best scale factor#scale[0.5]{ }#alpha_{reso}", 0.5, 2.2);
    xjjroot::setthgrstyle(h1s_alpha_best[type_data], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);
    h1s_chi2_best[type_data] = make_hvsy(Form("h1_chi2_%s_best", type_data.c_str()), "Best#scale[0.5]{ }#chi^{2} / ndf", 0., 3.5);
    xjjroot::setthgrstyle(h1s_chi2_best[type_data], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);

    // set bin content
    for (int i=0; i<ny; i++) {
      auto* dir = inf->GetDirectory(Form("dir__y-%d/dir_%s", i, type_data.c_str()));
      if (!dir) {
        __XJJLOG << "!! bad dir, abort." << std::endl;
        return 2;
      }
      auto fp = init(dir, Form("%s__y-%d", type_data.c_str(), i));
      if (fp.ibin_wp.at("best") < 0) 
        return 2;
      if (nsf != fp.h1s.at("chi2")->GetXaxis()->GetNbins()) {
        __XJJLOG << "!! inconsistent binning, abort." << std::endl;
        return 2;
      }
      xjjroot::setthgrstyle(fp.h1s.at("chi2"), ss.color, ss.mstyle, 1.7, ss.color, 1, 1);
      xjjroot::setthgrstyle(fp.grs["fprompt"], ss.color, ss.mstyle, 1.7, ss.color, 1, 1);
      xjjroot::setthgrstyle(fp.h1s.at("fprompt"), ss.color, ss.mstyle, 1.7, ss.color, 1, 1, -1, -1, -1, 0.2, 0.2);
      for (int k=0; k<nsf; k++) {
        const auto fprompt_k = fp.h1s.at("fprompt")->GetBinContent(k+1), fprompterr_k = fp.h1s.at("fprompt")->GetBinError(k+1);
        h1sfs_fprompt[type_data][k]->SetBinContent(i+1, fprompt_k);
        h1sfs_fprompt[type_data][k]->SetBinError(i+1, 0.001);
        // h1sfs_fprompt[type_data][k]->SetBinError(i+1, fprompterr_k);
        const auto x_k = h1sfs_fprompt[type_data][k]->GetBinCenter(i+1), ex_k = h1sfs_fprompt[type_data][k]->GetBinWidth(i+1)/2.;
        const auto eyl_k = fp.grs["fprompt"]->GetErrorYlow(k), eyh_k = fp.grs["fprompt"]->GetErrorYhigh(k);
        if (k == fp.ibin_wp.at("best")-1) {
          grs_fprompt_best[type_data]->SetPoint(i, x_k, fprompt_k);
          grs_fprompt_best[type_data]->SetPointError(i, ex_k, ex_k, eyl_k, eyh_k);
          h1s_alpha_best[type_data]->SetBinContent(i+1, h1_bins_sf->GetBinCenter(k+1));
          h1s_alpha_best[type_data]->SetBinError(i+1, 0);
          h1s_chi2_best[type_data]->SetBinContent(i+1, fp.h1s.at("chi2")->GetBinContent(k+1));
          h1s_chi2_best[type_data]->SetBinError(i+1, 0);
        }
        if (k == fp.ibin_wp.at("fix")-1) { 
          grs_fprompt_fix[type_data]->SetPoint(i, x_k, fprompt_k);
          grs_fprompt_fix[type_data]->SetPointError(i, ex_k, ex_k, eyl_k, eyh_k);
        }
      }
      fpys[type_data].push_back(fp);
    }
  }

  // auto var_tex = xjjc::str_erasestar(h1s_fprompt_best["data-sub"]->GetXaxis()->GetTitle(), " [*");
  auto* h1_dump_fprompt = (TH1D*)h1sfs_fprompt["data-sub"].front()->Clone("h1_dump_fprompt");

  xjjroot::setgstyle(1);
  gStyle->SetLineScalePS(1.5);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf", "c", xjjroot::mypdf::w_default*2, xjjroot::mypdf::h_default);
  // pdf->prepare();
  // h1_dump_fprompt->Draw("axis");
  // for (auto& [_, vhs] : h1sfs_fprompt) {
  //   for (auto& h : vhs)
  //     h->Draw("pe1 same");
  // }
  // for (auto& [_, h] : grs_fprompt_fix)
  //   h->Draw("pe1 same");
  // for (auto& [_, h] : grs_fprompt_best)
  //   h->Draw("pe1 same");
  // xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
  // xjjroot::drawtexgroup(0.24, 0.86, {
  //     tbins.label_pt(-1),
  //     infos.at("data").at("cut_tex"),
  //     "Fit on#scale[0.5]{ }#bf{" + var_tex + "}",
  //   }, 0.038, 13, 42, 1.2);
  // xjjroot::drawtexgroup(0.91, 0.86, {
  //     "Uncert. by toy smearing",
  //   }, 0.038, 33, 42, 1.2);
  // auto* leg1 = new TLegend(0.45, 0.35-2*0.038*1.2, 0.8, 0.35);
  // xjjroot::setleg(leg1, 0.038);
  // leg1->SetNColumns(2);
  // for (auto& [key, h] : grs_fprompt_fix)
  //   leg1->AddEntry(h, style_data(key).title.c_str(), "p");
  // for (auto& [key, h] : grs_fprompt_best)
  //   leg1->AddEntry(h, style_data(key).title.c_str(), "p");
  // leg1->Draw();
  // xjjroot::drawtex(0.45-0.01, 0.35-0.005, ss_fix.title.c_str(), 0.038, 33, 62, ss_fix.color);
  // xjjroot::drawtex(0.45-0.01, 0.35-0.005-0.038*1.2, ss_best.title.c_str(), 0.038, 33, 62, ss_best.color);
  // pdf->getc()->RedrawAxis();
  // pdf->write();

  auto* leg2 = new TLegend(0.24, 0.86-0.01-0.038*1.2*(3+2), 0.5, 0.86-0.01-0.038*1.2*3);
  xjjroot::setleg(leg2, 0.038);
  for (auto& [key, h] : h1s_alpha_best)
    leg2->AddEntry(h, style_data(key).title.c_str(), "p");

  auto draw_global = [&infos, &tbins, &var_tex](int ibin_y = -1) {
    xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
    std::vector<std::string> tlist = {
      tbins.label_pt(-1),
      infos.at("data").at("cut_tex"),
      "Fit on#scale[0.5]{ }#bf{" + var_tex + "}",
    };
    if (ibin_y >= 0) tlist.insert(tlist.begin() + 1, tbins.label_y(ibin_y));
    xjjroot::drawtexgroup(0.24, 0.86, tlist, 0.038, 13, 42, 1.2);
  };
  
  pdf->prepare();
  pdf->getc()->Divide(2, 1);
  pdf->getc()->cd(1);
  std::cout<<gPad->GetWNDC()<<", "<<gPad->GetHNDC()<<std::endl;
  h1s_alpha_best["data-sub"]->Draw("axis");
  xjjroot::drawbox(h1s_alpha_best["data-sub"]->GetXaxis()->GetXmin(), h1_bins_sf->GetBinCenter(1),
                   h1s_alpha_best["data-sub"]->GetXaxis()->GetXmax(), h1_bins_sf->GetBinCenter(h1_bins_sf->GetNbinsX()), kBlack, 0.05);
  for (auto& [_, h] : h1s_alpha_best)
    h->Draw("pl same");
  leg2->Draw();
  draw_global();
  pdf->getc()->RedrawAxis();

  pdf->getc()->cd(2);
  h1s_chi2_best["data-sub"]->Draw("axis");
  xjjroot::drawbox(h1s_chi2_best["data-sub"]->GetXaxis()->GetXmin(), 0.8,
                   h1s_chi2_best["data-sub"]->GetXaxis()->GetXmax(), 1.2, kBlack, 0.05);
  for (auto& [_, h] : h1s_chi2_best)
    h->Draw("pl same");
  draw_global();
  leg2->Draw();
  pdf->getc()->RedrawAxis();
  pdf->write();

  for (int i=0; i<ny; i++) {
    pdf->prepare();
    pdf->getc()->Divide(2, 1);
    std::vector<TH1D*> hlist;

    pdf->getc()->cd(1);
    hlist.clear();
    for (auto& [_, fp] : fpys) {
      hlist.push_back(fp[i].h1s.at("chi2"));
    }
    xjjana::sethsmin(hlist, 0);
    xjjana::sethsmax(hlist, 1.5);
    
    hlist.front()->Draw("axis");
    for (auto& [_, fp] : fpys) {
      xjjroot::drawline(fp[i].h1s.at("chi2")->GetBinCenter(fp[i].ibin_wp.at("best")), 0,
                        fp[i].h1s.at("chi2")->GetBinCenter(fp[i].ibin_wp.at("best")), fp[i].h1s.at("chi2")->GetBinContent(fp[i].ibin_wp.at("best")),
                        fp[i].h1s.at("chi2")->GetLineColor(), 2, 2);
    }
    for (auto& h : hlist) h->Draw("p same");
    draw_global(i);
    xjjroot::moveleg_n_draw(leg2, -1, 0.86-0.01-0.038*1.2*4);

    pdf->getc()->cd(2);
    hlist.clear();
    for (auto& [_, fp] : fpys) {
      hlist.push_back(fp[i].h1s.at("fprompt"));
    }
    xjjana::sethsmin(hlist, 0.);
    xjjana::sethsmax(hlist, 1.5);
    
    hlist.front()->Draw("axis");
    for (auto& [_, fp] : fpys) {
      xjjroot::drawline(fp[i].h1s.at("fprompt")->GetBinCenter(fp[i].ibin_wp.at("best")), 0,
                        fp[i].h1s.at("fprompt")->GetBinCenter(fp[i].ibin_wp.at("best")), fp[i].h1s.at("fprompt")->GetBinContent(fp[i].ibin_wp.at("best")),
                        fp[i].grs["fprompt"]->GetLineColor(), 2, 2);
    }
    
    for (auto& h : hlist) h->Draw("pe1 same");
    for (auto& [_, fp] : fpys) fp[i].grs["fprompt"]->Draw("pe1 same");
    draw_global(i);
    xjjroot::drawtexgroup(0.91, 0.86, {
        "Uncert. by toy smearing",
        "GetParError() as light lines"
      }, 0.038, 33, 42, 1.2, 1, 0.2, { kBlack, kGray });
    xjjroot::moveleg_n_draw(leg2, 0.6, 0.5);

    pdf->write();

    for (auto& [type_data, fp] : fpys) {
      for (auto q : { Qual::pull, Qual::ratio }) {
        pdf->prepare();
        pdf->getc()->Divide(2, 1);
        int kp = 1;
        for (const std::string wp : { "fix", "best" }) {
          auto* fitter = fp[i].frs.at(wp);
          auto pads = fitter->draw(static_cast<TPad*>(pdf->getc()->cd(kp)), q);
          pads[0]->cd();
          xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)", 1./fpfitter::pratio);
          xjjroot::drawtexgroup(fpfitter::xleft - 0.01, fpfitter::ytop - 0.005, {
              style_data(type_data).title.c_str(),
              tbins.label_pt(-1), tbins.label_y(i),
              "#alpha_{reso} = " + std::string(Form("%.2f", fp[i].h1s.at("chi2")->GetBinCenter(fp[i].ibin_wp.at(wp)))),
            }, fpfitter::tsize/fpfitter::pratio, 33);
          kp++;
        }
        pdf->write();
      }
    }
  }
  
  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h3_bins);
  xjjroot::writehist(h1_bins_sf);
  for (auto& [type, fps] : fpys) {
    outf->cd();
    auto* dir = outf->mkdir(Form("dir_%s", type.c_str()));
    dir->cd();
    xjjroot::writehist(grs_fprompt_best.at(type));
    xjjroot::writehist(grs_fprompt_fix.at(type));
    xjjroot::writehist(h1s_alpha_best.at(type));
    xjjroot::writehist(h1s_chi2_best.at(type));

    auto *dir_best = dir->mkdir("fit_best"), *dir_fix = dir->mkdir("fit_fix");
    for (auto& fp : fps) {
      xjjroot::writehist(fp.h1s.at("fprompt"));
      xjjroot::writehist(fp.h1s.at("chi2"));
      // dir_best->cd();
      // xjjroot::writehist(fp.frs.at("best"));
      // dir_fix->cd();
      // xjjroot::writehist(fp.frs.at("fix"));
      // for (auto& [key, h] : fp.h1s) {
      //   if (xjjc::str_contains(key, "best-")) {
      //     dir_best->cd();
      //     xjjroot::writehist(h);
      //   } else if (xjjc::str_contains(key, "fix-")) {
      //     dir_fix->cd();
      //     xjjroot::writehist(h);
      //   }
      // }
      // dir->cd();
    }
    dir->mkdir("dir_sfs")->cd();
    for (auto& h : h1sfs_fprompt.at(type))
      xjjroot::writehist(h);
  }
  outf->cd();

  outf->mkdir("info")->cd();
  for (auto& [iname, info] : infos) {
    auto* t_data = new TTree(iname.c_str(), "");
    for (auto& [key, content] : info) {
      t_data->Branch(key.c_str(), &content);
    }
    t_data->Fill();
    t_data->Write();
  }
  outf->cd();
  xjjroot::closefile(outf);
  inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
