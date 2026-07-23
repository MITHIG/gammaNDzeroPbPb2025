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

struct Input {
  draw::bintex tbins;
  std::map<std::string, xjjc::info> infos;
};

struct Fprompt {
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, TGraphAsymmErrors*> grs;
  std::map<std::string, std::vector<TH1D*>> h1ys;
  std::vector<TH1D*> h1s_fprompt_sf;
  std::string type_data = "", var = "";
  // std::map<std::string, std::vector<fpfitter*>> fitterys;
};

Fprompt init(TDirectory* dir) {
  Fprompt fp;
  for (const std::string key : {
      "fprompt_.+_best",
      "fprompt_.+_fix",
      "alpha_.+_best",
      "chi2_.+_best",
    }) {
    const auto name = xjjc::str_eraseall(key, "_.+");
    auto* h = xjjana::getobj_regexp_first<TH1D>(dir, "h1_.*" + key, "", false);
    if (h) fp.h1s[name] = h;
    else {
      auto* g = xjjana::getobj_regexp_first<TGraphAsymmErrors>(dir, "gr_.*" + key, "", false);
      if (g) fp.grs[name] = g;
    }
  }
  for (const std::string &key : { "fprompt", "chi2" }) {
    fp.h1ys[key] = xjjana::getobj_regexp<TH1D>(dir, "h1_.*" + key + ".*__y-.+");
    for (auto& h : fp.h1ys[key])
      xjjana::sethminmax(h, 0., 1.5);
  }
  fp.h1s_fprompt_sf = xjjana::getobj_regexp<TH1D>(dir->GetDirectory("dir_sfs"), ".+_sf-.+");
  for (auto& h : fp.h1s_fprompt_sf) {
    xjjana::sethabsminmax(h, 0., 1.5);
  }
  xjjroot::print_tab(fp.h1s, 0);
  xjjroot::print_tab(fp.grs, 0);

  const auto ny = fp.h1ys["chi2"].size();
  
  // for (const std::string type : { "best", "fix" }) {
  //   auto* dir_fit = dir->GetDirectory("dir_best");
  //   for (int i=0; i<ny; i++) {
  //     auto* hdata = xjjana::getobj_regexp_first<TH1D>(dir_fit, Form("h1_%s-data.+__y-%d.*", type.c_str(), i));
  //     auto* htotal = xjjana::getobj_regexp_first<TH1D>(dir_fit, Form("h1_%s-total.+__y-%d.*", type.c_str(), i));
  //     auto* hnonprompt = xjjana::getobj_regexp_first<TH1D>(dir_fit, Form("h1_%s-nonprompt.+__y-%d.*", type.c_str(), i));
  //     auto* fresult = xjjana::getobj_regexp_first<TFitResult>(dir_fit, ".+");
  //     auto* fitter = new fpfitter(hdata, htotal, hnonprompt, fresult);
  //     fp.fitterys.push_back(fitter);
  //   }
  // }
  return fp;
}

void style_fp(Fprompt& fp, const xjjroot::thgrstyle& ss) {
  for (auto& [_, h] : fp.h1s) {
    xjjroot::sethempty(h, 0, 0.1);
    xjjroot::setthgrstyle(h, ss);
  }
  for (auto& [_, hh] : fp.h1ys)
    for (auto& h : hh) {
      xjjroot::sethempty(h, 0, 0.1);
      xjjroot::setthgrstyle(h, ss);
    }
  for (auto& [_, h] : fp.grs)
    xjjroot::setthgrstyle(h, ss);
  for (auto& h : fp.h1s_fprompt_sf) {
    xjjroot::sethempty(h, 0, 0.1);
    xjjroot::setthgrstyle(h, xjjroot::color_alpha(ss.mcolor, 0.1), ss.mstyle, 0, xjjroot::color_alpha(ss.lcolor, 0.1), ss.lstyle);
  }
}

int macro(const std::vector<std::string>& inputnames, const std::string& outputname) {
  //
  std::map<std::string, std::vector<Fprompt>> fps;
  std::vector<Input> inputs;
  for (const auto &inputname : inputnames) {
    auto* inf = TFile::Open(inputname.c_str());
    auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
    Input in = { .tbins = draw::bintex(h3_bins, 0, 2) };
    for (const std::string name : { "data", "prompt", "nonprompt", "fit" } ) {
      in.infos[name] = xjjana::get_info(inf, Form("info/%s", name.c_str()));
    }
    inputs.push_back(in);

    const auto var = in.infos.at("fit").at("var");
    const auto ss_var = style_data(var, 1);

    for (const std::string &type_data : { "sigswap", "sub" }) {
      const auto ss_type = style_data(type_data);
      auto fp = init(inf->GetDirectory(Form("dir_data-%s", type_data.c_str())));
      fp.type_data = type_data;
      fp.var = var;
      style_fp(fp, { .mcolor = ss_var.color, .mstyle = ss_type.mstyle, .msize = 1.5, .lcolor = ss_var.color, .lstyle = 1, .lwidth = 1 });

      fps[type_data].push_back(fp);
    }
  }
  const auto nvar = inputs.size();
  // const auto var_tex = style_data(var, 1).title;
  // std::cout<<fps["sub"].front().h1s_fprompt_sf.size()<<std::endl;
  const auto& fp_dump = fps["sub"].front();
  auto* h1_dump_fprompt = (TH1D*)fp_dump.h1s_fprompt_sf.front()->Clone("h1_dump_fprompt");
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  auto* leg1 = new TLegend(0.46, 0.35-2*0.038*1.2, 0.8, 0.35);
  xjjroot::setleg(leg1, 0.038);
  leg1->SetNColumns(2);
  for (int i=0; i<nvar; i++) {
    for (auto& [key, fp] : fps)
      leg1->AddEntry(fp[i].h1s["alpha_best"], style_data(key).title.c_str(), "p");
  }
  leg1->Draw();

  for (const std::string& htype : { "best", "fix" }) {
    pdf->prepare();
    h1_dump_fprompt->Draw("axis");
    for (auto& [_, vfp] : fps) {
      for (auto& fp : vfp)
        for (auto& h : fp.h1s_fprompt_sf)
          h->Draw("pe1 same");
    }
    for (auto& [_, vfp] : fps) {
      for (auto& fp : vfp)
        fp.grs["fprompt_"+htype]->Draw("pe1 same");
    }
    xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
    xjjroot::drawtexgroup(0.24, 0.86, {
        inputs.front().tbins.label_pt(-1),
        inputs.front().infos.at("data").at("cut_tex"),
        "#bf{" + style_data(htype).title + "}",
      }, 0.038, 13, 42, 1.2);
    xjjroot::drawtexgroup(0.91, 0.86, {
        "Uncert. by toy smearing",
      }, 0.038, 33, 42, 1.2);
    leg1->Draw();
    xjjroot::drawtexgroup(0.46-0.01, 0.35-0.005, {
        style_data(fps["sub"][0].var, 1).title,
        style_data(fps["sub"][1].var, 1).title,      
      }, 0.038, 33, 62, 1.2, 1, 0, { style_data(fps["sub"][0].var, 1).color, style_data(fps["sub"][1].var, 1).color });
    pdf->getc()->RedrawAxis();
    pdf->write();
  }

  // auto* leg2 = new TLegend(0.24, 0.86-0.01-0.038*1.2*(3+2), 0.5, 0.86-0.01-0.038*1.2*3);
  // xjjroot::setleg(leg2, 0.038);
  // for (auto& [key, h] : h1s_alpha_best)
  //   leg2->AddEntry(h, style_data(key).title.c_str(), "p");
  
  // pdf->prepare();
  // h1s_alpha_best["data-sub"]->Draw("axis");
  // xjjroot::drawbox(h1s_alpha_best["data-sub"]->GetXaxis()->GetXmin(), h1_bins_sf->GetBinCenter(1),
  //                  h1s_alpha_best["data-sub"]->GetXaxis()->GetXmax(), h1_bins_sf->GetBinCenter(h1_bins_sf->GetNbinsX()), kBlack, 0.05);
  // for (auto& [_, h] : h1s_alpha_best)
  //   h->Draw("pl same");
  // xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
  // xjjroot::drawtexgroup(0.24, 0.86, {
  //     tbins.label_pt(-1),
  //     inputs.front().infos.at("data").at("cut_tex"),
  //     "#bf{Fit on " + var_tex + "}",
  //   }, 0.038, 13, 42, 1.2);
  // leg2->Draw();
  // pdf->getc()->RedrawAxis();
  // pdf->write();

  // pdf->prepare();
  // h1s_chi2_best["data-sub"]->Draw("axis");
  // xjjroot::drawbox(h1s_chi2_best["data-sub"]->GetXaxis()->GetXmin(), 0.8,
  //                  h1s_chi2_best["data-sub"]->GetXaxis()->GetXmax(), 1.2, kBlack, 0.05);
  // for (auto& [_, h] : h1s_chi2_best)
  //   h->Draw("pl same");
  // xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
  // xjjroot::drawtexgroup(0.24, 0.86, {
  //     tbins.label_pt(-1),
  //     inputs.front().infos.at("data").at("cut_tex"),
  //     "#bf{Fit on " + var_tex + "}",
  //   }, 0.038, 13, 42, 1.2);
  // leg2->Draw();
  // pdf->getc()->RedrawAxis();
  // pdf->write();

  // for (int i=0; i<ny; i++) {
  //   std::vector<TH1D*> hlist;
  //   for (auto& [_, fp] : fpys) {
  //     hlist.push_back(fp[i].h1s["chi2"]);
  //   }
  //   xjjana::sethsmin(hlist, 0);
  //   xjjana::sethsmax(hlist, 1.5);
    
  //   pdf->prepare();
  //   hlist.front()->Draw("axis");
  //   for (auto& [_, fp] : fpys) {
  //     xjjroot::drawline(fp[i].h1s["chi2"]->GetBinCenter(fp[i].ibin_best), 0,
  //                       fp[i].h1s["chi2"]->GetBinCenter(fp[i].ibin_best), fp[i].h1s["chi2"]->GetBinContent(fp[i].ibin_best),
  //                       fp[i].h1s["chi2"]->GetLineColor(), 2, 2);
  //   }
  //   for (auto& h : hlist) h->Draw("p same");
  //   xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
  //   xjjroot::drawtexgroup(0.24, 0.86, {
  //       tbins.label_pt(-1),
  //       tbins.label_y(i),
  //       inputs.front().infos.at("data").at("cut_tex"),
  //       "#bf{Fit on " + var_tex + "}",
  //     }, 0.038, 13, 42, 1.2);
  //   xjjroot::moveleg_n_draw(leg2, -1, 0.86-0.01-0.038*1.2*4);
  //   pdf->write();
  // }
  
  // for (int i=0; i<ny; i++) {
  //   std::vector<TH1D*> hlist;
  //   for (auto& [_, fp] : fpys) {
  //     hlist.push_back(fp[i].h1s["fprompt"]);
  //   }
  //   xjjana::sethsmin(hlist, 0.);
  //   xjjana::sethsmax(hlist, 1.5);
    
  //   pdf->prepare();
  //   hlist.front()->Draw("axis");
  //   for (auto& [_, fp] : fpys) {
  //     xjjroot::drawline(fp[i].h1s["fprompt"]->GetBinCenter(fp[i].ibin_best), 0,
  //                       fp[i].h1s["fprompt"]->GetBinCenter(fp[i].ibin_best), fp[i].h1s["fprompt"]->GetBinContent(fp[i].ibin_best),
  //                       fp[i].grs["fprompt"]->GetLineColor(), 2, 2);
  //   }
    
  //   for (auto& h : hlist) h->Draw("pe1 same");
  //   for (auto& [_, fp] : fpys) fp[i].grs["fprompt"]->Draw("pe1 same");
  //   xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
  //   xjjroot::drawtexgroup(0.24, 0.86, {
  //       tbins.label_pt(-1),
  //       tbins.label_y(i),
  //       inputs.front().infos.at("data").at("cut_tex"),
  //       "#bf{Fit on " + var_tex + "}",
  //     }, 0.038, 13, 42, 1.2);
  //   xjjroot::drawtexgroup(0.91, 0.86, {
  //       "Uncert. by toy smearing",
  //       "GetParError() as light lines"
  //     }, 0.038, 33, 42, 1.2, 1, 0.2, { kBlack, kGray });
  //   xjjroot::moveleg_n_draw(leg2, 0.6, 0.5);
  //   pdf->write();
  // }
  
  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    const auto inputs = xjjc::str_divide_trim(argv[1], ",");
    return macro(inputs, argv[2]);
  }
  return 1;
}
