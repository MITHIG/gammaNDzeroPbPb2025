#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/draw.h"

#define __COOK_NAME__
#define __DRAW_STYLE__
#include "style.h"

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
};

Fprompt init(TDirectory* dir) {
  Fprompt fp;
  for (const std::string key : {
      "fprompt.+_best",
      "fprompt.+_fix",
      "alpha.+_best",
      "chi2.+_best",
    }) {
    const auto name = xjjc::str_eraseall(key, ".+");
    auto* h = xjjana::getobj_regexp_first<TH1D>(dir, "h1_.*" + key, "", false);
    if (h) fp.h1s[name] = h;
    if (xjjc::str_contains(key, "fprompt")) {
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
  const auto& fp_dump = fps["sub"].front();
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  std::vector<std::string> titles_var; std::vector<Color_t> colors_var;
  auto* leg1 = new TLegend(0.47, 0.35-2*0.038*1.2, 0.85, 0.35);
  xjjroot::setleg(leg1, 0.038);
  leg1->SetNColumns(2);
  for (int i=0; i<nvar; i++) {
    auto var = fps.begin()->second.at(i).var;
    for (auto& [key, fp] : fps)
      leg1->AddEntry(fp[i].h1s.at("alpha_best"), style_data(key).title.c_str(), "p");
    titles_var.push_back(style_data(var, 1).title);
    colors_var.push_back(style_data(var, 1).color);
  }
  leg1->Draw();
  auto gtex1 = xjjroot::drawtexgroup(0.46-0.01, 0.35-0.005, titles_var, 0.038, 33, 62, 1.2, 1, 0, colors_var);

  for (const std::string& htype : { "best", "fix" }) {
    pdf->prepare();
    fp_dump.h1s_fprompt_sf.front()->Draw("axis");
    for (auto& [_, vfp] : fps) {
      for (auto& fp : vfp)
        for (auto& h : fp.h1s_fprompt_sf)
          h->Draw("pe1 same");
    }
    for (auto& [_, vfp] : fps) {
      for (auto& fp : vfp)
        fp.grs.at("fprompt_" + htype)->Draw("pe1 same");
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
    for (auto& t : gtex1) t->Draw();
    pdf->getc()->RedrawAxis();
    pdf->write();
  }

  for (const std::string& xvar : { "chi2", "alpha" }) {
    pdf->prepare();
    auto* h_dump = fp_dump.h1s.at(xvar + "_best");
    h_dump->Draw("axis");
    xjjroot::drawline(h_dump->GetXaxis()->GetXmin(), 1, h_dump->GetXaxis()->GetXmax(), 1, kGray+1, 2, 1);
    for (auto& [_, vfp] : fps) {
      for (auto& fp : vfp)
        fp.h1s.at(xvar + "_best")->Draw("pl same");
    }
    xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
    xjjroot::drawtexgroup(0.24, 0.86, {
        inputs.front().tbins.label_pt(-1),
        inputs.front().infos.at("data").at("cut_tex"),
        // "#bf{" + style_data("best").title + "}",
      }, 0.038, 13, 42, 1.2);
    xjjroot::moveleg_n_draw(leg1, -1, 0.75);
    xjjroot::movetexgroup_n_draw(gtex1, -1, 0.75-0.005);
    pdf->getc()->RedrawAxis();
    pdf->write();
  }

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
