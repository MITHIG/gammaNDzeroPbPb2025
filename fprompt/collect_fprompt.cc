#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "draw.h"

#define __COOK_NAME__
#define __DRAW_STYLE__
#include "style.h"

#define __BINS_PTY_PLACEHOLDER__
#include "bins.h"

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
    const auto name = xjjc::str_eraseall(key, { ".+" });
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

void envelope_to_hist(const std::vector<TGraphAsymmErrors*>& gs, TH1D* h);

int macro(const std::vector<std::string>& inputnames, const std::string& outputname,
          const std::string& newbintag, int save_png) {
  bins::print();
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

  auto* h1_fprompt = (TH1D*)fp_dump.h1s_fprompt_sf.front()->Clone("h1_fprompt");
  h1_fprompt->Reset("ICES");
  std::vector<TGraphAsymmErrors*> h1s_best;
  for (auto& [_, vfp] : fps) {
    for (auto& fp : vfp)
      h1s_best.push_back(fp.grs.at("fprompt_best"));
  }
  envelope_to_hist(h1s_best, h1_fprompt);
  auto* h1_fprompt_rebin = new TH1D("h1_fprompt_rebin", "", bins::ybins.size()-1, bins::ybins.data());
  xjjroot::setthgrstyle(h1_fprompt_rebin, kBlack, 47, 1.4, kGray+2, 1, 1, kGray, 0.5, 1001);
  auto* h2_fprompt_rebin = new TH2D("h2_fprompt_rebin", ";y;#it{p}_{T} (GeV)",
                                    bins::ybins.size()-1, bins::ybins.data(),
                                    bins::ptbins.size()-1, bins::ptbins.data());
  for (int i = 0; i < h2_fprompt_rebin->GetNbinsX(); i++) {
    const double x = h1_fprompt_rebin->GetBinCenter(i+1);
    const int k = h1_fprompt->FindBin(x);
    const auto content = h1_fprompt->GetBinContent(k),
      error = h1_fprompt->GetBinError(k);
    h1_fprompt_rebin->SetBinContent(i+1, content);
    h1_fprompt_rebin->SetBinError  (i+1, error);
    for (int j=0; j < h2_fprompt_rebin->GetNbinsY(); j++) {
      h2_fprompt_rebin->SetBinContent(i+1, j+1, content);
      h2_fprompt_rebin->SetBinError  (i+1, j+1, error);
    }
  }
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = draw::png_name(pdf);
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
        inputs.front().infos.at("data").at("ecut_tex"),
        "#bf{" + style_data(htype).title + "}",
      }, 0.038, 13, 42, 1.2);
    xjjroot::drawtexgroup(0.91, 0.86, {
        "Uncert. by toy smearing",
      }, 0.038, 33, 42, 1.2);
    leg1->Draw();
    for (auto& t : gtex1) t->Draw();
    pdf->getc()->RedrawAxis();
    pdf->write(Form("%s/fprompt_best.pdf", name_png.c_str()), save_png && htype == "best" ? "" : "X");
  }

  pdf->prepare();
  fp_dump.h1s_fprompt_sf.front()->Draw("axis");
  h1_fprompt_rebin->Draw("pe2 same");
  xjjroot::drawCMS(xjjroot::CMS::internal, inputs.front().infos.at("data").at("input_tex") + " (5.36 TeV)");
  xjjroot::drawtexgroup(0.24, 0.86, {
      inputs.front().tbins.label_pt(-1),
      inputs.front().infos.at("data").at("ecut_tex"),
    }, 0.038, 13, 42, 1.2);
  pdf->getc()->RedrawAxis();
  pdf->write(Form("%s/fprompt_syst.pdf", name_png.c_str()), save_png ? "" : "X");
  
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
        inputs.front().infos.at("data").at("ecut_tex"),
        // "#bf{" + style_data("best").title + "}",
      }, 0.038, 13, 42, 1.2);
    xjjroot::moveleg_n_draw(leg1, -1, 0.75);
    xjjroot::movetexgroup_n_draw(gtex1, -1, 0.75-0.005);
    pdf->getc()->RedrawAxis();
    pdf->write(Form("%s/%s.pdf", name_png.c_str(), xvar.c_str()), save_png ? "" : "X");
  }

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + "/" + newbintag.c_str() + ".root");
  xjjroot::writehist(h2_fprompt_rebin);
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 7) {
    bins::ybins = xjjc::str_convert_vector<double>(argv[4], ",");
    bins::ptbins = xjjc::str_convert_vector<double>(argv[5], ",");
    const auto inputs = xjjc::str_divide_trim(argv[1], ",");
    return macro(inputs, argv[2], argv[3], std::atoi(argv[6]));
  }
  // if (argc == 3) {
  //   const auto inputs = xjjc::str_divide_trim(argv[1], ",");
  //   return macro(inputs, argv[2]);
  // }
  return 1;
}

void envelope_to_hist(const std::vector<TGraphAsymmErrors*>& gs, TH1D* h) {
  for (int i = 0; i < h->GetNbinsX(); ++i) {
    double ymin =  std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    for (const auto* g : gs) {
      if (!g || i >= g->GetN())
        continue;

      const double y = g->GetPointY(i);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    if (!std::isfinite(ymin))
      continue;
    h->SetBinContent(i + 1, 0.5 * (ymax + ymin));
    h->SetBinError  (i + 1, 0.5 * (ymax - ymin));
  }
}
