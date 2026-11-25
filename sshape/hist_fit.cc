#include <TH2D.h>

#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"

#include "dfitter.h"
#include "draw.h"

struct fitfunc {
  std::string name, formula, title;
  std::vector<double> init;
  Color_t color;
};
std::string dumpfunc(std::string formula, const TF1* f);

int macro(const std::string& input_data, const std::string& input_template, const std::string& outputname,
          const std::string& fit_opt = "3G-Peaky", int save_png = 0) {

  std::cout<<std::endl;
  // double xfitmin = -1.5, double xfitmax = 1.5;
  std::vector<fitfunc> funcs = {
    { "pol0", "[0]", "Constant", { 0.95 }, xjjroot::mycolor_middle["blue"] },
    { "mshape2", "[0]+[1]*(exp(-0.5*((x-0.8)/[2])^2)+exp(-0.5*((x+0.8)/[2])^2))", "Reflected two Gaus", { 0.93, 0.03, 0.3 }, xjjroot::mycolor_middle["red"] },
    // { "mshape3", "[0]+[1]*(exp(-0.5*((x-[3])/[2])^2)+exp(-0.5*((x+[3])/[2])^2))", "[0] + [1]#upoint2Gaus(#pm[3], [2])", { 0.93, 0.03, 0.3, 0.8 }, xjjroot::mycolor_middle["green"] },
  };
  const auto fitopt = xjjroot::parse_input(fit_opt);
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

  std::map<std::string, TH2D*> h2s;
  std::map<std::string, std::vector<TH1D*>> h1vys;
  std::map<std::string, TH1D*> h1ys;

  auto read_hists = [&h2s, &h1vys](TFile* inf) {
    for (auto* h2 : xjjana::getobj_regexp<TH2D>(inf, "h2_.+")) {
      auto name = xjjc::str_eraseall(h2->GetName(), { "h2_" });
      if (h2s.find(name) != h2s.end()) {
        __XJJLOG << "!! error: name " << name << " already in the map." << std::endl;
        xjjroot::print_tab(h2s, 0);
        return 3;
      }

      std::cout<<name<<std::endl;
      auto ny = h2->GetXaxis()->GetNbins();
      std::cout<<ny<<std::endl;
      for (int i=0; i<h2->GetXaxis()->GetNbins(); i++) {
        auto* h1_mass = h2->ProjectionY(Form("h1_mass_%s__y-%d", name.c_str(), i),
                                        i+1, i+1, // y bin
                                        "e");
        h1vys[name].push_back(h1_mass);
      }
      h2s[name] = h2;
    }
    return 0;
  };

  if (read_hists(inf_data)) { return 3; }
  if (read_hists(inf_template)) { return 3; }

  auto* h2_bins = static_cast<TH2D*>(h2s.at("data")->Clone("h2_bins")); h2_bins->Reset();
  const draw::bintex tbins(h2_bins, 0, 0);

  // prepare TH1
  auto make_h1_y = [&tbins, &h1ys](const std::string& name, const std::string& ytitle) {
    auto* h = tbins.make_h1_y<TH1D>(Form("h1_y_%s", name.c_str()));
    h->GetYaxis()->SetTitle(ytitle.c_str());
    h1ys[name] = h;
  };
  make_h1_y("yield", "Raw Yield");
  make_h1_y("width68mc", "Signal Effective#scale[0.4]{ }#sigma in MC [GeV]");
  make_h1_y("walpha", "#alpha =#scale[0.4]{ }#frac{#sigma_{data}}{#sigma_{MC}}");

  xjjroot::print_tab(h2s, 0);
  xjjroot::print_tab(h1ys, 0);

  xjjroot::setgstyle(1);
  auto dfs = std::vector<xjjroot::dfitter*>(tbins.ny(), nullptr);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });
  for (int i=0; i<tbins.ny(); i++) {
    pdf->prepare();
    const auto &h = h1vys.at("data")[i],
      &hmc = h1vys.at("mc-match")[i], &hmcswap = h1vys.at("mc-swap")[i],
      &hmckk = h1vys.at("mc-kk")[i], &hmcpipi = h1vys.at("mc-pipi")[i];

    auto* df = new xjjroot::dfitter(fitopt.content.c_str());
    df->fit(h, hmc, hmcswap, hmckk, hmcpipi);
    xjjroot::drawtexgroup(0.24, 0.86, {
        tbins.label_y(i),
        // "#bf{" + info_data.at("cut_tex") + "}"
      }, 0.035, 13);
    df->draw_leg();
    df->draw_result(0.24, 0.86-1*0.035*1.15, 0.035);
    xjjroot::drawCMS(xjjroot::CMS::internal, info_data.at("input_tex"));

    h1ys.at("yield")->SetBinContent(i+1, df->yield());
    h1ys.at("yield")->SetBinError(i+1, df->yieldErr());
    const auto w68mc = df->width_mc_match(xjjana::frac_1sigma);
    h1ys.at("width68mc")->SetBinContent(i+1, w68mc.first);
    h1ys.at("width68mc")->SetBinError(i+1, w68mc.second);
    const auto walpha = df->alpha();
    h1ys.at("walpha")->SetBinContent(i+1, 1+walpha.first);
    h1ys.at("walpha")->SetBinError(i+1, walpha.second);

    dfs[i] = df;
    gPad->RedrawAxis();
    pdf->write();
  }

  for (int i=0; i<tbins.ny(); i++) {
    pdf->prepare();
    auto* df = dfs[i];

    for (const std::string s : { "match", "swap", "kk", "pipi" }) {
      if (!df->haskkpipi() && (s == "kk" || s == "pipi")) continue;
      auto* hmc = h1vys.at("mc-" + s)[i];
      df->set_hist(hmc);
      if (xjjroot::fstyle.find(s) != xjjroot::fstyle.end() && s != "match")
        xjjroot::setthgrstyle(hmc, xjjroot::fstyle.at(s).lcolor, -1, -1, xjjroot::fstyle.at(s).lcolor, -1, -1, -1, -1, -1, 0.6, 0.6);
      hmc->Draw(s == "match" ? "pe1" : "pe1 same");
    }
    df->draw_fmc();
    df->draw_legmc();
    xjjroot::drawtexgroup(0.24, 0.86, { tbins.label_y(i) }, 0.035, 13);
    df->draw_params(0.24, 0.86-0.035*1.15, 0.035);
    // xjjroot::drawtexgroup(0.88, 0.86, { info_fit.at("fitopt_tex") }, 0.035, 33);
    xjjroot::drawCMS(xjjroot::CMS::simulation, info_template.at("input_tex"));
    gPad->RedrawAxis();
    pdf->write();
  }
  
  for (const std::string& p : { "walpha", "yield", "width68mc" }) {
    auto* h = h1ys.at(p);
    pdf->prepare();
    const auto ymin = xjjana::gethminimum(h), ymax = xjjana::gethmaximum(h);
    auto yaxismin = ymin >= 0 ? 0 : 1.6*ymin;
    auto yaxismax = (ymax >= 0 ? ymax : (0-ymin)*0.5) * 1.6;
    if (p == "walpha") {
      yaxismin = ymin * 0.8;
      yaxismax = ymax * 1.1;
    }
    // __XJJLOG << ">> " << p << " [" << ymin << " - " << ymax << "]"
    //          << " -> [" << yaxismin << " - " << yaxismax << "]" << std::endl;
    // if (ymax <= 0 && ymin < 0)
    //   ymax = 0 - ymin;
    auto* hempty = new TH2F("hempty", Form(";y;%s", h->GetYaxis()->GetTitle()),
                            10, -2, 2,
                            10, yaxismin, yaxismax);
    xjjroot::sethempty(hempty, 0, 0.3);

    const auto cc = kBlack;
    xjjroot::setthgrstyle(h, kBlack, 21, 1.5, kBlack, 1, 1);
    hempty->Draw("axis");
    h->Draw("pe1 same");
    // xjjroot::drawtexgroup(0.24, 0.85, {
    //     info_data.at("cut_tex"),
    //     info_fit.at("fitopt_tex"),
    //   }, 0.035, 13, 42, 1.25);
    xjjroot::drawCMS(xjjc::str_contains(p, "mc") ? xjjroot::CMS::simulation : xjjroot::CMS::internal,
                     xjjc::str_contains(p, "mc") ? info_template.at("input_tex") : info_data.at("input_tex"));

    if (p == "walpha") {
      xjjroot::drawline_horizon(1, hempty, kBlack, 2, 2);
      const auto fmin = hempty->GetXaxis()->GetXmin(), fmax = hempty->GetXaxis()->GetXmax(),
        fitmin = h->GetXaxis()->GetXmin(), fitmax = h->GetXaxis()->GetXmax();
      xjjroot::drawbox(fmin, hempty->GetYaxis()->GetXmin(), fitmin, hempty->GetYaxis()->GetXmax(), kBlack, 0.05);
      xjjroot::drawbox(fmax, hempty->GetYaxis()->GetXmin(), fitmax, hempty->GetYaxis()->GetXmax(), kBlack, 0.05);
      auto* leg = new TLegend(0.40, 0.48-0.035*1.25*funcs.size(), 0.70, 0.48);
      xjjroot::setleg(leg, 0.035);
      for (auto& ff : funcs) {
        auto* f = new TF1(Form("f_walpha_%s", ff.name.c_str()), ff.formula.c_str(), hempty->GetXaxis()->GetXmin(), hempty->GetXaxis()->GetXmax());
        xjjroot::settfstyle(f, ff.color, 1, 3);
        f->SetNpx(400);
        for (size_t i = 0; i < ff.init.size(); i++) f->SetParameter(i, ff.init[i]);
        auto r = h->Fit(f, "QNS0R", "", fitmin, fitmax);
        // auto r = h->Fit(f, "S", "", fitmin, fitmax);
        // f->SetLineStyle(2);
        // f->SetLineColorAlpha(ff.color, 0.4);
        f->Draw("same");
        leg->AddEntry(f, ff.title.c_str(), "l");
        xjjroot::drawtexgroup(xjjroot::get_pad_center('X'), 0.85, {
            "Signal resolution difference",
            "(Data / MC)"
          }, 0.035, 23);
        std::cout << dumpfunc(ff.formula, f) << std::endl;
      }
      leg->Draw();
      // gROOT->GetListOfFunctions()->Print();
    }
    auto save_indi = save_png && p == "walpha";
    pdf->write(Form("%s_%s.pdf", name_png.c_str(), p.c_str()), save_indi ? "" : "X");

    delete hempty;
  }

  pdf->draw_cover( {
      "#bf{Data} " + info_data.at("input"),
      "#bf{Cut} " + info_data.at("cut"),
      "#bf{Template} " + info_template.at("input"),
      "#bf{Cut} " + info_template.at("cut"),
    }, 0.03);

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (const auto& [_, h] : h1ys)
    xjjroot::writehist(h);
  // for (const auto& [_, hh] : h1vys)
  //   for (const auto& h : hh)
  //     xjjroot::writehist(h);
  xjjroot::writehist(h2_bins);
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 6)
    return macro(argv[1], argv[2], argv[3], argv[4], std::atoi(argv[5]));
  return 1;
}

std::string dumpfunc(std::string formula, const TF1* f) {
  for (int i = f->GetNpar() - 1; i >= 0; i--) {
    const std::string key = "[" + std::to_string(i) + "]";
    // const std::string val = Form("(%.2e)", f->GetParameter(i));
    const std::string val = Form("%.3f", f->GetParameter(i));
    for (size_t pos = formula.find(key); pos != std::string::npos; pos = formula.find(key, pos + val.size()))
      formula.replace(pos, key.size(), val);
  }
  return formula;
}
