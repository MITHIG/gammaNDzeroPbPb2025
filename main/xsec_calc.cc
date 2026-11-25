#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "xjjstruct.h"

#include "measurements.h"
#include "draw.h"

namespace global {
  float BR_DtoKpi = 0.03936, err_BR_DtoKpi = 0.030*1.e-2;
}

int xsec_calc(const std::string& inputname_raw, const std::string& inputname_effd,
          const std::string& inputname_effevent, // const std::string inputname_fprompt,
          float lumi, const std::string& outputdir, float f_emdpu, int save_png = 0) {
  xjjc::whats_this("calculate cross sections", "xsec_calc.cc");
  __XJJLOG << ">> lumi: " << lumi << " nb-1" << std::endl;
  __XJJLOG << ">> f_emdpu: " << f_emdpu << std::endl;
  std::map<std::string, std::vector<TH1D*>> h1pts;
  // std::map<std::string, std::vector<TH1D*>> h1ys;
  std::map<std::string, std::map<std::string, std::string>> infos;
  auto* h3_bins = xjjana::getobj<TH3D>(inputname_raw + "::h3_bins");
  if (!h3_bins)
    return 2;
  draw::bintex tbins(h3_bins, 0, 2);
  //
  auto get_h1pts = [&h1pts, &infos](const std::string &inputname, const std::string &category,
                                    const std::vector<std::string>& h1names,
                                    const std::vector<std::string>& infots) {
    __XJJLOG << "++ [" << category << "] \e[2m" << inputname << "\e[0m" << std::endl;
    auto* inf = xjjroot::openfile(inputname.c_str());
    if (!inf) {
      __XJJLOG << "?? no " << category << " input file: " << inputname << ", skip." << std::endl;
      return 2;
    }
    // TH1D
    std::vector<Color_t> colors = { xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["red"] };
    for (auto& name : h1names) {
      // auto vh = xjjana::getobj_regexp<TH1D>(inf, name + ".*__pt-.+");
      auto vh = xjjana::getobj_regexp<TH1D>(inf, name + "__pt-.+");
      if (vh.empty()) {
        __XJJLOG << "?? TH1D " << name << " not found, skip." << std::endl;
        continue;
      }
      for (int i=0; i<vh.size(); i++) {
        xjjroot::sethempty(vh[i], 0, 0.4);
        auto cc = vh.size() > 1 ? colors[i] : kBlack;
        xjjroot::setthgrstyle(vh[i], cc, 20, 1.6, cc, 1, 1);
      }
      h1pts[xjjc::str_eraseall(name, { "h1_" })] = vh;
    }
    // info
    for (auto& tr : infots) {
      auto info = xjjroot::read_info(inf, tr);
      auto name = category;
      if (tr != "info")
        name += ("_" + xjjc::str_eraseall(tr, { "/info" }));
      infos[name] = info;
      __XJJLOG << ">> [" << category << "] info : " << name << std::endl;
    }
    return 0;
  };

  get_h1pts(inputname_raw, "raw", { "h1_y_yield" }, { "data/info", "template/info" });
  get_h1pts(inputname_effd, "effd", { "h1_y_eff__rebin" }, { "info" });
  get_h1pts(inputname_effevent, "effevent", { "h1_y_evteff" }, { "info" });

  Event event_is = xjjc::str_contains(infos.at("raw_data").at("cut_tag"), "gammaN") ? Event::gammaN :
    (xjjc::str_contains(infos.at("raw_data").at("cut_tag"), "Ngamma") ? Event::Ngamma : Event::Other);
  __XJJLOG << ">> event_is = " << static_cast<int>(event_is) << ", accordint to " << infos.at("raw_data").at("cut_tag") << std::endl;

  xjjroot::print_tab(h1pts, 0);
  for (int j=0; j<tbins.npt(); j++) {
    const auto* h_yield = h1pts.at("y_yield")[j];

    auto* h_yieldlumi = (TH1D*)h_yield->Clone(xjjc::str_replaceall(h_yield->GetName(), "yield", "yield-lumi").c_str());
    h_yieldlumi->Scale(1./(lumi*1.e6));
    h_yieldlumi->GetYaxis()->SetTitle("Raw Yield / Luminosity [mb]");
    h1pts["y_yield-lumi"].push_back(h_yieldlumi);

    auto* h_corr = (TH1D*)h_yield->Clone(xjjc::str_replaceall(h_yield->GetName(), "yield", "corr").c_str());
    h_corr->Divide(h1pts.at("y_eff__rebin")[j]);
    h_corr->GetYaxis()->SetTitle("Corrected Yield");
    h1pts["y_corr"].push_back(h_corr);

    auto* h_xsec_nofemdpu = (TH1D*)h_corr->Clone(xjjc::str_replaceall(h_corr->GetName(), "corr", "xsec-nofemdpu").c_str());
    h_xsec_nofemdpu->Divide(h1pts.at("y_evteff")[j]);
    h_xsec_nofemdpu->Scale(0.5/global::BR_DtoKpi/(lumi*1.e6)/tbins.binwidth_pt(j)/*dpt*/, "width");
    h_xsec_nofemdpu->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{y}d#it{p}_{T}} [mb/GeV] (no EMD PU corr)");
    xjjroot::sethempty(h_xsec_nofemdpu, 0, 0.2); // closer to axis than other hists
    h1pts["y_xsec-nofemdpu"].push_back(h_xsec_nofemdpu);

    auto* h_xsec = (TH1D*)h_xsec_nofemdpu->Clone(xjjc::str_replaceall(h_xsec_nofemdpu->GetName(), "xsec-nofemdpu", "xsec").c_str());
    h_xsec->Scale(1./f_emdpu);
    h_xsec->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#it{y}d#it{p}_{T}} [mb/GeV]");
    xjjroot::sethempty(h_xsec, 0, 0.2); // closer to axis than other hists
    h1pts["y_xsec"].push_back(h_xsec);

    if (tbins.npt() == 1) {
      xjjroot::setthgrstyle(h_xsec_nofemdpu, xjjroot::mycolor_middle["red"], -1, -1, xjjroot::mycolor_middle["red"]);
      xjjroot::setthgrstyle(h_xsec, xjjroot::mycolor_middle["red"], -1, -1, xjjroot::mycolor_middle["red"]);
    }
  }

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputdir + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/", "figs/" }, { ".pdf", "" }});

  auto draw_global = [&infos, &lumi](bool divide2) {
    xjjroot::drawCMS(xjjroot::CMS::internal, xjjc::str_replaceall(infos.at("raw_data").at("input_tex"), ")", Form(", %.1f#scale[0.3]{ }#mub^{-1})", lumi*1.e3)));
    xjjroot::drawtexgroup(0.23, 0.85, { infos.at("raw_data").at("cut_tex") }, 0.04, 13, 42, 1.25);
    xjjroot::drawtexgroup(0.88, 0.86, { divide2 ? xjjroot::CMS::DzDzbar2 : xjjroot::CMS::DznDzbar }, 0.043, 33, 42, 1.25);
  };
  
  xjjana::sethsmin(h1pts.at("y_xsec"), 0.);
  auto ymax = xjjana::sethsmax(h1pts.at("y_xsec"), 1.8);
  if (ymax < 2.5) xjjana::sethsabsmax(h1pts.at("y_xsec"), 3.2);
  else if (ymax >= 2.5 && ymax < 5.) xjjana::sethsabsmax(h1pts.at("y_xsec"), 6.6);

  const auto x1 = (event_is==Event::Ngamma ? 0.222 : 0.54), y1 = 0.755,
    tsize = 0.038, lspace = 1.25, tlsize = tsize*lspace;
  auto* legvs23 = new TLegend(x1, y1-tlsize*0.9 - tlsize*(tbins.npt()+2), x1+0.3, y1-tlsize*0.9);
  xjjroot::setleg(legvs23, tsize);
  for (int i=0; i<tbins.npt(); i++)
    legvs23->AddEntry(h1pts.at("y_xsec")[i], tbins.label_pt(i).c_str(), "p");
  legvs23->AddEntry((TObject*)0, "", NULL);
  legvs23->AddEntry(measurement::get_style(), xjjroot::str_fixspace(xjjc::number_range_string(float(2), float(5), "#it{p}_{T}") + " GeV").c_str(), "pf");

  pdf->prepare();
  h1pts.at("y_xsec").front()->Draw("axis");
  measurement::draw_HIN_25_002(event_is);
  for (auto& h : h1pts.at("y_xsec"))
    h->Draw("pe1 same");
  draw_global(true);
  legvs23->Draw();
  xjjroot::drawtex(legvs23->GetX1() + 0.008, legvs23->GetY2()+tlsize*0.5 - 0.005, "This analysis", tsize, 12);
  xjjroot::drawtex(legvs23->GetX1() + 0.008, legvs23->GetY1()+tlsize*1.5 - 0.005, "2023 PbPb (HIN-25-002)", tsize, 12);
  pdf->write();

  // for (auto& h : h1pts.at("y_xsec"))
  //   xjjroot::print_th(h);
  
  auto* leg = new TLegend(x1, y1-tlsize*tbins.npt(), x1+0.3, y1);
  xjjroot::setleg(leg, tsize);
  for (int i=0; i<tbins.npt(); i++)
    leg->AddEntry(h1pts.at("y_yield")[i], tbins.label_pt(i).c_str(), "p");
  leg->Draw();

  for (auto& t : { "y_yield", "y_yield-lumi", "y_eff__rebin", "y_evteff", "y_corr", "y_xsec-nofemdpu", "y_xsec" }) {
    pdf->prepare();
    xjjana::sethsmin(h1pts.at(t), 0.);
    xjjana::sethsmax(h1pts.at(t), 1.8);
    h1pts.at(t).front()->Draw("axis");
    for (auto& h : h1pts.at(t)) {
      if (h1pts.size() == 1)
        xjjroot::setthgrstyle(h, kBlack, -1, -1, kBlack, -1, -1);
      h->Draw("pe1 same");
    }
    leg->Draw();
    draw_global(xjjc::str_contains(t, "xsec"));
    pdf->write();
    // pdf->write(name_png + xjjc::str_replaceall(t, "y_", "_") + ".pdf", save_png ? "" : "X");
  }

  pdf->draw_cover({
      "#bf{Data}",
      "#bf{Input} " + infos.at("raw_data")["input"],
      "#bf{Cut} " + infos.at("raw_data")["cut"],
    }, 0.038);
  pdf->draw_cover({
      "#bf{Mass template}",
      "#bf{Input} " + infos.at("raw_template")["input"],
      "#bf{Cut} " + infos.at("raw_template")["cut"],
    }, 0.038);
  pdf->draw_cover({
      "#bf{D efficiency}",
      "#bf{Input} " + infos.at("effd")["inputmc"],
      "#bf{Numerator} " + infos.at("effd")["cut_eff_num"],
      "#bf{Denominator} " + infos.at("effd")["cut_eff_den"],
    }, 0.038);
  pdf->draw_cover({
      "#bf{Event selection efficiency}",
      "#bf{D cut} " + infos.at("effevent")["cut_d"],
      "#bf{Event selection} " + infos.at("effevent")["cut_evt"],
    }, 0.038);

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputdir + ".root");
  // for (auto& [_, h] : h2s) xjjroot::writehist(h);
  for (auto& [_, hh] : h1pts)
    for (auto& h : hh)
      xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);

  xjjroot::Writeinfo tinfo("info");
  tinfo.cast_branch("lumi", lumi);
  tinfo.cast_branch("f_emdpu", f_emdpu);
  tinfo.cast_branch("event_is", static_cast<int>(event_is));
  tinfo.cast_branch("data_input", infos.at("raw_data")["input"]);
  tinfo.cast_branch("data_input_tex", infos.at("raw_data")["input_tex"]);
  tinfo.cast_branch("data_cut", infos.at("raw_data")["cut"]);
  tinfo.cast_branch("data_cut_tex", infos.at("raw_data")["cut_tex"]);
  tinfo.cast_branch("template_input", infos.at("raw_template")["input"]);
  tinfo.cast_branch("template_input_tex", infos.at("raw_template")["input_tex"]);
  tinfo.cast_branch("template_cut", infos.at("raw_template")["cut"]);
  tinfo.cast_branch("template_cut_tex", infos.at("raw_template")["cut_tex"]);
  tinfo.cast_branch("effd_input", infos.at("effd")["inputmc"]);
  tinfo.cast_branch("effd_input_tex", infos.at("effd")["inputmc_tex"]);
  tinfo.cast_branch("deff_cutd", infos.at("effd")["cutd"]);
  tinfo.cast_branch("deff_cutd_tex", infos.at("effd")["cutd_tex"]);
  tinfo.cast_branch("deff_cutevt", infos.at("effd")["cutevt"]);
  tinfo.cast_branch("deff_cutevt_tex", infos.at("effd")["cutevt_tex"]);
  tinfo.cast_branch("deff_cuteffnum", infos.at("effd")["cut_eff_num"]);
  tinfo.cast_branch("deff_cuteffden", infos.at("effd")["cut_eff_den"]);
  tinfo.cast_branch("evteff_cutd", infos.at("effevent")["cut_d"]);
  tinfo.cast_branch("evteff_cutd_tex", infos.at("effevent")["cut_d_tex"]);
  tinfo.cast_branch("evteff_cutevt", infos.at("effevent")["cut_evt"]);
  tinfo.cast_branch("evteff_cutevt_tex", infos.at("effevent")["cut_evt_tex"]);
  tinfo.close();

  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 8) {
    return xsec_calc(argv[1], argv[2], argv[3], atof(argv[4]), argv[5], std::atof(argv[6]), std::atoi(argv[7]));
  }
  if (argc == 7) {
    return xsec_calc(argv[1], argv[2], argv[3], atof(argv[4]), argv[5], std::atof(argv[6]));
  }
  return 1;
}

