#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"

#include "draw.h"

struct Property {
  std::string tex;
  int rorder = 0;
};

const float ytop = 0.84, lspace = 1.2, tsize = 0.036;
float ypos(float i = 0, float margin = 0) { return ytop+margin-i*(tsize*lspace); }

int syst_bdt(const std::string& inputstr, const std::string& outputname, int save_png = 0) {
  xjjc::whats_this("extrapolate BDT syst", "syst_bdt.cc");  
  const auto pinput = xjjroot::parse_input(inputstr);
  auto* inf = xjjroot::openfile(pinput.content);
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins");
  draw::bintex tbins(h3_bins, 0, 2);
  const auto colors_y = xjjroot::grayscales_color(tbins.ny()),
    colors_f = xjjroot::grayscales_color(tbins.ny(), xjjroot::mycolor_dark["red"]);
  const auto is_ptdiff = xjjc::str_contains(inputstr, "b-ptdiff");
  
  const std::map<std::string, Property> texs = {
    { "sum", { .tex = "Xn0n + 0nXn (#it{y} #rightarrow -#it{y})" } },
    { "gammaN", { .tex = "Xn0n (#gammaN)" } },
    { "Ngamma", { .tex = "Xn0n (N#gamma)", .rorder = 1 } },
    { "NgammaRef", { .tex = "0nXn (N#gamma) (#it{y} #rightarrow -#it{y})" } },
  };

  // const std::vector<float> nocut_point(tbins.ny(), -0.2);
  const std::vector<float> nocut_point = { -0.22, -0.22, -0.22, -0.22, -0.2, -0.2, -0.14, -0.14 };
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), { { "figspdf/" , "figs/" }, { ".pdf", "" } });
  for (int j=0; j<tbins.npt(); j++) {
    auto* dir = inf->GetDirectory(Form("dir__pt-%d", j));
    auto* dir_out = outf->mkdir(Form("dir__pt-%d", j));
    dir_out->cd();

    auto* dir_scan = inf->GetDirectory(Form("dir__pt-%d/scan", j));
    for (auto& [key, prop] : texs) {
      auto vg = xjjana::getobj_regexp<TGraphErrors>(dir_scan, Form("gr-scan_y_xsec_%s__pt-%d__y-[0-9]+", key.c_str(), j));
      const auto ymax = xjjana::gethsmaximum(vg);
      for (int i=0; i<tbins.ny(); i++)
        xjjroot::setthgrstyle(vg[i], colors_y[i], xjjroot::mstylelist_solid(i), 1.5, colors_y[i], 1, 1);
      auto* hempty = new TH2F(Form("hempty-scan_xsec_%s__pt-%d", key.c_str(), j), Form(";%s;%s", pinput.tex.c_str(), draw::title_xsec.c_str()),
                              10, -0.24, 0.06,
                              10, 0, ymax * 1.8);
      xjjroot::sethempty(hempty, 0, 0.2);
      hempty->GetXaxis()->SetNdivisions(505);

      auto* h1 = xjjana::getobj<TH1D>(dir, "h1_y_xsec_" + key);
      auto* gr_err = new TGraphAsymmErrors(tbins.ny());
      gr_err->SetName(("gr-err_y_xsec_" + key).c_str());
      xjjroot::setthgrstyle(gr_err, kBlack, 21, 1.5, kBlack, 1, 1, kBlack, 0.1, 1001);
      auto* gr_rel = new TGraphAsymmErrors(tbins.ny());
      gr_rel->SetName(("gr-rel_y_xsec_" + key).c_str());
      xjjroot::setthgrstyle(gr_rel, kBlack, 21, 1.5, kBlack, 1, 1, kBlack, 0.1, 1001);
      
      auto *leg_y = new TLegend(0.21, ypos(0.7 + std::ceil(tbins.ny()/3.)), 0.91, ypos(0.7));
      xjjroot::setleg(leg_y, tsize);
      leg_y->SetNColumns(3);
      for (int k=0; k<tbins.ny(); k++)
        leg_y->AddEntry(vg[k], tbins.label_y(k).c_str(), "p");
        
      pdf->prepare();
      hempty->Draw("axis");
      for (int i=0; i<tbins.ny(); i++) {
        if (is_ptdiff && i==0) continue;
        
        const auto ax = h1->GetBinCenter(i+1), axerr = h1->GetBinWidth(i+1)/2., ay = h1->GetBinContent(i+1);
        vg[i]->Draw("pe same");
        auto* f_fit = new TF1(Form("f_fit__pt-%d__y-%d", j, i), "[0]+[1]*x", nocut_point[i], 0.1);
        xjjroot::settfstyle(f_fit, colors_f[i], 1, 2);
        f_fit->SetParameters(ay, 0);
        double fit_min = -0.14, fit_max = 0.02;
        if (is_ptdiff && i==2) { fit_min = -0.1; fit_max = 0; }
        if (i < 2) fit_max = 0;
        
        auto r = vg[i]->Fit(f_fit, "SQ", "", fit_min, fit_max);
        // f_fit->SetLineStyle(2);
        f_fit->Draw("l same");

        const auto val_nocut = f_fit->Eval(nocut_point[i]), err_nocut = val_nocut - ay;
        // const auto err_max = std::max(val_nocut, xjjana::gethmaximum(vg[i])) - ay, err_min = std::min(val_nocut, xjjana::gethminimum(vg[i])) - ay;
        const auto err_max = val_nocut - ay, err_min = val_nocut - ay;
        const auto err_low = err_min < 0 ? std::abs(err_min) : 0;
        const auto err_high = err_max > 0 ? std::abs(err_max) : 0;
        gr_err->SetPoint(i, ax, ay);
        gr_err->SetPointError(i, axerr, axerr, err_low, err_high);
        gr_rel->SetPoint(i, ax, 0);
        gr_rel->SetPointError(i, axerr, axerr, err_low/ay, err_high/ay);
      }
      leg_y->Draw();
      // (prop.rorder ? leg_y_rorder : leg_y)->Draw();
      // xjjroot::drawtexgroup(0.23, ypos(-0.2), xjjc::str_divide_trim(tags, ","), tsize, 13);
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      const auto to_savepng = save_png && key=="sum";
      if (!prop.rorder)
        pdf->write(Form("%s_extrapolation_%s_pt-%d.pdf", name_png.c_str(), key.c_str(), j), to_savepng ? "" : "X");

      pdf->prepare();
      h1->Draw("pe1");
      gr_err->Draw("2 same");
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex, "Uncert. from BDT extrapolation" }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      if (!prop.rorder)
        pdf->write(Form("%s_err_%s_pt-%d.pdf", name_png.c_str(), key.c_str(), j), to_savepng ? "" : "X");

      const auto ymax_rel = std::max(xjjana::gethwerrmaximum(gr_rel)*2, 0.2),
        ymin_rel = std::min(xjjana::gethwerrminimum(gr_rel)*1.1, -0.1);
      auto* hempty_rel = new TH2D("hempty_rel", Form(";%s;Relative Uncertainty", h1->GetXaxis()->GetTitle()),
                                  10, h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),
                                  10, ymin_rel, ymax_rel);
      xjjroot::sethempty(hempty_rel);
      pdf->prepare();
      hempty_rel->Draw("axis");
      xjjroot::drawline_horizon(0, hempty_rel, kBlack, 2);
      gr_rel->Draw("p2 same");
      xjjroot::drawtexgroup(0.23, ypos(-0.2), { prop.tex, "Uncert. from BDT extrapolation" }, tsize, 13);
      xjjroot::drawtexgroup(0.90, ypos(-0.2), { tbins.label_pt(j) }, tsize, 33);
      xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
      if (!prop.rorder)
        pdf->write(Form("%s_rel_%s_pt-%d.pdf", name_png.c_str(), key.c_str(), j), to_savepng ? "" : "X");

      xjjroot::writehist(h1);
      xjjroot::writehist(gr_err);
      xjjroot::writehist(gr_rel);
    }
  }

  pdf->close();
  xjjroot::writehist(h3_bins);
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4)
    return syst_bdt(argv[1], argv[2], std::atoi(argv[3]));
  if (argc == 3)
    return syst_bdt(argv[1], argv[2]);
}
