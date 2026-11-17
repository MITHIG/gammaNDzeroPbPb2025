#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"

#include "draw.h"

int macro(const std::string inputname, const std::string& outputname) {
  auto* inf = xjjroot::openfile(inputname);
  if (!inf) return 2;

  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
  const draw::bintex tbins(h3_bins, 0, 2);
  
  auto* dir = (TDirectory*)inf->GetDirectory("dir_Dmva_BDT");
  std::map<std::string, std::vector<TH1D*>> h1s;
  std::vector<TF1*> fs;
  h1s["sub"] = xjjana::getobj_regexp<TH1D>(dir, "h1_Dmva_BDT_data-sub__y-.+");
  h1s["sigswap"] = xjjana::getobj_regexp<TH1D>(dir, "h1_Dmva_BDT_data-sigswap__y-.+");
  h1s["sideband"] = xjjana::getobj_regexp<TH1D>(dir, "h1_Dmva_BDT_data-sideband__y-.+");
  h1s["match"] = xjjana::getobj_regexp<TH1D>(dir, "h1_Dmva_BDT_mc-match__y-.+");

  const float fdiv = 1.9/3, tsize = 0.038, lspace = 1.2, tsize_up = tsize/fdiv*0.95;
  auto y_pos_up = [&tsize_up, &lspace, &fdiv](float n = 0) -> float {
    const float up = 1 - xjjroot::margin_pad_top - 1.5*tsize_up;
    // __XJJLOG << ">> up: " << up << ", this y: " << up - n*tsize_up*lspace << std::endl;
    return up - n*tsize_up*lspace;
  };
  auto* leg = new TLegend(0.52, y_pos_up(3), 0.70, y_pos_up());
  xjjroot::setleg(leg, tsize_up);
  for (auto& [key, vhs] : h1s) {
    for (auto& h : vhs) {
      xjjroot::sethempty(h, 0, 0);
      h->Scale(1./h->Integral(), "width");
      h->GetYaxis()->SetTitle("Self-normalized");
      xjjana::sethminmax(h, 0., 2.5);
      if (key == "sideband")
        xjjroot::setthgrstyle(h, kGray, 20, 1.5, kGray, 1, 1, kGray, 0.2, 1001);
    }
  }
  for (int i=0; i<tbins.ny(); i++) {
    auto* h_sub = (TH1D*)h1s.at("sub")[i]->Clone(xjjc::str_replaceall(h1s.at("sub")[i]->GetName(), {{ "data-sub", "ratio-sub" }}).c_str());
    h_sub->Divide(h1s.at("match")[i]);
    xjjana::sethminmax(h_sub, 0., 1.5);
    if (h_sub->GetMaximum() > 5)
      h_sub->SetMaximum(2.9);
    h_sub->GetYaxis()->SetTitle("Data / MC");
    h1s["ratio-sub"].push_back(h_sub);
    // fit ratio with a smoothed step; skip bins where the MC shape is negligible
    auto* h_fit = (TH1D*)h_sub->Clone(Form("%s_forfit", h_sub->GetName()));
    for (int b=1; b<=h_fit->GetNbinsX(); b++) {
      if (h1s.at("match")[i]->GetBinContent(b) < 0.05) {
        h_fit->SetBinContent(b, 0);
        h_fit->SetBinError(b, 0);
      }
    }
    // y-0, y-2: step with fixed width; y-1: step with free width; y-3: wide parabola
    const double xmin = h_sub->GetXaxis()->GetXmin(), xmax = h_sub->GetXaxis()->GetXmax();
    TF1* f = nullptr;
    if (i == 3) {
      // vertex form: maximum [0] at x = [1], reaches 0 at x = [1] +- [2]
      f = new TF1(Form("f_ratio__y-%d", i), "[0]*(1-((x-[1])/[2])^2)", xmin, xmax);
      f->SetParameters(1.1, 0.08, 0.5);
      f->FixParameter(1, 0.08);
      f->SetParLimits(2, 0.3, 5);
    } else if (i == 1) {
      f = new TF1(Form("f_ratio__y-%d", i), "[0]+[1]*0.5*TMath::Erfc((x-[2])/[3])", xmin, xmax);
      f->SetParameters(0.8, 0.8, 0.05, 0.05);
      f->SetParLimits(2, -0.15, 0.3);
      f->SetParLimits(3, 0.005, 0.5);
    } else {
      f = new TF1(Form("f_ratio__y-%d", i), "[0]+[1]*0.5*TMath::Erfc((x-[2])/0.03)", xmin, xmax);
      f->SetParameters(0.8, 0.8, 0.05);
      f->SetParLimits(2, -0.15, 0.3);
    }
    h_fit->Fit(f, "Q0N");
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    std::cout << f->GetName() << ": chi2/ndf = " << f->GetChisquare() << "/" << f->GetNDF();
    for (int p=0; p<f->GetNpar(); p++) std::cout << ", p" << p << " = " << f->GetParameter(p) << " +- " << f->GetParError(p);
    std::cout << std::endl;
    // formula with fitted values in place of [0], [1], ...
    std::string formula = xjjc::str_replaceall(f->GetTitle(), {{ "(x-", "(Dmva_BDT-" }});
    for (int p=0; p<f->GetNpar(); p++) {
      const double v = f->GetParameter(p);
      formula = xjjc::str_replaceall(formula, {{ Form("[%d]", p), Form(v < 0 ? "(%.3g)" : "%.3g", v) }});
    }
    std::cout << f->GetName() << ": " << formula << std::endl;
    fs.push_back(f);
  }
  
  leg->AddEntry(h1s["match"].front(), "Signal in MC", "f");
  leg->AddEntry(h1s["sub"].front(), "Signal in data", "pe");
  leg->AddEntry(h1s["sideband"].front(), "Sideband in data", "f");

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (int i=0; i<tbins.ny(); i++) {
    pdf->prepare();
    auto cs = xjjroot::twopads(pdf->getc(), h1s.at("sub")[i], h1s.at("ratio-sub")[i], fdiv);
    cs.front()->cd();
    h1s.at("sideband")[i]->Draw("hist same");
    h1s.at("match")[i]->Draw("hist same");
    h1s.at("sub")[i]->Draw("pe1 same");
    xjjroot::drawtexgroup(0.25, y_pos_up(), { tbins.label_y(i), tbins.label_pt() }, tsize_up, 13);
    xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)", 1./fdiv);
    leg->Draw();
    cs[1]->cd();
    xjjroot::drawline_horizon(1, h1s.at("ratio-sub")[i], h1s.at("match")[i]->GetLineColor(), 2);   
    fs[i]->Draw("same");
    h1s.at("ratio-sub")[i]->Draw("pe1 same");
    pdf->getc()->cd();
    pdf->write();
  }

  for (auto& [_, vhs] : h1s) {
    for (auto& h : vhs) {
      xjjroot::sethempty(h, 0, 0);
    }
  }
  for (int i=0; i<tbins.ny(); i++) {
    pdf->prepare();
    h1s.at("ratio-sub")[i]->Draw("pe1");
    xjjroot::drawline_horizon(1, h1s.at("ratio-sub")[i], h1s.at("match")[i]->GetLineColor(), 2);
    fs[i]->Draw("same");
    h1s.at("ratio-sub")[i]->Draw("pe1 same");
    xjjroot::drawtexgroup(0.25, 0.86, { tbins.label_y(i), tbins.label_pt(),
        Form("#chi^{2}/ndf = %.1f/%d", fs[i]->GetChisquare(), fs[i]->GetNDF()) }, 0.038, 13);
    xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
    pdf->write();
  }
  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3)
    return macro(argv[1], argv[2]);

  return 1;
}
