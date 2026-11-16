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
    // h1s.at("ratio-sub")[i]->Draw("pe1 same");
    xjjroot::drawtexgroup(0.25, 0.86, { tbins.label_y(i), tbins.label_pt() }, 0.038, 13);
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
