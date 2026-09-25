#include "../../include/xjjanauti.h"
#include "../../include/xjjmypdf.h"
#include "../../include/util.h"
#include "../../include/draw.h"

int macro(const std::string& inputfile, const std::string& outputname) {
  auto* inf = TFile::Open(inputfile.c_str());
  if (xjjroot::failfile(inf)) return 2;

  const auto info = util::read_info(inf);
  
  const std::vector<std::string> labels = { "Plus", "Minus" },
    texs = { "HF+", "HF-" };
  const std::vector<Color_t> colors = { xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"] };

  auto* leg = new TLegend(0.24, 0.85-0.038*1.2*2, 0.24+0.2, 0.85);
  xjjroot::setleg(leg, 0.038);
  std::vector<TGraph*> grs;
  double xmin = 1.e3, xmax = -1;
  for (int i=0; i<labels.size(); i++) {
    auto* gr = xjjana::getobj<TGraph>(inf, "gr_rgap_percent_HFEMax" + labels[i]);
    xjjroot::setthgrstyle(gr, colors[i], xjjroot::mstylelist_solid(i), 1.5, colors[i], 1, 1);
    leg->AddEntry(gr, texs[i].c_str(), "p");
    xmin = std::min(xmin, xjjana::gethminimumX(gr));
    xmax = std::max(xmax, xjjana::gethmaximumX(gr));
    grs.push_back(gr);
  }

  auto* hempty = new TH2D("hempty", Form(";Gap Threshold in HF %s;Fraction kept of EmpyBX events", info.at("hfvar_tex").c_str()),
                          10, xmin-1, xmax+1, 10, 0.97, 1);
  xjjroot::sethempty(hempty, 0, 0.2);

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto name_png = draw::png_name(pdf);
  pdf->prepare();
  hempty->Draw("axis");
  for (auto g : grs) g->Draw("pl same");
  leg->Draw();
  xjjroot::drawCMS(xjjroot::CMS::internal, info.at("input_tex") + " (PbPb 5.36 TeV)");
  pdf->write(name_png + ".pdf");
  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) return macro(argv[1], argv[2]);
  return 1;
}
