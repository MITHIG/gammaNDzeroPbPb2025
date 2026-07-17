#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "fpfitter.h"

int macro(const std::string& inputname, const std::string& outputname = "") {

  auto* inf = TFile::Open(inputname.c_str());

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (const std::string var : { "Dip3D", "Dip3Dsig" }) {
    std::map<std::string, std::vector<TH1D*>> h1ys;
    auto* dir = inf->GetDirectory(Form("dir_%s", var.c_str()));
    for (const std::string type : { "mc-prompt", "mc-nonprompt",
                                   "data-sub", "data-sigswap" }) {
      h1ys[type] = xjjana::getobj_regexp<TH1D>(dir, ".+_" + var + ".*_" + type + ".*__y-.+");
    }
    auto fit_type = [&h1ys, &pdf, &var](std::string type, std::string title) {
      pdf->draw_cover({ "#bf{Variable} " + var, "#bf{Signal extraction} " + title });
      for (int i=0; i<h1ys.at("data-" + type).size(); i++) {
        auto* fitter = new fpfitter(h1ys.at("data-"+type)[i],
                                    h1ys.at("mc-prompt")[i],
                                    h1ys.at("mc-nonprompt")[i], type);
        if (fitter->status() > 0)
          continue;
        fitter->fit();
        fitter->print_fitresult();
        pdf->prepare();
        auto pads = fitter->draw(pdf->getc());
        pads.front()->cd();
        xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)", 1./fpfitter::pratio);
        xjjroot::drawtex(0.46, 0.835, title.c_str(), fpfitter::tsize/fpfitter::pratio, 33);
        pdf->write();
      }
    };
    fit_type("sub", "Sideband sub");
    fit_type("sigswap", "RooFit sPlot");
  }

  pdf->close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 2) return macro(argv[1]);
  if (argc == 3) return macro(argv[1], argv[2]);
  return 1;
}
