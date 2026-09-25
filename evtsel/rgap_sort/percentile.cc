#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "../../include/xjjanauti.h"
#include "../../include/xjjstruct.h"
#include "../../include/xjjmypdf.h"
#include "../../include/util.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

int macro(const std::string& inputnamestr, const std::string& hfvarstr, const std::string& outputname,
          double pmin = 0.975, double pmax = 0.995, double pstep = 0.001)
{
  const auto pinput = xjjroot::parse_input(inputnamestr);
  const auto phfvar = xjjroot::parse_input(hfvarstr);
  
  auto* inf = TFile::Open(pinput.content.c_str());
  if (xjjroot::failfile(inf)) return 2;

  auto* tree = dynamic_cast<TTree*>(inf->Get("Tree"));
  if (!tree) {
    std::cerr << "!! failed to get Tree, abort." << std::endl;
    xjjroot::closefile(inf);
    return 1;
  }
  auto t_HFEMaxPlus = xjjc::str_replaceall(phfvar.content, "HFEMax", "HFEMaxPlus"),
    t_HFEMaxMinus = xjjc::str_replaceall(phfvar.content, "HFEMax", "HFEMaxMinus");

  auto set_branch = [&tree]<typename T>(const std::string& t_br, T& br, bool required = true) -> bool {
    if (tree->GetBranch(t_br.c_str())) {
      tree->SetBranchStatus(t_br.c_str(), 1);
      tree->SetBranchAddress(t_br.c_str(), &br);
    } else {
      __XJJLOG << "!! no branch: " << t_br << std::endl;
      if (required) return false;
    }
    return true;
  };
  float HFEMaxPlus; if (!set_branch(t_HFEMaxPlus, HFEMaxPlus)) return 2; 
  float HFEMaxMinus; if (!set_branch(t_HFEMaxMinus, HFEMaxMinus)) return 2; 

  bool isNotBptxOR; if (!set_branch("isNotBptxOR", isNotBptxOR)) return 2; 

  std::vector<float> l_HFEMaxPlus, l_HFEMaxMinus;
  l_HFEMaxPlus.reserve(tree->GetEntries("isNotBptxOR"));
  l_HFEMaxMinus.reserve(tree->GetEntries("isNotBptxOR"));
  
  auto nentries = tree->GetEntries();
  for (long long int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries);
    tree->GetEntry(i);

    if (!isNotBptxOR) continue;
    l_HFEMaxPlus.push_back(HFEMaxPlus);
    l_HFEMaxMinus.push_back(HFEMaxMinus);
  }
  xjjc::progressbar_summary(nentries);
  
  if (l_HFEMaxPlus.empty() || l_HFEMaxMinus.empty()) {
    __XJJLOG << "!! no events pass isNotBptxOR, abort." << std::endl;
    xjjroot::closefile(inf);
    return 2;
  }

  std::sort(l_HFEMaxPlus.begin(), l_HFEMaxPlus.end());
  std::sort(l_HFEMaxMinus.begin(), l_HFEMaxMinus.end());

  auto percent = [&pmin, &pmax, &pstep](const std::vector<float>& values, const std::string& name) -> TGraph* {
    TGraph* gr = nullptr;
    
    std::vector<std::pair<double, double>> points;
    const auto npass = values.size();
    for (double target_fraction = pmin; target_fraction <= pmax; target_fraction += pstep) {
      // const double target_fraction = ipercent / 1000.;
      const auto nkeep = static_cast<std::size_t>(std::ceil(target_fraction * npass));
      const auto index = std::min(nkeep, npass) - 1;
      const float percentile_value = values[index];
      const float strict_cut = std::nextafter(percentile_value, std::numeric_limits<float>::infinity());
      const auto actual_events = static_cast<std::size_t>(std::lower_bound(values.begin(), values.end(), strict_cut) - values.begin());
      const double actual_fraction = static_cast<double>(actual_events) / npass;

      points.emplace_back(actual_fraction, strict_cut);
      // std::cout << target_fraction << " "
      //           << nkeep << " "
      //           << percentile_value << " "
      //           << strict_cut << " "
      //           << actual_fraction << " "
      //           << actual_events << "\n";
    }
    xjjc::print_tab(points, 0);
    gr = new TGraph(points.size());
    gr->SetName(Form("gr_rgap_percent_%s", name.c_str()));
    for (std::size_t i = 0; i < points.size(); ++i)
      gr->SetPoint(i, points[i].second, points[i].first);
    xjjroot::setthgrstyle(gr, kBlack, 20, 1.4, kBlack, 1, 1);

    return gr;
  };

  auto* gr_HFEMaxPlus = percent(l_HFEMaxPlus, "HFEMaxPlus");
  auto* gr_HFEMaxMinus = percent(l_HFEMaxMinus, "HFEMaxMinus");

  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  pdf->prepare();
  gr_HFEMaxPlus->Draw("ALP");
  pdf->write();
  pdf->prepare();
  gr_HFEMaxMinus->Draw("ALP");
  pdf->write();

  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(gr_HFEMaxPlus);
  xjjroot::writehist(gr_HFEMaxMinus);
  util::Writeinfo tinfo("info");
  tinfo.cast_branch("input", pinput.content);
  tinfo.cast_branch("input_tex", pinput.tex);
  tinfo.cast_branch("hfvar", phfvar.content);
  tinfo.cast_branch("hfvar_tex", phfvar.tex);
  tinfo.close();
  xjjroot::closefile(outf);

  inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4)
    return macro(argv[1], argv[2], argv[3]);

  return 1;
}
