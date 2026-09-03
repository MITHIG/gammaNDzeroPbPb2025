#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

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

int macro(const std::string& inputname = "HIEmptyBX.root",
          const std::string& outputname = "percentile.root")
{
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf || inf->IsZombie()) {
    std::cerr << "failed to open input file: " << inputname << std::endl;
    return 1;
  }

  auto* tree = dynamic_cast<TTree*>(inf->Get("Tree"));
  if (!tree) {
    std::cerr << "failed to get TTree named Tree from " << inputname << std::endl;
    inf->Close();
    return 1;
  }

  TTreeReader reader(tree);
  TTreeReaderValue<Float_t> hfemax(reader, "HFEMaxPlus_eta5");
  TTreeReaderValue<Bool_t> isNotBptxOR(reader, "isNotBptxOR");

  std::vector<float> values;
  values.reserve(tree->GetEntries());

  Long64_t nread = 0;
  while (reader.Next()) {
    ++nread;
    if (!(*isNotBptxOR)) continue;
    values.push_back(*hfemax);
  }

  if (values.empty()) {
    std::cerr << "no events pass isNotBptxOR" << std::endl;
    inf->Close();
    return 1;
  }

  std::sort(values.begin(), values.end());

  const auto npass = values.size();
  std::vector<std::pair<double, double>> points;

  std::cout << std::setprecision(10);
  std::cout << "input: " << inputname << "\n";
  std::cout << "tree entries read: " << nread << "\n";
  std::cout << "events passing isNotBptxOR: " << npass << "\n";
  std::cout << "target_fraction target_events order_value strict_cut actual_fraction actual_events\n";

  for (int ipercent = 980; ipercent <= 995; ++ipercent) {
    const double target_fraction = ipercent / 1000.;
    const auto nkeep = static_cast<std::size_t>(std::ceil(target_fraction * npass));
    const auto index = std::min(nkeep, npass) - 1;
    const float percentile_value = values[index];
    const float strict_cut = std::nextafter(percentile_value, std::numeric_limits<float>::infinity());
    const auto actual_events = static_cast<std::size_t>(
      std::lower_bound(values.begin(), values.end(), strict_cut) - values.begin());
    const double actual_fraction = static_cast<double>(actual_events) / npass;

    points.emplace_back(strict_cut, actual_fraction);
    std::cout << target_fraction << " "
              << nkeep << " "
              << percentile_value << " "
              << strict_cut << " "
              << actual_fraction << " "
              << actual_events << "\n";
  }

  auto* outf = TFile::Open(outputname.c_str(), "recreate");
  if (!outf || outf->IsZombie()) {
    std::cerr << "failed to create output file: " << outputname << std::endl;
    inf->Close();
    return 1;
  }

  auto* graph = new TGraph(points.size());
  graph->SetName("gr_fraction_vs_HFEMaxPlus_eta5");
  graph->SetTitle(";HFEMaxPlus_eta5 cut XX;Fraction kept with HFEMaxPlus_eta5 < XX");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.9);
  graph->SetLineWidth(2);
  for (std::size_t i = 0; i < points.size(); ++i)
    graph->SetPoint(i, points[i].first, points[i].second);
  graph->Write();

  TCanvas canvas("c_fraction_vs_HFEMaxPlus_eta5", "", 700, 600);
  graph->Draw("ALP");
  canvas.SaveAs("percentile.pdf");

  outf->Close();
  std::cout << "wrote graph: " << outputname << ":gr_fraction_vs_HFEMaxPlus_eta5\n";
  std::cout << "wrote plot: percentile.pdf\n";

  inf->Close();
  return 0;
}

int main(int argc, char* argv[])
{
  const std::string inputname = argc > 1 ? argv[1] : "HIEmptyBX.root";
  const std::string outputname = argc > 2 ? argv[2] : "percentile.root";
  return macro(inputname, outputname);
}
