#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TTree.h>

struct Stats {
  double mean = 0;
  double mean_err = 0;
  double rms = 0;
  double rms_err = 0;
  double gaus_mean = 0;
  double gaus_mean_err = 0;
  double sigma = 0;
  double sigma_err = 0;
  int n = 0;
};

struct Dataset {
  std::string label;
  Color_t color;
  std::vector<std::vector<Stats>> stats;
  bool valid = false;
};

bool exists(const std::string& path) {
  std::ifstream file(path);
  return file.good();
}

Stats get_stats(const std::vector<double>& values, const std::string& name) {
  Stats result;
  result.n = values.size();
  if (values.empty()) return result;
  double sum = 0, sum2 = 0;
  for (const auto value : values) {
    sum += value;
    sum2 += value * value;
  }
  result.mean = sum / values.size();
  result.rms = values.size() > 1
    ? std::sqrt(std::max(0., (sum2 - values.size() * result.mean * result.mean) /
                            (values.size() - 1)))
    : 0;
  result.mean_err = values.size() > 0 ? result.rms / std::sqrt(values.size()) : 0;
  result.rms_err = values.size() > 1 ? result.rms / std::sqrt(2. * (values.size() - 1)) : 0;

  TH1D h((name + "_h").c_str(), "", 80, -5, 5);
  for (const auto value : values) h.Fill(value);
  TF1 gaus((name + "_gaus").c_str(), "gaus", -2.5, 2.5);
  gaus.SetParameters(h.GetMaximum(), result.mean, result.rms);
  h.Fit(&gaus, "Q0R");
  result.gaus_mean = gaus.GetParameter(1);
  result.gaus_mean_err = gaus.GetParError(1);
  result.sigma = std::abs(gaus.GetParameter(2));
  result.sigma_err = gaus.GetParError(2);
  if (!std::isfinite(result.sigma) || result.sigma <= 0) result.sigma = result.rms;
  return result;
}

Dataset read_dataset(const std::string& path, bool closure, int npt, int ny) {
  Dataset result;
  result.label = closure ? "Fitted-model closure" : "Poisson bootstrap";
  result.color = closure ? kRed + 1 : kBlue + 1;
  if (!exists(path)) return result;
  auto* file = TFile::Open(path.c_str(), "READ");
  if (!file || file->IsZombie()) return result;
  auto* tree = static_cast<TTree*>(file->Get("toys"));
  if (!tree) return result;
  if (closure && !tree->GetBranch("truth_yield")) return result;
  if (!closure && !tree->GetBranch("nominal_yield")) return result;

  int pt = 0, y = 0;
  double yield = 0, yield_err = 0, reference = 0, reference_err = 0;
  int fit_ok = 1;
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("yield", &yield);
  tree->SetBranchAddress("yield_err", &yield_err);
  if (closure) {
    tree->SetBranchAddress("truth_yield", &reference);
    if (tree->GetBranch("fit_ok")) tree->SetBranchAddress("fit_ok", &fit_ok);
  } else {
    tree->SetBranchAddress("nominal_yield", &reference);
    tree->SetBranchAddress("nominal_yield_err", &reference_err);
  }

  std::vector<std::vector<std::vector<double>>> pulls(
      npt, std::vector<std::vector<double>>(ny));
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    if (pt < 0 || pt >= npt || y < 0 || y >= ny || !fit_ok || yield < 0) continue;
    const double denominator = closure ? yield_err : reference_err;
    if (denominator <= 0) continue;
    pulls[pt][y].push_back((yield - reference) / denominator);
  }
  result.stats.assign(npt, std::vector<Stats>(ny));
  for (int ipt = 0; ipt < npt; ++ipt)
    for (int iy = 0; iy < ny; ++iy)
      result.stats[ipt][iy] = get_stats(pulls[ipt][iy],
                                        Form("pull_%d_%d_%d", closure, ipt, iy));
  result.valid = true;
  return result;
}

int main(int argc, char* argv[]) {
  if (argc < 3 || argc > 5) return 1;
  gROOT->SetBatch(true);

  const std::string toy_file = argv[1];
  const std::string closure_file = argc >= 4 ? argv[2] : "";
  const std::string output = argc >= 4 ? argv[3] : argv[2];
  std::string fit_output = argc == 5 ? argv[4] : output;
  if (argc != 5) {
    const auto dot = fit_output.rfind('.');
    fit_output = (dot == std::string::npos ? fit_output : fit_output.substr(0, dot)) + ".root";
  }

  std::string binning_file = exists(toy_file) ? toy_file : closure_file;
  if (!exists(binning_file)) return 2;
  auto* binfile = TFile::Open(binning_file.c_str(), "READ");
  auto* h3_bins = binfile ? static_cast<TH3D*>(binfile->Get("h3_bins")) : nullptr;
  if (!h3_bins) return 3;
  const int npt = h3_bins->GetZaxis()->GetNbins();
  const int ny = h3_bins->GetXaxis()->GetNbins();

  auto bootstrap = read_dataset(toy_file, false, npt, ny);
  auto closure = read_dataset(closure_file, true, npt, ny);
  if (!bootstrap.valid && !closure.valid) return 4;

  auto* fitfile = TFile::Open(fit_output.c_str(), "RECREATE");
  auto* fit_tree = new TTree("pull_fit", "Pull-distribution fit results");
  std::string dataset;
  int pt_index = -1, y_index = -1, n = 0;
  double mean = 0, mean_err = 0, rms = 0, rms_err = 0;
  double gaus_mean = 0, gaus_mean_err = 0, gaus_sigma = 0, gaus_sigma_err = 0;
  fit_tree->Branch("dataset", &dataset);
  fit_tree->Branch("pt", &pt_index);
  fit_tree->Branch("y", &y_index);
  fit_tree->Branch("n", &n);
  fit_tree->Branch("mean", &mean);
  fit_tree->Branch("mean_err", &mean_err);
  fit_tree->Branch("rms", &rms);
  fit_tree->Branch("rms_err", &rms_err);
  fit_tree->Branch("gaus_mean", &gaus_mean);
  fit_tree->Branch("gaus_mean_err", &gaus_mean_err);
  fit_tree->Branch("gaus_sigma", &gaus_sigma);
  fit_tree->Branch("gaus_sigma_err", &gaus_sigma_err);
  auto save_stats = [&](const Dataset& data) {
    if (!data.valid) return;
    dataset = data.label;
    for (pt_index = 0; pt_index < npt; ++pt_index) {
      for (y_index = 0; y_index < ny; ++y_index) {
        const auto& s = data.stats[pt_index][y_index];
        n = s.n; mean = s.mean; mean_err = s.mean_err; rms = s.rms; rms_err = s.rms_err;
        gaus_mean = s.gaus_mean; gaus_mean_err = s.gaus_mean_err;
        gaus_sigma = s.sigma; gaus_sigma_err = s.sigma_err;
        fit_tree->Fill();
      }
    }
  };
  save_stats(bootstrap);
  save_stats(closure);

  const auto* yaxis = h3_bins->GetXaxis();
  std::vector<double> yedges(ny + 1);
  for (int i = 0; i < ny; ++i) yedges[i] = yaxis->GetBinLowEdge(i + 1);
  yedges[ny] = yaxis->GetBinUpEdge(ny);

  TCanvas canvas("toy_pull", "toy_pull", 1000, 700);
  canvas.Print((output + "[").c_str());
  for (int ipt = 0; ipt < npt; ++ipt) {
    canvas.Clear();
    TH2D frame(Form("pull_frame_%d", ipt),
               Form("Pull distribution, p_{T} bin %d; y; pull", ipt),
               ny, yedges.data(), 1, -3.5, 3.5);
    frame.SetStats(false);
    frame.Draw();

    TGraphErrors bootstrap_gaus, closure_gaus;
    bootstrap_gaus.SetMarkerStyle(20);
    bootstrap_gaus.SetMarkerColor(kBlue + 1);
    bootstrap_gaus.SetLineColor(kBlue + 1);
    closure_gaus.SetMarkerStyle(20);
    closure_gaus.SetMarkerColor(kRed + 1);
    closure_gaus.SetLineColor(kRed + 1);

    for (int iy = 0; iy < ny; ++iy) {
      const double xlow = yedges[iy], xhigh = yedges[iy + 1];
      const double xcenter = 0.5 * (xlow + xhigh);
      const double box_width = 0.16 * (xhigh - xlow);
      const double xshift = box_width;
      if (bootstrap.valid) {
        const auto& s = bootstrap.stats[ipt][iy];
        const int point = bootstrap_gaus.GetN();
        const double x = xcenter - xshift;
        bootstrap_gaus.SetPoint(point, x, s.gaus_mean);
        bootstrap_gaus.SetPointError(point, 0, s.sigma);
        auto* box = new TBox(x - box_width, s.mean - s.rms,
                             x + box_width, s.mean + s.rms);
        box->SetFillColorAlpha(kBlue + 1, 0.20);
        box->SetLineColor(kBlue + 1);
        box->SetLineStyle(2);
        box->Draw();
      }
      if (closure.valid) {
        const auto& s = closure.stats[ipt][iy];
        const int point = closure_gaus.GetN();
        const double x = xcenter + xshift;
        closure_gaus.SetPoint(point, x, s.gaus_mean);
        closure_gaus.SetPointError(point, 0, s.sigma);
        auto* box = new TBox(x - box_width, s.mean - s.rms,
                             x + box_width, s.mean + s.rms);
        box->SetFillColorAlpha(kRed + 1, 0.20);
        box->SetLineColor(kRed + 1);
        box->SetLineStyle(2);
        box->Draw();
      }
    }

    if (bootstrap.valid) bootstrap_gaus.Draw("PZ same");
    if (closure.valid) closure_gaus.Draw("PZ same");
    TLine zero(yedges.front(), 0, yedges.back(), 0);
    zero.SetLineStyle(2);
    zero.Draw();
    TLine plus(yedges.front(), 1, yedges.back(), 1);
    plus.SetLineStyle(3);
    plus.Draw();
    TLine minus(yedges.front(), -1, yedges.back(), -1);
    minus.SetLineStyle(3);
    minus.Draw();

    auto* leg = new TLegend(0.62, 0.72, 0.90, 0.90);
    if (bootstrap.valid) leg->AddEntry(&bootstrap_gaus, "Bootstrap Gaussian mean #pm #sigma", "pe");
    if (closure.valid) leg->AddEntry(&closure_gaus, "Closure Gaussian mean #pm #sigma", "pe");
    leg->AddEntry((TObject*)nullptr, "Dashed boxes: sample mean #pm RMS", "");
    leg->Draw();
    canvas.Print(output.c_str());
  }
  canvas.Print((output + "]").c_str());
  fit_tree->Write();
  fitfile->Close();
  return 0;
}
