#include "xjjanauti.h"
#include "xjjstruct.h"
#include "draw.h"
#include "util.h"

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
  TH1D* h = nullptr;
  TF1* gaus = nullptr;
};

struct Dataset {
  int toytype;
  std::string label;
  std::vector<std::vector<Stats>> stats;
  bool valid = false;
  int ntoys;
  // xjjc::info info;
};

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

  auto* h = new TH1D(("h1_pull_" + name).c_str(), ";Pull;", 80, -5, 5);
  xjjroot::sethempty(h);
  for (const auto value : values) h->Fill(value);
  auto* gaus = new TF1(("f1_gaus_" + name).c_str(), "gaus", -2.5, 2.5);
  gaus->SetParameters(h->GetMaximum(), result.mean, result.rms);
  h->Fit(gaus, "Q0R");
  result.gaus_mean = gaus->GetParameter(1);
  result.gaus_mean_err = gaus->GetParError(1);
  result.sigma = std::abs(gaus->GetParameter(2));
  result.sigma_err = gaus->GetParError(2);
  if (!std::isfinite(result.sigma) || result.sigma <= 0) result.sigma = result.rms;
  result.h = h;
  result.gaus = gaus;
  
  return result;
}

Dataset read_dataset(const std::string& path, const int npt, const int ny) {
  Dataset result;
  auto* inf = xjjroot::readfile(path);
  if (!inf) return result;
  auto* tree = static_cast<TTree*>(inf->Get("toys"));
  if (!tree) return result;
  const auto info = util::read_info(inf, "info");
  result.toytype = std::atoi(info.at("toytype").c_str());
  result.label = info.at("toytype_tex");
  result.ntoys = std::atoi(info.at("ntoys").c_str());
  
  // int pt = 0, y = 0;
  // double yield = 0, yield_err = 0, reference = 0, reference_err = 0;
  // int fit_ok = 1;
  int pt = 0; tree->SetBranchAddress("ipt", &pt);
  int y = 0; tree->SetBranchAddress("iy", &y);
  double yield = 0; tree->SetBranchAddress("yield", &yield);
  double yield_err = 0; tree->SetBranchAddress("yield_err", &yield_err);
  double reference = 0; tree->SetBranchAddress("nominal_yield", &reference);
  double reference_err = 0; tree->SetBranchAddress("nominal_yield_err", &reference_err);

  auto pulls = xjjc::array2d<std::vector<double>>(npt, ny);
  for (Long64_t ientry = 0; ientry < tree->GetEntries(); ++ientry) {
    tree->GetEntry(ientry);
    if (pt < 0 || pt >= npt || y < 0 || y >= ny || yield < 0) continue;
    const double denominator = result.toytype == 1 ? yield_err : reference_err;
    if (denominator <= 0) continue;
    pulls[pt][y].push_back((yield - reference) / denominator);
  }
  // result.stats.assign(npt, std::vector<Stats>(ny));
  result.stats = xjjc::array2d<Stats>(npt, ny);
  for (int ipt = 0; ipt < npt; ++ipt)
    for (int iy = 0; iy < ny; ++iy)
      result.stats[ipt][iy] = get_stats(pulls[ipt][iy],
                                        Form("pull_%d_%d_%d", result.toytype, ipt, iy));
  result.valid = true;
  return result;
}

int macro(const std::vector<std::string>& inputnames, const std::string& outputname,
          int draw_rms = 0) {
  
  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", inputnames.front().c_str()));
  if (!h3_bins) return 2;
  const draw::bintex tbins(h3_bins, 0, 2);
  const auto npt = tbins.npt(), ny = tbins.ny();
  auto* hempty = new TH2D("hempty", ";y;Pull Mean#scale[0.4]{ }#pm#scale[0.4]{ }#sigma", 10, tbins.edgelow_y(), tbins.edgeup_y(), 10, -2.5, 3.5);
  xjjroot::sethempty(hempty, 0, 0.1);

  std::vector<Dataset> datasets;
  for (const auto& input : inputnames) {
    auto data = read_dataset(input, npt, ny);
    if (!data.valid) continue;
    datasets.push_back(data);
  }

  const float tsize = 0.04, lspace = 1.2, y_top = 0.85;
  auto y_pos = [&tsize, &lspace, &y_top](float i = 0) -> float {
    return y_top - i*lspace*tsize;
  };
  xjjroot::setgstyle(1);
  gStyle->SetEndErrorSize(8);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (int ipt = 0; ipt < npt; ++ipt) {    

    pdf->prepare();
    hempty->Draw("axis");
    if (!draw_rms)
      xjjroot::drawbox(hempty->GetXaxis()->GetXmin(), -1, hempty->GetXaxis()->GetXmax(), 1, kBlack, 0.03);
    for (const auto p : std::vector<std::pair<float, Style_t>>{ {0, 2}, {-1, 3}, {1, 3} })
      xjjroot::drawline_horizon(p.first, hempty, kBlack, p.second);
    const auto ndata = datasets.size();
    auto* leg = new TLegend(0.6, y_pos(ndata), 0.85, y_pos());
    xjjroot::setleg(leg, tsize);

    for (int idata = 0; idata < ndata; idata++) {
      const auto data = datasets[idata];
      const auto itype = data.toytype;

      auto* gr_rms = new TGraphErrors(ny);
      gr_rms->SetName(Form("gr_rms_%d__pt-%d", data.toytype, ipt));
      const auto color = draw::colors[data.toytype];
      xjjroot::setthgrstyle(gr_rms, color, xjjroot::mstylelist_solid(itype), 1.5,
                            -1, -1, -1,
                            color, 0.1, 1001);
      auto* gr_gaus = new TGraphErrors(ny);
      gr_gaus->SetName(Form("gr_gaus_%d__pt-%d", data.toytype, ipt));
      xjjroot::setthgrstyle(gr_gaus, color, xjjroot::mstylelist_solid(itype), 1.5,
                            color, 1, 1);

      leg->AddEntry(gr_gaus, data.label.c_str(), "pe");
    
      for (int iy = 0; iy < ny; ++iy) {
        const auto& s = data.stats[ipt][iy]; // 
        const auto y_width = std::min(0.15, tbins.binwidth_y(iy)*1./(ndata+1)),
          y_center = tbins.edgelow_y(iy) + (tbins.binwidth_y(iy)-y_width*ndata)/2. + (idata+0.5)*y_width;
        // __XJJLOG << "y_width = " << y_width << ", ycenter = " << y_center << std::endl;

        gr_gaus->SetPoint(iy, y_center, s.gaus_mean);
        gr_gaus->SetPointError(iy, 0, s.sigma);
        gr_rms->SetPoint(iy, y_center, s.mean);
        gr_rms->SetPointError(iy, y_width/2., s.rms);
      }
      if (draw_rms)
      gr_rms->Draw("2 same");
      gr_gaus->Draw("pe1 same");
    } // for (int idata = 0; idata < ndata; idata++) {
    leg->Draw();
    xjjroot::drawtexgroup(0.25, y_pos(0.1), {
        tbins.label_pt(ipt),
        Form("N_{toys} = %d", datasets.front().ntoys),
      }, tsize, 13, 42, lspace);
    xjjroot::drawCMS(xjjroot::CMS::internal, "PbPb (5.36 TeV)");
    pdf->write(); 
  }
  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4)
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], std::atoi(argv[3]));
  if (argc == 3)
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2]);

  return 1;
}
