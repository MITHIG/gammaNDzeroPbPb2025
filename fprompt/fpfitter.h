#pragma once

#include "TF1.h"
#include "TFitResult.h"
#include "TRandom.h"

#define __COOK_NAME__
#include "style.h"

#include "xjjanauti.h"

enum class Qual { pull = 0, ratio = 1 };
std::vector<double> yref = { 0, 1 };
std::vector<std::string> drawopt = { "hist", "pe1" };

class fpfitter {
public:
  fpfitter(TH1D* hdata, TH1D* hprompt, TH1D* hnonprompt, std::string name_data = "", std::string name_mc = "");
  fpfitter(TH1D* hdata, TH1D* htotal_fitted, TH1D* hnonprompt_fitted, TFitResult* fresult, TH1D* hpull = nullptr, TH1D* hratio = nullptr);
  void fit(double xmin = 0, double xmax = 0);  
  std::vector<TPad*> draw(TPad* c, Qual q = Qual::pull);
  void print_fitresult();
  int status() const { return status_; }
  double fit_xmin() const { return hdata_->GetXaxis()->GetBinLowEdge(ibin_fit_min_); }
  double fit_xmax() const { return hdata_->GetXaxis()->GetBinLowEdge(ibin_fit_max_) + hdata_->GetXaxis()->GetBinWidth(ibin_fit_max_); }
  double fprompt() const { return fitresult_->Parameter(0); }
  double fprompt_err_par() const { return fitresult_->ParError(0); }
  double fprompt_err_low() const { return fprompt_err_low_; }
  double fprompt_err_high() const { return fprompt_err_high_; }
  int ndf() const { return (ibin_fit_max_ - ibin_fit_min_ + 1)/*nbins fitted*/ - 1/*npars*/; }
  double chi2() const { return chi2_; }

  void write_to_file();
  
  // style
  static constexpr float pratio = 2./3, tsize = 0.036,
    xleft = 0.48, ytop = 0.84, ybottom = ytop - tsize/pratio*1.1*4 - 0.01 - tsize/pratio*1.2*2;
  
private:
  TH1D* hdata_, *hprompt_, *hnonprompt_;
  const std::string name_data_, name_mc_;
  int status_;
  double area_hprompt_, area_hnonprompt_;
  TF1 *func_, *func_smear_;
  TFitResultPtr fitresult_;
  int ibin_fit_min_, ibin_fit_max_, nbin_fit_;
  double chi2_;
  TH1D *hprompt_fitted_, *hnonprompt_fitted_, *htotal_fitted_, *hpull_, *hratio_;
  double fprompt_err_low_, fprompt_err_high_;

  double template_model(double* x, double* par) {
    const double fpr = par[0];
    const int bin = hprompt_->FindFixBin(x[0]);
    if (bin < 1 || bin > hprompt_->GetNbinsX()) return 0;
    const double p = hprompt_->GetBinContent(bin) / area_hprompt_;
    const double n = hnonprompt_->GetBinContent(bin) / area_hnonprompt_;
    return fpr * p + (1. - fpr) * n;
  }
  void calc_pull_chi2(int ibin_min, int ibin_max, bool fill_hists = true);
  static bool same_axis(const TH1D* a, const TH1D* b);
  TH1D* clone_h1(TH1D* h, const std::string& suffix);
  void random_smear(TH1D* h0, TH1D* h) {
    for(int i = 0; i < h0->GetXaxis()->GetNbins(); i++) {
      h->SetBinContent(i+1, gRandom->Gaus(h0->GetBinContent(i+1), h0->GetBinError(i+1)));
    }
  }
};

fpfitter::fpfitter(TH1D* hdata, TH1D* hprompt, TH1D* hnonprompt,
                   std::string name_data, std::string name_mc)
  : hdata_(hdata), hprompt_(hprompt), hnonprompt_(hnonprompt),
    name_data_(name_data.empty() ? xjjc::unique_str() : name_data),
    name_mc_(name_mc.empty() ? xjjc::unique_str() : name_mc),
    status_(-1), area_hprompt_(0), area_hnonprompt_(0),
    func_(nullptr), func_smear_(nullptr), fitresult_(nullptr), chi2_(0),
    hprompt_fitted_(nullptr), hnonprompt_fitted_(nullptr), htotal_fitted_(nullptr),
    fprompt_err_low_(-1), fprompt_err_high_(-1)
{
  for (auto* h : { hdata, hprompt, hnonprompt } ) {
    if (!h) {
      __XJJLOG << "!! bad histogram, abort." << std::endl;
      status_ = 2;
      return;
    }
    if (h->Integral("width") <= 0) {
      __XJJLOG << "!! bad histogram area, abort." << std::endl;
      status_ = 2;
      return;
    }
  }
  if (!same_axis(hdata_, hprompt_) || !same_axis(hdata_, hnonprompt_)) {
    __XJJLOG << "!! histograms do not have matching binning" << std::endl;
    status_ = 2;
    return;
  }

  // !! how to check if the distribution has been scaled by widths?
  // used in fitting functions
  area_hprompt_ = hprompt_->Integral("width");
  area_hnonprompt_ = hnonprompt_->Integral("width");

  ibin_fit_min_ = 1;
  ibin_fit_max_ = hdata_->GetXaxis()->GetNbins();
  nbin_fit_ = 0;
  for (int i=0; i<hdata_->GetXaxis()->GetNbins(); i++) {
    // if (hdata_->GetBinError(i+1) > 0)
    nbin_fit_++;
  }
  func_ = new TF1(Form("func_%s", add_suffix__y(xjjc::str_eraseall(hdata_->GetName(), std::vector<std::string>{ "h1_", "norm_" }), name_mc_).c_str()),
                  [this](double* x, double* par) { return template_model(x, par); },
                  hdata_->GetXaxis()->GetXmin(), hdata_->GetXaxis()->GetXmax(), 1);
  func_->SetParName(0, "fprompt");
  func_->SetNpx(1000);
  func_smear_ = new TF1("func_smear",
                        [this](double* x, double* par) { return template_model(x, par); },
                        hdata_->GetXaxis()->GetXmin(), hdata_->GetXaxis()->GetXmax(), 1);
}

fpfitter::fpfitter(TH1D* hdata, TH1D* htotal_fitted, TH1D* hnonprompt_fitted, TFitResult* fresult, TH1D* hpull, TH1D* hratio) :
  hdata_(hdata), htotal_fitted_(htotal_fitted), hnonprompt_fitted_(hnonprompt_fitted),
  hpull_(hpull), hratio_(hratio), fitresult_(fresult),
  status_(-1), area_hprompt_(0), area_hnonprompt_(0),
  fprompt_err_low_(-1), fprompt_err_high_(-1)
{
  bool fill_hpull = !hpull_ || !hratio_;
  ibin_fit_min_ = 1;
  ibin_fit_max_ = hdata_->GetXaxis()->GetNbins(); // !! need to find a way to tell the fitting range
  calc_pull_chi2(ibin_fit_min_, ibin_fit_max_, fill_hpull);
}

void fpfitter::fit(double xmin, double xmax) {
  if (!xjjc::almost_eq(xmin, xmax)) {
    ibin_fit_min_ = 0;
    ibin_fit_max_ = 0;
    nbin_fit_ = 0;
    for (int i=0; i<hdata_->GetXaxis()->GetNbins(); i++) {
      if (xjjc::almost_eq(hdata_->GetBinLowEdge(i+1), xmin)) {
        ibin_fit_min_ = i+1;
      }
      if (xjjc::almost_eq(hdata_->GetBinLowEdge(i+1) + hdata_->GetBinWidth(i+1), xmax)) {
        ibin_fit_max_ = i+1;
      }
    }
  }
  // __XJJLOG << ">> fitting range: [ " << fit_xmin() << " - " << fit_xmax() << " ]" << std::endl;
  if (fit_xmin() >= fit_xmax()) {
    __XJJLOG << "!! bad fitting range, abort." << std::endl;
    return;
  }
  
  func_->SetParameter(0, 0.9);
  func_->SetParLimits(0, 0., 1.);
  // "I" uses the function integral over each bin, which is the correct
  // comparison for variable-width, bin-width-normalized histograms.
  // https://root.cern.ch/doc/v632/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
  fitresult_ = hdata_->Fit(func_, "S R I Q N", "", fit_xmin(), fit_xmax()); // S: TFitResultPtr behaves as a smart pointer to the TFitResult object
  fitresult_->SetName(xjjc::str_replaceall(func_->GetName(), "func_", "fresult_").c_str());
  status_ = int(fitresult_);
  if (int(fitresult_) > 0) {
    __XJJLOG << "++ bad fitting for " << func_->GetName() << std::endl;
    __XJJLOG << ">> fitting result code : " << int(fitresult_) << std::endl;
  }
  
  hprompt_fitted_ = clone_h1(hprompt_, "_fitted" + name_data_);
  hprompt_fitted_->Scale(fprompt() / area_hprompt_);
  hnonprompt_fitted_ = clone_h1(hnonprompt_, "_fitted" + name_data_);
  hnonprompt_fitted_->Scale((1-fprompt()) / area_hnonprompt_);
  htotal_fitted_ = clone_h1(hprompt_fitted_, "_total");
  htotal_fitted_->Add(hnonprompt_fitted_);

  auto* hdata_smear = (TH1D*)hdata_->Clone("hdata_smear");
  std::vector<double> toys;
  const int nSmearData = 1000;
  for (int j=0; j<nSmearData; j++) {
    random_smear(hdata_, hdata_smear);
    func_smear_->SetParameter(0, 0.9);
    func_smear_->SetParLimits(0, 0., 1.);
    auto r = hdata_smear->Fit(func_smear_, "SRIQN", "", fit_xmin(), fit_xmax());
    if (r.Get() && r->Status() == 0)
      toys.push_back(func_smear_->GetParameter(0));
  }
  delete hdata_smear;

  std::sort(toys.begin(), toys.end());
  auto quantile = [&toys](double q) {
    const double x = q * (toys.size() - 1);
    const auto i = static_cast<size_t>(x);
    const auto a = x - i;
    if (i + 1 >= toys.size()) return toys.back();
    return toys[i] * (1 - a) + toys[i + 1] * a;
  };
  if (toys.size() < 100) {
    __XJJLOG << "!! too few successful toys" << std::endl;
  } else {
    const double q16 = quantile(0.1587);
    const double q84 = quantile(0.8413);
    fprompt_err_low_  = fprompt() - q16;
    fprompt_err_high_ = q84 - fprompt();
  }

  calc_pull_chi2(ibin_fit_min_, ibin_fit_max_, true /*fill hpull_ and hratio_*/);
}

void fpfitter::calc_pull_chi2(int ibin_min, int ibin_max, bool fill_hists) {
  if (fill_hists) {
    hpull_ = clone_h1(hdata_, "_pull" + name_mc_);
    hpull_->Reset("ICESM");
    hpull_->GetYaxis()->SetTitle("Pull");
    hratio_ = clone_h1(hdata_, "_ratio" + name_mc_);
    hratio_->Reset("ICESM");
    hratio_->GetYaxis()->SetTitle("Data / Fit");
  }
  
  chi2_ = 0;
  for (int i = ibin_min; i <= ibin_max; ++i) {
    const double d = hdata_->GetBinContent(i);
    double e = hdata_->GetBinError(i);
    bool emptybin = e <= 0;
    if (e <= 0) e = std::sqrt(std::max(std::abs(d), 1.));
    const double ftotal = htotal_fitted_->GetBinContent(i); 
    const double r = (d - ftotal) / e;
    chi2_ += r * r;

    if (fill_hists) {
      hpull_->SetBinContent(i, emptybin ? 0 : r);
      hpull_->SetBinError(i, 1.e-5);
      hratio_->SetBinContent(i, ftotal==0 ? 0 : d/ftotal);
      hratio_->SetBinError(i, ftotal==0 ? 0 : hdata_->GetBinError(i)/ftotal);
    }
  }
}

std::vector<TPad*> fpfitter::draw(TPad* c, Qual q) {
  xjjroot::setthgrstyle(hdata_, kBlack, 20, 1.7, kBlack, 1, 1);
  xjjroot::setthgrstyle(htotal_fitted_, xjjroot::mycolor_middle["red"], 21, 0, xjjroot::mycolor_middle["red"], 1, 2, xjjroot::mycolor_middle["red"], 0.4, 1001);
  xjjroot::setthgrstyle(hnonprompt_fitted_, xjjroot::mycolor_middle["blue"], 21, 0, xjjroot::mycolor_middle["blue"], 2, 2, xjjroot::mycolor_middle["blue"], 0.45, 1001);
  xjjroot::setthgrstyle(hpull_, xjjroot::mycolor_middle["red"], 21, 0, xjjroot::mycolor_middle["red"], 2, 2, 0);
  xjjroot::setthgrstyle(hratio_, kBlack, 20, 1.7, kBlack, 1, 1, 0);
  auto* hpull_temp = (TH1D*)hpull_->Clone(xjjc::unique_str().c_str());
  xjjroot::setthgrstyle(hpull_temp, xjjroot::mycolor_middle["red"], 21, 0, xjjroot::mycolor_middle["red"], 1, 4, 0);
  auto* hnonprompt_fitted_temp = (TH1D*)hnonprompt_fitted_->Clone(xjjc::unique_str().c_str());
  xjjroot::setthgrstyle(hnonprompt_fitted_temp, 0, 0, 0, 0, 0, 0, 10, 1, 1001);

  std::vector<TH1D*> hlist = { hdata_, htotal_fitted_ };
  xjjana::sethsnonzeromin(hlist, 0.5);
  xjjana::sethsmax(hlist, 10.);
  xjjana::sethabsminmax(hpull_, -3.95, 3.95);
  xjjana::sethabsminmax(hratio_, 0.01, 1.99);
 
  auto* h_low = (q == Qual::pull ? hpull_ : hratio_);
  const auto iq = static_cast<size_t>(q);
  auto yref_q = yref[iq];
  auto drawopt_q = drawopt[iq];
  
  float tsizes = tsize/pratio;
  auto* leg = new TLegend(xleft, ytop-tsizes*1.1*4, xleft+tsizes*5, ytop);
  xjjroot::setleg(leg, tsizes);
  leg->AddEntry(hdata_, "Data (signal extracted)", "pe");
  leg->AddEntry(htotal_fitted_, "Fit", "l");
  xjjroot::addentrybystyle(leg, "Prompt", "f", { .lcolor = 0, .lstyle = 0, .lwidth = 0, .fcolor = htotal_fitted_->GetFillColor(), .falpha = 0.3 });
  leg->AddEntry(hnonprompt_fitted_, "Nonprompt", "f");
  
  auto pads = xjjroot::twopads(c, hdata_, h_low, pratio);
  pads[0]->SetLogy(1);
  pads[0]->cd();
  htotal_fitted_->Draw("hist same");
  hnonprompt_fitted_temp->Draw("hist same");
  hnonprompt_fitted_->Draw("hist same");
  hdata_->Draw("pe1 same");
  leg->Draw();
  xjjroot::drawtexgroup(xleft + 0.01, ytop - tsizes*1.1*4 - 0.01, {
      Form("#it{f}_{prompt} = %.2f%s (#pm%.2f)", fprompt(), ((fprompt_err_low_>=0 && fprompt_err_high_>=0) ? Form("^{ +%.2f}_{ -%.2f}", fprompt_err_high_, fprompt_err_low_) : ""), fprompt_err_par()),
      Form("#chi^{2} / ndf = %.2f (%.2f) / %d", chi2(), fitresult_->Chi2(), ndf())
    }, tsizes, 13, 42, 1.2);
  pads[0]->RedrawAxis();
  pads[1]->cd();
  xjjroot::drawline(h_low->GetXaxis()->GetXmin(), yref_q, h_low->GetXaxis()->GetXmax(), yref_q, kGray+3, 2, 2);
  h_low->Draw(Form("%s same", drawopt_q.c_str()));
  if (q == Qual::pull) hpull_temp->Draw("pe same");
  c->cd();
  return pads;
}

TH1D* fpfitter::clone_h1(TH1D* h, const std::string& suffix) {
  auto* hnew = (TH1D*)h->Clone(add_suffix__y(h->GetName(), suffix).c_str());
  return hnew;
}

bool fpfitter::same_axis(const TH1D* a, const TH1D* b) {
  if (a->GetNbinsX() != b->GetNbinsX()) return false;
  for (int i = 1; i <= a->GetNbinsX() + 1; ++i) {
    if (!xjjc::almost_eq(a->GetBinLowEdge(i), b->GetBinLowEdge(i))) return false;
  }
  return true;
}

void fpfitter::print_fitresult() {
  __XJJLOG << "++ fitting result" << std::endl;
  xjjc::print_tab<std::string>({ { "f_prompt", Form("%f +%f -%f (+/- %f)", fprompt(), fprompt_err_high_, fprompt_err_low_, fprompt_err_par()) },
                                 { "f_nonprompt", Form("%f +/- %f", 1-fprompt(), fprompt_err_par()) },
                                 { "chi2 (manual)", Form("%f", chi2_) },
                                 { "chi2 (fitresult)", Form("%f", fitresult_->Chi2()) },
                                 { "ndf", Form("%d", ndf()) },
    }, 0);
}

void fpfitter::write_to_file() {
  xjjroot::writehist(hdata_);
  xjjroot::writehist(htotal_fitted_);
  xjjroot::writehist(hnonprompt_fitted_);
  xjjroot::writehist(hpull_);
  xjjroot::writehist(hratio_);
  xjjroot::writehist(fitresult_.Get());
}
