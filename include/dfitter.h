#ifndef _XJJROOT_DFITTER_H_
#define _XJJROOT_DFITTER_H_

/*******************************************************************************************
 * Class : xjjroot::dfitter                                                                * 
 * This class provides tools fitting Dzero invariant mass spectra.                         * 
 * The object to be used can be declared by                                                * 
 *                                                                                         * 
 *    xjjroot::dfitter obj(options);                                                       * 
 *                                                                                         * 
 * Options supported are listed below                                                      * 
 *                                                                                         * 
 *   "3"  : Using 3-Gaussian function to model signal (default is 2-Gaussian function)     * 
 *   "P"  : Have peaky background                                                          *
 *   "S"  : Draw significance info and lines at signal region                              * 
 *   "V"  : Switch off Quiet mode of fitting                                               * 
 *                                                                                         * 
 * The core function of this class is                                                      * 
 *                                                                                         * 
 *    TF1* fit(TH1*, TH1*, TH1*, std::string, std::string, std::vector<std::string>)       * 
 *                                                                                         * 
 ******************************************************************************************/

#include <TString.h>
#include <TMath.h>
#include <TFitResult.h>
#include <Math/DistFunc.h>

namespace xjjroot {
  const std::map<std::string, thgrstyle> fstyle = {
    { "h", thgrstyle(kBlack, 20, 1.4, kBlack, 1, 1) },
    { "f", thgrstyle(-1, -1, -1, kRed+1, 1, 3) },
    { "match", thgrstyle(kBlack, 20, 1.4, kOrange-3, 2, 3, kOrange-3, 0.2, 1001) },
    { "matchtot", thgrstyle(kBlack, 20, 1.4, kOrange-3, 1, 3, kOrange-3, 0.2, 1001) },
    { "swap", thgrstyle(kCyan+3, 20, 1.4, kCyan+3, 1, 3, kCyan+3, 1, 3004, -1, 0.6) },
    { "kk", thgrstyle(kMagenta-5, 20, 1.4, kMagenta-5, 1, 3, kMagenta-5, 1, 3005, -1, 0.6) },
    { "pipi", thgrstyle(kBlue-5, 20, 1.4, kBlue-5, 1, 3, kBlue-5, 1, 3006, -1, 0.6) },
    { "background", thgrstyle(kAzure-6, 20, 1.4, kAzure-6, 2, 3) },
    { "notmatch", thgrstyle(kGray, -1, -1, kGray, 2, 3) },
  };

  class dfitter
  {
  public:
    dfitter(Option_t* option = "");
    ~dfitter() {};

    void fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped,
             const TH1* hmassMCKK, const TH1* hmassMCPiPi);
    bool fitted() const { return fitted_; }
    void reset();

    static void set_hist(TH1* h);

    std::vector<std::string> draw_result(float x = 0.25, float y = 0.86, float tsize = 0.035, float lspacescale = 1.15) const;
    void draw_params(float x = 0.25, float y = 0.86, float tsize = 0.035, float lspacescale = 1.15) const;
    void draw_fmc() const;
    void draw_leg(float x1 = 0.65, float y2 = 0.88) { xjjroot::moveleg_n_draw(leg_, x1, y2); }
    void draw_legmc(float x1 = 0.69, float y2 = 0.70) { xjjroot::moveleg_n_draw(legmc_, x1, y2); }
    
    double yield() const { return yield_; }
    double yieldErr() const { return yieldErr_; }
    double S() const { return S_; }
    double B() const { return B_; }
    double chi2() const { return 2.*r_->MinFcnValue(); }
    double ndf() const { return fun_f_->GetNDF(); }
    double chi2prob() const { return TMath::Prob(chi2(), ndf()); }
    std::pair<double, double> width_mc_match(double frac = 0.682689) const;

    bool haskkpipi() const { return opt_haskkpipi_; }

    TF1* f_f(const std::string& name) const;
    TF1* f_match(const std::string& name = "") const;
    TF1* f_swap(const std::string& name = "") const;
    TF1* f_kk(const std::string& name = "") const;
    TF1* f_pipi(const std::string& name = "") const;
    TF1* f_background(const std::string& name = "") const;
    TF1* f_notmatch(const std::string& name = "") const;

  private:
    double S_;
    double B_;
    double yield_;
    double yieldErr_;

    TF1* fun_f_;
    TF1* fun_mc_match_; // need this member to save the parameters before changing by fitting on data
    TF1* fun_mc_swap_;
    TF1* fun_mc_kk_;
    TF1* fun_mc_pipi_;
    std::vector<TF1*> vfun_mc_match_;
    TLegend *leg_, *legmc_;
    
    TFitResultPtr r_;

    std::string option_;
    bool opt_3gaus_;
    bool opt_haskkpipi_;
    bool opt_sig_;
    bool opt_verbose_;

    // status
    bool fitted_;

    const float tsize_ = 0.038;
    double xmin_ = 0, xmax_ = 0, binwidth_ = 0;
    double signal_region_l_ = 1.8649 - 0.045;
    double signal_region_h_ = 1.8649 + 0.045;

    void make_legs();
    void calculate_SnB();
    void parse_opt();
    void parse_fmc();

    TF1* clone_fun(const TF1* fun, const std::string& fun_name) const;
  };
}

xjjroot::dfitter::dfitter(Option_t* option) :
  option_(option), fun_f_(nullptr), fun_mc_match_(nullptr), fun_mc_swap_(nullptr),
  fun_mc_kk_(nullptr), fun_mc_pipi_(nullptr) {
  vfun_mc_match_.clear();
  parse_opt();
  reset();
  make_legs();
}

void xjjroot::dfitter::parse_opt() {
  opt_3gaus_ = xjjc::str_contains(option_, "3");
  opt_haskkpipi_ = xjjc::str_contains(option_, "P");
  opt_sig_ = xjjc::str_contains(option_, "S");
  opt_verbose_ = xjjc::str_contains(option_, "V");
}

void xjjroot::dfitter::reset() {
  fitted_ = false;
  // for (auto& s : { S_, B_, yield_, yieldErr_ }) s = -1;
  // for (auto& f : { fun_f_, fun_mc_match_, fun_mc_swap_, fun_mc_kk_, fun_mc_pipi_ }) f = nullptr;
  S_ = -1;
  B_ = -1;
  yield_ = -1;
  yieldErr_ = -1;
  delete fun_f_; fun_f_ = nullptr;
  delete fun_mc_match_; fun_mc_match_ = nullptr;
  delete fun_mc_swap_; fun_mc_swap_ = nullptr;
  delete fun_mc_kk_; fun_mc_kk_ = nullptr;
  delete fun_mc_pipi_; fun_mc_pipi_ = nullptr;
  for (auto& f : vfun_mc_match_) delete f;
  vfun_mc_match_.clear();
}

void xjjroot::dfitter::make_legs() {
  leg_ = new TLegend(0.63, 0.865-tsize_*1.25*(opt_haskkpipi_ ? 7 : 5), 0.63+0.25, 0.865);
  xjjroot::setleg(leg_, tsize_);
  xjjroot::addentrybystyle(leg_, "Data", "pe", fstyle.at("h"));
  xjjroot::addentrybystyle(leg_, "Fit", "l", fstyle.at("f"));
  xjjroot::addentrybystyle(leg_, xjjroot::CMS::Dz + " +#scale[0.4]{ }" + xjjroot::CMS::Dzbar + " Signal", "f", fstyle.at("match"));
  xjjroot::addentrybystyle(leg_, "K-#pi swapped", "f", fstyle.at("swap"));
  if (opt_haskkpipi_) {
    xjjroot::addentrybystyle(leg_, xjjroot::CMS::Dz + "(" + xjjroot::CMS::Dzbar + ")#scale[0.4]{ }#rightarrow KK", "f", fstyle.at("kk"));
    xjjroot::addentrybystyle(leg_, xjjroot::CMS::Dz + "(" + xjjroot::CMS::Dzbar + ")#scale[0.4]{ }#rightarrow#scale[0.4]{ }#pi#pi", "f", fstyle.at("pipi"));
  }
  xjjroot::addentrybystyle(leg_, "Combinatorial", "l", fstyle.at("background"));
  leg_->Draw();

  legmc_ = new TLegend(0.66, 0.63-tsize_*1.25*(opt_haskkpipi_ ? 4 : 2), 0.66+0.2, 0.63);
  xjjroot::setleg(legmc_, tsize_);
  xjjroot::addentrybystyle(legmc_, "Signal", "pl", fstyle.at("matchtot"));
  xjjroot::addentrybystyle(legmc_, "K-#pi swapped", "pl", fstyle.at("swap"));
  if (opt_haskkpipi_) {
    xjjroot::addentrybystyle(legmc_, xjjroot::CMS::Dz + "(" + xjjroot::CMS::Dzbar + ")#scale[0.4]{ }#rightarrow KK", "pl", fstyle.at("kk"));
    xjjroot::addentrybystyle(legmc_, xjjroot::CMS::Dz + "(" + xjjroot::CMS::Dzbar + ")#scale[0.4]{ }#rightarrow#scale[0.4]{ }#pi#pi", "pl", fstyle.at("pipi"));
  }
  legmc_->Draw();
}

void xjjroot::dfitter::fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped,
                           const TH1* hmassMCKK, const TH1* hmassMCPiPi) {

  reset();

  if (!hmass || !hmassMCSignal || !hmassMCSwapped) {
    __XJJLOG << "!! bad histograms" << std::endl;
    return;
  }
  if (!hmassMCKK || !hmassMCPiPi) opt_haskkpipi_ = false;
  xmin_ = hmass->GetXaxis()->GetXmin(); //
  xmax_ = hmass->GetXaxis()->GetXmax(); //
  binwidth_ = (xmax_ - xmin_) / hmass->GetXaxis()->GetNbins();
  if (fabs(binwidth_ - hmass->GetBinWidth(1)) > 1.e-7) {
    __XJJLOG << "!! bad bin width: xmin = " << xmin_
             << ", xmax = " << xmax_
             << ", binwidth = (xmax_ - xmin_) / hmass->GetXaxis()->GetNbins() = \e[1m" << binwidth_ << "\e[0m"
             << "vs hmass->GetBinWidth(1) = \e[1m" << hmass->GetBinWidth(1) << "\e[0m"
             << std::endl;
    if (binwidth_ == 0) binwidth_ = hmass->GetBinWidth(1);
  }

  fitted_ = true;
  
  std::string str_fun_f = 
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*([12]*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11]))+(1-[12])*TMath::Gaus(x,[1],[13]*(1+[11]))/(sqrt(2*3.14159)*[13]*(1+[11]))))+[15]*TMath::Gaus(x,[14],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11]))+[16]*ROOT::Math::crystalball_pdf(x,[19],[20],[18],[17])+(1-[7]-[15]-[16])*ROOT::Math::crystalball_pdf(-x,[23],[24],[22],-[21]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x";

  fun_f_ = new TF1(Form("f_%s", xjjc::unique_str().c_str()), str_fun_f.c_str(), xmin_, xmax_);
  fun_f_->SetNpx(2000);
  xjjroot::setthgrstyle(fun_f_, fstyle.at("f"));
  
  auto* h = (TH1F*)hmass->Clone(Form("h_%s", xjjc::unique_str().c_str()));
  set_hist(h);
  auto* hMCSignal = (TH1F*)hmassMCSignal->Clone(Form("hMCSignal_%s", xjjc::unique_str().c_str()));
  set_hist(hMCSignal);
  auto* hMCSwapped = (TH1F*)hmassMCSwapped->Clone(Form("hMCSwapped_%s", xjjc::unique_str().c_str()));
  set_hist(hMCSwapped);
  auto* hMCKK = hmassMCKK ? (TH1F*)hmassMCKK->Clone(Form("hMCKK_%s", xjjc::unique_str().c_str())) : nullptr;
  set_hist(hMCKK);
  auto* hMCPiPi = hmassMCPiPi ? (TH1F*)hmassMCPiPi->Clone(Form("hMCPiPi_%s", xjjc::unique_str().c_str())) : nullptr;
  set_hist(hMCPiPi);
  const auto peak_mc_kk = hMCKK ? hMCKK->GetBinCenter(hMCKK->GetMaximumBin()) : (double)0;
  const auto peak_mc_pipi = hMCPiPi ? hMCPiPi->GetBinCenter(hMCPiPi->GetMaximumBin()) : (double)0;

  const char* fitopt = opt_verbose_?"L m":"L m q";
  
  const double param_init_0 = 100.,
    param_init_1 = 1.8649,
    param_init_14 = 1.8649,
    param_init_2 = 0.03,
    param_init_10 = 0.005,
    param_init_13 = 0.002,
    param_init_8 = 0.1,
    param_init_9 = 0.1,
    param_init_12 = 0.5,
    param_init_17 = 1.75,
    param_init_18 = 0.04,
    param_init_19 = 1.5,
    param_init_20 = 3.,
    param_init_21 = 1.98,
    param_init_22 = 0.04,
    param_init_23 = 1.5,
    param_init_24 = 3.;

  //
  fun_f_->FixParameter(8,  param_init_8);
  fun_f_->FixParameter(17, param_init_17);
  fun_f_->FixParameter(18, param_init_18);
  fun_f_->FixParameter(19, param_init_19);
  fun_f_->FixParameter(20, param_init_20);
  fun_f_->FixParameter(21, param_init_21);
  fun_f_->FixParameter(22, param_init_22);
  fun_f_->FixParameter(23, param_init_23);
  fun_f_->FixParameter(24, param_init_24);

  // -- fit MC
  fun_f_->FixParameter(3, 0);
  fun_f_->FixParameter(4, 0);
  fun_f_->FixParameter(5, 0);
  fun_f_->FixParameter(6, 0);
  fun_f_->FixParameter(11, 0);

  //  - fit signal
  // std::cout<<"fit signal"<<std::endl;
  fun_f_->FixParameter(7,  1);
  fun_f_->FixParameter(15, 0);
  fun_f_->FixParameter(16, 0);
  //
  fun_f_->FixParameter(1,  param_init_1);
  fun_f_->FixParameter(14, param_init_14);
  fun_f_->SetParameter(0,  param_init_0);
  fun_f_->SetParameter(2,  param_init_2);
  fun_f_->SetParLimits(2,  0.01,  0.5);
  fun_f_->SetParameter(10, param_init_10);
  fun_f_->SetParLimits(10, 0.001, 0.05);
  fun_f_->SetParameter(13, param_init_13); // 3gaus
  fun_f_->SetParLimits(13, 0.002, 0.1); // 3gaus
  fun_f_->SetParameter(9,  param_init_9);
  fun_f_->SetParLimits(9, 0, 1);
  if (opt_3gaus_) {
    fun_f_->SetParameter(12, param_init_12);
    fun_f_->SetParLimits(12, 0, 1);
  } else {
    fun_f_->FixParameter(12, 1);
  }
  hMCSignal->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  fun_f_->ReleaseParameter(1);
  fun_f_->SetParLimits(1, 1.85, 1.9);
  hMCSignal->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSignal->Fit(fun_f_->GetName(), fitopt, "", xmin_, xmax_);

  fun_mc_match_ = f_match(Form("%s_match_mc", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_match_, -1, -1, -1, fstyle.at("match").lcolor, 1, fstyle.at("match").lwidth, 0, 0, 0);

  const auto norm_mc_match = fun_f_->GetParameter(0);
  fun_f_->FixParameter(1, fun_f_->GetParameter(1));
  fun_f_->FixParameter(2, fun_f_->GetParameter(2));
  fun_f_->FixParameter(10, fun_f_->GetParameter(10));
  fun_f_->FixParameter(13, fun_f_->GetParameter(13)); // 3gaus
  fun_f_->FixParameter(9, fun_f_->GetParameter(9));
  if (opt_3gaus_) fun_f_->FixParameter(12, fun_f_->GetParameter(12));

  //   - fit swapped
  // std::cout<<"fit swap"<<std::endl;
  fun_f_->FixParameter(7, 0);
  fun_f_->FixParameter(15, 1);
  fun_f_->FixParameter(16, 0);
  fun_f_->ReleaseParameter(8);
  fun_f_->SetParameter(8, param_init_8);
  fun_f_->SetParLimits(8,  0.02,  0.2);
  fun_f_->FixParameter(14, fun_f_->GetParameter(1));
  //  
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  fun_f_->ReleaseParameter(14);
  fun_f_->SetParLimits(14, 1.85, 1.9);
  hMCSwapped->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  hMCSwapped->Fit(fun_f_->GetName(), fitopt,"", xmin_, xmax_);
  
  fun_mc_swap_ = f_swap(Form("%s_mc_swap", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_swap_, -1, -1, -1, fstyle.at("swap").lcolor, 1, fstyle.at("swap").lwidth, 0, 0, 0);

  const auto norm_mc_swap = fun_f_->GetParameter(0);
  fun_f_->FixParameter(14, fun_f_->GetParameter(14)); // swap mean
  fun_f_->FixParameter(8, fun_f_->GetParameter(8)); // swap width

  //   - fit KK
  if (opt_haskkpipi_) {
    // std::cout<<"fit KK"<<std::endl;
    fun_f_->FixParameter(7, 0);
    fun_f_->FixParameter(15, 0);
    fun_f_->FixParameter(16, 1);

    fun_f_->ReleaseParameter(17);
    fun_f_->ReleaseParameter(18);
    fun_f_->ReleaseParameter(19);
    fun_f_->ReleaseParameter(20);
    fun_f_->SetParameter(17, peak_mc_kk);
    fun_f_->SetParameter(18, param_init_18);
    fun_f_->SetParameter(19, param_init_19);
    fun_f_->SetParameter(20, param_init_20);
    fun_f_->SetParLimits(17, std::max(xmin_, peak_mc_kk - 0.05), std::min(xmax_, peak_mc_kk + 0.05));
    fun_f_->SetParLimits(18, 0.002, 0.2);
    fun_f_->SetParLimits(19, 0.1,   10.);
    fun_f_->SetParLimits(20, 1.01,  100.);

    hMCKK->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCKK->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCKK->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCKK->Fit(fun_f_->GetName(), fitopt,"", xmin_, xmax_);
  } else {
    fun_f_->SetParameter(0, 0);
  }
  fun_mc_kk_ = f_kk(Form("%s_mc_kk", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_kk_, fstyle.at("kk"));
  const auto norm_mc_kk = fun_f_->GetParameter(0);
  fun_f_->FixParameter(17, fun_f_->GetParameter(17));
  fun_f_->FixParameter(18, fun_f_->GetParameter(18));
  fun_f_->FixParameter(19, fun_f_->GetParameter(19));
  fun_f_->FixParameter(20, fun_f_->GetParameter(20));

  //   - fit pipi
  if (opt_haskkpipi_) {
    // std::cout<<"fit pipi"<<std::endl;
    fun_f_->FixParameter(7, 0);
    fun_f_->FixParameter(15, 0);
    fun_f_->FixParameter(16, 0);

    fun_f_->ReleaseParameter(21);
    fun_f_->ReleaseParameter(22);
    fun_f_->ReleaseParameter(23);
    fun_f_->ReleaseParameter(24);
    fun_f_->SetParameter(21, peak_mc_pipi);
    fun_f_->SetParameter(22, param_init_22);
    fun_f_->SetParameter(23, param_init_23);
    fun_f_->SetParameter(24, param_init_24);
    fun_f_->SetParLimits(21, std::max(xmin_, peak_mc_pipi - 0.05), std::min(xmax_, peak_mc_pipi + 0.05));
    fun_f_->SetParLimits(22, 0.002, 0.2);
    fun_f_->SetParLimits(23, 0.1,   10.);
    fun_f_->SetParLimits(24, 1.01,  100.);

    hMCPiPi->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCPiPi->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCPiPi->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
    hMCPiPi->Fit(fun_f_->GetName(), fitopt,"", xmin_, xmax_);
  } else {
    fun_f_->SetParameter(0, 0);
  }
  fun_mc_pipi_ = f_pipi(Form("%s_mc_pipi", fun_f_->GetName()));
  xjjroot::setthgrstyle(fun_mc_pipi_, fstyle.at("pipi"));
  const auto norm_mc_pipi = fun_f_->GetParameter(0);
  fun_f_->FixParameter(21, fun_f_->GetParameter(21));
  fun_f_->FixParameter(22, fun_f_->GetParameter(22));
  fun_f_->FixParameter(23, fun_f_->GetParameter(23));
  fun_f_->FixParameter(24, fun_f_->GetParameter(24));

  parse_fmc();
  
  //  -- fit data
  const auto norm_mc_total = norm_mc_match + norm_mc_swap + norm_mc_kk + norm_mc_pipi;
  fun_f_->FixParameter(7, norm_mc_match/norm_mc_total);
  fun_f_->FixParameter(15, norm_mc_swap/norm_mc_total);
  fun_f_->FixParameter(16, norm_mc_kk/norm_mc_total);
  fun_f_->ReleaseParameter(3);
  fun_f_->ReleaseParameter(4);
  fun_f_->ReleaseParameter(5);
  fun_f_->ReleaseParameter(6);
  
  h->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "q", "", xmin_, xmax_);
  fun_f_->ReleaseParameter(1);
  fun_f_->SetParLimits(1, 1.86, 1.87);
  // fun_f_->ReleaseParameter(11);
  // fun_f_->SetParLimits(11, -0.5, 0.5);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  h->Fit(fun_f_->GetName(), "L q", "", xmin_, xmax_);
  r_ = h->Fit(fun_f_->GetName(), Form("%s S", fitopt),"", xmin_, xmax_);

  auto* fun_background = f_background();
  auto* fun_match = f_match();
  auto* fun_swap = f_swap();
  auto* fun_kk = f_kk();
  auto* fun_pipi = f_pipi();
  auto* fun_notmatch = f_notmatch();

  yield_ = fun_match->Integral(xmin_, xmax_)/binwidth_;
  yieldErr_ = fun_match->Integral(xmin_, xmax_)/binwidth_*fun_match->GetParError(0)/fun_match->GetParameter(0);

  calculate_SnB();
  
  h->Draw("pe1");
  fun_match->Draw("same");
  fun_background->Draw("same");
  fun_swap->Draw("same");
  if (opt_haskkpipi_) {
    fun_kk->Draw("same");
    fun_pipi->Draw("same");
  }
  if (opt_sig_) {
    fun_notmatch->SetRange(signal_region_l_, signal_region_h_);
    fun_notmatch->Draw("same");
    xjjroot::drawline(signal_region_l_, 0, signal_region_l_, fun_f_->Eval(signal_region_l_), fun_notmatch->GetLineColor(), fun_notmatch->GetLineStyle(), fun_notmatch->GetLineWidth());
    xjjroot::drawline(signal_region_h_, 0, signal_region_h_, fun_f_->Eval(signal_region_h_), fun_notmatch->GetLineColor(), fun_notmatch->GetLineStyle(), fun_notmatch->GetLineWidth());
  }
  fun_f_->Draw("same");
}

std::vector<std::string> xjjroot::dfitter::draw_result(float x, float y, float tsize, float lspacescale) const {
  std::vector<std::string> text = {
    Form("N = %.0f#scale[0.5]{ }#pm %.0f", yield_, yieldErr_),
    Form("#chi^{2} / ndf = %.1f / %.0f", chi2(), ndf()),
    Form("Prob = %.2f%s", chi2prob()*100, "%"),
    Form("#bar{m} = %.3f%s", fun_f_->GetParameter(1), (fun_f_->GetParError(1)==0 ? " (fixed)" : Form("#scale[0.4]{ }#pm %.3f", fun_f_->GetParError(1)))),
  };
  std::vector<Color_t> color(text.size(), kBlack);
  if (opt_sig_) {
    text.push_back(Form("S = %.0f, B = %.0f", S_, B_));
    text.push_back(Form("S/#sqrt{S+B} = %.1f", S_/TMath::Sqrt(S_ + B_)));
    xjjc::vec_append(color, std::vector<Color_t>(2, fstyle.at("notmatch").lcolor));
  }
  xjjroot::drawtexgroup(x, y, text, tsize, 13, 42, lspacescale,
                        1, 0.2, color);

  // if (fdrawdetail) {
  //   drawtex(texxpos, texypos=(texypos-texdypos), Form("N#scale[0.6]{#lower[0.7]{sig}}/(N#scale[0.6]{#lower[0.7]{sig}}+N#scale[0.6]{#lower[0.7]{swap}}) = %.2f",fun_f_->GetParameter(7)));
  // }
  return text;
}

void xjjroot::dfitter::draw_fmc() const {
  fun_mc_swap_->Draw("same");
  fun_mc_match_->Draw("same");
  for (const auto& f : vfun_mc_match_)
    f->Draw("same");
  if (opt_haskkpipi_) {
    fun_mc_kk_->Draw("same");
    fun_mc_pipi_->Draw("same");
  }
}


void xjjroot::dfitter::calculate_SnB() {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return;
  }
  auto* fun_match = f_match(Form("fun_match_%s", xjjc::unique_str().c_str()));
  S_ = fun_match->Integral(signal_region_l_, signal_region_h_)/binwidth_;
  delete fun_match;
  auto* fun_notmatch = f_notmatch(Form("fun_notmatch_%s", xjjc::unique_str().c_str()));
  B_ = fun_notmatch->Integral(signal_region_l_, signal_region_h_)/binwidth_;
  delete fun_notmatch;
}

TF1* xjjroot::dfitter::f_match(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_match", fun_f_->GetName()) : name;
  std::string str_fun_match = "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))";
  auto* fun = new TF1(fname.c_str(), str_fun_match.c_str(), fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 1 }, { 2, 2 }, { 3, 7 }, { 4, 9 }, { 5, 10 }, { 6, 11 }, { 7, 12 }, { 8, 13 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  fun->SetNpx(2000);
  xjjroot::setthgrstyle(fun, fstyle.at("match"));

  return fun;
}

TF1* xjjroot::dfitter::f_swap(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_swap", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]*[2]*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 14 }, { 2, 15 }, { 3, 8 }, { 4, 11 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("swap"));

  return fun;
}

TF1* xjjroot::dfitter::f_kk(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_kk", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]*[1]*ROOT::Math::crystalball_pdf(x,[4],[5],[3],[2])", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 16 }, { 2, 17 }, { 3, 18 }, { 4, 19 }, { 5, 20 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("kk"));

  return fun;
}

TF1* xjjroot::dfitter::f_pipi(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_pipi", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]*(1-[1]-[2]-[3])*ROOT::Math::crystalball_pdf(-x,[6],[7],[5],-[4])", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 7 }, { 2, 15 }, { 3, 16 }, { 4, 21 }, { 5, 22 }, { 6, 23 }, { 7, 24 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("pipi"));

  return fun;
}

TF1* xjjroot::dfitter::f_background(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_background", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]+[1]*x+[2]*x*x+[3]*x*x*x", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 3 }, { 1, 4 }, { 2, 5 }, { 3, 6 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("background"));

  return fun;
}

TF1* xjjroot::dfitter::f_notmatch(const std::string& name) const {
  if (!fitted_) {
    __XJJLOG << "!! not fitted yet" << std::endl;
    return nullptr;
  }
  std::string fname = name.empty() ? Form("%s_notmatch", fun_f_->GetName()) : name;
  auto* fun = new TF1(fname.c_str(), "[0]*([2]*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))+[9]*ROOT::Math::crystalball_pdf(x,[12],[13],[11],[10])+(1-[14]-[2]-[9])*ROOT::Math::crystalball_pdf(-x,[17],[18],[16],-[15]))+[5]+[6]*x+[7]*x*x+[8]*x*x*x", fun_f_->GetXmin(), fun_f_->GetXmax());
  std::map<int, int> params = {
    { 0, 0 }, { 1, 14 }, { 2, 15 }, { 3, 8 }, { 4, 11 },
    { 5, 3 }, { 6, 4 }, { 7, 5 }, { 8, 6 },
    { 9, 16 }, { 10, 17 }, { 11, 18 }, { 12, 19 }, { 13, 20 },
    { 14, 7 }, { 15, 21 }, { 16, 22 }, { 17, 23 }, { 18, 24 }
  };
  for (const auto& p : params) {
    fun->SetParameter(p.first, fun_f_->GetParameter(p.second));
    fun->SetParError(p.first, fun_f_->GetParError(p.second));
  }
  xjjroot::setthgrstyle(fun, fstyle.at("notmatch"));

  return fun;
}

std::pair<double, double> xjjroot::dfitter::width_mc_match(double frac) const {
  auto width = xjjana::tf_width(fun_mc_match_, fun_mc_match_->GetParameter(1), frac, { 0, 1, 3, 4, 6, 7 });
  return width;
}

TF1* xjjroot::dfitter::clone_fun(const TF1* fun, const std::string& name) const {
  auto* newfun = new TF1(*fun);
  newfun->SetName(name.c_str());
  return newfun;
}

void xjjroot::dfitter::set_hist(TH1* h) {
  if (!h) return;
  h->SetXTitle("m_{#piK} [GeV/c^{2}]");
  h->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", h->GetBinWidth(1)*1.e+3));
  xjjroot::sethempty(h, 0., 0.5);
  xjjroot::setthgrstyle(h, fstyle.at("h"));
  h->SetMaximum(-1111);
  h->SetAxisRange(0, h->GetMaximum()*1.4*1.2, "Y");
}

void xjjroot::dfitter::parse_fmc() {
  std::string str_fun_match = "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))";
  std::vector<std::vector<double>> params = {
    { fun_mc_match_->GetParameter(0)*fun_mc_match_->GetParameter(3)*fun_mc_match_->GetParameter(4), fun_mc_match_->GetParameter(1), fun_mc_match_->GetParameter(2) },
    { fun_mc_match_->GetParameter(0)*fun_mc_match_->GetParameter(3)*(1-fun_mc_match_->GetParameter(4))*fun_mc_match_->GetParameter(7), fun_mc_match_->GetParameter(1), fun_mc_match_->GetParameter(5) },
  };
  std::vector<std::vector<double>> parerrors = {
    { fun_mc_match_->GetParError(0)*fun_mc_match_->GetParameter(3)*fun_mc_match_->GetParameter(4), fun_mc_match_->GetParError(1), fun_mc_match_->GetParError(2) },
    { fun_mc_match_->GetParError(0)*fun_mc_match_->GetParameter(3)*(1-fun_mc_match_->GetParameter(4))*fun_mc_match_->GetParameter(7), fun_mc_match_->GetParError(1), fun_mc_match_->GetParError(5) },    
  };
  if (opt_3gaus_) {
    params.push_back( { fun_mc_match_->GetParameter(0)*fun_mc_match_->GetParameter(3)*(1-fun_mc_match_->GetParameter(4))*(1-fun_mc_match_->GetParameter(7)), fun_mc_match_->GetParameter(1), fun_mc_match_->GetParameter(8) } );
    parerrors.push_back( { fun_mc_match_->GetParError(0)*fun_mc_match_->GetParameter(3)*(1-fun_mc_match_->GetParameter(4))*(1-fun_mc_match_->GetParameter(7)), fun_mc_match_->GetParError(1), fun_mc_match_->GetParError(8) } );
  }
  for (int i=0; i<params.size(); i++) {
    auto* f1 = new TF1(Form("%s-%d", fun_mc_match_->GetName(), i), "[0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])", xmin_, xmax_);
    for (int j=0; j<params.at(i).size(); j++) {
      f1->SetParameter(j, params.at(i).at(j));
      f1->SetParError(j, parerrors.at(i).at(j));
    }
    xjjroot::setlinestyle(f1, fun_mc_match_->GetLineColor(), 2, fun_mc_match_->GetLineWidth());
    f1->SetNpx(2000);
    vfun_mc_match_.push_back(f1);
  }
}

void xjjroot::dfitter::draw_params(float x, float y, float tsize, float lspacescale) const {
  float norm = 0;
  std::vector<std::pair<double, double>> params_match;
  for (const auto& f : vfun_mc_match_) {
    params_match.push_back( { f->GetParameter(2), f->GetParameter(0) } );
    norm += f->GetParameter(0);
  }
  std::sort(params_match.begin(), params_match.end(), [](std::pair<double, double> a, std::pair<double, double> b) {
    return a.first < b.first; // width
  });
  std::string sigma_match, frac_match;
  for (auto& p : params_match) {
    sigma_match += std::string(Form("%s%.3f", (sigma_match.empty() ? "" : ", "), p.first));
    p.second /= norm;
    frac_match += std::string(Form("%s%.2f", (frac_match.empty() ? "" : ", "), p.second));
  }
  const auto norm_masswin = fun_mc_match_->Integral(xmin_, xmax_) + fun_mc_swap_->Integral(xmin_, xmax_) + fun_mc_kk_->Integral(xmin_, xmax_) + fun_mc_pipi_->Integral(xmin_, xmax_);
  std::vector<std::string> rtex = {
    // "Mean (data) = " + std::string(Form("%.3f%s", fun_f_->GetParameter(1), (fun_f_->GetParError(1)==0 ? " (fixed)" : Form("#scale[0.4]{ }#pm %.3f", fun_f_->GetParError(1))))),
    "Mean = " + std::string(Form("%.3f (signal), %.3f (swap)", fun_mc_match_->GetParameter(1), fun_mc_swap_->GetParameter(1))),
    "Signal#scale[0.5]{ }#sigma in MC = " + sigma_match,
    "Fraction of each gaus = " + frac_match,
    "#sigma_{data}/#sigma_{MC} - 1 = " + std::string(Form("%.2f %s", fun_f_->GetParameter(11), (fun_f_->GetParError(11)==0 ? "(fixed)" : Form("#pm %.2f", fun_f_->GetParError(11))))),
    Form("N#scale[0.7]{sig} / (N#scale[0.7]{sig} + N#scale[0.7]{swap}) = %.2f", fun_f_->GetParameter(7)/(fun_f_->GetParameter(7)+fun_f_->GetParameter(15))),
  };
  if (fun_f_->GetParError(11)==0) rtex.erase(rtex.begin() + 3);
  if (opt_haskkpipi_) {
    xjjc::vec_append(rtex, {
        "In mass window:",
        Form("#rightarrow Signal %.2f", fun_mc_match_->Integral(xmin_, xmax_) / norm_masswin),
        Form("#rightarrow Swap %.2f", fun_mc_swap_->Integral(xmin_, xmax_) / norm_masswin),
        Form("#rightarrow KK %.2f", fun_mc_kk_->Integral(xmin_, xmax_) / norm_masswin),
        Form("#rightarrow#scale[0.4]{ }#pi#pi %.2f", fun_mc_pipi_->Integral(xmin_, xmax_) / norm_masswin),
        // Form("Fractions: sig %.2f, swap %.2f, KK %.2f, #pi#pi %.2f",
        //      fun_f_->GetParameter(7), fun_f_->GetParameter(14), fun_f_->GetParameter(15),
        //      1. - fun_f_->GetParameter(7) - fun_f_->GetParameter(14) - fun_f_->GetParameter(15)),
      });
  }
  xjjroot::drawtexgroup(x, y, rtex, tsize, 13, 42, lspacescale);
}

#endif
