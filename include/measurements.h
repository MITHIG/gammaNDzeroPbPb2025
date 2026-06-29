
namespace HIN_25_002_gammaN {
  std::vector<double> xbins = { -2, -1, 0, 1, 2 };
  std::vector<double> x = { -1.5, -0.5, 0.5, 1.5 };
  std::vector<double> xerr = { 0.5, 0.5, 0.5, 0.5 };
  std::vector<double> y = { 1.8097112860892386, 1.258530183727034, 0.6649505350292748, 0.1844336765596606 };
  std::vector<double> ystat = { 0.28619, 0.137795, 0.09363, 0.07950 };
  std::vector<double> ysyst = { 0.43988, 0.300323, 0.24909, 0.08480 };

  TGraphErrors* draw() {
    auto* gstat = new TGraphErrors(xbins.size()-1, x.data(), y.data(), xerr.data(), ystat.data());
    xjjroot::setthgrstyle(gstat, kBlack, 20, 1.6, kBlack, 1, 1);
    auto* gsyst = new TGraphErrors(xbins.size()-1, x.data(), y.data(), xerr.data(), ysyst.data());
    xjjroot::setthgrstyle(gsyst, kBlack, 20, 1.6, 0, 0, 0, kBlack, 0.9, 3004);

    gsyst->Draw("2 same");
    gstat->Draw("pe1 same");
    
    return gsyst;
  }
}
