enum Event { gammaN, Ngamma, Other };

namespace measurement {
  TGraphErrors* draw(const std::vector<double>& xbins, const std::vector<double>& y,
                     const std::vector<double>& ystat, const std::vector<double>& ysyst) {
    std::vector<double> x, xerr;
    for (int i=0; i<xbins.size()-1; i++) {
      x.push_back((xbins[i+1] + xbins[i]) / 2.);
      xerr.push_back((xbins[i+1] - xbins[i]) / 2.);
    }
    auto color = kGray+1;
    auto* gstat = new TGraphErrors(xbins.size()-1, x.data(), y.data(), xerr.data(), ystat.data());
    xjjroot::setthgrstyle(gstat, color, 20, 1.6, color, 1, 1);
    auto* gsyst = new TGraphErrors(xbins.size()-1, x.data(), y.data(), xerr.data(), ysyst.data());
    xjjroot::setthgrstyle(gsyst, color, 20, 1.6, 0, 0, 0, color, 0.8, 3004);

    gsyst->Draw("2 same");
    gstat->Draw("pe1 same");
    
    return gsyst;
  }

  TGraphErrors* draw_HIN_25_002(Event e);
}


namespace HIN_25_002_gammaN {
  const std::vector<double> xbins = { -2, -1, 0, 1, 2 };
  const std::vector<double> y = { 1.80971, 1.25853, 0.66495, 0.18443 };
  const std::vector<double> ystat = { 0.28619, 0.137795, 0.09363, 0.07950 };
  const std::vector<double> ysyst = { 0.43988, 0.300323, 0.24909, 0.08480 };
  TGraphErrors* draw() { return measurement::draw(xbins, y, ystat, ysyst); }
}

namespace HIN_25_002_Ngamma {
  const std::vector<double> xbins = { -2, -1, 0, 1, 2 };
  const std::vector<double> y = { 0.32143, 0.615659, 1.36978, 1.29313 };
  const std::vector<double> ystat = { 0.08901, 0.09148, 0.15824, 0.27198};
  const std::vector<double> ysyst = { 0.13104, 0.17308, 0.38571, 0.50440 };
  TGraphErrors* draw() { return measurement::draw(xbins, y, ystat, ysyst); }
}

TGraphErrors* measurement::draw_HIN_25_002(Event e) {
  if (e == Event::gammaN)
    return HIN_25_002_gammaN::draw();
  else if (e == Event::Ngamma)
    return HIN_25_002_Ngamma::draw();
  return nullptr;
}
