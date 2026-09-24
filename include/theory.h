
namespace fonll {
  struct Model {
    std::string name = "";
    Color_t color = kBlack;
    std::string tex = "";
    std::map<std::string, TGraphAsymmErrors*> gs;
  };

  class DrawSets {
  public:
    DrawSets() { ; }
    DrawSets(const std::string&, double ptmin, double ptmax);
    TLegend* draw();
    bool valid() const { return valid_; }
    size_t n() const { return models.size(); }
    double ymaximum();
  private:
    bool valid_;
    std::vector<Model> models;
  };

  const std::vector<Model> configs = {
    { .name = "CT18ANLO_mc13", .color = static_cast<Color_t>(TColor::GetColor("#ed9418")), .tex = "CT18ANLO (pPDF, m_{c}=1.3 GeV)" },
    { .name = "CT18ANLO_mc15", .color = static_cast<Color_t>(TColor::GetColor("#e44760")), .tex = "CT18ANLO (pPDF, m_{c}=1.5 GeV)" },
    { .name = "EPPS21Pb_mc13", .color = static_cast<Color_t>(TColor::GetColor("#94529a")), .tex = "EPPS21 (nPDF, m_{c}=1.3 GeV)" },
    { .name = "EPPS21Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#3279df")), .tex = "EPPS21 (nPDF, m_{c}=1.5 GeV)" },
    { .name = "nNNPDF30Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#26a4ac")), .tex = "nNNPDF3.0 (nPDF, m_{c}=1.5 GeV)" },
    { .name = "EPPS16_mc13", .color = static_cast<Color_t>(TColor::GetColor("#8a8645")), .tex = "EPPS16 (nPDF, m_{c}=1.3 GeV)" }
  };

}

fonll::DrawSets::DrawSets(const std::string& inputname, double ptmin, double ptmax)
  : valid_(false) {
  auto* inf = TFile::Open(inputname.c_str());
  if (!inf) {
    __XJJLOG << "!! bad input file, abort." << std::endl;
    return;
  }
  std::cout<< Form("pt_%.0f_%.0f", ptmin, ptmax) << std::endl;
  auto* dir = (TDirectory*)inf->Get(Form("pt_%.0f_%.0f", ptmin, ptmax));
  if (!dir) {
    __XJJLOG << "!! no dir, abort." << std::endl;
    return;
  }
  for (const auto& m : configs) {
    auto* mdir = (TDirectory*)dir->Get(m.name.c_str());
    if (!mdir) continue;
    std::map<std::string, TGraphAsymmErrors*> gs;
    for (const std::string t : { "Central", "Scale7", "PDF" }) {
      auto* g = xjjana::getobj<TGraphAsymmErrors>(mdir, "g" + t);
      if (!g) break;
      gs[t] = g; // 
    }
    if (gs.size() == 0) continue;

    models.push_back({ .name = m.name, .color = m.color, .tex = m.tex, .gs = gs });
  }
  if (models.empty()) return;

  valid_ = true;
};

double fonll::DrawSets::ymaximum() {
  double ymax = -1.;
  for (auto& m : models) {
    auto* g = m.gs.at("Scale7");
    for (int j=0; j<g->GetN(); j++) {
      const auto ym = g->GetPointY(j) + g->GetErrorYhigh(j);
      ymax = std::max(ymax, ym);
    }
  }
  return ymax;
}

TLegend* fonll::DrawSets::draw() {
  auto draw_one = [this](const std::string& name, int i, Color_t color) {
    if (models[i].gs.find(name) == models[i].gs.end()) return;
    auto* g = models[i].gs.at(name);
    if (!g) return;
    for (int j=0; j<g->GetN(); j++) {
      const auto y = g->GetPointY(j), ylow = g->GetErrorYlow(j), yhigh = g->GetErrorYhigh(j);
      const auto x = g->GetPointX(j), xlow = g->GetErrorXlow(j), xhigh = g->GetErrorXhigh(j);
      const auto width = (xlow + xhigh)/(models.size() + 1), xleft = x - xlow + (0.5+i)*width;
      xjjroot::drawbox(xleft, y-ylow, xleft+width, y+yhigh, color, 1, 1001, color, 1, 3);
    }
  };

  auto* leg = new TLegend(0.50, 0.76-0.028*1.2*models.size(), 0.71, 0.76);
  xjjroot::setleg(leg, 0.028);
  for (int i=0; i<models.size(); i++) {
    const auto& m = models[i];
    draw_one("Scale7", i, xjjroot::color_alpha(models[i].color, 0.25));
    draw_one("PDF", i, xjjroot::color_alpha(models[i].color, 0.7));
    // draw_one("Central", i, xjjroot::color_alpha_black(models[i].color, 0));
    draw_one("Central", i, models[i].color);
    xjjroot::addentrybystyle(leg, m.tex, (m.gs.find("PDF")!=m.gs.end()) ? "f" : "l",
                             0, 0, 0,
                             models[i].color, 1, 2,
                             xjjroot::color_alpha(models[i].color, 0.7), 1, 1001);
  }
  return leg;
}
