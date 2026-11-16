
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
    void draw_base(const std::vector<std::string>& dmodels, float ltsize);
    void draw();
    void draw_mc15();
    void draw_mc13();
    bool valid() const { return valid_; }
    size_t n() const { return models.size(); }
    double ymaximum();
  private:
    bool valid_;
    std::vector<Model> models;
    TLegend *leg_n, *leg_p;
    float tsize = 0.031, lspace = 1.2, lheight = tsize*lspace;
  };

  const std::vector<Model> configs = {
    { .name = "CT18ANLO_mc13", .color = static_cast<Color_t>(TColor::GetColor("#015AA6")), .tex = "CT18ANLO (pPDF, m_{c}=1.3 GeV)", },
    { .name = "CT18ANLO_mc15", .color = static_cast<Color_t>(TColor::GetColor("#458984")), .tex = "CT18ANLO (pPDF, m_{c}=1.5 GeV)" },
    { .name = "EPPS21Pb_mc13", .color = static_cast<Color_t>(TColor::GetColor("#F2B445")), .tex = "EPPS21 (nPDF, m_{c}=1.3 GeV)" },
    { .name = "EPPS21Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#A6442E")), .tex = "EPPS21 (nPDF, m_{c}=1.5 GeV)" },
    { .name = "nNNPDF30Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#D95829")), .tex = "nNNPDF3.0 (nPDF, m_{c}=1.5 GeV)" },
    { .name = "EPPS16_mc13", .color = static_cast<Color_t>(TColor::GetColor("#734886")), .tex = "EPPS16 (nPDF, m_{c}=1.3 GeV)" }

    // { .name = "CT18ANLO_mc13", .color = static_cast<Color_t>(TColor::GetColor("#ed9418")), .tex = "CT18ANLO (pPDF, m_{c}=1.3 GeV)", },
    // { .name = "CT18ANLO_mc15", .color = static_cast<Color_t>(TColor::GetColor("#e44760")), .tex = "CT18ANLO (pPDF, m_{c}=1.5 GeV)" },
    // { .name = "EPPS21Pb_mc13", .color = static_cast<Color_t>(TColor::GetColor("#94529a")), .tex = "EPPS21 (nPDF, m_{c}=1.3 GeV)" },
    // { .name = "EPPS21Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#3279df")), .tex = "EPPS21 (nPDF, m_{c}=1.5 GeV)" },
    // { .name = "nNNPDF30Pb_mc15", .color = static_cast<Color_t>(TColor::GetColor("#26a4ac")), .tex = "nNNPDF3.0 (nPDF, m_{c}=1.5 GeV)" },
    // { .name = "EPPS16_mc13", .color = static_cast<Color_t>(TColor::GetColor("#8a8645")), .tex = "EPPS16 (nPDF, m_{c}=1.3 GeV)" }
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

void fonll::DrawSets::draw_base(const std::vector<std::string>& draws, float ltsize) {
  auto draw_one = [this](const std::string& name, int i, Color_t color, int nmodels) {
    if (models[i].gs.find(name) == models[i].gs.end()) return;
    auto* g = models[i].gs.at(name);
    if (!g) return;
    for (int j=0; j<g->GetN(); j++) {
      const auto y = g->GetPointY(j), eylow = g->GetErrorYlow(j), eyhigh = g->GetErrorYhigh(j);
      const auto x = g->GetPointX(j), exlow = g->GetErrorXlow(j), exhigh = g->GetErrorXhigh(j);
      const auto width = (exlow + exhigh)/(nmodels + 1), xleft = x - exlow + (0.5+i)*width;
      if (eylow + eyhigh > 0) 
        xjjroot::drawbox(xleft, y-eylow, xleft+width, y+eyhigh, color, 1, 1001, color, 1, 1);
      else
        xjjroot::drawline(xleft, y, xleft+width, y, color, 1, 2);
    }
  };

  int im = 0, nm_p = 0, nm_n = 0;
  int has_mc13 = 0, has_mc15 = 0;
  for (const auto& m : models) {
    if (std::find(draws.begin(), draws.end(), m.name) == draws.end()) continue;
    if (xjjc::str_contains(m.tex, "nPDF")) nm_n++;
    if (xjjc::str_contains(m.tex, "pPDF")) nm_p++;
    if (xjjc::str_contains(m.name, "mc13")) has_mc13++;
    if (xjjc::str_contains(m.name, "mc15")) has_mc15++;

    draw_one("Scale7", im, xjjroot::color_alpha(m.color, 0.15), draws.size());
    draw_one("PDF", im, xjjroot::color_alpha(m.color, 0.7), draws.size());
    draw_one("Central", im, m.color, draws.size());

    im++;
  }

  const std::string tcommon = has_mc13 && !has_mc15 ? "(m_{c}=1.3 GeV)" : (has_mc15 && !has_mc13 ? "(m_{c}=1.5 GeV)" : "");

  float x1 = 0.23, y2 = 0.27;
  leg_p = new TLegend(x1, y2 - lheight*nm_p, x1+0.2, y2);
  xjjroot::setleg(leg_p, ltsize);
  x1 = 0.55; y2 = 0.72;
  leg_n = new TLegend(x1, y2 - lheight*nm_n, x1+0.2, y2);
  xjjroot::setleg(leg_n, ltsize);
  for (const auto& m : models) {
    if (std::find(draws.begin(), draws.end(), m.name) == draws.end()) continue;
    auto* leg = xjjc::str_contains(m.tex, "nPDF") ? leg_n : leg_p;
    auto tleg = xjjc::str_replaceall_regex(m.tex, ".PDF, ", "");
    if (!tcommon.empty()) tleg = xjjc::str_eraseall(tleg, { " " + tcommon });
    xjjroot::addentrybystyle(leg, tleg, (m.gs.find("PDF") != m.gs.end()) ? "f" : "l",
                             0, 0, 0,
                             m.color, 1, 2,
                             xjjroot::color_alpha(m.color, 0.7), 1, 1001);
  }
  // xjjroot::autoleg_n_draw(leg_p, 0., 0.20, 0.03, 1.2);
  // xjjroot::autoleg_n_draw(leg_n, 0.50, 0.76, 0.03, 1.2);
  leg_p->Draw();
  leg_n->Draw();
  
  // float x1 = leg_p->GetX1NDC() > 0 ? leg_p->GetX1NDC() : leg_p->GetX1();
  // float y2 = leg_p->GetX1NDC() > 0 ? leg_p->GetY2NDC() : leg_p->GetY2();
  // __XJJLOG << leg_p->GetX1() << ", " << leg_p->GetY2() << " | " << leg_p->GetX1NDC() << ", " << leg_p->GetY2NDC() << std::endl;

  auto x1pos = [](TLegend* leg) { return leg->GetX1NDC() > 0 ? leg->GetX1NDC() : leg->GetX1(); };
  auto y2pos = [](TLegend* leg) { return leg->GetX1NDC() > 0 ? leg->GetY2NDC() : leg->GetY2(); };
  auto y1pos = [](TLegend* leg) { return leg->GetX1NDC() > 0 ? leg->GetY1NDC() : leg->GetY1(); };

  if (nm_p) {
    std::vector<std::string> texs = { "G#gammaA-FONLL + pPDF" };
    xjjroot::drawtexgroup(x1pos(leg_p) + 0.005, y2pos(leg_p) + 0.01, texs, ltsize, 11);
  }
  if (nm_n) {
    std::vector<std::string> texs = { "G#gammaA-FONLL + nPDF" };
    xjjroot::drawtexgroup(x1pos(leg_n) + 0.005, y2pos(leg_n) + 0.01, texs, ltsize, 11);
  }
  xjjroot::drawtex(x1pos(leg_n) + 0.005, y1pos(leg_n) - 0.01, tcommon.c_str(), ltsize, 13);
}

void fonll::DrawSets::draw() {
  std::vector<std::string> dmodels;
  for (auto& m : models) {
    dmodels.push_back(m.name);
  }
  draw_base(dmodels, tsize);
}

void fonll::DrawSets::draw_mc15() {
  draw_base({
      "CT18ANLO_mc15", "EPPS21Pb_mc15", "nNNPDF30Pb_mc15"
    }, tsize*1.1);
}

void fonll::DrawSets::draw_mc13() {
  draw_base({
      "CT18ANLO_mc13", "EPPS21Pb_mc13"
    }, tsize*1.1);
}
