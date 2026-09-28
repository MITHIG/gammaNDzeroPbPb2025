
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
    void draw();
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
    { .name = "CT18ANLO_mc13", .color = static_cast<Color_t>(TColor::GetColor("#ed9418")), .tex = "CT18ANLO (pPDF, m_{c}=1.3 GeV)", },
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

  int nm_p = 0, nm_n = 0;
  for (const auto& m : models) {
    if (xjjc::str_contains(m.tex, "nPDF")) nm_n++;
    if (xjjc::str_contains(m.tex, "pPDF")) nm_p++;
  }

  float x1 = 0.23, y2 = 0.27;
  leg_p = new TLegend(x1, y2 - lheight*nm_p, x1+0.2, y2);
  xjjroot::setleg(leg_p, tsize);
  x1 = 0.55; y2 = 0.72;
  leg_n = new TLegend(x1, y2 - lheight*nm_n, x1+0.2, y2);
  xjjroot::setleg(leg_n, tsize);
  for (int i=0; i<models.size(); i++) {
    const auto& m = models[i];
    auto* leg = xjjc::str_contains(m.tex, "nPDF") ? leg_n : leg_p;
    auto tleg = xjjc::str_replaceall_regex(m.tex, ".PDF, ", "");
    xjjroot::addentrybystyle(leg, tleg, (m.gs.find("PDF") != m.gs.end()) ? "f" : "l",
                             0, 0, 0,
                             models[i].color, 1, 2,
                             xjjroot::color_alpha(models[i].color, 0.7), 1, 1001);
  }
  // xjjroot::autoleg_n_draw(leg_p, 0., 0.20, 0.03, 1.2);
  // xjjroot::autoleg_n_draw(leg_n, 0.50, 0.76, 0.03, 1.2);
  leg_p->Draw();
  leg_n->Draw();
  
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

void fonll::DrawSets::draw() {
  auto draw_one = [this](const std::string& name, int i, Color_t color) {
    if (models[i].gs.find(name) == models[i].gs.end()) return;
    auto* g = models[i].gs.at(name);
    if (!g) return;
    for (int j=0; j<g->GetN(); j++) {
      const auto y = g->GetPointY(j), ylow = g->GetErrorYlow(j), yhigh = g->GetErrorYhigh(j);
      const auto x = g->GetPointX(j), xlow = g->GetErrorXlow(j), xhigh = g->GetErrorXhigh(j);
      const auto width = (xlow + xhigh)/(models.size() + 1), xleft = x - xlow + (0.5+i)*width;
      if (ylow + yhigh > 0) 
        xjjroot::drawbox(xleft, y-ylow, xleft+width, y+yhigh, color, 1, 1001, color, 1, 1);
      else
        xjjroot::drawline(xleft, y, xleft+width, y, color, 1, 2);
    }
  };

  for (int i=0; i<models.size(); i++) {
    const auto& m = models[i];
    draw_one("Scale7", i, xjjroot::color_alpha(models[i].color, 0.25));
    draw_one("PDF", i, xjjroot::color_alpha(models[i].color, 0.7));
    draw_one("Central", i, models[i].color);
  }

  leg_p->Draw();
  leg_n->Draw();

  float x1 = leg_p->GetX1NDC() > 0 ? leg_p->GetX1NDC() : leg_p->GetX1();
  float y2 = leg_p->GetX1NDC() > 0 ? leg_p->GetY2NDC() : leg_p->GetY2();
  // __XJJLOG << leg_p->GetX1() << ", " << leg_p->GetY2() << " | " << leg_p->GetX1NDC() << ", " << leg_p->GetY2NDC() << std::endl;

  auto x1pos = [](TLegend* leg) { return leg->GetX1NDC() > 0 ? leg->GetX1NDC() : leg->GetX1(); };
  auto y2pos = [](TLegend* leg) { return leg->GetX1NDC() > 0 ? leg->GetY2NDC() : leg->GetY2(); };
  
  xjjroot::drawtex(x1pos(leg_p) + 0.005, y2pos(leg_p) + 0.01, "G#gammaA-FONLL + pPDF", tsize, 11);
  xjjroot::drawtex(x1pos(leg_n) + 0.005, y2pos(leg_n) + 0.01, "G#gammaA-FONLL + nPDF", tsize, 11);
}
