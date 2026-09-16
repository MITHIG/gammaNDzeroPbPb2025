#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "draw.h"
#include "util.h"

int macro(const std::vector<std::string>& inputnames, const std::string& outputname,
          const std::string& tags = "", const std::string& title = "",
          const std::string& excepth = "",
          float legdx = 0, float legdy = 0) {
  __XJJLOG << ">> " << inputnames.size() << std::endl;
  const std::vector<Color_t> colors = { kBlack, xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["magenta"], xjjroot::mycolor_middle["cyan"], xjjroot::mycolor_middle["violet"] };
  const float cr = 2./3, ytop = 0.84, lspace = 1.2, tsize = 0.038, tsize_up = 0.038/cr;

  std::vector<std::regex> excepts;
  for (const auto e : xjjc::str_divide_trim(excepth, ",")) {
    std::regex re(e);
    excepts.push_back(re);
  }
  auto* h3_bins = xjjana::getobj<TH3D>(Form("%s::h3_bins", util::parse_input(inputnames.front()).content.c_str()));
  TH2D* h2_bins = nullptr;
  if (!h3_bins) xjjana::getobj<TH2D>(Form("%s::h2_bins", util::parse_input(inputnames.front()).content.c_str()));
  draw::bintex tbins;
  if (h3_bins) {
    tbins.seth(h3_bins, 0, 2);
  } else if (h2_bins) {
    tbins.seth(h2_bins, 0, 1);
  }
  
  std::map<std::string, std::vector<TH1D*>> hs, hratios;
  TLegend* leg = nullptr;
  for (int i=0; i<inputnames.size(); i++) {
    const auto inputp = util::parse_input(inputnames[i]);
    auto* inf = TFile::Open(inputp.content.c_str());
    auto vec_to_map = [&inputp, &excepts](std::vector<TH1D*> vh1s){
      std::map<std::string, TH1D*> ih1s;
      for (auto& h : vh1s) {
        bool skip = false;
        for (auto re : excepts) {
          if (std::regex_match(h->GetName(), re)) {
            skip = true;
            break;
          }
        }
        if (!skip) {
          ih1s[h->GetName()] = h;
          if (!inputp.tag.empty())
            h->SetName(Form("%s_%s", h->GetName(), inputp.tag.c_str()));
        }
      }
      return ih1s;
    };
    auto ih1s = vec_to_map(xjjana::getobj_regexp<TH1D>(inf));
    auto cc = colors[i%colors.size()];
    for (auto& [_, h] : ih1s) {
      xjjroot::setthgrstyle(h, cc, xjjroot::markerlist_solid[i%xjjroot::markerlist_solid.size()], 1.5, cc, 1, 1);
    }
    if (hs.empty()) { // push based on ih1s key
      for (auto& [key, h] : ih1s)
        hs[key].push_back(h);
    } else { // push based on hs key
      for (auto& [key, vh] : hs) {
        if (ih1s.find(key) != ih1s.end()) {
          vh.push_back(ih1s.at(key));
        } else {
          vh.push_back(nullptr);
        }
      }
    }
    if (!leg) {
      leg = new TLegend(0.60+legdx, ytop+legdy-tsize_up*lspace*inputnames.size(), 0.85+legdx, ytop+legdy);
      xjjroot::setleg(leg, tsize_up);
    }
    leg->AddEntry(ih1s.begin()->second, inputp.tex.c_str(), "p");
  }
  // clean which hists to draw
  for (const auto& [key, vh] : hs) {
    int bad = 0;
    for (auto& h : vh)
      if (!h) {
        bad++;
        break;
      }
    if (bad) {
      __XJJLOG << "!! not all input files have: " << key << ", erase it." << std::endl;
      hs.erase(key);
    }
  }
  xjjroot::print_tab(hs, 0);

  if (hs.empty())
    return 2;

  for (auto& [key, vh] : hs) {
    for (auto& h : vh) {
      auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      hratio->Divide(vh.front());
      hratio->GetYaxis()->SetTitle("Ratio");
      // xjjana::sethabsminmax(hratio, 0.71, 1.29);
      hratios[key].push_back(hratio);
    }
    xjjana::sethsabsmin(vh, 1.e-2);
    xjjana::sethsmax(vh, 1.8);
    xjjana::sethsmin(hratios[key], 0.71);
    xjjana::sethsmax(hratios[key], 1.29);
  }
  xjjroot::print_tab(hratios, 0);

  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (auto& [key, vh] : hs) {
    pdf->prepare();
    auto pads = xjjroot::twopads(pdf->getc(), vh.front(), hratios[key].front(), cr, 0.3);
    pads.front()->cd();
    vh.front()->Draw("axis");
    for (auto& h : vh)
      h->Draw("pe1 same");
    leg->Draw();
    auto labels = xjjc::str_divide_trim(tags, ",");
    if (tbins.valid()) {
      const auto index_pt = xjjc::str_extract_index(key, "__pt-");
      labels.push_back(tbins.label_pt(index_pt));
      const auto index_y = xjjc::str_extract_index(key, "__y-");
      labels.push_back(tbins.label_y(index_y));
    }
    xjjroot::drawtexgroup(0.24, ytop-(lspace-1)*tsize_up/2, labels, tsize_up, 13, 42, lspace);
    xjjroot::drawCMS(xjjroot::CMS::internal, title, 1./cr);

    pads.back()->cd();
    hratios[key].front()->Draw("axis");
    xjjroot::drawline(hratios[key].front()->GetXaxis()->GetXmin(), 1, hratios[key].front()->GetXaxis()->GetXmax(), 1,
                      hratios[key].front()->GetLineColor(), 2, 1);
    xjjroot::drawbox(hratios[key].front()->GetXaxis()->GetXmin(), 1.05, hratios[key].front()->GetXaxis()->GetXmax(), 0.95,
                     hratios[key].front()->GetLineColor(), 0.05, 1001);
    for (auto& h : hratios[key])
      if (h != hratios[key].front())
        h->Draw("pe1 same");
    pdf->getc()->cd();
    pdf->write();
  }

  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  __XJJLOG << ">> argc" << argc << std::endl;
  if (argc == 8) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], argv[4], argv[5], std::atof(argv[6]), std::atof(argv[7]));
  }
  if (argc == 6) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], argv[4], argv[5]);
  }
  if (argc == 5) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], argv[3], argv[4]);
  }
  if (argc == 3) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2]);
  }
  return 1;
}
