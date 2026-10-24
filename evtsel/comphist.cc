#include <TFile.h>
#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjstruct.h"
#include "xjjmypdf.h"
#include "util.h"
#include "draw.h"

#include "variables.h"

int macro(const std::vector<std::string>& inputnames, const std::string& outputname,
          int no_ratio = 1, const std::string& drawh = ".+",
          float legdx = 0, float legdy = 0, int save_png = 0) {
  if (inputnames.empty()) {
    __XJJLOG << "!! bad inputs, abort." << std::endl;
    return 2;
  }
  const float pratio = (no_ratio ? 1. : 2./3), ytop = 0.87, lspace = 1.2, tsize = 0.038, tsize_up = 0.038/pratio;
  auto y_pos = [&ytop, &lspace, &tsize_up](float n = 0.) -> float {
    return ytop - n*lspace*tsize_up;
  };
  std::vector<std::regex> draws;
  for (const auto e : xjjc::str_divide_trim(drawh, ",")) {
    std::regex re(e);
    draws.push_back(re);
  }
  auto* inf0 = xjjroot::openfile(xjjroot::parse_input(inputnames.front()).content);
  if (!inf0) {
    __XJJLOG << "!! bad reference file, abort." << std::endl;
  }
  const auto info = util::read_info(inf0);
  auto this_is = [&info](const std::string& item) {
    return static_cast<bool>(std::atoi(info.at(item).c_str()));
  };
  auto* h3_bins = xjjana::getobj<TH3D>(inf0, "h3_bins");
  if (!h3_bins) {
    __XJJLOG << "!! cannot find h3_bins, abort." << std::endl;
    xjjroot::closefile(inf0);
    return 2;
  }
  draw::bintex tbins(h3_bins, 0, 2);

  auto varname_key = [](const std::string& tkey) -> std::string {
    auto name = xjjc::str_extract_regex(tkey, "h1_(.+)_pt-.+_y-.+");
    if (name.empty()) return "";
    return name.front();
  };

  std::map<std::string, std::vector<TH1D*>> hs, hs_norm, hratios, hratios_norm;
  TLegend* leg = nullptr;
  for (int i=0; i<inputnames.size(); i++) {
    const auto inputp = xjjroot::parse_input(inputnames[i]);
    auto* inf = xjjroot::openfile(inputp.content);
    if (!inf) {
      __XJJLOG << "!! bad input, abort." << std::endl;
      return 2;
    }
    auto vec_to_map = [&inputp, &draws, &varname_key](std::vector<TH1D*> vh1s){
      std::map<std::string, TH1D*> ih1s;
      for (auto& h : vh1s) {
        bool keep = false;
        for (auto re : draws) {
          if (std::regex_match(h->GetName(), re)) {
            keep = true;
            break;
          }
        }
        if (!keep) continue;
        auto it_var = std::find_if(variables.begin(), variables.end(), [&h, &varname_key](const global::variable& v) {
          return v.varname == varname_key(h->GetName());
        });
        if (it_var == variables.end()) {
          __XJJLOG << "?? unknown varname from: " << h->GetName() << ", skip." << std::endl;
          continue;
        }

        ih1s[h->GetName()] = h;
        // if (!inputp.tag.empty())
        //   h->SetName(Form("%s_%s", h->GetName(), inputp.tag.c_str()));
      }
      return ih1s;
    };
    auto ih1s = vec_to_map(xjjana::getobj_regexp<TH1D>(inf, "h1_.+_pt-.+_y-.+", "", false));
    const auto cc = draw::icolor(i);
    for (auto& [_, h] : ih1s) {
      xjjroot::setthgrstyle(h, cc, xjjroot::mstylelist_solid(i), 1.5, cc, 1, 1);
    }
    if (hs.empty()) { // push based on ih1s key
      for (auto& [key, h] : ih1s) {
        hs[key].push_back(h);
      }
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
      const size_t nlabelr = this_is("has_dreq") ? 5 : 1;
      leg = new TLegend(0.60+legdx, y_pos(nlabelr + inputnames.size())+legdy, 0.85+legdx, y_pos(nlabelr)+legdy);
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

  //
  for (auto& [key, vh] : hs) {
    auto it_var = std::find_if(variables.begin(), variables.end(), [&key, &varname_key](const global::variable& v) {
      return v.varname == varname_key(key);
    });
    for (auto& h : vh) {
      // std::cout << h->GetName() << " -> " << (*it_var).varname << std::endl;      
      auto* h_norm = (TH1D*)h->Clone(Form("%s_norm", h->GetName()));
      h_norm->Scale(1./h_norm->Integral(), "width");
      h_norm->GetYaxis()->SetTitle("Self-normalized");
      hs_norm[key].push_back(h_norm);
      
      auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      hratio->Divide(vh.front());
      hratio->GetYaxis()->SetTitle("Ratio");
      // xjjana::sethabsminmax(hratio, 0.71, 1.29);
      hratios[key].push_back(hratio);

      auto* hratio_norm = (TH1D*)h->Clone(Form("%s_norm_ratio", h->GetName()));
      hratio_norm->Divide(hs_norm[key].front());
      hratio_norm->GetYaxis()->SetTitle("Ratio");
      // xjjana::sethabsminmax(hratio_norm, 0.71, 1.29);
      hratios_norm[key].push_back(hratio_norm);
    }
    if ((*it_var).logy) {
      xjjana::sethsnonzeromin(hs_norm[key], 0.9);
      xjjana::sethsnonzeromin(hs[key], 0.9);
      xjjana::sethsmax(hs_norm[key], 20);
      xjjana::sethsmax(hs[key], 20);
    } else {
      xjjana::sethsabsmin(hs_norm[key], no_ratio ? 0 : 1.e-3);
      xjjana::sethsabsmin(hs[key], no_ratio ? 0 : 1.e-3);
      xjjana::sethsmax(hs_norm[key], 1.5);
      xjjana::sethsmax(hs[key], 1.5);
    }
    xjjana::sethsmin(hratios[key], 0.71);
    xjjana::sethsmax(hratios[key], 1.29);
    xjjana::sethsmin(hratios_norm[key], 0.71);
    xjjana::sethsmax(hratios_norm[key], 1.29);
  }

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  const auto png_name = draw::png_name(pdf);

  auto draw_hists = [&pdf, &pratio, &tbins, &info](const std::vector<TH1D*>& vh, const std::vector<TH1D*>& vhratios, const std::string& key, int logy, int has_dreq) {
    auto pads = xjjroot::twopads(pdf->getc(), vh.front(), vhratios.front(), pratio, 0.3);
    pads.front()->SetLogy(logy);
    pads.front()->cd(); //
    vh.front()->Draw("axis");
    for (auto& h : vh)
      h->Draw("histe1 same");
    xjjroot::drawCMS(xjjroot::CMS::internal, draw::title_ph, 1./pratio);

    if (pads.size() > 1) {
      pads.back()->cd(); //
      vhratios.front()->Draw("axis");
      xjjroot::drawline(vhratios.front()->GetXaxis()->GetXmin(), 1, vhratios.front()->GetXaxis()->GetXmax(), 1,
                        vhratios.front()->GetLineColor(), 2, 1);
      xjjroot::drawbox(vhratios.front()->GetXaxis()->GetXmin(), 1.05, vhratios.front()->GetXaxis()->GetXmax(), 0.95,
                       vhratios.front()->GetLineColor(), 0.05, 1001);
      for (auto& h : vhratios)
        if (h != vhratios.front())
          h->Draw("pe1 same");
    }
    std::vector<std::string> labels_right;
    if (has_dreq) {
      if (tbins.valid()) {
        const auto index_pt = xjjc::str_extract_index(key, "_pt-");
        labels_right.push_back(tbins.label_pt(index_pt));
        const auto index_y = xjjc::str_extract_index(key, "_y-");
        if (index_y >= 0)
          labels_right.push_back(tbins.label_y(index_y));
      }
      if (!labels_right.empty()) {
        labels_right.insert(labels_right.begin(), std::string("At least a D cand in"));
        labels_right.push_back(info.at("dcut_tex"));
      }
    }
    pads.front()->cd(); //
    return labels_right;
  };

  auto treat_tex = [&this_is, &info](const std::string& tex_tochange) {
    return xjjc::str_replaceall(tex_tochange, {
        { "1nXOR", this_is("is_gammaN") ? "Xn0n" : "0nXn" },
      });
  };

  
  for (auto& [key, _] : hs) {
    auto it_var = std::find_if(variables.begin(), variables.end(), [&key, &varname_key](const global::variable& v) {
      return v.varname == varname_key(key);
    });
    auto varnote = xjjc::str_divide_trim(treat_tex((*it_var).note), ",");
    auto texs_left = xjjc::vec_append_into(xjjc::str_divide_trim(info.at("str_evtsel"), ","), varnote);
       
    pdf->prepare();
    auto texs_right = draw_hists(hs_norm[key], hratios_norm[key], key, (*it_var).logy, this_is("has_dreq"));
    leg->Draw();
    xjjroot::drawtexgroup(0.24, y_pos(), texs_left, tsize_up, 13, 42, lspace);
    xjjroot::drawtexgroup(0.55, y_pos(1), texs_right, tsize_up, 13, 42, lspace);
    pdf->getc()->cd();
    pdf->write(png_name + "/" + key + "_norm.pdf", save_png ? "" : "X");

    // pdf->prepare();
    // draw_hists(hs[key], hratios[key], key, (*it_var).logy, this_is("has_dreq"));
    // leg->Draw();
    // xjjroot::drawtexgroup(0.24, y_pos(), texs_left, tsize_up, 13, 42, lspace);
    // xjjroot::drawtexgroup(0.55, y_pos(1), texs_right, tsize_up, 13, 42, lspace);
    // pdf->getc()->cd();
    // pdf->write(png_name + "/" + key + ".pdf", save_png ? "" : "X");
  }

  for (int j=0; j<tbins.npt(); j++) {
    for (int i=0; i<tbins.ny(); i++) {
      auto key_pt_y = [](const std::string& var) -> std::string {
        return Form("");
      };
    }
  }

  pdf->close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(xjjc::str_divide_trim(argv[1], ","), argv[2], std::atoi(argv[3]), argv[4]);
  }
  return 1;
}
