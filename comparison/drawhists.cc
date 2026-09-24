#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "variables.h"

#include "draw.h"

std::vector<Color_t> colors = { kBlack, xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["blue"], xjjroot::mycolor_middle["green"] };

int macro(const std::string& inputnames, const std::string& output,
          int save_png = 0, int no_ratio = 0,
          const std::string& title = "PbPb (5.36 TeV)")
{
  std::string varname, cut_tex;
  int isDGvar = -1;
  TH3D* h3_bins = nullptr;

  auto inputs = xjjc::str_divide_trim(inputnames, ",");
  
  float pratio = no_ratio ? 1. : 2./3;
  float tsize = 0.033/pratio;
  // auto* leg = new TLegend(0.86-tsize*5, 0.80-tsize*1.1*inputs.size(), 0.86, 0.80);
  auto* leg = new TLegend(0.75-tsize*4, 0.80-tsize*1.1*inputs.size(), 0.75, 0.80);
  xjjroot::setleg(leg, tsize);

  std::map<std::string, std::vector<TH1D*>> h1s, h1s_norm, h1s_lumi;
  std::map<std::string, xjjc::array2D<TH1D*>> h1ys, h1ys_norm, h1ys_lumi;
  
  for (int i=0 ;i<inputs.size(); i++) {
    const xjjroot::thgrstyle style = i ?
      xjjroot::thgrstyle{ colors[i%colors.size()], xjjroot::mstylelist_solid(i), 1.4, colors[i%colors.size()], 1, 1, 0, 0, 0, 1, 1 } :
      xjjroot::thgrstyle{ 0, 0, 0, kBlack, 1, 2, kBlack, 0.05, 1001, 0.3, 0.3 };

    auto* inf = TFile::Open(inputs[i].c_str());
    if (!inf) {
      __XJJLOG << "!! warning: bad input " << inputs[i] << ", skip."<< std::endl;
      continue;
    }
    __XJJLOG << "++ " << inputs[i] << std::endl;
    const auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
    __XJJLOG << ">> info" << std::endl;
    xjjc::print_tab(info, -1);
    auto check_consistency = [&info](const std::string& key, std::string& holder) {
      if (holder.empty()) {
        holder = info.at(key);
      } else if (holder != info.at(key)) {
        __XJJLOG << "!! error: different " << key << ", abort." << std::endl;
        return 2;
      }
      return 0;
    };
    if (check_consistency("varname", varname)) return 2;
    check_consistency("cut_tex", cut_tex);
    //
    const auto lumi = std::atof(info.at("lumi").c_str());
    const auto lumi_is_nevt = std::atoi(info.at("lumi_is_nevt").c_str());
    __XJJLOG << ">> lumi: " << lumi << std::endl;
    
    auto* h3 = xjjana::getobj_regexp_first<TH3D>(inf, "h3_.+_.+_.+");
    if (!h3) {
      __XJJLOG << "!! error: no h3_*_*_*, skip." << std::endl;
      continue;
    }
    if (!h3_bins) {
      h3_bins = (TH3D*)h3->Clone("h3_bins");
      h3_bins->Reset("ICEM");
    }
    if (isDGvar < 0) isDGvar = static_cast<int>(xjjc::str_contains(h3->GetName(), "y_var_pt"));

    for (auto& h : xjjana::getobj_regexp<TH1D>(inf)) {
      h->GetYaxis()->SetTitle("Entries");
      xjjroot::sethempty(h);
      xjjroot::setthgrstyle(h, style);
      h->GetXaxis()->SetNdivisions(505);
 
      auto* hnorm = (TH1D*)h->Clone(Form("%s__norm", h->GetName()));
      hnorm->Scale(1./hnorm->Integral(), "width");
      hnorm->GetYaxis()->SetTitle("Self normalized");

      auto* hlumi = (TH1D*)h->Clone(Form("%s__lumi", h->GetName()));
      hlumi->Scale(1./lumi);
      hlumi->GetYaxis()->SetTitle(lumi_is_nevt ? "Entries / Number of Events" : "Entries / \"Lumi\"");

      auto name = xjjc::str_eraseall(h->GetName(), { "h1_" });
      const auto t_suffix = xjjc::str_extract_regex(name, "(__[a-z]+-[0-9]+)").front();
  
      if (t_suffix.empty()) {
        h1s[name].push_back(h);      
        h1s_norm[name].push_back(hnorm);
        h1s_lumi[name].push_back(hlumi);
      } else {
        name = xjjc::str_eraseall(name, { t_suffix });
        if (h1ys.find(name) == h1ys.end()) {
          h1ys[name].resize(h3->GetNbinsX());
          h1ys_norm[name].resize(h3->GetNbinsX());
          h1ys_lumi[name].resize(h3->GetNbinsX());
        }
        const int j = std::atoi(xjjc::str_extract_regex(t_suffix, "([0-9]+)").front().c_str());
        if (j < 0 || j >= h1ys[name].size()) {
          __XJJLOG << "!! cannot extract index j : " << j << ", skip." << std::endl;
          continue;
        }
        h1ys[name][j].push_back(h);
        h1ys_norm[name][j].push_back(hnorm);
        h1ys_lumi[name][j].push_back(hlumi);
      }
    }
    xjjroot::addentrybystyle(leg, info.at("input_tex").c_str(), style.fstyle>0?"l":"pe", style);
  }
  xjjroot::print_tab(h1s, 0);
  for (const auto& [_, hs] : h1ys)
    xjjroot::print_tab(hs, 0);

  auto it_var = std::find_if(vars.begin(), vars.end(), [&varname](const xjjana::variable& v) {
    return v.varname == varname;
  });
  if (it_var == vars.end()) {
    __XJJLOG << "!! bad varname: " << varname << std::endl;
    return 2;
  }
  const auto& the_var = *it_var;

  draw::bintex tbin(h3_bins, 0, 2);
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + output + "/" + varname + ".pdf");
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), {{ "figspdf", "figs" }, { "/"+varname, "_"+varname }, { ".pdf", "" }});

  auto draw_hist_list = [&inputs, &pdf, &pratio, &tsize, &cut_tex, &title, &leg](std::vector<TH1D*> result, int logy) {
    if (logy) xjjana::sethsnonzeromin(result, 0.5);
    else xjjana::sethsmin(result, 0.1);
    xjjana::sethsmax(result, logy ? 1.e2 : 1.6);

    std::vector<TH1D*> result_ratio;
    for (auto& h : result) {
      auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      hratio->Divide(result.front());
      hratio->SetFillStyle(0);
      hratio->GetYaxis()->SetTitle("Ratio");
      if (result_ratio.size() == 0) {
        xjjana::sethunivalue(hratio, 1., 0.);
        hratio->SetLineStyle(2);
      }
      result_ratio.push_back(hratio);
    }
    xjjana::sethsnonzeromin(result_ratio, 0.95);
    xjjana::sethsmax(result_ratio, 1.05);    

    auto pads = xjjroot::twopads(pdf->getc(), result.front(), result_ratio.front(), pratio);
    pads[0]->SetLogy(logy);
    pads[0]->cd();
    for (int i=0; i<result.size(); i++)
      result[i]->Draw(i?"pe1 same":"hist same");
    leg->Draw();
    xjjroot::drawCMS(xjjroot::CMS::internal, title.c_str(), 1./pratio);
    xjjroot::drawtexgroup(0.23, 0.86, xjjc::str_divide_trim(cut_tex, "%%"), tsize, 13, 42, 1.1);
    if (pads.size() > 1) {
      pads[1]->cd();
      for (int i=0; i<result_ratio.size(); i++)
        result_ratio[i]->Draw(i?"pe1 same":"hist same");
    }
    pads[0]->cd();
    return pads;
  };

  for (const auto scale : std::map<int, std::string>{
      { 0, "Linear scale" },
      { 1, "Log scale" },
    }) {
    pdf->prepare();
    draw_hist_list(h1s["var"], scale.first);
    if (isDGvar) {
      xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*xjjc::str_divide_trim(cut_tex, "%%").size(), {
          tbin.label_y(), tbin.label_pt()
        }, tsize, 13);
    }
    pdf->getc()->cd();
    if (save_png && (scale.first == the_var.logy)) {
      pdf->write(name_png + ".pdf");
    } else pdf->write();

    if (isDGvar) {
      for (int j=0; j<h1ys["var"].size(); j++) {
        pdf->prepare();
        draw_hist_list(h1ys["var"][j], scale.first);
        xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
            tbin.label_y(j), tbin.label_pt()
          }, tsize, 13);
        pdf->getc()->cd();
        pdf->write();
      }
    }

    pdf->draw_cover({ "#bf{Self-normalized}", scale.second });
  
    pdf->prepare();
    draw_hist_list(h1s_norm["var"], scale.first);
    if (isDGvar) {
      xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
          tbin.label_y(), tbin.label_pt()
        }, tsize, 13);
    }
    pdf->getc()->cd();
    if (save_png && (scale.first == the_var.logy)) {
      pdf->write(name_png + "_norm.pdf");
    } else pdf->write();

    if (isDGvar) {
      for (int j=0; j<h1ys_norm["var"].size(); j++) {
        pdf->prepare();
        draw_hist_list(h1ys_norm["var"][j], scale.first);
        xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
            tbin.label_y(j), tbin.label_pt()
          }, tsize, 13);
        pdf->getc()->cd();
        pdf->write();
      }
    }

    pdf->draw_cover({ "#bf{\"Lumi-normalized\"}", scale.second });
  
    pdf->prepare();
    draw_hist_list(h1s_lumi["var"], scale.first);
    if (isDGvar) {
      xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
          tbin.label_y(), tbin.label_pt()
        }, tsize, 13);
    }
    pdf->getc()->cd();
    if (save_png && (scale.first == the_var.logy)) {
      pdf->write(name_png + "_lumi.pdf");
    } else pdf->write();

    if (isDGvar) {
      for (int j=0; j<h1ys_lumi["var"].size(); j++) {
        pdf->prepare();
        draw_hist_list(h1ys_lumi["var"][j], scale.first);
        xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
            tbin.label_y(j), tbin.label_pt()
          }, tsize, 13);
        pdf->getc()->cd();
        pdf->write();
      }
    }

    // pdf->prepare();
    // draw_hist_list(h1s_lumi["var"], scale.first);
    // if (isDGvar) {
    //   xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
    //       tbin.label_y(), tbin.label_pt()
    //     }, tsize, 13);
    // }
    // pdf->getc()->cd();
    // if (save_png && (scale.first == the_var.logy)) {
    //   pdf->write(name_png + ".pdf");
    // } else pdf->write();

    // if (isDGvar) {
    //   for (int j=0; j<h1ys["var"].size(); j++) {
    //     pdf->prepare();
    //     draw_hist_list(h1ys["var"][j], scale.first);
    //     xjjroot::drawtexgroup(0.23, 0.86-tsize*1.1*2, {
    //         tbin.label_y(j), tbin.label_pt()
    //       }, tsize, 13);
    //     pdf->getc()->cd();
    //     pdf->write();
    //   }
    // }    
  }
  
  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(argv[1], argv[2], std::atoi(argv[3]), std::atoi(argv[4]));
  }
  if (argc == 4) {
    return macro(argv[1], argv[2], atoi(argv[3]));
  }
  return 1;
}
