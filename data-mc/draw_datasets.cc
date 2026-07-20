#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/util.h"
#include "../include/draw.h"
#define __VARIABLES_ROOSPLOT__
#include "variables.h"

class MassRange {
public:
  MassRange(TTree* t) { initialize(t); }
  double mass_signal_low;
  double mass_signal_high;
  double deltamass_sideband_low;
  double deltamass_sideband_high;

  double sideband_scale() { return (mass_signal_high-mass_signal_low)*0.5 / (deltamass_sideband_high-deltamass_sideband_low); }
  std::string str_signal() { return tex_signal; }
  std::string str_sideband() { return tex_sideband; }
private:
  std::string tex_signal, tex_sideband;
  void initialize(TTree* t);
};

struct hName {
  std::string varname = "";
  std::string type = "";
  int index_y = -1;
};
hName parse_hname(std::string hname, bool verbose = false);

// double scale_sideband = 0.5;
int macro(std::string inputname, std::string outputname) {
  const auto inputfile = util::parse_input(inputname).content;
  auto* inf = TFile::Open(inputfile.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! bad file: " << inputfile << ", abort." << std::endl;
    return 2;
  }
  std::map<std::string, xjjc::info> infos;
  for (std::string tr : { "data/info", "template/info" }) {
    auto info = xjjana::getval_regexp((TTree*)inf->Get(tr.c_str()));
    infos[xjjc::str_eraseall(tr, "/info")] = info;
  }
  for (auto& [key, info] : infos) {
    __XJJLOG << "++ infos [" << key << "]" << std::endl;
    xjjc::print_tab<std::string, std::string>(info, -1);
  }

  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
  draw::bintex btex(h3_bins, 0, 2);
  auto ny = h3_bins->GetXaxis()->GetNbins();

  //
  std::vector<MassRange*> mranges(ny, nullptr);
  for (int i=0; i<ny; i++) {
    auto* dir = (TDirectory*)inf->GetDirectory(Form("dir__y-%d", i));
    if (!dir) {
      __XJJLOG << "!! bad directory, skip." << std::endl;
      continue;
    }
    mranges[i] = new MassRange((TTree*)dir->Get("mass_range"));
  }

  // 
  auto vh1s_all = xjjana::getobj_regexp_recur<TH1D>(inf, ".+", "", false);
  auto collect_h1ys = [&vh1s_all, &ny](std::string type, xjjroot::thgrstyle style = {}, bool verbose = false) {
    std::map<std::string, std::vector<TH1D*>> h1ys;
    for (auto& h : vh1s_all) {
      auto hpar = parse_hname(h->GetName());
      if (hpar.type != type) continue;
      if (hpar.index_y < 0 || hpar.index_y >= ny) {
        __XJJLOG << "!! bad index_y : " << hpar.index_y << std::endl;
        h1ys.clear();
        return h1ys;
      }
      if (h1ys.find(hpar.varname) == h1ys.end()) {
        h1ys[hpar.varname].resize(ny, nullptr);
      }
      xjjroot::sethempty(h, 0, 0.4);
      xjjroot::setthgrstyle(h, style);
      h1ys[hpar.varname][hpar.index_y] = h;
    }
    if (h1ys.empty()) {
      __XJJLOG << "!! no matching histograms : " << type << std::endl;
    } else if (verbose) {
      xjjroot::print_tab(h1ys, 0);
    }
    return h1ys;
  };

  auto make_h1ys_norm = [](std::map<std::string, std::vector<TH1D*>>& h1ys, xjjroot::thgrstyle style = {}) {
    std::map<std::string, std::vector<TH1D*>> h1ys_norm;
    for (auto& [name, vh] : h1ys) {
      // h1ys[hpar.varname].resize(ny, nullptr);
      for (auto& h : vh) {
        auto* h_norm = h ? (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "__y", "_norm__y").c_str()) : nullptr; // !! without y
        if (h_norm) {
          h_norm->Scale(1./h->Integral(), "width");
          h_norm->GetYaxis()->SetTitle("Self-normalized / bin width");
          xjjroot::setthgrstyle(h_norm, style);
        }
        h1ys_norm[name].push_back(h_norm);
      }
    }
    return h1ys_norm;
  };
  
  auto h1ys_data_main = collect_h1ys("data-main", { kGray+3, 24, 1.5, kGray+3, 1, 1, kGray+3, 0.1, 1001 }, true);
  auto h1ys_norm_data_main = make_h1ys_norm(h1ys_data_main, { kGray, 20, 1.5, kGray, 1, 1, kGray, 0.2, 1001 });
  if (h1ys_norm_data_main.empty()) return 2;
  auto h1ys_data_sigswap = collect_h1ys("data-sigswap", { xjjroot::mycolor_middle["blue"], 25, 1.5, xjjroot::mycolor_middle["blue"], 1, 1 }, true);
  auto h1ys_norm_data_sigswap = make_h1ys_norm(h1ys_data_sigswap);
  if (h1ys_norm_data_sigswap.empty()) return 2;
  // auto h1ys_data_sideband = collect_h1ys("data-sideband", { xjjroot::mycolor_middle["cyan"], 24, 1.5, xjjroot::mycolor_middle["cyan"], 1, 1, xjjroot::mycolor_middle["cyan"], 1, 3004 }, true);
  auto h1ys_data_sideband = collect_h1ys("data-sideband", { kGray, 20, 1.5, kGray, 1, 1, kGray, 0.2, 1001 }, true);
  if (h1ys_data_sideband.empty()) return 2;
  auto h1ys_norm_data_sideband = make_h1ys_norm(h1ys_data_sideband);
  auto h1ys_mc_match = collect_h1ys("mc-match", { xjjroot::mycolor_middle["red"], 21, 0, xjjroot::mycolor_middle["red"], 1, 1, xjjroot::mycolor_middle["red"], 1, 3004 }, true);
  auto h1ys_norm_mc_match = make_h1ys_norm(h1ys_mc_match);
  
  auto make_h1ys_scale = [&mranges](std::map<std::string, std::vector<TH1D*>>& h1ys, xjjroot::thgrstyle style = {},
                                    double scale_force = -1) {
    std::map<std::string, std::vector<TH1D*>> h1ys_scale;
    for (auto& [name, vh] : h1ys) {
      for (auto& h : vh) {
        auto* h_scale = h ? (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "__y", "_scale__y").c_str()) : nullptr;
        if (h_scale) {
          auto scale = scale_force > 0 ? scale_force : mranges[parse_hname(h->GetName()).index_y]->sideband_scale();
          // __XJJLOG << ">> scale [" << parse_hname(h->GetName()).index_y << " ] : " << scale << std::endl;
          h_scale->Scale(scale);
          xjjroot::setthgrstyle(h_scale, style);
        }
        h1ys_scale[name].push_back(h_scale);
      }
    }
    return h1ys_scale;
  };  
  auto make_h1ys_sub = [](const std::map<std::string, std::vector<TH1D*>>& h1ys_main,
                          const std::map<std::string, std::vector<TH1D*>>& h1ys_sideband_scaled,
                          xjjroot::thgrstyle style = {}) {
    std::map<std::string, std::vector<TH1D*>> h1ys_sub;
    for (auto& [name, vh] : h1ys_main) {
      for (auto& h : vh) {
        auto* h_sub = h ? (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "-main", "-sub").c_str()) : nullptr;
        if (h_sub) {
          auto index_y = parse_hname(h->GetName()).index_y;
          if (h1ys_sideband_scaled.at(name)[index_y]) {
            h_sub->Add(h1ys_sideband_scaled.at(name)[index_y], -1);
            xjjroot::setthgrstyle(h_sub, style);
          }
        }
        h1ys_sub[name].push_back(h_sub);
      }
    }
    return h1ys_sub;
  };
  
  auto h1ys_data_sideband_scaled = make_h1ys_scale(h1ys_data_sideband, { xjjroot::mycolor_middle["cyan"], 24, 1.5, xjjroot::mycolor_middle["cyan"], 1, 1, xjjroot::mycolor_middle["cyan"], 1, 3004 });
  auto h1ys_data_sub = make_h1ys_sub(h1ys_data_main, h1ys_data_sideband_scaled,
                                     { xjjroot::mycolor_middle["blue"], 20, 1.7, xjjroot::mycolor_middle["blue"], 1, 1 });
  auto h1ys_norm_data_sub = make_h1ys_norm(h1ys_data_sub);

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  
  for (auto& [name, _] : h1ys_data_main) {
    pdf->draw_cover({ "#bf{" + name + "}" }, 0.05);
    const auto the_var = var_by_name(xjjc::str_eraseall(name, "_norm"));
    if (the_var.varname.empty()) continue;

    // preparations
    auto set_hsminmax = [&the_var](std::vector<TH1D*> hlist) {
      if (the_var.logy) {
        xjjana::sethsnonzeromin(hlist, 0.5);
        xjjana::sethsmax(hlist, 50.);
      } else {
        xjjana::sethsmin(hlist, 0.);
        xjjana::sethsmax(hlist, 1.5);
      }
    };
    auto init_leg = [](int nrow) {
      auto* leg = new TLegend(0.52, 0.86-0.040*nrow, 0.70, 0.86);
      xjjroot::setleg(leg, 0.035);
      return leg;
    };

    for (int i=0; i<ny; i++) {
      // norm
      auto* leg_norm = init_leg(5);
      leg_norm->AddEntry(h1ys_norm_data_sideband.at(name)[i], "Data in sideband", "f");
      leg_norm->AddEntry((TObject*)0, mranges[i]->str_sideband().c_str(), "");
      leg_norm->AddEntry(h1ys_norm_mc_match.at(name)[i], "Gen-matched MC", "f");
      leg_norm->AddEntry(h1ys_norm_data_sigswap.at(name)[i], "sPlot signal + swap", "pe");
      leg_norm->AddEntry(h1ys_norm_data_sub.at(name)[i], "Sideband subtracted", "pe");

      set_hsminmax({ h1ys_norm_data_sigswap.at(name)[i], h1ys_norm_mc_match.at(name)[i], h1ys_norm_data_sideband.at(name)[i], h1ys_norm_data_sub.at(name)[i] });
      pdf->prepare();
      pdf->getc()->SetLogy(0);
      if (the_var.logy)
        pdf->getc()->SetLogy();
      h1ys_norm_data_sideband.at(name)[i]->Draw("hist");
      h1ys_norm_mc_match.at(name)[i]->Draw("hist e1 same");
      h1ys_norm_data_sigswap.at(name)[i]->Draw("pe1 same");
      h1ys_norm_data_sub.at(name)[i]->Draw("pe1 same");
      leg_norm->Draw();
      xjjroot::drawCMS(xjjroot::CMS::internal, infos["data"]["input_tex"] + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt(), infos["data"]["cut_tex"] }, 0.038, 13);
      pdf->write();
    }
    
    if (name == "Dy") { // overlay Dy
      auto* leg_norm_Dy = init_leg(4);
      leg_norm_Dy->AddEntry(h1ys_norm_data_sideband.at(name).front(), "Data in sideband", "f");
      leg_norm_Dy->AddEntry(h1ys_norm_mc_match.at(name).front(), "Gen-matched MC", "f");
      leg_norm_Dy->AddEntry(h1ys_norm_data_sigswap.at(name).front(), "sPlot signal + swap", "pe");
      leg_norm_Dy->AddEntry(h1ys_norm_data_sub.at(name).front(), "Sideband subtracted", "pe");
      pdf->prepare();
      h1ys_norm_data_sideband.at(name).back()->Draw("axis");
      for (int i=0; i<ny; i++) {
        h1ys_norm_data_sideband.at(name)[i]->Draw("hist same");
        h1ys_norm_mc_match.at(name)[i]->Draw("hist e1 same");
        h1ys_norm_data_sigswap.at(name)[i]->Draw("pe1 same");
        h1ys_norm_data_sub.at(name)[i]->Draw("pe1 same");
      }
      leg_norm_Dy->Draw();
      xjjroot::drawCMS(xjjroot::CMS::internal, infos["data"]["input_tex"] + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(-1), btex.label_pt(), infos["data"]["cut_tex"] }, 0.038, 13);
      pdf->write();
    }

    for (int i=0; i<ny; i++) {
      // sub
      auto* leg_sub = init_leg(5);
      leg_sub->AddEntry(h1ys_data_main.at(name)[i], "Data in signal region", "pf");
      leg_sub->AddEntry((TObject*)0, mranges[i]->str_signal().c_str(), "");
      leg_sub->AddEntry(h1ys_data_sideband_scaled.at(name)[i], "Sideband scaled", "pf");
      leg_sub->AddEntry((TObject*)0, mranges[i]->str_sideband().c_str(), "");
      leg_sub->AddEntry(h1ys_data_sub.at(name)[i], "Sideband subtracted", "pe");

      set_hsminmax({ h1ys_data_main.at(name)[i], h1ys_data_sideband_scaled.at(name)[i], h1ys_data_sub.at(name)[i] });
      pdf->prepare();
      pdf->getc()->SetLogy(0);
      if (the_var.logy)
        pdf->getc()->SetLogy();
      h1ys_data_main.at(name)[i]->Draw("hist e1");
      h1ys_data_sideband_scaled.at(name)[i]->Draw("hist e1 same");
      h1ys_data_sub.at(name)[i]->Draw("pe1 same");
      leg_sub->Draw();
      xjjroot::drawCMS(xjjroot::CMS::internal, infos["data"]["input_tex"] + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt(), infos["data"]["cut_tex"] }, 0.038, 13);
      pdf->write();
    }
    
  }

  pdf->draw_cover({
      "#bf{data input} " + infos["data"]["input"],
      "#bf{MC input} " + infos["template"]["input"],
    });
  
  pdf->close();

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (auto& [name, _] : h1ys_data_main) {
    outf->mkdir(Form("dir_%s", name.c_str()))->cd();
    for (auto* h1ys : { &h1ys_data_sigswap, &h1ys_data_sub,
                       &h1ys_data_main, &h1ys_data_sideband, &h1ys_data_sideband_scaled,
                       &h1ys_mc_match } ) {
      for (auto h : h1ys->at(name)) {
        if (h) xjjroot::writehist(h);
      }
    }
  }
  outf->cd();
  xjjroot::writehist(h3_bins);
  for (auto& [iname, info] : infos) {
    outf->mkdir(iname.c_str())->cd();
    auto* t_data = new TTree("info", "");
    for (auto& [key, content] : info) {
      t_data->Branch(key.c_str(), &content);
    }
    t_data->Fill();
    t_data->Write();
    outf->cd();
  }
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}

void MassRange::initialize(TTree* t) {
#define GET_VAL(q)                              \
  double q##_ = 0;                              \
  t->SetBranchAddress( #q, &q##_);

  GET_VAL(mass_signal_low);
  GET_VAL(mass_signal_high);
  GET_VAL(deltamass_sideband_low);
  GET_VAL(deltamass_sideband_high);
  
  t->GetEntry(0);

#define SET_VAL(q)                                      \
  q = q##_;                                             \
  __XJJLOG << ">> " #q " -> " << q << std::endl;
  
  SET_VAL(mass_signal_low);
  SET_VAL(mass_signal_high);
  SET_VAL(deltamass_sideband_low);
  SET_VAL(deltamass_sideband_high);

  tex_signal = Form("%.3f < m < %.3f GeV", mass_signal_low, mass_signal_high);
  __XJJLOG << ">> signal range -> " << tex_signal << std::endl;
  tex_sideband = Form("%.3f <#scale[0.5]{ }#Deltam < %.3f GeV", deltamass_sideband_low, deltamass_sideband_high);
  __XJJLOG << ">> sideband range -> " << tex_sideband << std::endl;
  __XJJLOG << ">> sideband scale -> " << sideband_scale() << std::endl;
}

hName parse_hname(std::string hname, bool verbose) {
  hName result;
  result.index_y = std::atoi(xjjc::str_erasestar(hname, "*__y-").c_str());
  auto t_noy_noh1 = xjjc::str_erasestar(hname, "__y-*");
  t_noy_noh1 = xjjc::str_eraseall(t_noy_noh1, "h1_");
  auto strs = xjjc::str_divide_trim(t_noy_noh1, "_");
  if (strs.size() > 0) {
    result.type = strs.back();
    result.varname = xjjc::str_eraseall(t_noy_noh1, "_"+result.type);
  }
  // h1_DsvpvDistance_2D_data-sideband__y-3
  if (verbose)
    __XJJLOG << ">> " << hname << " : [" << result.varname << ", " << result.type << ", " << result.index_y << " ]" << std::endl;

  return result;
}
