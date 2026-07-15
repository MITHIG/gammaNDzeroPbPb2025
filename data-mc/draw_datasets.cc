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

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  for (int i=0; i<h3_bins->GetXaxis()->GetNbins(); i++) {
    auto* dir = (TDirectory*)inf->GetDirectory(Form("dir__y-%d", i));
    if (!dir) {
      __XJJLOG << "!! bad directory, skip." << std::endl;
      continue;
    }
    auto make_h1s = [&dir, &i](std::string type, xjjroot::thgrstyle style = {}) {
      const auto suffix = std::string(Form("%s__y-%d", type.c_str(), i));
      auto vh1s = xjjana::getobj_regexp<TH1D>(dir, "h1_.+" + suffix);
      std::map<std::string, TH1D*> h1s;
      for (auto& h : vh1s) {
        xjjroot::sethempty(h, 0, 0.4);
        xjjroot::setthgrstyle(h, style);
        auto name = xjjc::str_eraseall(h->GetName(), { "h1_", suffix });
        h1s[name] = h;
      }
      if (h1s.empty()) {
        __XJJLOG << "!! no matching histograms : " << suffix << std::endl;
      }
      return h1s;
    };
    auto make_h1s_norm = [](std::map<std::string, TH1D*>& h1s, xjjroot::thgrstyle style = {}) {
      std::map<std::string, TH1D*> h1s_norm;
      for (auto& [name, h] : h1s) {
        auto* hnorm = (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "h1_", "h1_norm_").c_str());
        hnorm->Scale(1./h->Integral(), "width");
        hnorm->GetYaxis()->SetTitle("Distribution (self-normalized)");
        xjjroot::sethempty(hnorm, 0, 0.4);
        xjjroot::setthgrstyle(hnorm, style);
        h1s_norm[name] = hnorm;
      }
      return h1s_norm;
    };

    auto* mrange = new MassRange((TTree*)dir->Get("mass_range"));

    auto make_h1s_scale = [](const std::map<std::string, TH1D*>& h1s, double scale) {
      std::map<std::string, TH1D*> h1s_scaled;
      for (auto& [name, h] : h1s) {
        auto* h_scaled = (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "h1_", "h1_sub_").c_str());
        // __XJJLOG << ">> scale = " << scale << std::endl;
        h_scaled->Scale(scale);
        xjjroot::sethempty(h_scaled, 0, 0.4);
        h1s_scaled[name] = h_scaled;
      }
      return h1s_scaled;
    };
    auto make_h1s_sub = [](const std::map<std::string, TH1D*>& h1s_main,
                           const std::map<std::string, TH1D*>& h1s_sideband_scaled,
                           xjjroot::thgrstyle style = {}) {
      std::map<std::string, TH1D*> h1s_sub;
      for (auto& [name, h] : h1s_main) {
        auto* hsub = (TH1D*)h->Clone(xjjc::str_replaceall(h->GetName(), "h1_", "h1_sub_").c_str());
        hsub->Add(h1s_sideband_scaled.at(name), -1);
        xjjroot::sethempty(hsub, 0, 0.4);
        xjjroot::setthgrstyle(hsub, style);
        h1s_sub[name] = hsub;
      }
      return h1s_sub;
    };
    
    auto h1s_data_sig = make_h1s("_data_sig", { xjjroot::mycolor_middle["green"], 20, 1.5, xjjroot::mycolor_middle["green"], 1, 1 });
    auto h1s_norm_data_sig = make_h1s_norm(h1s_data_sig);
    if (h1s_norm_data_sig.empty()) continue;
    auto h1s_data_sigswap = make_h1s("_data_sigswap", { xjjroot::mycolor_middle["blue"], 25, 1.5, xjjroot::mycolor_middle["blue"], 1, 1 });
    auto h1s_norm_data_sigswap = make_h1s_norm(h1s_data_sigswap);
    if (h1s_norm_data_sigswap.empty()) continue;
    auto h1s_mc_match = make_h1s("_mc_match", { xjjroot::mycolor_middle["red"], 21, 0, xjjroot::mycolor_middle["red"], 1, 1, xjjroot::mycolor_middle["red"], 1, 3004 });
    auto h1s_norm_mc_match = make_h1s_norm(h1s_mc_match);
    if (h1s_norm_mc_match.empty()) continue;

    auto h1s_data_main = make_h1s("_data_main", { kGray+2, 24, 1.5, kGray+2, 1, 1, kGray+2, 0.1, 1001 });
    if (h1s_data_main.empty()) continue;
    auto h1s_data_sideband = make_h1s("_data_sideband", { xjjroot::mycolor_middle["cyan"], 24, 1.5, xjjroot::mycolor_middle["cyan"], 1, 1, xjjroot::mycolor_middle["cyan"], 1, 3004 });
    if (h1s_data_sideband.empty()) continue;
    auto h1s_data_sideband_scaled = make_h1s_scale(h1s_data_sideband, mrange->sideband_scale());
    auto h1s_data_sub = make_h1s_sub(h1s_data_main, h1s_data_sideband_scaled,
                                     { xjjroot::mycolor_middle["blue"], 20, 1.7, xjjroot::mycolor_middle["blue"], 1, 1 });
    auto h1s_norm_data_sub = make_h1s_norm(h1s_data_sub);
    auto h1s_norm_data_main = make_h1s_norm(h1s_data_main, { kGray, 20, 1.5, kGray, 1, 1, kGray, 0.2, 1001 });
    
    pdf->draw_cover({ "#bf{" + btex.label_y(i) + "}" }, 0.05);

    for (auto& [name, _] : h1s_norm_data_sig) {
      const auto the_var = var_by_name(xjjc::str_eraseall(name, "_norm"));
      if (the_var.varname.empty()) continue;

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
      
      auto* leg_sub = init_leg(5);
      leg_sub->AddEntry(h1s_data_main.at(name), "Data in signal region", "pf");
      leg_sub->AddEntry((TObject*)0, mrange->str_signal().c_str(), "");
      leg_sub->AddEntry(h1s_data_sideband_scaled.at(name), "Sideband scaled", "pf");
      leg_sub->AddEntry((TObject*)0, mrange->str_sideband().c_str(), "");
      leg_sub->AddEntry(h1s_data_sub.at(name), "Sideband subtracted", "pe");
      
      set_hsminmax({ h1s_data_main.at(name), h1s_data_sideband_scaled.at(name), h1s_data_sub.at(name) });
      pdf->prepare();
      pdf->getc()->SetLogy(0);
      if (the_var.logy)
        pdf->getc()->SetLogy();
      h1s_data_main.at(name)->Draw("hist e1");
      h1s_data_sideband_scaled.at(name)->Draw("hist e1 same");
      h1s_data_sub.at(name)->Draw("pe1 same");
      leg_sub->Draw();
      xjjroot::drawCMS(xjjroot::CMS::internal, infos["data"]["input_tex"] + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt(), infos["data"]["cut_tex"] }, 0.038, 13);
      pdf->write();

      auto* leg_norm = init_leg(6);
      leg_norm->AddEntry(h1s_norm_data_main.at(name), "Data in signal region", "f");
      leg_norm->AddEntry((TObject*)0, mrange->str_signal().c_str(), "");
      leg_norm->AddEntry(h1s_norm_mc_match.at(name), "Gen-matched MC", "f");
      leg_norm->AddEntry(h1s_norm_data_sigswap.at(name), "sPlot signal+swap", "pe");
      // leg_norm->AddEntry(h1s_norm_data_sig.at(name), "sPlot signal", "pe");
      leg_norm->AddEntry(h1s_norm_data_sub.at(name), "Sideband subtracted", "pe");
      leg_norm->AddEntry((TObject*)0, mrange->str_sideband().c_str(), "");

      auto hlist_norm = std::vector<TH1D*>{ h1s_norm_data_sig.at(name), h1s_norm_data_sigswap.at(name), h1s_norm_mc_match.at(name), h1s_norm_data_main.at(name) };
      if (name != "Dmass") hlist_norm.push_back(h1s_norm_data_sub.at(name));
      set_hsminmax(hlist_norm);
      pdf->prepare();
      pdf->getc()->SetLogy(0);
      if (the_var.logy)
        pdf->getc()->SetLogy();
      h1s_norm_data_main.at(name)->Draw("hist");
      h1s_norm_mc_match.at(name)->Draw("hist e1 same");
      h1s_norm_data_sigswap.at(name)->Draw("pe1 same");
      // h1s_norm_data_sig.at(name)->Draw("pe1 same");
      if (name != "Dmass") h1s_norm_data_sub.at(name)->Draw("pe1 same");
      leg_norm->Draw();
      xjjroot::drawCMS(xjjroot::CMS::internal, infos["data"]["input_tex"] + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.25, 0.86, { btex.label_y(i), btex.label_pt(), infos["data"]["cut_tex"] }, 0.038, 13);
      pdf->write();
    }

  } // loop y 
  
  pdf->close();

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
