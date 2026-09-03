#include <TFile.h>
#include <TH3D.h>
#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "util.h"
#include "draw.h"

#include "variables.h"

int macro(const std::string& inputname, const std::string& outputname) {
  const auto inputp = util::parse_input(inputname);
  auto* inf = TFile::Open(inputp.content.c_str());
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins");
  const auto ny = h3_bins->GetXaxis()->GetNbins(), npt = h3_bins->GetZaxis()->GetNbins();
  const draw::bintex tbins(h3_bins, 0, 2);
  auto info = xjjana::getval_regexp((TTree*)inf->Get("info"));
  __XJJLOG << "++ info" << std::endl;
  xjjc::print_tab(info, -1);
  const bool has_dreq = !xjjc::str_contains(info.at("dcut"), "Dnoreq");
  auto this_is = [&info](const std::string& item) {
    return static_cast<bool>(std::atoi(info.at(item).c_str()));
  };
  const auto colors = xjjroot::grayscales_color(ny, xjjroot::mycolor_dark["red"], 1, 0.3);

  std::map<std::string, TH3D*> h3s;
  std::map<std::string, xjjc::array2D<TH1D*>> h1ptys, h1ptys_refy;
  for (auto* h3 : xjjana::getobj_regexp<TH3D>(inf, "h3_y_.+_pt")) {
    h3->Sumw2();
    const auto name = xjjc::str_eraseall(h3->GetName(), std::vector<std::string>{ "h3_y_", "_pt" });
    auto it_var = std::find_if(variables.begin(), variables.end(), [&name](const global::variable& v) {
      return v.varname == name;
    });
    if (it_var == variables.end()) {
      __XJJLOG << "?? unknown varname: " << name << ", skip." << std::endl;
      continue;
    }
    const auto logy = (*it_var).logy;
    
    h3s[name] = h3;
    h1ptys[name] = xjjc::array2d<TH1D*>(npt, ny);
    h1ptys_refy[name] = xjjc::array2d<TH1D*>(npt, ny);
    for (int i=0; i<npt; i++) {
      for (int j=0; j<ny; j++) {
        const int refj = this_is("dir_Ngamma") ? ny-j-1 : j;
        auto* h1 = h3->ProjectionY(Form("h1_%s_pt-%d_y-%d", name.c_str(), i, j),
                                   j+1, j+1, // y bin
                                   i+1, i+1, // pt bin
                                   "e");
        h1->GetYaxis()->SetTitle("Entries");
        h1->GetXaxis()->SetNdivisions(505);
        xjjroot::sethempty(h1, 0, 0);
        xjjroot::setthgrstyle(h1, colors[refj], 20, 1.3, colors[refj], 1, 1);
        h1ptys.at(name)[i][refj] = h1;
        h1ptys_refy.at(name)[i][refj] = (TH1D*)h1->Clone(xjjc::str_replaceall(h1->GetName(), Form("_y-%d", j), Form("_refy-%d", refj)).c_str());
      }
      if (logy) {
        xjjana::sethsnonzeromin(h1ptys.at(name)[i], 0.9);
        xjjana::sethsmax(h1ptys.at(name)[i], 20);
      } else {
        xjjana::sethsabsmin(h1ptys.at(name)[i], 1.e-3);
        xjjana::sethsmax(h1ptys.at(name)[i], 1.5);
      }
    }
  }

  auto pos_y = [](int n = 0) { return 0.87-0.042*n; };
  float legx1 = 0.55, legy2 = pos_y(3), legy1 = pos_y(3+ny);
  auto* leg = new TLegend(legx1, legy1, 0.75, legy2);
  xjjroot::setleg(leg, 0.038);
  for (int j=0; j<ny; j++) {
    auto* h = h1ptys.begin()->second.front()[j];
    leg->AddEntry(h, tbins.label_y(xjjc::str_extract_index(h->GetName(), "_y-")).c_str(), "p");
  }
  std::vector<std::string> ts_common = {
    // std::string(this_is("l1_ZDCOr") ? "L1_ZDCOr + " : (this_is("l1_ZeroBias") ? "L1_ZeroBias + " : "")) + "PV filter + cscTightHalo",
    std::string(this_is("l1_zdcor") ? "L1_ZDCOr + " : (this_is("l1_zerobias") ? "L1_ZeroBias + " : "")) + "PV filter + cscTightHalo",
    (this_is("dir_gammaN") ? "#gamma #rightarrow Z+ (#gammaN)" : "#gamma #rightarrow Z- (N#gamma)"),
    // info["evt_tex"],
  };

  auto treat_tex = [&this_is, &info](const std::string& tex_tochange) {
    return xjjc::str_replaceall(tex_tochange, {
        { "1nXOR", this_is("dir_gammaN") ? "Xn0n" : "0nXn" },
      });
  };
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  for (auto& [name, vh1s] : h1ptys) {
    auto var = var_by_name(name);
    auto varnote = xjjc::str_divide_trim(treat_tex(var.note), ",");
    auto texs_left = xjjc::vec_append_into(ts_common, varnote);
    for (int i=0; i<npt; i++) {
      pdf->prepare();
      if (var.logy) pdf->getc()->SetLogy();
      else pdf->getc()->SetLogy(0);
      vh1s[i].front()->Draw("axis");
      for (int j=0; j<ny; j++) {
        vh1s[i][j]->Draw("pe same");
      }
      xjjroot::drawCMS(xjjroot::CMS::internal, info["input_tex"] + " (5.36 TeV)");
      if (has_dreq) {
        leg->Draw();
        xjjroot::drawtexgroup(legx1+0.005, legy2+0.042*2, { "At least a D cand in", tbins.label_pt(i) }, 0.038, 13, 42, 1.1);
        xjjroot::drawtexgroup(legx1+0.005, legy1, xjjc::str_divide_trim(info["dcut_tex"], ","), 0.038, 13, 42, 1.1);
      }
      xjjroot::drawtexgroup(0.24, pos_y(), texs_left, 0.038, 13, 42, 1.1);
      gPad->RedrawAxis();
      pdf->write();
    }
  }
  pdf->close();
  
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (const auto& [_, h] : h3s) xjjroot::writehist(h);
  for (const auto& [_, hhh] : h1ptys_refy)
    for (const auto& hh : hhh)
      for (const auto& h : hh)
        xjjroot::writehist(h);
  for (const auto& h : xjjana::getobj_regexp<TH1D>(inf, "h1_.+"))
    xjjroot::writehist(h);
  xjjroot::writehist(h3_bins);
  auto* t = new TTree("info", "");
  std::map<std::string, std::string> t_container; // for lifetime of the intermediate string
  auto cast_branch = [&t, &t_container]<typename T>(const std::string& name, const T& x) {
    t_container[name] = xjjc::to_string(x);
    t->Branch(name.c_str(), &(t_container[name]));
  };
  for (const auto& [key, val] : info) {
    cast_branch(key, val);
  }
  t->Fill();
  t->Write();
  outf->cd();
  xjjroot::closefile(outf);
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3)
    return macro(argv[1], argv[2]);
}
