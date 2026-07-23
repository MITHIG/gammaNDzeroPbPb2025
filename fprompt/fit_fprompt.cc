#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "TGraphAsymmErrors.h"

#define __COOK_NAME__
#include "style.h"
#include "../include/draw.h"

#include "fpfitter.h"

int macro(const std::string& inputname, const std::string& outputname) {
  //
  auto* inf = TFile::Open(inputname.c_str());
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
  draw::bintex tbins(h3_bins, 0, 2);
  auto ny = h3_bins->GetXaxis()->GetNbins();
  auto* h1_bins_sf = xjjana::getobj<TH1D>(inf, "h1_bins_sf");
  auto nsf = h1_bins_sf->GetXaxis()->GetNbins();

  auto* dir_var = xjjana::getobj_regexp_first<TDirectory>(inf, "dir_.+");
  if (!dir_var) {
    __XJJLOG << "!! no directory dir_.+, abort." << std::endl;
    return 2;
  }
  const auto var = xjjc::str_eraseall(dir_var->GetName(), "dir_");
  
  std::map<std::string, xjjc::info> infos;
  for (const std::string name : { "data", "prompt", "nonprompt" } ) {
    infos[name] = xjjana::get_info(inf, Form("info/%s", name.c_str()));
  }
  infos["fit"]["var"] = var;

  for (auto& [name, info] : infos) {
    __XJJLOG << ">> " << name << std::endl;
    xjjc::print_tab(info, -1);
  }

  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  
  //
  for (int i=0; i < ny; i++) {
    auto* dir_y_output = outf->mkdir(Form("dir__y-%d", i));
    auto* dir_y = dir_var->GetDirectory(Form("dir__y-%d", i));

    // prepare MC histograms
    std::map<std::string, std::vector<TH1D*>> h1sfs_mc;
    for (const std::string type_mc : { "mc-prompt", "mc-nonprompt" }) {
      h1sfs_mc[type_mc] = xjjana::getobj_regexp<TH1D>(dir_y, ".+" + type_mc + ".*_sf-.+__y-" + xjjc::to_string(i), "", false);
      std::sort(h1sfs_mc[type_mc].begin(), h1sfs_mc[type_mc].end(),
                [](const TH1D* ha, const TH1D* hb) {
                  return index_sf(ha) < index_sf(hb);
                });
      xjjroot::print_vec_v(h1sfs_mc[type_mc], 0);
    }
    if (nsf != h1sfs_mc.at("mc-prompt").size()) {
      __XJJLOG << "!! inconsistent number of scale factors, abort." << std::endl;
      return 2;
    }

    for (const std::string type_data : { "data-sub", "data-sigswap" }) {
      auto* dir_type_output = dir_y_output->mkdir(Form("dir_%s", type_data.c_str()));

      auto* h1_data = xjjana::getobj_regexp_first<TH1D>(dir_y, ".+" + type_data + ".*__y-" + xjjc::to_string(i), "", true);
      if (!h1_data)
        return 2;

      auto make_h1_sf = [&h1_bins_sf, &var, &type_data, &i](std::string yvar, std::string ytitle) {
        auto* h1_sf = (TH1D*)h1_bins_sf->Clone(Form("h1_%s-sf_%s_%s__y-%d", yvar.c_str(), var.c_str(), type_data.c_str(), i)); //
        h1_sf->GetYaxis()->SetTitle(ytitle.c_str());
        xjjroot::sethempty(h1_sf, 0, 0.3);
        xjjroot::setthgrstyle(h1_sf, kBlack, 20, 1.5, kBlack, 1, 1);
        return h1_sf;
      };
      auto* h1_chi2 = make_h1_sf("chi2", "#chi^{2} / ndf");
      auto* h1_fprompt = make_h1_sf("fprompt", "#it{f}_{prompt}");
      xjjroot::setthgrstyle(h1_fprompt, kGray, -1, -1, kGray);

      auto* gr_fprompt = new TGraphAsymmErrors();
      gr_fprompt->SetName(xjjc::str_replaceall(h1_fprompt->GetName(), "h1_", "gr_").c_str());
      xjjroot::setthgrstyle(gr_fprompt, kBlack, 21, 1.3, kBlack, 1, 1);
     
      // fit 
      pdf->draw_cover({ "#bf{Variable} " + var, "#bf{Signal extraction} " + style_data(type_data).title, tbins.label_y(i) });

      int ngr_fprompt = 0;
      std::vector<fpfitter*> results(nsf, nullptr);
      for (int k=0; k<nsf; k++) {
        auto *h1_mc_prompt = h1sfs_mc.at("mc-prompt")[k], *h1_mc_nonprompt = h1sfs_mc.at("mc-nonprompt")[k];
        xjjc::progressslide(k, nsf, 2);

        auto* fitter = new fpfitter(h1_data, h1_mc_prompt, h1_mc_nonprompt, "-" + xjjc::str_eraseall(type_data, "data-"), "-sf-" + xjjc::to_string(k));
        if (fitter->status() > 0)
          continue;
        fitter->fit();
        // fitter->print_fitresult();
        h1_fprompt->SetBinContent(k+1, fitter->fprompt());
        h1_fprompt->SetBinError(k+1, fitter->fprompt_err_par());
        h1_chi2->SetBinContent(k+1, fitter->chi2() / fitter->ndf());
        h1_chi2->SetBinError(k+1, 0);

        const double x = h1_bins_sf->GetBinCenter(k+1);
        const double ex = 0.5 * h1_bins_sf->GetBinWidth(k+1);
        gr_fprompt->SetPoint(ngr_fprompt, x, fitter->fprompt());
        gr_fprompt->SetPointError(ngr_fprompt, ex, ex,
                                  fitter->fprompt_err_low(),
                                  fitter->fprompt_err_high());
        ngr_fprompt++;
        
        // draw 
        pdf->prepare();
        auto pads = fitter->draw(pdf->getc());
        pads.front()->cd();
        xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)", 1./fpfitter::pratio);
        xjjroot::drawtexgroup(fpfitter::xleft - 0.01, fpfitter::ytop - 0.005, {
            style_data(type_data).title.c_str(),
            tbins.label_pt(-1), tbins.label_y(i),
            "#alpha_{reso} = " + std::string(Form("%.2f", h1_bins_sf->GetBinCenter(k+1))),
          }, fpfitter::tsize/fpfitter::pratio, 33);
        pdf->getc()->cd();
        pdf->write();

        results[k] = fitter;
      } // for (int k=0; k<nsf; k++) {
      
      auto ibin_best = h1_chi2->GetMinimumBin(), ibin_fix = h1_chi2->FindBin(1.);
      if (ibin_fix <= 0 || ibin_fix >= h1_chi2->GetNbinsX()) {
        __XJJLOG << "!! can not find ibin_fix, abort." << std::endl;
        return 3;
      }

      xjjana::sethminmax(h1_chi2, 0, 1.3);
      pdf->prepare();
      h1_chi2->Draw("p");
      xjjroot::drawline(h1_chi2->GetBinCenter(ibin_best), h1_chi2->GetMinimum(), h1_chi2->GetBinCenter(ibin_best), h1_chi2->GetBinContent(ibin_best), kGray+3, 2, 2);
      xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.55, fpfitter::ytop - 0.005, {
          style_data(type_data).title.c_str(),
          tbins.label_pt(-1), tbins.label_y(i),
        }, fpfitter::tsize, 13);
      pdf->write();

      xjjana::sethabsminmax(h1_fprompt, 0, 1.5);
      pdf->prepare();
      h1_fprompt->Draw("pe1");
      gr_fprompt->Draw("pe1 same");
      xjjroot::drawline(h1_fprompt->GetBinCenter(ibin_best), h1_fprompt->GetMinimum(), h1_fprompt->GetBinCenter(ibin_best), h1_fprompt->GetBinContent(ibin_best), kGray+3, 2, 2);
      xjjroot::drawCMS(xjjroot::CMS::internal, infos.at("data").at("input_tex") + " (5.36 TeV)");
      xjjroot::drawtexgroup(0.55, fpfitter::ytop - 0.005, {
          style_data(type_data).title.c_str(),
          tbins.label_pt(-1), tbins.label_y(i),
        }, fpfitter::tsize, 13);
      pdf->write();

      dir_type_output->cd();
      xjjroot::writehist(h1_chi2);
      xjjroot::writehist(h1_fprompt);
      xjjroot::writehist(gr_fprompt);
      dir_type_output->mkdir("fit_best")->cd();
      results[ibin_best-1]->write_to_file();
      dir_type_output->mkdir("fit_fix")->cd();
      results[ibin_fix-1]->write_to_file();
      dir_type_output->cd();
    } // for (auto& [type_data, h1_data] : h1s_data) {
    
  } // for (int i=0; i < ny; i++) {

  pdf->close();
  outf->cd();
  xjjroot::writehist(h3_bins);
  xjjroot::writehist(h1_bins_sf);
  outf->mkdir("info")->cd();
  for (auto& [iname, info] : infos) {
    auto* t_data = new TTree(iname.c_str(), "");
    for (auto& [key, content] : info)
      t_data->Branch(key.c_str(), &content);
    t_data->Fill();
    t_data->Write();
  }
  outf->cd();
  xjjroot::closefile(outf);
  inf->Close();
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
