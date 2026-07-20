#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"

#include "xjjanauti.h"
#include "xjjmypdf.h"

#include "../include/util.h"
#define __VARIABLES_ROOSPLOT__
#include "variables.h"

double signal_effective_sigma(RooWorkspace *ws, double area_frac);
std::vector<std::pair<RooRealVar*, bool>> fix_shape_parameters(RooAbsPdf* pdf, RooAbsData* data, const std::vector<std::string> &floating_yields);
void restore_parameter_states(const std::vector<std::pair<RooRealVar*, bool>> &old_states);

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

  auto wsys = xjjana::getobj_regexp<RooWorkspace>(inf, "ws__y-.+");
  auto* h3_bins = xjjana::getobj<TH3D>(inf, "h3_bins_y-mass-pt");
  if (wsys.size() != h3_bins->GetXaxis()->GetNbins()) {
    __XJJLOG << "!! inconsistent bin number: " << "ws__y-* vs h3_bins, abort." << std::endl;
    return 2;
  }
  // draw::bintex btex(h3_bins, 0, 2);
  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  xjjroot::writehist(h3_bins);
  
  // xjjroot::setgstyle(1);
  // auto* pdf = new xjjroot::mypdf("figspdf/" + outputname + ".pdf");

  // std::map<std::string, std::vector<TH1D*>> h1ys;
  std::vector<std::map<std::string, TH1D*>> h1ys; 
  for (int i=0; i<wsys.size(); i++) {
    // main pdf and dataset
    auto* ws = wsys.at(i);
    auto* pdf_total = ws->pdf("pdf_total"); // RooAbsPdf*
    auto* ds_data = (RooDataSet*)ws->data(Form("data_main__y-%d", i));
    auto* ds_mc_match = (RooDataSet*)ws->data(Form("mc_match__y-%d", i));
    if (!pdf_total || !ds_data) {
      __XJJLOG << "!! missing pdf_total or data_main__y-" << i << ", skip." << std::endl;
      continue;
    }
    if (!ds_mc_match) {
      __XJJLOG << "?? missing mc_match__y-" << i << "; skip matched-MC histograms." << std::endl;
    }

    // yield variables
    auto* Dmass = ws->var("Dmass");
    RooArgSet mass_obs(*Dmass); // used for get pdf_sig, pdf_swap value
    auto* n_sigswap = ws->var("n_sigswap"); // RooRealVar*
    auto* n_bkg = ws->var("n_bkg");

    // to calculate splitting of sig and swap
    auto* pdf_sig = ws->pdf("pdf_sig");
    auto* pdf_swap = ws->pdf("pdf_swap");
    const auto val_frac_sig = ws->var("par_sigswap_frac")->getVal();

    // calculate signal region
    const auto val_mean = ws->var("par_mean")->getVal();
    const auto eff_sigma = signal_effective_sigma(ws, xjjana::frac_1sigma);
    const auto width_signal = eff_sigma*2;
    double mass_signal_low = val_mean - width_signal,
      mass_signal_high = val_mean + width_signal;
    double deltamass_sideband_low = width_signal;
    __XJJLOG << "-- try to set sideband to 2sigma - 7*sigma" << std::endl;
    const auto width_try = 7*eff_sigma;
    double deltamass_sideband_high = std::min(width_try, val_mean - Dmass->getMin());
    deltamass_sideband_high = std::min(deltamass_sideband_high, Dmass->getMax() - val_mean);
    if (val_mean-width_try < Dmass->getMin() || val_mean+width_try > Dmass->getMax()) {
      __XJJLOG << "?? side band range too wide, adjusted to Dmass range: " << std::endl;
    }
    std::cout << "    [ " << val_mean-deltamass_sideband_high << " - " << val_mean-deltamass_sideband_low << " ] , "
              << "[ " << val_mean+deltamass_sideband_low << " - " << val_mean+deltamass_sideband_high << " ]"
              << std::endl;

    // make splot
    auto* ds_data_sigswap = dynamic_cast<RooDataSet*>(ds_data->Clone(Form("%s_splot_sigswap", ds_data->GetName())));
    auto old_states = fix_shape_parameters(pdf_total, ds_data_sigswap, { "n_sigswap", "n_bkg" });
    __XJJLOG << "++ make splot" << std::endl;
    auto* splot_sigswap = new RooStats::SPlot("splot_sigswap", "sPlot for signal+swap",
                                              *ds_data_sigswap, pdf_total,
                                              RooArgList(*n_sigswap, *n_bkg));
    restore_parameter_states(old_states);
    
    __XJJLOG << "++ loop variables" << std::endl;
    std::map<std::string, TH1D*> h1s_data_main, h1s_data_sigswap, h1s_data_sig, h1s_data_sideband, h1s_mc_match;
    const RooArgSet* columns = ds_data->get();
    for (RooAbsArg* arg : *columns) {
      auto *var = dynamic_cast<RooRealVar*>(arg);
      if (!var) continue;
      const std::string name = var->GetName();
      if (name.empty() || name[0] != 'D') continue;
      const auto the_var = var_by_name(name);
      if (the_var.varname.empty()) continue;
      
      const auto* hname = Form("h1_%s_data-main__y-%d", name.c_str(), i);
      if (!the_var.bins.empty()) {
        h1s_data_main[name] = new TH1D(hname, Form(";%s;Entries", the_var.vartex.c_str()), the_var.bins.size()-1, the_var.bins.data());
      } else {
        h1s_data_main[name] = new TH1D(hname, Form(";%s;Entries", the_var.vartex.c_str()), the_var.nbin, the_var.varmin, the_var.varmax);
      }
      h1s_data_main[name]->Sumw2();
      h1s_data_sigswap[name] = (TH1D*)h1s_data_main[name]->Clone(xjjc::str_replaceall(h1s_data_main[name]->GetName(), "data-main", "data-sigswap").c_str());
      h1s_data_sig[name] = (TH1D*)h1s_data_main[name]->Clone(xjjc::str_replaceall(h1s_data_main[name]->GetName(), "data-main", "data-sig").c_str());
      h1s_data_sideband[name] = (TH1D*)h1s_data_main[name]->Clone(xjjc::str_replaceall(h1s_data_main[name]->GetName(), "data-main", "data-sideband").c_str());
      h1s_mc_match[name] = (TH1D*)h1s_data_main[name]->Clone(xjjc::str_replaceall(h1s_data_main[name]->GetName(), "data-main", "mc-match").c_str());
      __XJJLOG << "     >> " << name << std::endl;
    }
    
    double sum_weight_sigswap = 0, sum_weight_sig = 0;
    const auto nentries = ds_data_sigswap->numEntries();
    for (int i = 0; i < nentries; ++i) {
      xjjc::progressslide(i, nentries);

      const RooArgSet* row_data = ds_data_sigswap->get(i);
      const auto mass = row_data->getRealValue("Dmass");

      // get weights(mass)
      const auto weight_sigswap = row_data->getRealValue("n_sigswap_sw");
      Dmass->setVal(mass);
      const double density_sig = pdf_sig->getVal(&mass_obs) * val_frac_sig,
        density_swap = pdf_swap->getVal(&mass_obs) * (1 - val_frac_sig),
        density_sigswap = density_sig + density_swap,
        absfrac_sig = (density_sigswap > 0.) ? density_sig/density_sigswap : 0.;
      const double weight_sig = weight_sigswap * absfrac_sig; // I think this way is wrong, only weight_sigswap is reliable
      sum_weight_sigswap += weight_sigswap;
      sum_weight_sig += weight_sig;

      for (auto& [name, _] : h1s_data_main) {
        const auto value = row_data->getRealValue(name.c_str());
        h1s_data_sigswap.at(name)->Fill(value, weight_sigswap);
        h1s_data_sig.at(name)->Fill(value, weight_sig);

        // fill signal region
        if (mass > mass_signal_low && mass < mass_signal_high) {
          h1s_data_main.at(name)->Fill(value);
        }
        // fill side band
        if (std::abs(mass - val_mean) > deltamass_sideband_low &&
            std::abs(mass - val_mean) < deltamass_sideband_high) {
          h1s_data_sideband.at(name)->Fill(value);
        }
      }
    }
    xjjc::progressbar_summary(nentries);
    __XJJLOG << ">> sPlot yield check" << std::endl
             << "   [ n_sigswap = " << n_sigswap->getVal()
             << ", sum(n_sigswap_sw) = " << sum_weight_sigswap
             << ", derived sum(signal sw) = " << sum_weight_sig
             << " ]" << std::endl;

    const auto nentries_mc = ds_mc_match ? ds_mc_match->numEntries() : 0;
    for (int i = 0; i < nentries_mc; ++i) {
      xjjc::progressslide(i, nentries_mc);

      const RooArgSet* row_mc = ds_mc_match->get(i);

      for (auto& [name, h1] : h1s_mc_match) {
        const auto value = row_mc->getRealValue(name.c_str());
        h1s_mc_match.at(name)->Fill(value);
      }      
    }
    if (ds_mc_match) xjjc::progressbar_summary(nentries_mc);

    outf->cd();
    outf->mkdir(Form("dir__y-%d", i))->cd();
    auto* t_mass = new TTree("mass_range", "");
    t_mass->Branch("mass_signal_low", &mass_signal_low, "mass_signal_low/D");
    t_mass->Branch("mass_signal_high", &mass_signal_high, "mass_signal_high/D");
    t_mass->Branch("deltamass_sideband_low", &deltamass_sideband_low, "deltamass_sideband_low/D");
    t_mass->Branch("deltamass_sideband_high", &deltamass_sideband_high, "deltamass_sideband_high/D");
    t_mass->Fill();
    t_mass->Write();
    for (auto& [_, h] : h1s_data_main) xjjroot::writehist(h);
    for (auto& [_, h] : h1s_data_sigswap) xjjroot::writehist(h);
    for (auto& [_, h] : h1s_data_sideband) xjjroot::writehist(h);
    for (auto& [_, h] : h1s_data_sig) xjjroot::writehist(h);
    for (auto& [_, h] : h1s_mc_match) xjjroot::writehist(h);

    outf->cd();
  } // loop y 
  
  // pdf->close();
  outf->cd();
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

double signal_effective_sigma(RooWorkspace *ws, double area_frac) {
  auto *pdf_sig = ws->pdf("pdf_sig");
  auto *Dmass = ws->var("Dmass");
  auto *mean = ws->var("par_mean");
  if (!pdf_sig || !Dmass || !mean) return -1.;
  if (area_frac <= 0. || area_frac >= 1.) return -1.;

  RooArgSet mass_obs(*Dmass);
  const double original_mass = Dmass->getVal();
  const double mass_min = Dmass->getMin();
  const double mass_max = Dmass->getMax();
  const double mass_mean = mean->getVal();

  auto integrate = [&](double xmin, double xmax) {
    if (xmax <= xmin) return 0.;

    constexpr int nsteps = 4000;
    const double step = (xmax - xmin) / nsteps;
    double area = 0.;

    for (int i = 0; i <= nsteps; ++i) {
      const double x = xmin + i * step;
      Dmass->setVal(x);
      const double y = pdf_sig->getVal(&mass_obs);
      const double coeff = (i == 0 || i == nsteps) ? 0.5 : 1.0;
      area += coeff * y;
    }

    return area * step;
  };

  const double total_area = integrate(mass_min, mass_max);
  if (total_area <= 0.) {
    Dmass->setVal(original_mass);
    return -1.;
  }

  double low = 0.;
  double high = std::min(mass_mean - mass_min, mass_max - mass_mean);

  for (int iter = 0; iter < 80; ++iter) {
    const double mid = 0.5 * (low + high);
    const double area = integrate(mass_mean - mid, mass_mean + mid);
    const double frac = area / total_area;

    if (frac < area_frac) {
      low = mid;
    } else {
      high = mid;
    }
  }

  Dmass->setVal(original_mass);
  return 0.5 * (low + high);
}

std::vector<std::pair<RooRealVar*, bool>> fix_shape_parameters(RooAbsPdf* pdf, RooAbsData* data, const std::vector<std::string> &floating_yields) {
  std::vector<std::pair<RooRealVar*, bool>> old_states;
  std::unique_ptr<RooArgSet> pars(pdf->getParameters(*data)); // why need a dataset?
  for (RooAbsArg *arg : *pars) {
    auto *var = dynamic_cast<RooRealVar*>(arg);
    if (!var) continue;
    
    bool is_yield = false;
    for (const auto &name : floating_yields) {
      if (name == var->GetName()) {
        is_yield = true;
        break;
      }
    }
    if (is_yield) continue;

    __XJJLOG << ">> " << var->GetName() << (var->isConstant() ? " : Constant" : "") << std::endl;
    old_states.emplace_back(var, var->isConstant());
    var->setConstant(true);
  }
  return old_states;
}

void restore_parameter_states(const std::vector<std::pair<RooRealVar*, bool>> &old_states) {
  for (const auto& [var, was_constant] : old_states) {
    var->setConstant(was_constant);
  }
}
