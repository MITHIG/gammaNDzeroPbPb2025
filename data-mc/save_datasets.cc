#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "xjjanauti.h"
#include "xjjstruct.h"

#define __VARIABLES_ROOSPLOT__
#include "variables.h"
#define __BINS_MASS__
#include "../include/bins.h"

struct Flatten {
  std::vector<float> *br;
  RooRealVar *roov;
};

enum class EcutPreset { none = 0, gammaN = 1, Ngamma = 2, twoDirs = 3 };
std::vector<std::string> ecut_name = { "none", "gammaN", "Ngamma", "gammaN + Ngamma" }; // only for print out
enum class GCutPreset { none = 0, match = 1, swap = 2 };
std::vector<std::string> gcut_name = { "none", "match", "swap" };
enum class DcutPreset { none = 0, BDT = 1, Loose = 2 };
std::vector<std::string> dcut_name = { "none", "Analysis BDT cut", "Loose cut" }; // only for print out

std::unique_ptr<RooDataSet> make_dataset(TTree* tree, const std::string& name, EcutPreset ecut, DcutPreset dcut, GCutPreset gcut = GCutPreset::none) {
  __XJJLOG << ">>                     name: " << name << std::endl;
  __XJJLOG << ">> event selection category: " << ecut_name[static_cast<int>(ecut)] << std::endl;
  __XJJLOG << ">> event selection category: " << dcut_name[static_cast<int>(dcut)] << std::endl;
  __XJJLOG << ">> gen-match       category: " << gcut_name[static_cast<int>(gcut)] << std::endl;  

  tree->SetBranchStatus("*", 0);

  __XJJLOG << "++ register variables" << std::endl;
  
  RooArgSet observables; // a set of RooRealVar
  std::map<std::string, Flatten> vars;
  for (auto& v : variables) {
    if (v.isbranch < 0) continue;

    __XJJLOG << "   >> " << v.varname << " // " << v.var << (v.isbranch ? "" : " \e[33m(no branch, need to calculate later)\e[0m") << std::endl;
    // create roorealvar
    vars[v.varname].roov = new RooRealVar(v.varname.c_str(), v.vartex.c_str(), v.varmin, v.varmax);
    observables.add(*(vars[v.varname].roov));

    // set branch address
    vars[v.varname].br = nullptr;
    if (!v.isbranch) continue;
    tree->SetBranchStatus(v.var.c_str(), 1);
    tree->SetBranchAddress(v.var.c_str(), &(vars[v.varname].br));
  }

  // branches to use 
#define SET_BRANCH(q, t, d)                     \
  t q = d;                                      \
  if (tree->GetBranch( #q )) {                  \
    tree->SetBranchStatus( #q , 1);             \
    tree->SetBranchAddress( #q , &q);           \
  }

  SET_BRANCH(Run, int, 0);
  SET_BRANCH(isL1ZDCOr, bool, true);
  SET_BRANCH(cscTightHalo2015Filter, bool, true);
  SET_BRANCH(selectedVtxFilter, bool, true);
  SET_BRANCH(ZDCgammaN, bool, true);
  SET_BRANCH(ZDCNgamma, bool, true);
  SET_BRANCH(HFEMaxPlus_eta5, float, 0.);
  SET_BRANCH(HFEMaxMinus_eta5, float, 0.);
  SET_BRANCH(Dsize, int, 0);

#define SET_BRANCH_VECTOR(q, t)                 \
  std::vector<t>* q = nullptr;                  \
  if (tree->GetBranch( #q )) {                  \
    tree->SetBranchStatus( #q , 1);             \
    tree->SetBranchAddress( #q , &q);           \
  }

  SET_BRANCH_VECTOR(Dgen, int);

  auto data = std::make_unique<RooDataSet>(name, "", observables);
  auto nentries = tree->GetEntries();
  for (long long int i=0; i<nentries; i++) {
    xjjc::progressslide(i, nentries, 100000);
    tree->GetEntry(i);

    if (!selectedVtxFilter) continue;
    if (Run > 10 && !(isL1ZDCOr && cscTightHalo2015Filter)) continue;
    bool evt_pass_gammaN = ZDCgammaN && HFEMaxPlus_eta5 < 16;
    bool evt_pass_Ngamma = ZDCNgamma && HFEMaxMinus_eta5 < 16;
    if (ecut == EcutPreset::gammaN && !evt_pass_gammaN) continue;
    if (ecut == EcutPreset::Ngamma && !evt_pass_Ngamma) continue;
    if (ecut == EcutPreset::twoDirs && !(evt_pass_gammaN || evt_pass_Ngamma)) continue;

    // std::cout<<Dsize<<std::endl;
    
#define VAL(q) vars[ #q ].br->at(j)
    for (int j=0; j<Dsize; j++) {
      // cut
      if (VAL(Dpt) < 2. || VAL(Dpt) > 5. || VAL(Dy) < -2. || VAL(Dy) > 2.) continue;
      
      if (gcut == GCutPreset::match && Dgen->at(j) != 23333) continue;
      if (gcut == GCutPreset::swap && Dgen->at(j) != 23344) continue;

      if (VAL(Dmass) < bins::minmass || VAL(Dmass) > bins::maxmass) continue; // can be reduced?

      if (!( std::abs(VAL(Dtrk1PtErr)/VAL(Dtrk1Pt)) < 0.1 && std::abs(VAL(Dtrk2PtErr)/VAL(Dtrk2Pt)) < 0.1 &&
             std::abs(VAL(Dtrk1Eta)) < 2.4 && std::abs(VAL(Dtrk2Eta)) < 2.4 &&
             VAL(Dtrk1Pt) > 0.5 && VAL(Dtrk2Pt) > 0.5 &&
             VAL(Dchi2cl) > 0.05 && (VAL(DsvpvDistance)/VAL(DsvpvDisErr)) > 1. )) continue;

      bool d_pass_gammaN = true;
      if (dcut == DcutPreset::BDT) d_pass_gammaN = ((VAL(Dy)<-1 && VAL(Dmva_BDT)>0.143) || (VAL(Dy)>=-1 && VAL(Dy)<0 && VAL(Dmva_BDT)>0.142) || (VAL(Dy)>=0 && VAL(Dy)<1 && VAL(Dmva_BDT)>0.123) || (VAL(Dy)>=1 && VAL(Dmva_BDT)>0.098));
      else if (dcut == DcutPreset::Loose) d_pass_gammaN = VAL(Dmva_BDT) > 0.;
      bool d_pass_Ngamma = true;
      if (dcut == DcutPreset::BDT) d_pass_Ngamma = ((VAL(Dy)>=1 && VAL(Dmva_BDT)>0.143) || (VAL(Dy)<1 && VAL(Dy)>=0 && VAL(Dmva_BDT)>0.142) || (VAL(Dy)<0 && VAL(Dy)>=-1 && VAL(Dmva_BDT)>0.123) || (VAL(Dy)<-1 && VAL(Dmva_BDT)>0.098));
      else if (dcut == DcutPreset::Loose) d_pass_Ngamma = VAL(Dmva_BDT) > 0.;

      const bool all_pass_gammaN = evt_pass_gammaN && d_pass_gammaN;
      const bool all_pass_Ngamma = evt_pass_Ngamma && d_pass_Ngamma;
      
      if (ecut == EcutPreset::gammaN && !all_pass_gammaN) continue;
      if (ecut == EcutPreset::Ngamma && !all_pass_Ngamma) continue;
      if (ecut == EcutPreset::twoDirs && !(all_pass_gammaN || all_pass_Ngamma)) continue;

      // set dataset values
      for (auto& [_, v] : vars) {
        if (v.br) {
          v.roov->setVal( v.br->at(j) );
          if (ecut == EcutPreset::twoDirs && evt_pass_Ngamma) {
            bool need_refl = false;
            for (const std::string& str_vref : { "Dy", "Eta" }) {
              if (xjjc::str_contains(v.roov->GetName(), str_vref)) {
                need_refl = true;
                break;
              }
            }
            if (need_refl)
              v.roov->setVal( 0. - v.br->at(j) ); // reflection
          }
        }
      }

      // !! add complicated variables
      vars.at("Ddls").roov->setVal(VAL(DsvpvDistance) / VAL(DsvpvDisErr));
      vars.at("Ddls_2D").roov->setVal(VAL(DsvpvDistance_2D) / VAL(DsvpvDisErr_2D));
      vars.at("Dip3Dsig").roov->setVal(VAL(Dip3D) / VAL(Dip3derr));

      data->add(observables);
    }
  }
  xjjc::progressbar_summary(nentries);

  __XJJLOG << ">> " << data->GetName() << " numEntries: " << data->numEntries() << std::endl;

  return data;
}

int macro(const std::string& inputstr, const std::string& outputname, const std::string& ecutstr, const std::string& dcutstr, int ismcref) {
  // parse cut
  auto ecuts = xjjroot::parse_input(ecutstr);
  auto ecut = static_cast<EcutPreset>(std::atoi(ecuts.content.c_str()));
  auto dcuts = xjjroot::parse_input(dcutstr);
  auto dcut = static_cast<DcutPreset>(std::atoi(dcuts.content.c_str()));

  // parse inputs
  auto inputs = xjjroot::parse_input(inputstr);
  const auto infname = inputs.content;
  auto* inf = TFile::Open(infname.c_str());
  if (!inf || inf->IsZombie()) {
    __XJJLOG << "!! bad file: " << infname << ", abort." << std::endl;
    return 2;
  }
  auto* tree = dynamic_cast<TTree*>(inf->Get("Tree"));
  if (!tree) {
    __XJJLOG << "!! bad tree: Tree, abort." << std::endl;
    return 2;
  }

  std::vector<std::unique_ptr<RooDataSet>> datasets;
  if (ismcref) {
    datasets.push_back( make_dataset(tree, "mc_match", ecut, dcut, GCutPreset::match) );
    datasets.push_back( make_dataset(tree, "mc_swap", ecut, dcut, GCutPreset::swap) );
  } else {
    datasets.push_back( make_dataset(tree, "data_main", ecut, dcut) );
  }

  auto* outf = xjjroot::newfile("rootfiles/" + outputname + ".root");
  for (auto& d : datasets)
    d->Write(d->GetName());
  auto* t = new TTree("info", "");
  t->Branch("input", &inputs.content);
  t->Branch("input_tex", &inputs.tex);
  t->Branch("input_tag", &inputs.tag);
  t->Branch("ecut", &ecuts.content);
  t->Branch("ecut_tex", &ecuts.tex);
  t->Branch("ecut_tag", &ecuts.tag);
  t->Branch("dcut", &dcuts.content);
  t->Branch("dcut_tex", &dcuts.tex);
  t->Branch("dcut_tag", &dcuts.tag);
  t->Fill();
  t->Write();
  xjjroot::closefile(outf);

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 6) {
    return macro(argv[1], argv[2], argv[3], argv[4], std::atoi(argv[5]));
  }
}
