#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "xjjanauti.h"

#define __VARIABLES_ROOSPLOT__
#include "variables.h"
#define __BINS_MASS__
#include "../include/bins.h"
#include "../include/util.h"

struct Flatten {
  std::vector<float> *br;
  RooRealVar *roov;
};

enum class ECutPreset { none = 0, gammaN = 1, Ngamma = 2 };
std::vector<std::string> ecut_name = { "none", "gammaN", "Ngamma" };
enum class GCutPreset { none = 0, match = 1, swap = 2 };
std::vector<std::string> gcut_name = { "none", "match", "swap" };

std::unique_ptr<RooDataSet> make_dataset(TTree* tree, std::string name, ECutPreset ecut, GCutPreset gcut = GCutPreset::none) {
  __XJJLOG << ">>                     name: " << name << std::endl;
  __XJJLOG << ">> event selection category: " << ecut_name[static_cast<int>(ecut)] << std::endl;
  __XJJLOG << ">> gen-match       category: " << gcut_name[static_cast<int>(gcut)] << std::endl;  

  tree->SetBranchStatus("*", 0);

  __XJJLOG << "++ register variables" << std::endl;
  
  RooArgSet observables; // a set of RooRealVar
  std::map<std::string, Flatten> vars;
  for (auto& v : variables) {
    if (v.isbranch < 0) continue;

    __XJJLOG << "   >> " << v.varname << " // " << v.var << (v.isbranch ? " (to set branch)" : "") << std::endl;
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
    xjjc::progressslide(i, nentries, 10000);
    tree->GetEntry(i);

    if (!selectedVtxFilter) continue;
    if (Run > 10 && !(isL1ZDCOr && cscTightHalo2015Filter)) continue;
    if (ecut == ECutPreset::gammaN && !(ZDCgammaN && HFEMaxPlus_eta5 < 16)) continue;
    if (ecut == ECutPreset::Ngamma && !(ZDCNgamma && HFEMaxMinus_eta5 < 16)) continue;

    // std::cout<<Dsize<<std::endl;
    
#define VAL(q) vars[ #q ].br->at(j)
    for (int j=0; j<Dsize; j++) {
      // cut
      // for (auto& [_, v] : vars) {
      //   if (v.br) {
      //     std::cout<<"   "<<j<<"  "<<v.roov->GetName() << " (" << v.br->size() << ")" << std::endl;
      //   }
      // }

      if (VAL(Dpt) < 2. || VAL(Dpt) > 5. || VAL(Dy) < -2. || VAL(Dy) > 2.) continue;
      
      if (gcut == GCutPreset::match && Dgen->at(j) != 23333) continue;
      if (gcut == GCutPreset::swap && Dgen->at(j) != 23344) continue;

      if (VAL(Dmass) < bins::minmass || VAL(Dmass) > bins::maxmass) continue; // can be reduced

      if (!( std::abs(VAL(Dtrk1PtErr)/VAL(Dtrk1Pt)) < 0.1 && std::abs(VAL(Dtrk2PtErr)/VAL(Dtrk2Pt)) < 0.1 &&
             std::abs(VAL(Dtrk1Eta)) < 2.4 && std::abs(VAL(Dtrk2Eta)) < 2.4 &&
             VAL(Dtrk1Pt) > 0.5 && VAL(Dtrk2Pt) > 0.5 &&
             VAL(Dchi2cl) > 0.05 && (VAL(DsvpvDistance)/VAL(DsvpvDisErr)) > 1. )) continue;

      if (ecut == ECutPreset::gammaN &&
          !((VAL(Dy)<-1 && VAL(Dmva_BDT)>0.143) || (VAL(Dy)>=-1 && VAL(Dy)<0 && VAL(Dmva_BDT)>0.142) || (VAL(Dy)>=0 && VAL(Dy)<1 && VAL(Dmva_BDT)>0.123) || (VAL(Dy)>=1 && VAL(Dmva_BDT)>0.098))) continue;
      if (ecut == ECutPreset::Ngamma &&
          !((VAL(Dy)>=1 && VAL(Dmva_BDT)>0.143) || (VAL(Dy)<1 && VAL(Dy)>=0 && VAL(Dmva_BDT)>0.142) || (VAL(Dy)<0 && VAL(Dy)>=-1 && VAL(Dmva_BDT)>0.123) || (VAL(Dy)<-1 && VAL(Dmva_BDT)>0.098))) continue;

      // set dataset values
      for (auto& [_, v] : vars) {
        if (v.br) {
          // std::cout<<"   "<<j<<"  "<<v.roov->GetName() << " (" << v.br->size() << ")" << std::endl;
          v.roov->setVal( v.br->at(j) );
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

int macro(std::string inputstr, std::string outputname, std::string ecutstr, int ismcref) {
  // parse cut
  auto ecuts = util::parse_input(ecutstr);
  auto ecut = static_cast<ECutPreset>(std::atoi(ecuts.content.c_str()));

  // parse inputs
  auto inputs = util::parse_input(inputstr);
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
    datasets.push_back( make_dataset(tree, "mc_match", ecut, GCutPreset::match) );
    datasets.push_back( make_dataset(tree, "mc_swap", ecut, GCutPreset::swap) );
  } else {
    datasets.push_back( make_dataset(tree, "data_main", ecut) );
  }

  auto* outf = xjjroot::newfile(outputname + ".root");
  for (auto& d : datasets)
    d->Write(d->GetName());
  auto* t = new TTree("info", "");
  t->Branch("input", &inputs.content);
  t->Branch("input_tex", &inputs.tex);
  t->Branch("input_tag", &inputs.tag);
  t->Branch("cut", &ecuts.content);
  t->Branch("cut_tex", &ecuts.tex);
  t->Branch("cut_tag", &ecuts.tag);
  t->Fill();
  t->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 5) {
    return macro(argv[1], argv[2], argv[3], std::atoi(argv[4]));
  }
}
