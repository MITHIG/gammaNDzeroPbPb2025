// Flatten the D-candidate vectors in skim_HiForest_2025PbPbUPC_1.root.
//
// The output contains one entry per D candidate.  Scalar branches are copied
// for every candidate, D* vector branches are written as scalar branches, and
// vector branches whose names do not start with D are omitted.
//
// Run from ROOT with:
//   root -l -b -q 'flatten.C()'
// or choose input/output files explicitly:
//   root -l -b -q 'flatten.C("input.root", "flattened.root")'

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"

#include <iostream>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../include/xjjcuti.h"

namespace {

  struct VectorColumn {
    enum class Kind { Float, Double, Int, Bool };

    std::string name;
    Kind kind;
    std::vector<float>* vf = nullptr;
    std::vector<double>* vd = nullptr;
    std::vector<int>* vi = nullptr;
    std::vector<bool>* vb = nullptr;
    float f = 0.f;
    double d = 0.;
    int i = 0;
    bool b = false;

    void bind(TTree* tree) {
      switch (kind) {
      case Kind::Float:  tree->SetBranchAddress(name.c_str(), &vf); break;
      case Kind::Double: tree->SetBranchAddress(name.c_str(), &vd); break;
      case Kind::Int:    tree->SetBranchAddress(name.c_str(), &vi); break;
      case Kind::Bool:   tree->SetBranchAddress(name.c_str(), &vb); break;
      }
    }

    void makeBranch(TTree* tree) {
      switch (kind) {
      case Kind::Float:  tree->Branch(name.c_str(), &f, (name + "/F").c_str()); break;
      case Kind::Double: tree->Branch(name.c_str(), &d, (name + "/D").c_str()); break;
      case Kind::Int:    tree->Branch(name.c_str(), &i, (name + "/I").c_str()); break;
      case Kind::Bool:   tree->Branch(name.c_str(), &b, (name + "/O").c_str()); break;
      }
    }

    std::size_t size() const {
      switch (kind) {
      case Kind::Float:  return vf ? vf->size() : 0;
      case Kind::Double: return vd ? vd->size() : 0;
      case Kind::Int:    return vi ? vi->size() : 0;
      case Kind::Bool:   return vb ? vb->size() : 0;
      }
      return 0;
    }

    void set(std::size_t index) {
      switch (kind) {
      case Kind::Float:  f = (*vf)[index]; break;
      case Kind::Double: d = (*vd)[index]; break;
      case Kind::Int:    i = (*vi)[index]; break;
      case Kind::Bool:   b = (*vb)[index]; break;
      }
    }
  };

  std::string normalizeType(std::string type) {
    const std::string prefix = "std::";
    while (type.compare(0, prefix.size(), prefix) == 0)
      type.erase(0, prefix.size());
    type.erase(0, type.find_first_not_of(" \t"));
    type.erase(type.find_last_not_of(" \t") + 1);
    return type;
  }

  bool isVectorType(const std::string& type) {
    const std::string t = normalizeType(type);
    return t.rfind("vector<", 0) == 0;
  }

} // namespace

void flatten(const char* input = "skim_HiForest_2025PbPbUPC_1.root",
             const char* output = "skim_HiForest_2025PbPbUPC_1_flattened.root") {
  std::unique_ptr<TFile> in(new TFile(input, "READ"));
  if (!in || in->IsZombie()) {
    std::cerr << "Cannot open input file: " << input << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(in->Get("Tree"));
  if (!tree) {
    std::cerr << "Input file does not contain a TTree named Tree" << std::endl;
    return;
  }

  // Keep a separate reader tree.  CloneTree(0) can alter branch state on its
  // source tree, so using a second handle makes the input addresses stable.
  std::unique_ptr<TFile> readFile(TFile::Open(input, "READ"));
  // std::unique_ptr<TFile> readFile(new TFile(input, "READ"));
  auto* readTree = readFile ? dynamic_cast<TTree*>(readFile->Get("Tree")) : nullptr;
  if (!readTree) {
    std::cerr << "Cannot create a second reader for Tree" << std::endl;
    return;
  }
  readTree->SetBranchStatus("*", 1);

  // Clone all non-vector branches.  Vector branches are added below only if
  // their names begin with D.
  tree->SetBranchStatus("*", 1);
  TObjArray* branches = tree->GetListOfBranches();
  std::vector<std::unique_ptr<VectorColumn>> columns;
  columns.reserve(branches->GetEntries());

  for (int ib = 0; ib < branches->GetEntries(); ++ib) {
    auto* branch = static_cast<TBranch*>(branches->At(ib));
    const std::string name = branch->GetName();
    const std::string type = normalizeType(branch->GetClassName());
    if (!isVectorType(type))
      continue;

    tree->SetBranchStatus(name.c_str(), 0);
    if (name.empty() || name[0] != 'D')
      continue;

    auto column = std::make_unique<VectorColumn>();
    column->name = name;
    if (type == "vector<float>")
      column->kind = VectorColumn::Kind::Float;
    else if (type == "vector<double>")
      column->kind = VectorColumn::Kind::Double;
    else if (type == "vector<int>")
      column->kind = VectorColumn::Kind::Int;
    else if (type == "vector<bool>")
      column->kind = VectorColumn::Kind::Bool;
    else {
      std::cerr << "Unsupported D vector branch " << name << " of type "
                << type << std::endl;
      continue;
    }
    columns.push_back(std::move(column));
  }

  // First scan the complete input. A D* branch is retained only if its
  // vector length is exactly Dsize for every event.
  int dsize = 0;
  if (!readTree->GetBranch("Dsize")) {
    std::cerr << "Input tree does not contain the required Dsize branch"
              << std::endl;
    return;
  }
  readTree->SetBranchAddress("Dsize", &dsize);
  for (const auto& column : columns)
    column->bind(readTree);

  __XJJLOG << "++ find bad D* branches whose size are not Dsize" << std::endl;
  std::vector<bool> bad(columns.size(), false);
  for (Long64_t entry = 0; entry < readTree->GetEntries(); ++entry) {
    xjjc::progressslide(entry, readTree->GetEntries(), 100000);
    readTree->GetEntry(entry);
    for (std::size_t icolumn = 0; icolumn < columns.size(); ++icolumn) {
      if (static_cast<Long64_t>(columns[icolumn]->size()) != dsize)
        bad[icolumn] = true;
    }
  }
  xjjc::progressbar_summary(readTree->GetEntries());

  std::vector<std::string> deletedBranches;
  for (std::size_t icolumn = 0; icolumn < columns.size(); ++icolumn) {
    if (bad[icolumn])
      deletedBranches.push_back(columns[icolumn]->name);
  }
  std::size_t badIndex = 0;
  columns.erase(std::remove_if(columns.begin(), columns.end(),
                               [&bad, &badIndex](const auto&) {
                                 return bad[badIndex++];
                               }),
                columns.end());
  xjjc::print_vec_v(deletedBranches, 0);

  std::unique_ptr<TFile> out(TFile::Open(output, "RECREATE"));
  // std::unique_ptr<TFile> out(new TFile(output, "RECREATE"));
  if (!out || out->IsZombie()) {
    std::cerr << "Cannot create output file: " << output << std::endl;
    return;
  }

  auto* flat = tree->CloneTree(0);
  flat->SetName("Tree");
  flat->SetTitle("Flattened D-candidate tree");
  for (const auto& column : columns) {
    readTree->SetBranchStatus(column->name.c_str(), 1);
    column->bind(readTree);
    column->makeBranch(flat);
  }

  const Long64_t entries = tree->GetEntries();
  Long64_t outputEntries = 0;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    // The cloned scalar branches remain attached to tree, while D vectors
    // are attached to readTree.  Advance both so each output candidate gets
    // the event-level values from its source event.
    xjjc::progressslide(entry, entries, 10000);
    tree->GetEntry(entry);
    readTree->GetEntry(entry);
    const std::size_t nCandidates = dsize < 0 ? 0 : dsize;

    for (std::size_t candidate = 0; candidate < nCandidates; ++candidate) {
      for (const auto& column : columns) {
        column->set(candidate);
      }
      flat->Fill();
      ++outputEntries;
    }
  }
  xjjc::progressbar_summary(entries);

  // Preserve the auxiliary metadata tree, if present.
  if (auto* info = dynamic_cast<TTree*>(in->Get("InfoTree")))
    info->CloneTree(-1, "fast")->Write();

  flat->Write();
  out->Write();
  out->Close();
  
  std::cout << "Wrote " << outputEntries << " flattened entries to "
            << output << std::endl;
  std::cout << "Deleted D* branches with a size mismatch:" << std::endl;
  if (deletedBranches.empty())
    std::cout << "  (none)" << std::endl;
  else
    for (const auto& name : deletedBranches)
      std::cout << "  " << name << std::endl;

  
}


int main(int argc, char* argv[]) {
  if (argc == 3) {
    flatten(argv[1], argv[2]);
    return 0;
  }
  return 1;
}
