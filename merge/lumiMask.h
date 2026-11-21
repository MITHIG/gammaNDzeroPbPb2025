#ifndef __LUMIMASK_H__
#define __LUMIMASK_H__

// Minimal CMS lumi-mask (golden JSON) reader, no external dependencies.
// Format: {"run": [[ls_lo, ls_hi], [ls_lo, ls_hi], ...], ...}

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>
#include <cctype>

class lumiMask {
public:
  lumiMask(const std::string& fjson) { parse(fjson); }

  // true if (run, lumi) is inside one of the certified ranges (inclusive)
  bool pass(int run, int lumi) const {
    auto it = mask_.find(run);
    if (it == mask_.end()) return false;
    for (const auto& r : it->second)
      if (lumi >= r.first && lumi <= r.second) return true;
    return false;
  }

  const std::map<int, std::vector<std::pair<int, int>>>& mask() const { return mask_; }

  void print() const {
    for (const auto& [run, ranges] : mask_) {
      std::cout << run << ":";
      for (const auto& r : ranges) std::cout << " [" << r.first << ", " << r.second << "]";
      std::cout << std::endl;
    }
  }

private:
  std::map<int, std::vector<std::pair<int, int>>> mask_;

  void parse(const std::string& fjson) {
    std::ifstream fin(fjson);
    if (!fin.is_open()) {
      std::cerr << "error: lumiMask: cannot open " << fjson << std::endl;
      return;
    }
    std::stringstream ss; ss << fin.rdbuf();
    const std::string s = ss.str();

    int depth = 0, run = -1;
    std::vector<int> nums; // numbers collected inside the current [lo, hi]
    for (size_t i = 0; i < s.size(); i++) {
      char c = s[i];
      if (c == '"') { // run number key
        size_t j = s.find('"', i + 1);
        run = std::stoi(s.substr(i + 1, j - i - 1));
        i = j;
      } else if (c == '[') {
        depth++;
        if (depth == 2) nums.clear();
      } else if (c == ']') {
        if (depth == 2 && nums.size() == 2)
          mask_[run].emplace_back(nums[0], nums[1]);
        else if (depth == 2)
          std::cerr << "warning: lumiMask: malformed range in run " << run << std::endl;
        depth--;
      } else if (std::isdigit(c) && depth == 2) {
        size_t j = i;
        while (j < s.size() && std::isdigit(s[j])) j++;
        nums.push_back(std::stoi(s.substr(i, j - i)));
        i = j - 1;
      }
    }
  }
};

#endif
