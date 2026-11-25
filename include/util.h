#pragma once

namespace util {
  int mirrorswap_hist(TH1 *h) {
    const auto nbins = h->GetNbinsX();
    if (nbins%2 != 0) {
      __XJJLOG << "!! the histogram can't be mirror swapped, abort." << std::endl;
      return 1;
    }
    for (int i = 1; i <= nbins / 2; ++i) {
      int j = nbins + 1 - i;
      const auto content_i = h->GetBinContent(i);
      const auto error_i   = h->GetBinError(i);
      const auto content_j = h->GetBinContent(j);
      const auto error_j   = h->GetBinError(j);
      h->SetBinContent(i, content_j);
      h->SetBinError(i, error_j);
      h->SetBinContent(j, content_i);
      h->SetBinError(j, error_i);
    }
    return 0;
  }
}

