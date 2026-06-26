#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "variables.h"

struct Input {
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, std::string> info;
};

std::vector<Color_t> colors = { kBlack, xjjroot::mycolor_satmiddle["red"], xjjroot::mycolor_satmiddle["azure"], xjjroot::mycolor_satmiddle["green"] };

int macro(std::string inputnames, std::string output)
{
  std::vector<Input> inputs;
  std::string varname, cut_tex;
  for (const auto& input : xjjc::str_divide_trim(inputnames, ",")) {
    auto* inf = TFile::Open(input.c_str());
    if (!inf) {
      __XJJLOG << "!! warning: bad input " << input << ", skip."<< std::endl;
      continue;
    }
    __XJJLOG << "++ " << input << std::endl;
    Input ii;
    // 
    ii.info = xjjana::getval_regexp((TTree*)inf->Get("info"));
    __XJJLOG << ">> info" << std::endl;
    xjjc::print_tab(ii.info, -1);
    auto check_consistency = [&ii](const std::string& key, std::string& holder) {
      if (holder.empty()) {
        holder = ii.info.at(key);
      } else if (holder != ii.info.at(key)) {
        __XJJLOG << "!! error: different " << key << ", abort." << std::endl;
        return 2;
      }
      return 0;
    };
    if (check_consistency("varname", varname)) return 2;
    if (check_consistency("cut_tex", cut_tex)) return 2;
    //
    auto hists = xjjana::getobj_regexp<TH1D>(inf);
    for (auto& h : hists) {
      xjjroot::sethempty(h);
      const xjjroot::thgrstyle style = inputs.size() ?
        xjjroot::thgrstyle{ colors[inputs.size()%colors.size()], 20, 1, colors[inputs.size()%colors.size()], 1, 2, 0, 0, 0, 1, 1 } :
        xjjroot::thgrstyle{ kBlack, 20, 1, kBlack, 1, 2, kBlack, 0.2, 1001, 0.3, 0.3 };
      xjjroot::setthgrstyle(h, style);
      ii.h1s[xjjc::str_eraseall(h->GetName(), "h1_")] = h;
    }
    xjjroot::print_tab(ii.h1s, 0);

    inputs.push_back(ii);
  }

  auto it_var = std::find_if(vars.begin(), vars.end(), [&varname](const xjjana::variable& v) {
    return v.varname == varname;
  });
  if (it_var == vars.end()) {
    __XJJLOG << "!! bad varname: " << varname << std::endl;
    return 2;
  }
  const auto& the_var = *it_var;
  
  xjjroot::setgstyle(1);
  auto* pdf = new xjjroot::mypdf("figspdf/" + output + "/" + varname + ".pdf");

  auto draw_hist_list = [&inputs, &pdf, &cut_tex](std::string key, int logy, float pratio) {
    std::vector<TH1D*> result;
    pratio = std::min(pratio, (float)1.);
    float tsize = 0.037/pratio;
    auto* leg = new TLegend(0.86-tsize*5, 0.80-tsize*1.1*inputs.size(), 0.86, 0.80);
    xjjroot::setleg(leg, tsize);
    for (auto& ii : inputs) {
      result.push_back(ii.h1s[key]);
      leg->AddEntry(ii.h1s[key], ii.info["input_tex"].c_str(), "f");
    }
    auto ymin = xjjana::sethsmin(result, logy ? 0.5 : 0);
    if (ymin == 0 && logy) {
      for (auto& h : result) {
        h->SetMinimum(1);
      }
    }
    xjjana::sethsmax(result, logy ? 5. : 1.4);

    std::vector<TH1D*> result_ratio;
    for (auto& h : result) {
      auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      hratio->Divide(result.front());
      hratio->SetFillStyle(0);
      hratio->GetYaxis()->SetTitle("Ratio");
      result_ratio.push_back(hratio);
    }
    xjjana::sethsmin(result_ratio, 0.9);
    xjjana::sethsmax(result_ratio, 1.05);    
    
    auto pads = xjjroot::twopads(pdf->getc(), result.front(), result_ratio.front(), pratio);
    pads[0]->SetLogy(logy);
    pads[0]->cd();
    for (auto& h : result) h->Draw("hist same");
    leg->Draw();
    xjjroot::drawCMS(xjjroot::CMS::internal, "2023 PbPb (5.36 TeV)", 1./pratio);
    xjjroot::drawtex(0.89, 0.82, cut_tex.c_str(), tsize, 31);
    pads[1]->cd();
    for (auto& h : result_ratio) h->Draw("hist same");

    pdf->getc()->cd();
    return pads;
  };

  pdf->prepare();
  draw_hist_list("var", the_var.logy, 2./3);
  pdf->getc()->cd();
  pdf->write();
  
  pdf->prepare();
  draw_hist_list("run", 0, 2./3);
  pdf->getc()->cd();
  pdf->write();
  
  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 3) {
    return macro(argv[1], argv[2]);
  }
  return 1;
}
