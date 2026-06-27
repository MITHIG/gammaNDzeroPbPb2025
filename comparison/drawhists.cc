#include "xjjanauti.h"
#include "xjjmypdf.h"
#include "variables.h"

struct Input {
  std::map<std::string, TH1D*> h1s;
  std::map<std::string, std::string> info;
};

std::vector<Color_t> colors = { kBlack, xjjroot::mycolor_satmiddle["red"], xjjroot::mycolor_satmiddle["azure"], xjjroot::mycolor_satmiddle["green"] };

int macro(std::string inputnames, std::string output, int save_png = 0, std::string title = "PbPb (5.36 TeV)")
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
    check_consistency("cut_tex", cut_tex);
    //
    auto hists = xjjana::getobj_regexp<TH1D>(inf);
    for (auto& h : hists) {
      const xjjroot::thgrstyle style = inputs.size() ?
        xjjroot::thgrstyle{ colors[inputs.size()%colors.size()], 20, 1, colors[inputs.size()%colors.size()], 1, 2, 0, 0, 0, 1, 1 } :
        xjjroot::thgrstyle{ kBlack, 20, 1, kBlack, 1, 2, kBlack, 0.05, 1001, 0.3, 0.3 };
      auto name = xjjc::str_eraseall(h->GetName(), "h1_");

      h->GetYaxis()->SetTitle("Entries");
      xjjroot::sethempty(h);
      xjjroot::setthgrstyle(h, style);
      ii.h1s[name] = h;

      auto* hnorm = (TH1D*)h->Clone(Form("%s__norm", h->GetName()));
      hnorm->Scale(1./hnorm->Integral(), "width");
      hnorm->GetYaxis()->SetTitle("Self normalized");
      xjjroot::sethempty(h);
      xjjroot::setthgrstyle(h, style);
      ii.h1s[name + "__norm"] = hnorm;
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
  auto name_png = xjjc::str_replaceall(pdf->getfilename(), {{ "figspdf", "figs" }, { "/"+varname, "_"+varname }, { ".pdf", "" }});

  auto draw_hist_list = [&inputs, &pdf, &cut_tex, &title](std::string key, int logy, float pratio) {
    std::vector<TH1D*> result;
    pratio = std::min(pratio, (float)1.);
    float tsize = 0.033/pratio;
    auto* leg = new TLegend(0.86-tsize*5, 0.80-tsize*1.1*inputs.size(), 0.86, 0.80);
    xjjroot::setleg(leg, tsize);
    for (auto& ii : inputs) {
      result.push_back(ii.h1s.at(key));
      leg->AddEntry(ii.h1s.at(key), ii.info["input_tex"].c_str(), ii.h1s.at(key)->GetFillStyle()>0?"f":"l");
    }
    if (logy) xjjana::sethsnonzeromin(result, 0.5);
    else xjjana::sethsmin(result, 0.1);
    xjjana::sethsmax(result, logy ? 50. : 1.4);

    std::vector<TH1D*> result_ratio;
    for (auto& h : result) {
      auto* hratio = (TH1D*)h->Clone(Form("%s_ratio", h->GetName()));
      hratio->Divide(result.front());
      hratio->SetFillStyle(0);
      hratio->GetYaxis()->SetTitle("Ratio");
      if (result_ratio.size() == 0) {
        xjjana::sethunivalue(hratio, 1., 0.);
        hratio->SetLineStyle(2);
      }
      result_ratio.push_back(hratio);
    }
    xjjana::sethsnonzeromin(result_ratio, 0.95);
    xjjana::sethsmax(result_ratio, 1.05);    
    
    auto pads = xjjroot::twopads(pdf->getc(), result.front(), result_ratio.front(), pratio);
    pads[0]->SetLogy(logy);
    pads[0]->cd();
    for (auto& h : result) h->Draw("hist same");
    leg->Draw();
    xjjroot::drawCMS(xjjroot::CMS::internal, title.c_str(), 1./pratio);
    xjjroot::drawtexgroup(0.23, 0.86, xjjc::str_divide_trim(cut_tex, "%%"), tsize, 13);
    pads[1]->cd();
    for (auto& h : result_ratio) h->Draw("hist same");

    pdf->getc()->cd();
    return pads;
  };

  pdf->prepare();
  draw_hist_list("var", the_var.logy, 2./3);
  pdf->getc()->cd();
  if (save_png) {
    pdf->write(name_png + ".pdf");
  } else pdf->write();
  
  pdf->prepare();
  draw_hist_list("var__norm", the_var.logy, 2./3);
  pdf->getc()->cd();
  if (save_png) {
    pdf->write(name_png + "__norm.pdf");
  } else pdf->write();
  
  pdf->close();
  
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc == 4) {
    return macro(argv[1], argv[2], atoi(argv[3]));
  }
  return 1;
}
