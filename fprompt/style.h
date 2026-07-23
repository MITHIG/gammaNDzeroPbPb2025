#pragma once

#ifdef __COOK_NAME__

std::string add_suffix__y(const std::string& name, const std::string& suffix) {
  auto tag = xjjc::str_erasestar(name, "__y-*");
  auto ypart = xjjc::str_eraseall(name, tag);
  std::string newname(Form("%s%s%s", tag.c_str(), suffix.c_str(), ypart.c_str()));
  // __XJJLOG << " >> " << name << " + " << suffix << " -> " << newname << std::endl;
  return newname;  
}

#include <stdexcept>
#include <cctype>

int extract_index(const std::string& name, const std::string& key = "sf-") {
  auto pos = name.find(key);
  if (pos == std::string::npos)
    throw std::runtime_error("Pattern <" + key + "> not found");

  pos += key.size();
  auto end = pos;
  while (end < name.size() &&
         std::isdigit(static_cast<unsigned char>(name[end]))) {
    ++end;
  }
  if (end == pos)
    throw std::runtime_error("No digits after <" + key + ">");
  return std::stoi(name.substr(pos, end - pos));
}

template<class T>
int index_sf(const T* h) {
  if (!h) { return -1; }
  return extract_index(h->GetName(), "sf-");
};

template<class T>
int index_y(const T* h) {
  if (!h) { return -1; }
  return extract_index(h->GetName(), "__y-");
};

#endif

struct Style {
  std::string title = "";
  Color_t color = kBlack;
  Style_t mstyle = 20;
  Style_t lstyle = 1;
};

const std::map<std::string, Style> m_style_data = {
  { "sub", Style{ .title = "Sideband sub", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "sigswap", Style{ .title = "sPlot", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
  { "fix", Style{ .title = "MC template", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "best", Style{ .title = "Best#scale[0.5]{ }#chi^{2}", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
  { "Dip3D", Style{ .title = "DCA", .color = xjjroot::mycolor_middle["blue"], .mstyle = 20, .lstyle = 1 } },
  { "Dip3Dsig", Style{ .title = "DCA /#scale[0.4]{ }#sigma_{DCA}", .color = xjjroot::mycolor_middle["red"], .mstyle = 21, .lstyle = 1 } },
};

Style style_data(const std::string& name, int exact = 0) {
  Style empty;
  for (const auto& [key, style] : m_style_data) {
    if (exact && name == key)
      return style;
    if (!exact && xjjc::str_contains(name, key))
      return style;
  }
  return empty;
}
