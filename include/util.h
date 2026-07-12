#include "xjjcuti.h"

namespace util {
  struct Inputpar {
    std::string content;
    std::string tex;
    std::string tag;
    std::vector<std::string> parse;
  };

  Inputpar parse_input(std::string inputname);
}

util::Inputpar util::parse_input(std::string inputname) {
  Inputpar p = { .content = "", .tex = "", .tag = "",
    .parse = xjjc::str_divide_trim(inputname, ";") };
  if (p.parse.size() > 0)
    p.content = p.parse[0];
  if (p.parse.size() > 1)
    p.tex = p.parse[1];
  if (p.parse.size() > 2)
    p.tag = p.parse[2];

  return p;
}

