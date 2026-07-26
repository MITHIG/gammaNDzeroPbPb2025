#include "xjjcuti.h"

namespace util {
  struct Inputpar {
    std::string content;
    std::string tex;
    std::string tag;
    // std::vector<std::string> parse;
  };

  Inputpar parse_input(std::string inputname);
}

util::Inputpar util::parse_input(std::string inputname) {
  auto parse = xjjc::str_divide_trim(inputname, ";");
  Inputpar p = { .content = "", .tex = "", .tag = "" };
  if (parse.size() > 0)
    p.content = parse[0];
  if (parse.size() > 1)
    p.tex = parse[1];
  if (parse.size() > 2)
    p.tag = parse[2];

  return p;
}

