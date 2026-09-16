/**
   auto* t = new TTree("info", "");
   std::map<std::string, std::string> t_cont;
   auto cast_branch = [&t, &t_cont]<typename T>(const std::string& name, const T& x) {
   t_cont[name] = xjjc::to_string(x);
   t->Branch(name.c_str(), &(t_cont[name]));
   };
   cast_branch("input", pinput.content);
   t->Fill();
   t->Write();

   auto* t = new TTree("info", "");
   for (auto& [key, content] : info) {
   t->Branch(key.c_str(), &content);
   }
   t->Fill();
   t->Write();
   xjjroot::closefile(outf);
**/

/**
   auto info = xjjana::getval_regexp(static_cast<TTree*>(inf->Get("info")));
   __XJJLOG << "++ info" << std::endl;
   xjjc::print_tab(info, -1);
**/

