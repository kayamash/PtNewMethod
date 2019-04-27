#include "PtNewMethod.C"

const string inputfilelist = "/home/kayamash/efflist/20190416data18_physics_Main_Ztap.list";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/20190407/20190407data18_physics_Main_ZtapMU6.root";

void run(){
  gROOT->LoadMacro("PtNewMethod.C");
  //Add chain
  TChain *chain = new TChain("t_tap");
  std::ifstream ifs(inputfilelist.c_str());
  std::string str;
  while(getline(ifs,str))chain->Add(str.c_str());
  
  PtNewMethod m(chain);
  m.Loop();
}

