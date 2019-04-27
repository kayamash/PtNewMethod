#include "PtNewMethod.C"

const string inputfilelist = "/home/kayamash/efflist/20190416data18_physics_Main_Ztap.list";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/LUT/data18_physics_Main_ZtapMU4.root";
const string LUTname = "/gpfs/fs7001/kayamash/Mywork/efficiencyloopoutput/LUT/kayamashNewMethod.LUT";
const string trigger = "mu4";//JPsimumu
//const string trigger = "mu26ivm";//Zmumu
const Int_t proc = 1;//Jpsitap = 1,Ztap = 3
const bool EventFullScan = kTRUE;//kTRUE = All Event,kFALSE = 1000000 events test

void run(){
  gROOT->LoadMacro("PtNewMethod.C");
  //Add chain
  TChain *chain = new TChain("t_tap");
  std::ifstream ifs(inputfilelist.c_str());
  std::string str;
  while(getline(ifs,str))chain->Add(str.c_str());
  
  PtNewMethod m(chain,trigger,proc);
  const Int_t events = (EventFullScan) ? (chain->GetEntries()) : (1000000);
  for(Int_t i = 0;i < events;i++){
    m.Loop(i);
  }

  TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
  m.Finalize(output_file,LUTname);
}

