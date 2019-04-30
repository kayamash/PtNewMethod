#include <vector>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBranch.h>
#include "PtNewMethod.cpp"

const string inputfilelist = "/home/kayamash/efflist/20190416data18_physics_Main_Ztap.list";
const string outputfilename = "/gpfs/fs6001/kayamash/Mywork/efficiencyloopoutput/LUT/data18_physics_Main_JPZtap.root";
const string LUTnameAlpha = "/gpfs/fs7001/kayamash/Mywork/efficiencyloopoutput/LUT/kayamashNewMethodAlpha.LUT";
const string LUTnameBeta = "/gpfs/fs7001/kayamash/Mywork/efficiencyloopoutput/LUT/kayamashNewMethodBeta.LUT";
const string trigger1 = "mu4";//JPsimumu
const string trigger2 = "mu26ivm";//Zmumu
const Int_t proc1 = 1;//Jpsitap = 1,Ztap = 3
const Int_t proc2 = 3;//Jpsitap = 1,Ztap = 3
const bool EventFullScan = kTRUE;//kTRUE = All Event,kFALSE = 1000000 events test

void run(){
  cout<<"start!"<<endl;
  //Add chain
  TChain *chain = new TChain("t_tap");
  std::ifstream ifs(inputfilelist.c_str());
  std::string str;
  while(getline(ifs,str)){
  	chain->Add(str.c_str());
  }

  if(!chain)cout<<"tree failed!"<<endl;
  cout<<"Total Events are "<<chain->GetEntries()<<endl;
  PtNewMethod m(chain);
  const Int_t events = (EventFullScan) ? (chain->GetEntries()) : (1000000);
  cout<<"loop start!"<<endl;
  for(Int_t i = 0;i < events;i++){
  	if(i%1000000 == 0) cout<<"The event is "<<i<<endl;
    m.Loop(i,trigger1,proc1);
    m.Loop(i,trigger2,proc2);
  }

  TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
  cout<<"finalize"<<endl;
  m.Finalize(output_file,LUTnameAlpha,LUTnameBeta);
  delete  output_file;
  delete chain;
}

