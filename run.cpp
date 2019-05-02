#include <vector>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TBranch.h>
#include "PtNewMethod.cpp"

const bool tsakaiMethod = kTRUE;
const string inputfiledata18list = "/home/kayamash/LUTlist/20190416data18_physics_Main_JPZtap.list";
const string inputfiledata17list = "/home/kayamash/LUTlist/20190430data17_physics_Main_JPZtap.list";
string outputfilename = "/gpfs/fs6001/kayamash/Mywork/LUT/data18_physics_Main_JPZtap.root";
string LUTnameAlpha = "/gpfs/fs7001/kayamash/Mywork/LUT/NewMethodAlphaJPZ.LUT";
string LUTnameBeta = "/gpfs/fs7001/kayamash/Mywork/LUT/kayamashNewMethodBetaJPZ.LUT";
const string triggermu4 = "mu4";//JPsimumu
const string triggermu26 = "mu26ivm";//Zmumu
const Int_t procJpsi = 1;//Jpsitap = 1,Ztap = 3
const Int_t procZ = 3;//Jpsitap = 1,Ztap = 3
const bool EventFullScan = kTRUE;//kTRUE = All Event,kFALSE = 1000000 events test
const bool usingdata18 = kTRUE;
const bool usingdata17 = kTRUE;
const bool usingJpsi = kTRUE;
const bool usingZ = kTRUE;

void run(){
  cout<<"start!"<<endl;
  if(tsakaiMethod){
    outputfilename = "/gpfs/fs6001/kayamash/Mywork/LUT/tsakai/data18_physics_Main_JPZtap.root";
    LUTnameAlpha = "/gpfs/fs7001/kayamash/Mywork/LUT/tsakai/NewMethodAlphaJPZ.LUT";
    LUTnameBeta = "/gpfs/fs7001/kayamash/Mywork/LUT/tsakai/NewMethodBetaJPZ.LUT";
  }

  //Add chain
  TChain *chain = new TChain("t_tap");
  if(usingdata18){
    std::ifstream ifs(inputfiledata18list.c_str());
    std::string str;
    while(getline(ifs,str)){
  	  chain->Add(str.c_str());
    }
  }
  if(usingdata17){
    std::ifstream ifs(inputfiledata17list.c_str());
    std::string str;
    while(getline(ifs,str)){
      chain->Add(str.c_str());
    }
  }

  if(!chain)cout<<"tree failed!"<<endl;
  cout<<"Total Events are "<<chain->GetEntries()<<endl;
  PtNewMethod m(chain);
  if(tsakaiMethod){
    m.Init(16,15);
  }else{
    m.Init(30,30);
  }
  const Int_t events = (EventFullScan) ? (chain->GetEntries()) : (1000000);
  cout<<"loop start!"<<endl;
  for(Int_t i = 0;i < events;i++){
  	if(i%1000000 == 0) cout<<"The event is "<<i<<endl;
    if(usingJpsi)m.Loop(i,triggermu4,procJpsi);
    if(usingZ)m.Loop(i,triggermu26,procZ);
  }

  TFile *output_file = new TFile(outputfilename.c_str(),"RECREATE");
  cout<<"finalize"<<endl;
  m.Finalize(output_file,LUTnameAlpha,LUTnameBeta);
  delete  output_file;
  delete chain;
}

