#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <sstream>
#include <TF1.h>
#include <TTree.h>
#include <TStyle.h>
#include <TText.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <vector>
#include <TMath.h>
#include <TBranch.h>
#include <TObject.h>
#include <TProfile.h>
#include "kayamashForLUT.chh"

bool kayamashForLUT::getLUTparameter(Double_t address,Double_t charge,Double_t eta,Double_t phi,Int_t (&par)[4]){
	for(Int_t i = 0; i < 4; ++i)par[i] = -1;
	if(address == 0.)par[0] = 0;//Large
	if(address == 1. && phi < -1.5)par[0] = 1;//LS sector 11
	if(address == 1. && phi > -1.5)par[0] = 2;//LS sector 15
	if(address == 2.)par[0] = 3;//Small
	if(address == 3.)par[0] = 4;//SS

	if(charge == 1.)par[1] = 0;//positive
	if(charge == -1.)par[1] = 1;//negative
	
	Int_t tmp_eta = static_cast<Int_t>(std::fabs(eta)*15/1.05);
    par[2] = (eta > 0) ? (tmp_eta) : (tmp_eta + 15);//Eta　divide to 30

	Double_t tmp_phi = phi + TMath::Pi();
    if(address == 0 || address == 1){//Large
    	while(tmp_phi > 3*TMath::Pi()/8.)tmp_phi -= TMath::Pi()/4.;
    	while(tmp_phi < TMath::Pi()/8.)tmp_phi += TMath::Pi()/4.;
    	tmp_phi -= TMath::Pi()/8.;
    }else if(address == 2 || address == 3){//Small
    	while(tmp_phi > TMath::Pi()/4.)tmp_phi -= TMath::Pi()/4.;
    }
    par[3] = static_cast<Int_t>(tmp_phi*120./TMath::Pi()); //Phi divide to 30

    if(par[0] >= 0 && par[1] >= 0 && par[2] >= 0 && par[3] >= 0)return kTRUE;
    return kFALSE;
}

bool kayamashForLUT::WriteLUT(TProfile *prof,Int_t par1,Int_t par2,Int_t par3,Int_t par4,bool alpha,std::string filename){
	Double_t fOrder = 0;
	Double_t sOrder = 0;
	if(alpha){//alpha
		if(par2 == 1){//charge
			switch(par1){//Sector
				case 0:
				fOrder = -0.0184048;
				sOrder = -3.62884;
				break;
				case 1:
				fOrder = -0.201341;
				sOrder = -2.69654;
				break;
				case 2:
				fOrder = -0.118256;
				sOrder = -2.91538;
				break;
				case 3:
				fOrder = -0.373771;
				sOrder = -4.20018;
				break;
				case 4:
				fOrder = -0.276773;
				sOrder = -5.23416;
				break;
			}
		}else{
			switch(par1){//Sector
				case 0:
				fOrder = 0.147342;
				sOrder = 2.77986;
				break;
				case 1:
				fOrder = -0.0414529;
				sOrder = 4.40735;
				break;
				case 2:
				fOrder = 0.191687;
				sOrder = 2.6177;
				break;
				case 3:
				fOrder = 0.250872;
				sOrder = 3.74082;
				break;
				case 4:
				fOrder = 0.225341;
				sOrder = 3.82399;
				break;
			}
		}
	}else{//beta
		if(par2 == 1){//charge
			switch(par1){//Sector
				case 0:
				fOrder = -0.0225169;
				sOrder = 4.08316;
				break;
				case 1:
				fOrder = -0.0044969;
				sOrder = 1.75247;
				break;
				case 2:
				fOrder = 0.0835498;
				sOrder = 1.01082;
				break;
				case 3:
				fOrder = 0.207164;
				sOrder = 7.40479;
				break;
				case 4:
				fOrder = 0.275761;
				sOrder = 7.10336;
				break;
			}
		}else{
			switch(par1){//Sector
				case 0:
				fOrder = -0.156158;
				sOrder = -2.28754;
				break;
				case 1:
				fOrder = 0.319591;
				sOrder = -4.28552;
				break;
				case 2:
				fOrder = 0.149876;
				sOrder = -2.47795;
				break;
				case 3:
				fOrder = -0.489921;
				sOrder = -3.20891;
				break;
				case 4:
				fOrder = -0.575507;
				sOrder = -2.90169;
				break;
			}
		}
	}
	ofstream ofs;
	ofs.open(filename.c_str(),std::ios::app);
	prof->Draw();
	Double_t p0 = 0;
	Double_t p1 = 0;
	Double_t chi = 0;
	if(prof->GetEntries() != 0){
		TF1 *fitProf = new TF1("fitProf","[0]*x + [1]*x*x",0.,0.3);
		fitProf->FixParameter(0,0.);
		fitProf->SetParameter(1,fOrder);
		fitProf->SetParameter(2,sOrder);
		prof->Fit(fitProf,"I","",0.,0.3);
		p0 = fitProf->GetParameter(1);
		p1 = fitProf->GetParameter(2);
		chi = fitProf->GetChisquare();
	}

	ofs<<par1<<"   "<<par2<<"   "<<par3<<"   "<<par4<<"   "<<p0<<"   "<<p1<<"   "<<chi<<std::endl;
	ofs.close();
	return kTRUE;
}