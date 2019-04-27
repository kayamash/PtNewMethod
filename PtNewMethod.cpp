#include "PtNewMethod.chh"
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

bool PtNewMethod::CutAll(Int_t Tagpass,Int_t L1pass){
	if(0.08 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && Tagpass > -1 && m_tag_proc == m_proc && m_poff_charge*m_tag_charge == -1){
    //if(m_sumReqdRL1 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && Tagpass > -1 && m_tag_proc == m_proc && L1pass > -1 && m_poff_charge*m_tag_charge == -1){
		return kTRUE;
	}else{
		return kFALSE;
	}
}

bool PtNewMethod::BarrelDicision(Float_t eta){
	if(std::fabs(m_poff_eta) < 1.05){
		return kTRUE;
	}else{
		return kFALSE;
	}
}

void PtNewMethod::Loop(Int_t ev){
	cout<<ev<<endl;
	tChain->GetEntry(ev);
	cout<<ev<<endl;
	Double_t pextL1_dR = 1; 
	Double_t pextSA_dR = 1; 
	Double_t pL1_pt = -99999;
	Double_t pL1_eta = 0;
	Double_t pL1_phi = 0;
	Int_t pL1_pass = 0;
	Double_t pSA_pt = -99999;
	Double_t pSA_eta = 0;
	Double_t pSA_phi = 0;
	Double_t pSA_dR = 1;
	Int_t pSA_pass = 0;
	Double_t pSA_sAddress = -1;
	float pSA_roieta = -99999;
	float pSA_roiphi = -99999;
	Double_t pSA_superpointZ_BI = 0;
	Double_t pSA_superpointZ_BM = 0;
	Double_t pSA_superpointZ_BO = 0;
	Double_t pSA_superpointZ_BME = 0;
	Double_t pSA_superpointZ_BEE = 0;
	Double_t pSA_superpointR_BI = 0;
	Double_t pSA_superpointR_BM = 0;
	Double_t pSA_superpointR_BO = 0;
	Double_t pSA_superpointR_BME = 0;
	Double_t pSA_superpointR_BEE = 0;
	Double_t pSA_superpointSlope_BI = 0;
	Double_t pSA_superpointSlope_BM = 0;
	Int_t pEFTAG_pass = -1;
	for(Int_t method = 0;method < 25;method++){
		if(m_mes_name->at(method) == m_method_name){
			pL1_pt = m_pL1_pt->at(method);
			pSA_pt = m_pSA_pt->at(method);
			pL1_eta = m_pL1_eta->at(method);
			pSA_eta = m_pSA_eta->at(method);
			pL1_phi = m_pL1_phi->at(method);
			pSA_phi = m_pSA_phi->at(method);
			pL1_pass = m_pL1_pass->at(method);
			pSA_pass = m_pSA_pass->at(method);
			pSA_dR = m_pSA_dR->at(method);
			pEFTAG_pass = m_pEFTAG_pass->at(method);
			pSA_sAddress = m_pSA_sAddress->at(method);
			pSA_roieta = m_pSA_roieta->at(method);
			pSA_roiphi = m_pSA_roiphi->at(method);
			pSA_superpointZ_BI = m_pSA_superpointZ_BI->at(method);
			pSA_superpointZ_BM = m_pSA_superpointZ_BM->at(method);
			pSA_superpointZ_BO = m_pSA_superpointZ_BO->at(method);
			pSA_superpointZ_BME = m_pSA_superpointZ_BME->at(method);
			pSA_superpointZ_BEE = m_pSA_superpointZ_BEE->at(method);
			pSA_superpointR_BI = m_pSA_superpointR_BI->at(method);
			pSA_superpointR_BM = m_pSA_superpointR_BM->at(method);
			pSA_superpointR_BO = m_pSA_superpointR_BO->at(method);
			pSA_superpointR_BME = m_pSA_superpointR_BME->at(method);
			pSA_superpointR_BEE = m_pSA_superpointR_BEE->at(method);
			pSA_superpointSlope_BI = m_pSA_superpointSlope_BI->at(method);
			pSA_superpointSlope_BM = m_pSA_superpointSlope_BM->at(method);
		}
	}
	if( !CutAll(pEFTAG_pass,pL1_pass) )return;

    //barrel alpha
	if(pSA_superpointR_BM != 0 && BarrelDicision(pSA_roieta) == kTRUE){
        Double_t barrelalpha = atan(pSA_superpointZ_BM/pSA_superpointR_BM) - atan(1.0/pSA_superpointSlope_BM);//Reciprocal number?
        m_h_BarrelAlpha->Fill(barrelalpha);
        m_h_PtvsBarrelAlpha->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    }
    //barrel alpha end

    //barrel beta
    if(pSA_superpointR_BI != 0 && pSA_superpointR_BM != 0 && BarrelDicision(pSA_roieta) == kTRUE){
        Double_t barrelbeta = atan(1.0/pSA_superpointSlope_BI) - atan(1.0/pSA_superpointSlope_BM);//Reciprocal number?
        m_h_BarrelBeta->Fill(barrelbeta);
        m_h_PtvsBarrelBeta->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    }
    //barrel beta end

}

void PtNewMethod::Finalize(TFile *tf1,std::string filename){
	tf1->cd();
	m_h_BarrelAlpha->Write();
	m_h_PtvsBarrelAlpha->Write();
	m_h_BarrelBeta->Write();
	m_h_PtvsBarrelBeta->Write();
	cout<<"finish!!"<<endl;
}
