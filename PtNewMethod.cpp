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
#include <TChain.h>

bool PtNewMethod::CutAll(Int_t Tagpass,Int_t L1pass){
	if(m_sumReqdRL1 < m_tp_extdR && 0.2 < m_tp_extdR && m_sumReqdREF < m_tp_dR && Tagpass > -1 && m_tag_proc == m_proc && m_poff_charge*m_tag_charge == -1){
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
	cout<<"check0"<<endl;
	tChain->GetEntry(ev);
	cout<<"check1"<<endl;
	Double_t pextL1_dR = 1; 
	Double_t pextSA_dR = 1; 
	Double_t pextCB_dR = 1; 
	Double_t pextEF_dR = 1; 
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
	Double_t pSA_phims = -99999;
	Double_t pSA_phibe = -99999;
	float pSA_roieta = -99999;
	float pSA_roiphi = -99999;
	vector<float> *pSA_rpcX = 0;
	vector<float> *pSA_rpcY = 0;
	vector<float> *pSA_rpcZ = 0;
	vector<float> *pSA_rpcR = 0;
	vector<float> *pSA_mdtZ = 0;
	vector<float> *pSA_mdtR = 0;
	vector<float> *pSA_mdtPhi = 0;
	vector<Int_t> *pSA_mdthitChamber = 0;
	Double_t pSA_ptTGC = -99999;
	Double_t pSA_ptalpha = -99999;
	Double_t pSA_ptbeta = -99999;
	Double_t pSA_superpointZ_BI = 0;
	Double_t pSA_superpointZ_BM = 0;
	Double_t pSA_superpointZ_BO = 0;
	Double_t pSA_superpointZ_BME = 0;
	Double_t pSA_superpointZ_EI = 0;
	Double_t pSA_superpointZ_EM = 0;
	Double_t pSA_superpointZ_EO = 0;
	Double_t pSA_superpointZ_EE = 0;
	Double_t pSA_superpointZ_CSC = 0;
	Double_t pSA_superpointZ_BEE = 0;
	Double_t pSA_superpointR_BI = 0;
	Double_t pSA_superpointR_BM = 0;
	Double_t pSA_superpointR_BO = 0;
	Double_t pSA_superpointR_BME = 0;
	Double_t pSA_superpointR_EI = 0;
	Double_t pSA_superpointR_EM = 0;
	Double_t pSA_superpointR_EO = 0;
	Double_t pSA_superpointR_EE = 0;
	Double_t pSA_superpointR_CSC = 0;
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
			pSA_phims = m_pSA_phims->at(method);
			pSA_phibe = m_pSA_phibe->at(method);
			pSA_roieta = m_pSA_roieta->at(method);
			pSA_roiphi = m_pSA_roiphi->at(method);
            	/*
            	pSA_rpcX = &(m_pSA_rpcX->at(method));
            	pSA_rpcY = &(m_pSA_rpcY->at(method));
            	pSA_rpcZ = &(m_pSA_rpcZ->at(method));
            	pSA_rpcR = &(m_pSA_rpcR->at(method));
            	pSA_mdtZ = &(m_pSA_mdtZ->at(method));
            	pSA_mdtR = &(m_pSA_mdtR->at(method));
            	pSA_mdtPhi = &(m_pSA_mdtPhi->at(method));
            	pSA_mdthitChamber = &(m_pSA_mdthitChamber->at(method));
            	*/
			pSA_ptTGC = m_pSA_pttgc->at(method);
			pSA_ptalpha = m_pSA_ptalpha->at(method);
			pSA_ptbeta = m_pSA_ptbeta->at(method);
			pSA_superpointZ_BI = m_pSA_superpointZ_BI->at(method);
			pSA_superpointZ_BM = m_pSA_superpointZ_BM->at(method);
			pSA_superpointZ_BO = m_pSA_superpointZ_BO->at(method);
			pSA_superpointZ_BME = m_pSA_superpointZ_BME->at(method);
			pSA_superpointZ_EI = m_pSA_superpointZ_EI->at(method);
			pSA_superpointZ_EM = m_pSA_superpointZ_EM->at(method);
			pSA_superpointZ_EO = m_pSA_superpointZ_EO->at(method);
			pSA_superpointZ_EE = m_pSA_superpointZ_EE->at(method);
			pSA_superpointZ_CSC = m_pSA_superpointZ_CSC->at(method);
			pSA_superpointZ_BEE = m_pSA_superpointZ_BEE->at(method);
			pSA_superpointR_BI = m_pSA_superpointR_BI->at(method);
			pSA_superpointR_BM = m_pSA_superpointR_BM->at(method);
			pSA_superpointR_BO = m_pSA_superpointR_BO->at(method);
			pSA_superpointR_BME = m_pSA_superpointR_BME->at(method);
			pSA_superpointR_EI = m_pSA_superpointR_EI->at(method);
			pSA_superpointR_EM = m_pSA_superpointR_EM->at(method);
			pSA_superpointR_EO = m_pSA_superpointR_EO->at(method);
			pSA_superpointR_EE = m_pSA_superpointR_EE->at(method);
			pSA_superpointR_CSC = m_pSA_superpointR_CSC->at(method);
			pSA_superpointR_BEE = m_pSA_superpointR_BEE->at(method);
			pSA_superpointSlope_BI = m_pSA_superpointSlope_BI->at(method);
			pSA_superpointSlope_BM = m_pSA_superpointSlope_BM->at(method);
		}
	}
	cout<<"check2"<<endl;
	if( !CutAll(pEFTAG_pass,pL1_pass) )return;
	cout<<ev<<endl;
    //segment
	Double_t segmentBISlope = 0;
	Double_t segmentBMSlope = 0;
	bool BIsegmentcheck = kFALSE;
	bool BMsegmentcheck = kFALSE;
	for(Int_t segmentNumber = 0; segmentNumber < 10; ++segmentNumber){
		Double_t tmp_segmentR = sqrt(m_probe_segment_x[segmentNumber]*m_probe_segment_x[segmentNumber] + m_probe_segment_y[segmentNumber]*m_probe_segment_y[segmentNumber]);
		Double_t tmp_segmentPR = sqrt(m_probe_segment_px[segmentNumber]*m_probe_segment_px[segmentNumber] + m_probe_segment_py[segmentNumber]*m_probe_segment_py[segmentNumber]);
		if(tmp_segmentR > 4000 && tmp_segmentR < 6500 && BIsegmentcheck == kFALSE){
			BIsegmentcheck = kTRUE;
			segmentBISlope = tmp_segmentPR/m_probe_segment_pz[segmentNumber];
		}
		if(tmp_segmentR > 6000 && tmp_segmentR < 9000 && BMsegmentcheck == kFALSE){
			BMsegmentcheck = kTRUE;
			segmentBMSlope = tmp_segmentPR/m_probe_segment_pz[segmentNumber];
		}
	}
	if(BIsegmentcheck)m_h_DeltaThetaBI->Fill(atan(segmentBISlope) - atan(1.0/pSA_superpointSlope_BI) );
	if(BMsegmentcheck)m_h_DeltaThetaBM->Fill(atan(segmentBMSlope) - atan(1.0/pSA_superpointSlope_BM) );

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

    //For LUT
    Double_t PhiIntegral = 0;
    if(pSA_sAddress == 0 || pSA_sAddress == 1){
    	if(pSA_phi < -2.5)PhiIntegral = pSA_phi + TMath::Pi();
    	if(pSA_phi > -2.5 || pSA_phi < -1.0)PhiIntegral = pSA_phi + TMath::Pi()/2.;
    	if(pSA_phi > -1.0 || pSA_phi < 0.24)PhiIntegral = pSA_phi;
    	if(pSA_phi > -0.24 || pSA_phi < 1.03)PhiIntegral = pSA_phi - TMath::Pi()/4.;
    	if(pSA_phi > 1.03 || pSA_phi < 1.80)PhiIntegral = pSA_phi - TMath::Pi()/2.;
    	if(pSA_phi > 1.80 || pSA_phi < 2.60)PhiIntegral = pSA_phi - 3*TMath::Pi()/4.;
    	if(pSA_phi > 2.60)PhiIntegral = pSA_phi - TMath::Pi();
    	m_h_LargePhi->Fill(PhiIntegral);
    }else if(pSA_sAddress == 2 || pSA_sAddress == 3){
    	if(pSA_phi < -2.0)PhiIntegral = pSA_phi + TMath::Pi();
    	if(pSA_phi > -2.0 || pSA_phi < 0.)PhiIntegral = pSA_phi + TMath::Pi()/4.;
    	if(pSA_phi > 0. || pSA_phi < 0.8)PhiIntegral = pSA_phi;
    	if(pSA_phi > 0.8 || pSA_phi < 1.55)PhiIntegral = pSA_phi - TMath::Pi()/4.;
    	if(pSA_phi > 1.55 || pSA_phi < 2.40)PhiIntegral = pSA_phi - TMath::Pi()/2.;
    	if(pSA_phi > 2.40)PhiIntegral = pSA_phi - 3*TMath::Pi()/4.;
    	m_h_SmallPhi->Fill(PhiIntegral);
    }
}

void PtNewMethod::Finalize(TFile *tf1,std::string filename){
	tf1->cd();
	m_h_BarrelAlpha->Write();
	m_h_PtvsBarrelAlpha->Write();
	m_h_BarrelBeta->Write();
	m_h_PtvsBarrelBeta->Write();
	m_h_DeltaThetaBI->Write();
	m_h_DeltaThetaBM->Write();
	m_h_LargePhi->Write();
	m_h_SmallPhi->Write();
	cout<<"finish!!"<<endl;
}
