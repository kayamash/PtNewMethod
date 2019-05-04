#include "../PtNewMethod/PtNewMethod.chh"
#include "kayamashForLUT.cpp"
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

const Int_t etaMax = 30;
const Int_t phiMax = 30;

bool PtNewMethod::CutAll(Int_t Tagpass,Int_t L1pass){
	if(L1pass > -1 && m_tag_proc == m_proc){
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

void PtNewMethod::Loop(Int_t ev,std::string name,Int_t proc){
	tChain->GetEntry(ev);
	m_method_name = name;
	m_proc = proc;
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
	Double_t pSA_phims = -99999;
	Double_t pSA_phibe = -99999;
	float pSA_roieta = -99999;
	float pSA_roiphi = -99999;
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
	if( !CutAll(pEFTAG_pass,pL1_pass) )return;
	if(!BarrelDicision(pSA_roieta))return;
	m_h_offPt->Fill(std::fabs(m_poff_pt*0.001));
	Double_t barrelalpha = -99999;
	Double_t barrelbeta = -99999;
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
	if(pSA_superpointR_BM != 0)barrelalpha = atan(pSA_superpointZ_BM/pSA_superpointR_BM) - atan(pSA_superpointSlope_BM);//Reciprocal number?
    //barrel alpha end

    //barrel beta
    if(pSA_superpointR_BI != 0 && pSA_superpointR_BM != 0)barrelbeta = atan(1.0/pSA_superpointSlope_BI) - atan(1.0/pSA_superpointSlope_BM);//Reciprocal number?
    //barrel beta end

    //Line BI
    if(pSA_superpointR_BI != 0){
    	Double_t deltaThetaLineBI = atan(pSA_superpointR_BI/pSA_superpointZ_BI) - atan(1.0/pSA_superpointSlope_BI);
    	m_h_DeltaThetaLineBI->Fill(deltaThetaLineBI);
    	m_h_PtvsDeltaThetaLineBI->Fill(std::fabs(m_poff_pt*0.001),deltaThetaLineBI);
    }

    Int_t LUTparameter[5];//sector charge eta phi sectorNumber
    for(Int_t i = 0; i < 5; ++i){
    	LUTparameter[i] = 0;
    }
    Double_t phiInteg = 0;
    kayamashForLUT LUT(0.,0.);
    bool LUTcheck = LUT.getLUTparameter(pSA_sAddress,m_poff_charge,m_poff_eta,m_poff_phi,LUTparameter,phiInteg);
    if(LUTcheck)switch(LUTparameter[0]){
    	case 0:
    	m_h_LargePhiEta[LUTparameter[2]]->Fill(m_poff_phi);
    	break;
    	case 1:
    	m_h_Sector11PhiEta[LUTparameter[2]]->Fill(m_poff_phi);
    	break;
    	case 2:
    	m_h_Sector15PhiEta[LUTparameter[2]]->Fill(m_poff_phi);
    	break;
    	case 3:
    	m_h_SmallPhiEta[LUTparameter[2]]->Fill(m_poff_phi);
    	break;
    	case 4:
    	m_h_SmallSpecialPhiEta[LUTparameter[2]]->Fill(m_poff_phi);
    	break;
    }
    if(LUTcheck && barrelalpha != -99999)m_h_PtvsBarrelAlpha_SectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]]->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    if(LUTcheck && barrelbeta != -99999)m_h_PtvsBarrelBeta_SectorChargeEtaPhi[LUTparameter[0]][LUTparameter[1]][LUTparameter[2]][LUTparameter[3]]->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    if(LUTcheck && barrelalpha != -99999){
    	m_h_BarrelAlpha->Fill(barrelalpha);
        m_h_PtvsBarrelAlpha->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    	switch(static_cast<Int_t>(pSA_sAddress)){
    		case 0:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelAlpha_LargePositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelAlpha_LargeNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		break;
    		case 1:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelAlpha_LSPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelAlpha_LSNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		break;
    		case 2:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelAlpha_SmallPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelAlpha_SmallNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		break;
    		case 3:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelAlpha_SSPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelAlpha_SSNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelalpha);
    		break;
    	}
    }
    if(LUTcheck && barrelbeta != -99999){
    	m_h_BarrelBeta->Fill(barrelbeta);
        m_h_PtvsBarrelBeta->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    	switch(static_cast<Int_t>(pSA_sAddress)){
    		case 0:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelBeta_LargePositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelBeta_LargeNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		break;
    		case 1:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelBeta_LSPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelBeta_LSNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		break;
    		case 2:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelBeta_SmallPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelBeta_SmallNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		break;
    		case 3:
    		if(m_poff_charge == 1.)m_h_PtvsBarrelBeta_SSPositive->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		if(m_poff_charge == -1.)m_h_PtvsBarrelBeta_SSNegative->Fill(1.0/std::fabs(m_poff_pt*0.001),barrelbeta);
    		break;
    	}
    }
}

void PtNewMethod::Finalize(TFile *tf1,std::string filenameA,std::string filenameB){
	ofstream ofs;
	ofs.open(filenameA.c_str());
	ofs.close();
	ofs.open(filenameB.c_str());
	ofs.close();
	kayamashForLUT LUT(16.,15.);

	tf1->cd();
	m_h_offPt->Write();
	m_h_BarrelAlpha->Write();
	m_h_PtvsBarrelAlpha->Write();
	m_h_BarrelBeta->Write();
	m_h_PtvsBarrelBeta->Write();
	m_h_DeltaThetaBI->Write();
	m_h_DeltaThetaBM->Write();
	m_h_DeltaThetaLineBI->Write();
	m_h_PtvsDeltaThetaLineBI->Write();
	//For LUT
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Large");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/LargeSpecial");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Small");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/SmallSpecial");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Large");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/LargeSpecial");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Small");
	tf1->mkdir("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/SmallSpecial");
	for(Int_t charge = 0; charge < 2;++charge){//positive,negative
		for(Int_t sector = 0; sector < 5;++sector){//Large,LargeSpecial sector11,LargeSpecial sector 15,Small,SmallSpecial
			if(charge == 0 && sector == 0)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Large");
			if(charge == 0 && (sector == 1 || sector == 2))tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/LargeSpecial");
			if(charge == 0 && sector == 3)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Small");
			if(charge == 0 && sector == 4)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/SmallSpecial");
			if(charge == 1 && (sector == 1 || sector == 2))tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Large");
			if(charge == 1 && sector == 0)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/LargeSpecial");
			if(charge == 1 && sector == 3)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Small");
			if(charge == 1 && sector == 4)tf1->cd("h_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/SmallSpecial");
        	for(Int_t eta = 0; eta < etaMax;++eta){
        		for(Int_t phi = 0; phi < phiMax;++phi){
        			m_h_PtvsBarrelAlpha_SectorChargeEtaPhi[sector][charge][eta][phi]->Write();
        			m_prof_PtvsBarrelAlpha_SectorChargeEtaPhi[sector][charge][eta][phi] = m_h_PtvsBarrelAlpha_SectorChargeEtaPhi[sector][charge][eta][phi]->ProfileX();
        			LUT.WriteLUT(m_prof_PtvsBarrelAlpha_SectorChargeEtaPhi[sector][charge][eta][phi],sector,charge,eta,phi,kTRUE,filenameA);
        			m_h_PtvsBarrelBeta_SectorChargeEtaPhi[sector][charge][eta][phi]->Write();
        			m_prof_PtvsBarrelBeta_SectorChargeEtaPhi[sector][charge][eta][phi] = m_h_PtvsBarrelBeta_SectorChargeEtaPhi[sector][charge][eta][phi]->ProfileX();
        			LUT.WriteLUT(m_prof_PtvsBarrelBeta_SectorChargeEtaPhi[sector][charge][eta][phi],sector,charge,eta,phi,kFALSE,filenameB);
        		}
        	}
        	tf1->cd();
    	}
    }
    tf1->cd();

    //Profile, only use for save profile to other directry,unnecessary?
    tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Large");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/LargeSpecial");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Small");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/SmallSpecial");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Large");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/LargeSpecial");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Small");
	tf1->mkdir("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/SmallSpecial");
	for(Int_t charge = 0; charge < 2;++charge){//positive,negative
		for(Int_t sector = 0; sector < 5;++sector){//Large,LargeSpecial sector11,LargeSpecial sector 15,Small,SmallSpecial
			if(charge == 0 && sector == 0)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Large");
			if(charge == 0 && (sector == 1 || sector == 2))tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/LargeSpecial");
			if(charge == 0 && sector == 3)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/Small");
			if(charge == 0 && sector == 4)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/positive/SmallSpecial");
			if(charge == 1 && sector == 0)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Large");
			if(charge == 1 && (sector == 1 || sector == 2))tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/LargeSpecial");
			if(charge == 1 && sector == 3)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/Small");
			if(charge == 1 && sector == 4)tf1->cd("Profile_PtvsBarrelAlpha_SectorChargeEtaPhi/negative/SmallSpecial");
        	for(Int_t eta = 0; eta < etaMax;++eta){
        		for(Int_t phi = 0; phi < phiMax;++phi){
        			m_prof_PtvsBarrelAlpha_SectorChargeEtaPhi[sector][charge][eta][phi]->Write();
        			m_prof_PtvsBarrelBeta_SectorChargeEtaPhi[sector][charge][eta][phi]->Write();
        		}
        	}
        	tf1->cd();
    	}
    }

    m_h_PtvsBarrelAlpha_LargePositive->Write();
    m_h_PtvsBarrelAlpha_LargeNegative->Write();
    m_h_PtvsBarrelAlpha_SmallPositive->Write();
    m_h_PtvsBarrelAlpha_SmallNegative->Write();
    m_h_PtvsBarrelAlpha_LSPositive->Write();
    m_h_PtvsBarrelAlpha_LSNegative->Write();
    m_h_PtvsBarrelAlpha_SSPositive->Write();
    m_h_PtvsBarrelAlpha_SSNegative->Write();
    m_h_PtvsBarrelBeta_LargePositive->Write();
    m_h_PtvsBarrelBeta_LargeNegative->Write();
    m_h_PtvsBarrelBeta_SmallPositive->Write();
    m_h_PtvsBarrelBeta_SmallNegative->Write();
    m_h_PtvsBarrelBeta_LSPositive->Write();
    m_h_PtvsBarrelBeta_LSNegative->Write();
    m_h_PtvsBarrelBeta_SSPositive->Write();
    m_h_PtvsBarrelBeta_SSNegative->Write();
    for(Int_t i = 0; i < 30;++i){
    	m_h_LargePhiEta[i]->Write();
    	m_h_Sector11PhiEta[i]->Write();
    	m_h_Sector15PhiEta[i]->Write();
    	m_h_SmallPhiEta[i]->Write();
    	m_h_SmallSpecialPhiEta[i]->Write();
    }
	cout<<"finish!!"<<endl;
}
