/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 28 22:40:34 2019 by ROOT version 5.34/38
// from TTree T1/WPrimeNtuple
// found on file: root://se01.indiacms.res.in//store/user/chatterj/TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_crab_L5JERC_TTBar_Mtt_700to1000_Autumn18_JECV19/191113_143123/0000/rootuple_jerc_l5_106.root
//////////////////////////////////////////////////////////

#ifndef Anal_Leptop_PROOF_h
#define Anal_Leptop_PROOF_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TMath.h>

#include "TLorentzVector.h"
#include <TProofOutputFile.h>
#include <TProofServ.h>

#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace TMVA;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

int getbinid(double val, int nbmx, float* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double eta_to_theta(double eta){
  return(2*atan(exp(-2*eta)));
}

double PhiInRange(const double& phi) {
  double phiout = phi;
  //added newly by Debarati//
  if(TMath::IsNaN(phiout)){
    gROOT->Error("PhiInRange","function called with NaN");
    return phiout;
  }
  //add ended//
  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout < -M_PI) phiout += 2*M_PI;
  else if (phiout >=  M_PI) phiout -= 2*M_PI;
  
  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

double diff_func(double f1, double f2){
  double ff = pow(f1-f2,2)*1./pow(f1+f2,2);
  return ff;
}


bool Muon_TightID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst,
				  float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
bool tightid = false;
if(muonisGL && muonisPF){
	if(muonchi<10 && muonhit>0 && muonmst>1){
		if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
			if(muonpixhit>0 && muontrklay>5){
				tightid = true;
				}
			} 
		}
	}	
return tightid;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon

bool Muon_Iso_ID(float muonpfiso)
{
bool isoid = false;	
if(muonpfiso<0.15) { isoid = true; } //SR
//if(muonpfiso>0.15) { isoid = true; }  // CR
return isoid;
}


double BTag_SF(int flavor, string syst, double pt){

// scale factors taken from the csv file in : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X (for medium WP)

double x = pt;
if(x>1000) { x = 1000; }

if(syst ==  "noSyst") {
	if(abs(flavor)==5||abs(flavor)==4){
        if (pt >= 20  && pt < 1000) return (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))));
     }else{
			if (pt >= 20  && pt < 1000) return (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x) ;
		  }
 }

if(syst ==  "down") {
	
	if(abs(flavor)==5){
		if(pt >= 20  && pt < 30) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.19459584355354309) );
		if(pt >= 30  && pt < 50) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.04693598672747612) );
		if(pt >= 50  && pt < 70) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.041476961225271225) );
		if(pt >= 70  && pt < 100) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.037213429808616638) );
		if(pt >= 100  && pt < 140) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.033781636506319046) );
		if(pt >= 140  && pt < 200) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.035268638283014297) );
		if(pt >= 200  && pt < 300) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.043516248464584351) );
		if(pt >= 300  && pt < 600) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.10369165241718292) );
		if(pt >= 600  && pt < 1000) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.29925653338432312) );
	}
	if(abs(flavor)==4){
		if(pt >= 20  && pt < 30) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.064865283668041229) );
		if(pt >= 30  && pt < 50) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.015645328909158707) );
		if(pt >= 50  && pt < 70) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.013825654052197933) );
		if(pt >= 70  && pt < 100) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.012404476292431355) );
		if(pt >= 100  && pt < 140) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.011260545812547207) );
		if(pt >= 140  && pt < 200) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.011756212450563908) );
		if(pt >= 200  && pt < 300) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.01450541615486145) );
		if(pt >= 300  && pt < 600) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.034563884139060974) );
		if(pt >= 600  && pt < 1000) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.099752180278301239) );
	}
	if(abs(flavor)!=5 && abs(flavor)!=4){
		if (pt >= 20  && pt < 1000) return ( (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x)*(1-(0.142253+0.000227323*x+-2.71704e-07*x*x)) );
	}
	
}

if(syst ==  "up") {
	
	if(abs(flavor)==5){
		if(pt >= 20  && pt < 30) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.19459584355354309 );
		if(pt >= 30  && pt < 50) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.04693598672747612 );
		if(pt >= 50  && pt < 70) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.041476961225271225 );
		if(pt >= 70  && pt < 100) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.037213429808616638 );
		if(pt >= 100  && pt < 140) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.033781636506319046 );
		if(pt >= 140  && pt < 200) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.035268638283014297 );
		if(pt >= 200  && pt < 300) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.043516248464584351 );
		if(pt >= 300  && pt < 600) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.10369165241718292 );
		if(pt >= 600  && pt < 1000) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.29925653338432312 );
	}
	if(abs(flavor)==4){
		if(pt >= 20  && pt < 30) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.064865283668041229 );
		if(pt >= 30  && pt < 50) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.015645328909158707 );
		if(pt >= 50  && pt < 70) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.013825654052197933 );
		if(pt >= 70  && pt < 100) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.012404476292431355 );
		if(pt >= 100  && pt < 140) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.011260545812547207 );
		if(pt >= 140  && pt < 200) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.011756212450563908 );
		if(pt >= 200  && pt < 300) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.01450541615486145 );
		if(pt >= 300  && pt < 600) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.034563884139060974 );
		if(pt >= 600  && pt < 1000) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.099752180278301239 );
	}
	if(abs(flavor)!=5 && abs(flavor)!=4){
		if (pt >= 20  && pt < 1000) return ( (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x)*(1+(0.115686+0.000270835*x+-3.2078e-07*x*x)) );
	}
	
}

return 1.0;

}

double BTag_MCEfficiency_TT(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.816517 ; 
		if(pt>=50 && pt<70) return 0.842834 ; 
		if(pt>=70 && pt<100) return 0.855092 ; 
		if(pt>=100 && pt<140) return 0.861484 ; 
		if(pt>=140 && pt<200) return 0.864907 ; 
		if(pt>=200 && pt<300) return 0.860499 ; 
		if(pt>=300 && pt<600) return 0.841221 ; 
		if(pt>=600 && pt<10000) return 0.794771 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.799168 ; 
		if(pt>=50 && pt<70) return 0.827512 ; 
		if(pt>=70 && pt<100) return 0.841101 ; 
		if(pt>=100 && pt<140) return 0.848418 ; 
		if(pt>=140 && pt<200) return 0.852231 ; 
		if(pt>=200 && pt<300) return 0.84615 ; 
		if(pt>=300 && pt<600) return 0.823063 ; 
		if(pt>=600 && pt<10000) return 0.765142 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.739489 ; 
		if(pt>=50 && pt<70) return 0.770812 ; 
		if(pt>=70 && pt<100) return 0.784899 ; 
		if(pt>=100 && pt<140) return 0.793595 ; 
		if(pt>=140 && pt<200) return 0.800014 ; 
		if(pt>=200 && pt<300) return 0.791883 ; 
		if(pt>=300 && pt<600) return 0.761082 ; 
		if(pt>=600 && pt<10000) return 0.693082 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.186714 ; 
		if(pt>=50 && pt<70) return 0.160875 ; 
		if(pt>=70 && pt<100) return 0.152947 ; 
		if(pt>=100 && pt<140) return 0.151804 ; 
		if(pt>=140 && pt<200) return 0.158477 ; 
		if(pt>=200 && pt<300) return 0.179461 ; 
		if(pt>=300 && pt<600) return 0.220294 ; 
		if(pt>=600 && pt<10000) return 0.238481 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.184478 ; 
		if(pt>=50 && pt<70) return 0.162409 ; 
		if(pt>=70 && pt<100) return 0.155749 ; 
		if(pt>=100 && pt<140) return 0.154941 ; 
		if(pt>=140 && pt<200) return 0.161089 ; 
		if(pt>=200 && pt<300) return 0.179591 ; 
		if(pt>=300 && pt<600) return 0.213557 ; 
		if(pt>=600 && pt<10000) return 0.225232 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.175049 ; 
		if(pt>=50 && pt<70) return 0.155982 ; 
		if(pt>=70 && pt<100) return 0.148945 ; 
		if(pt>=100 && pt<140) return 0.148611 ; 
		if(pt>=140 && pt<200) return 0.157628 ; 
		if(pt>=200 && pt<300) return 0.17258 ; 
		if(pt>=300 && pt<600) return 0.198044 ; 
		if(pt>=600 && pt<10000) return 0.213212 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.016934 ; 
		if(pt>=50 && pt<70) return 0.00924123 ; 
		if(pt>=70 && pt<100) return 0.00748564 ; 
		if(pt>=100 && pt<140) return 0.0069227 ; 
		if(pt>=140 && pt<200) return 0.00732266 ; 
		if(pt>=200 && pt<300) return 0.00898066 ; 
		if(pt>=300 && pt<600) return 0.0149634 ; 
		if(pt>=600 && pt<10000) return 0.0302075 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0180451 ; 
		if(pt>=50 && pt<70) return 0.0102217 ; 
		if(pt>=70 && pt<100) return 0.00832923 ; 
		if(pt>=100 && pt<140) return 0.00770286 ; 
		if(pt>=140 && pt<200) return 0.00802331 ; 
		if(pt>=200 && pt<300) return 0.00966405 ; 
		if(pt>=300 && pt<600) return 0.0159261 ; 
		if(pt>=600 && pt<10000) return 0.0340812 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0240714 ; 
		if(pt>=50 && pt<70) return 0.0142675 ; 
		if(pt>=70 && pt<100) return 0.0119286 ; 
		if(pt>=100 && pt<140) return 0.0116722 ; 
		if(pt>=140 && pt<200) return 0.0132621 ; 
		if(pt>=200 && pt<300) return 0.0177477 ; 
		if(pt>=300 && pt<600) return 0.0306986 ; 
		if(pt>=600 && pt<10000) return 0.0619629 ; 
	}
}


return 1.0;

}


double BTag_MCEfficiency_ST(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.808507 ; 
		if(pt>=50 && pt<70) return 0.842147 ; 
		if(pt>=70 && pt<100) return 0.857556 ; 
		if(pt>=100 && pt<140) return 0.863114 ; 
		if(pt>=140 && pt<200) return 0.86546 ; 
		if(pt>=200 && pt<300) return 0.867303 ; 
		if(pt>=300 && pt<600) return 0.852594 ; 
		if(pt>=600 && pt<10000) return 0.795692 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.790807 ; 
		if(pt>=50 && pt<70) return 0.826757 ; 
		if(pt>=70 && pt<100) return 0.843925 ; 
		if(pt>=100 && pt<140) return 0.850409 ; 
		if(pt>=140 && pt<200) return 0.854573 ; 
		if(pt>=200 && pt<300) return 0.854686 ; 
		if(pt>=300 && pt<600) return 0.834977 ; 
		if(pt>=600 && pt<10000) return 0.763769 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.725811 ; 
		if(pt>=50 && pt<70) return 0.765861 ; 
		if(pt>=70 && pt<100) return 0.786233 ; 
		if(pt>=100 && pt<140) return 0.79589 ; 
		if(pt>=140 && pt<200) return 0.798712 ; 
		if(pt>=200 && pt<300) return 0.794785 ; 
		if(pt>=300 && pt<600) return 0.772344 ; 
		if(pt>=600 && pt<10000) return 0.690087 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.17957 ; 
		if(pt>=50 && pt<70) return 0.154237 ; 
		if(pt>=70 && pt<100) return 0.147169 ; 
		if(pt>=100 && pt<140) return 0.147091 ; 
		if(pt>=140 && pt<200) return 0.152557 ; 
		if(pt>=200 && pt<300) return 0.167154 ; 
		if(pt>=300 && pt<600) return 0.217309 ; 
		if(pt>=600 && pt<10000) return 0.236696 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.177132 ; 
		if(pt>=50 && pt<70) return 0.156133 ; 
		if(pt>=70 && pt<100) return 0.149773 ; 
		if(pt>=100 && pt<140) return 0.150236 ; 
		if(pt>=140 && pt<200) return 0.154026 ; 
		if(pt>=200 && pt<300) return 0.173306 ; 
		if(pt>=300 && pt<600) return 0.210182 ; 
		if(pt>=600 && pt<10000) return 0.245206 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.168265 ; 
		if(pt>=50 && pt<70) return 0.150387 ; 
		if(pt>=70 && pt<100) return 0.143042 ; 
		if(pt>=100 && pt<140) return 0.142441 ; 
		if(pt>=140 && pt<200) return 0.150273 ; 
		if(pt>=200 && pt<300) return 0.159623 ; 
		if(pt>=300 && pt<600) return 0.186812 ; 
		if(pt>=600 && pt<10000) return 0.245013 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0189811 ; 
		if(pt>=50 && pt<70) return 0.00925074 ; 
		if(pt>=70 && pt<100) return 0.00743353 ; 
		if(pt>=100 && pt<140) return 0.00695717 ; 
		if(pt>=140 && pt<200) return 0.00722728 ; 
		if(pt>=200 && pt<300) return 0.00901535 ; 
		if(pt>=300 && pt<600) return 0.0169166 ; 
		if(pt>=600 && pt<10000) return 0.0329311 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0200871 ; 
		if(pt>=50 && pt<70) return 0.0102896 ; 
		if(pt>=70 && pt<100) return 0.00807831 ; 
		if(pt>=100 && pt<140) return 0.00764342 ; 
		if(pt>=140 && pt<200) return 0.00821108 ; 
		if(pt>=200 && pt<300) return 0.00976979 ; 
		if(pt>=300 && pt<600) return 0.0159442 ; 
		if(pt>=600 && pt<10000) return 0.0330585 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0260421 ; 
		if(pt>=50 && pt<70) return 0.014171 ; 
		if(pt>=70 && pt<100) return 0.011781 ; 
		if(pt>=100 && pt<140) return 0.0119679 ; 
		if(pt>=140 && pt<200) return 0.0143276 ; 
		if(pt>=200 && pt<300) return 0.0193145 ; 
		if(pt>=300 && pt<600) return 0.0306877 ; 
		if(pt>=600 && pt<10000) return 0.0618204 ; 
	}
}
	
	
return 1.0;	
}

double BTag_MCEfficiency_DIB(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.804116 ; 
		if(pt>=50 && pt<70) return 0.833425 ; 
		if(pt>=70 && pt<100) return 0.845986 ; 
		if(pt>=100 && pt<140) return 0.85206 ; 
		if(pt>=140 && pt<200) return 0.853713 ; 
		if(pt>=200 && pt<300) return 0.85249 ; 
		if(pt>=300 && pt<600) return 0.8764 ; 
		if(pt>=600 && pt<10000) return 0.915479 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.785403 ; 
		if(pt>=50 && pt<70) return 0.821003 ; 
		if(pt>=70 && pt<100) return 0.831992 ; 
		if(pt>=100 && pt<140) return 0.839383 ; 
		if(pt>=140 && pt<200) return 0.838345 ; 
		if(pt>=200 && pt<300) return 0.835367 ; 
		if(pt>=300 && pt<600) return 0.849568 ; 
		if(pt>=600 && pt<10000) return 0.901798 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.716813 ; 
		if(pt>=50 && pt<70) return 0.756991 ; 
		if(pt>=70 && pt<100) return 0.772715 ; 
		if(pt>=100 && pt<140) return 0.781882 ; 
		if(pt>=140 && pt<200) return 0.791533 ; 
		if(pt>=200 && pt<300) return 0.790272 ; 
		if(pt>=300 && pt<600) return 0.793086 ; 
		if(pt>=600 && pt<10000) return 0.893946 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.181339 ; 
		if(pt>=50 && pt<70) return 0.157787 ; 
		if(pt>=70 && pt<100) return 0.149497 ; 
		if(pt>=100 && pt<140) return 0.147983 ; 
		if(pt>=140 && pt<200) return 0.151851 ; 
		if(pt>=200 && pt<300) return 0.178893 ; 
		if(pt>=300 && pt<600) return 0.232613 ; 
		if(pt>=600 && pt<10000) return 0.250121 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.178676 ; 
		if(pt>=50 && pt<70) return 0.158658 ; 
		if(pt>=70 && pt<100) return 0.152222 ; 
		if(pt>=100 && pt<140) return 0.146962 ; 
		if(pt>=140 && pt<200) return 0.158739 ; 
		if(pt>=200 && pt<300) return 0.178602 ; 
		if(pt>=300 && pt<600) return 0.214291 ; 
		if(pt>=600 && pt<10000) return 0.250597 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.169733 ; 
		if(pt>=50 && pt<70) return 0.152602 ; 
		if(pt>=70 && pt<100) return 0.144413 ; 
		if(pt>=100 && pt<140) return 0.14245 ; 
		if(pt>=140 && pt<200) return 0.151063 ; 
		if(pt>=200 && pt<300) return 0.16838 ; 
		if(pt>=300 && pt<600) return 0.217661 ; 
		if(pt>=600 && pt<10000) return 0.277058 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.014781 ; 
		if(pt>=50 && pt<70) return 0.00782419 ; 
		if(pt>=70 && pt<100) return 0.00667191 ; 
		if(pt>=100 && pt<140) return 0.00635735 ; 
		if(pt>=140 && pt<200) return 0.00633714 ; 
		if(pt>=200 && pt<300) return 0.00891228 ; 
		if(pt>=300 && pt<600) return 0.015369 ; 
		if(pt>=600 && pt<10000) return 0.0337295 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0157996 ; 
		if(pt>=50 && pt<70) return 0.00873308 ; 
		if(pt>=70 && pt<100) return 0.00731909 ; 
		if(pt>=100 && pt<140) return 0.00710014 ; 
		if(pt>=140 && pt<200) return 0.00770687 ; 
		if(pt>=200 && pt<300) return 0.00908425 ; 
		if(pt>=300 && pt<600) return 0.0164606 ; 
		if(pt>=600 && pt<10000) return 0.0318129 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.022855 ; 
		if(pt>=50 && pt<70) return 0.0131073 ; 
		if(pt>=70 && pt<100) return 0.0110448 ; 
		if(pt>=100 && pt<140) return 0.0109769 ; 
		if(pt>=140 && pt<200) return 0.0126815 ; 
		if(pt>=200 && pt<300) return 0.0165187 ; 
		if(pt>=300 && pt<600) return 0.0303801 ; 
		if(pt>=600 && pt<10000) return 0.0459584 ; 
	}
}

	
return 1.0;	
}

double BTag_MCEfficiency_WJ(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.773179 ; 
		if(pt>=50 && pt<70) return 0.802888 ; 
		if(pt>=70 && pt<100) return 0.811475 ; 
		if(pt>=100 && pt<140) return 0.810867 ; 
		if(pt>=140 && pt<200) return 0.815031 ; 
		if(pt>=200 && pt<300) return 0.817471 ; 
		if(pt>=300 && pt<600) return 0.817811 ; 
		if(pt>=600 && pt<10000) return 0.800319 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.755851 ; 
		if(pt>=50 && pt<70) return 0.782044 ; 
		if(pt>=70 && pt<100) return 0.791109 ; 
		if(pt>=100 && pt<140) return 0.791923 ; 
		if(pt>=140 && pt<200) return 0.793282 ; 
		if(pt>=200 && pt<300) return 0.799978 ; 
		if(pt>=300 && pt<600) return 0.793994 ; 
		if(pt>=600 && pt<10000) return 0.775056 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.700553 ; 
		if(pt>=50 && pt<70) return 0.73214 ; 
		if(pt>=70 && pt<100) return 0.734361 ; 
		if(pt>=100 && pt<140) return 0.732253 ; 
		if(pt>=140 && pt<200) return 0.73254 ; 
		if(pt>=200 && pt<300) return 0.732826 ; 
		if(pt>=300 && pt<600) return 0.723183 ; 
		if(pt>=600 && pt<10000) return 0.701779 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.182172 ; 
		if(pt>=50 && pt<70) return 0.165362 ; 
		if(pt>=70 && pt<100) return 0.16155 ; 
		if(pt>=100 && pt<140) return 0.15973 ; 
		if(pt>=140 && pt<200) return 0.165449 ; 
		if(pt>=200 && pt<300) return 0.18322 ; 
		if(pt>=300 && pt<600) return 0.208985 ; 
		if(pt>=600 && pt<10000) return 0.232441 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.178307 ; 
		if(pt>=50 && pt<70) return 0.163694 ; 
		if(pt>=70 && pt<100) return 0.160277 ; 
		if(pt>=100 && pt<140) return 0.160435 ; 
		if(pt>=140 && pt<200) return 0.164053 ; 
		if(pt>=200 && pt<300) return 0.179322 ; 
		if(pt>=300 && pt<600) return 0.199861 ; 
		if(pt>=600 && pt<10000) return 0.216861 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.163093 ; 
		if(pt>=50 && pt<70) return 0.149782 ; 
		if(pt>=70 && pt<100) return 0.144074 ; 
		if(pt>=100 && pt<140) return 0.14598 ; 
		if(pt>=140 && pt<200) return 0.154604 ; 
		if(pt>=200 && pt<300) return 0.167698 ; 
		if(pt>=300 && pt<600) return 0.179704 ; 
		if(pt>=600 && pt<10000) return 0.202716 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0211166 ; 
		if(pt>=50 && pt<70) return 0.0117541 ; 
		if(pt>=70 && pt<100) return 0.00898972 ; 
		if(pt>=100 && pt<140) return 0.00772426 ; 
		if(pt>=140 && pt<200) return 0.0076329 ; 
		if(pt>=200 && pt<300) return 0.00867429 ; 
		if(pt>=300 && pt<600) return 0.0138477 ; 
		if(pt>=600 && pt<10000) return 0.0314047 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0231075 ; 
		if(pt>=50 && pt<70) return 0.0133186 ; 
		if(pt>=70 && pt<100) return 0.0101603 ; 
		if(pt>=100 && pt<140) return 0.0088122 ; 
		if(pt>=140 && pt<200) return 0.00866231 ; 
		if(pt>=200 && pt<300) return 0.00976423 ; 
		if(pt>=300 && pt<600) return 0.0145544 ; 
		if(pt>=600 && pt<10000) return 0.0338159 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0326849 ; 
		if(pt>=50 && pt<70) return 0.0182817 ; 
		if(pt>=70 && pt<100) return 0.0140935 ; 
		if(pt>=100 && pt<140) return 0.0135605 ; 
		if(pt>=140 && pt<200) return 0.0155933 ; 
		if(pt>=200 && pt<300) return 0.0189553 ; 
		if(pt>=300 && pt<600) return 0.028047 ; 
		if(pt>=600 && pt<10000) return 0.0582431 ; 
	}
}

return 1.0;

}

double BTag_MCEfficiency_DY(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.773204 ; 
		if(pt>=50 && pt<70) return 0.804343 ; 
		if(pt>=70 && pt<100) return 0.816013 ; 
		if(pt>=100 && pt<140) return 0.826033 ; 
		if(pt>=140 && pt<200) return 0.834742 ; 
		if(pt>=200 && pt<300) return 0.842605 ; 
		if(pt>=300 && pt<600) return 0.83394 ; 
		if(pt>=600 && pt<10000) return 0.80125 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.762663 ; 
		if(pt>=50 && pt<70) return 0.78978 ; 
		if(pt>=70 && pt<100) return 0.800517 ; 
		if(pt>=100 && pt<140) return 0.808844 ; 
		if(pt>=140 && pt<200) return 0.821501 ; 
		if(pt>=200 && pt<300) return 0.824157 ; 
		if(pt>=300 && pt<600) return 0.812596 ; 
		if(pt>=600 && pt<10000) return 0.76391 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.705352 ; 
		if(pt>=50 && pt<70) return 0.73545 ; 
		if(pt>=70 && pt<100) return 0.737149 ; 
		if(pt>=100 && pt<140) return 0.748329 ; 
		if(pt>=140 && pt<200) return 0.757617 ; 
		if(pt>=200 && pt<300) return 0.758913 ; 
		if(pt>=300 && pt<600) return 0.741301 ; 
		if(pt>=600 && pt<10000) return 0.69709 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.180718 ; 
		if(pt>=50 && pt<70) return 0.165395 ; 
		if(pt>=70 && pt<100) return 0.161978 ; 
		if(pt>=100 && pt<140) return 0.162888 ; 
		if(pt>=140 && pt<200) return 0.169717 ; 
		if(pt>=200 && pt<300) return 0.188192 ; 
		if(pt>=300 && pt<600) return 0.21083 ; 
		if(pt>=600 && pt<10000) return 0.227585 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.176288 ; 
		if(pt>=50 && pt<70) return 0.160393 ; 
		if(pt>=70 && pt<100) return 0.158328 ; 
		if(pt>=100 && pt<140) return 0.161196 ; 
		if(pt>=140 && pt<200) return 0.168502 ; 
		if(pt>=200 && pt<300) return 0.182771 ; 
		if(pt>=300 && pt<600) return 0.199411 ; 
		if(pt>=600 && pt<10000) return 0.216308 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.160878 ; 
		if(pt>=50 && pt<70) return 0.144542 ; 
		if(pt>=70 && pt<100) return 0.139217 ; 
		if(pt>=100 && pt<140) return 0.142293 ; 
		if(pt>=140 && pt<200) return 0.15166 ; 
		if(pt>=200 && pt<300) return 0.164528 ; 
		if(pt>=300 && pt<600) return 0.178783 ; 
		if(pt>=600 && pt<10000) return 0.197509 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0287349 ; 
		if(pt>=50 && pt<70) return 0.0156881 ; 
		if(pt>=70 && pt<100) return 0.0113781 ; 
		if(pt>=100 && pt<140) return 0.00957575 ; 
		if(pt>=140 && pt<200) return 0.00916245 ; 
		if(pt>=200 && pt<300) return 0.0103327 ; 
		if(pt>=300 && pt<600) return 0.0169311 ; 
		if(pt>=600 && pt<10000) return 0.0324204 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0312084 ; 
		if(pt>=50 && pt<70) return 0.0174152 ; 
		if(pt>=70 && pt<100) return 0.0127746 ; 
		if(pt>=100 && pt<140) return 0.0107694 ; 
		if(pt>=140 && pt<200) return 0.0103107 ; 
		if(pt>=200 && pt<300) return 0.0115971 ; 
		if(pt>=300 && pt<600) return 0.0176692 ; 
		if(pt>=600 && pt<10000) return 0.0339471 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.040997 ; 
		if(pt>=50 && pt<70) return 0.0222552 ; 
		if(pt>=70 && pt<100) return 0.0160718 ; 
		if(pt>=100 && pt<140) return 0.0147675 ; 
		if(pt>=140 && pt<200) return 0.0163437 ; 
		if(pt>=200 && pt<300) return 0.0200816 ; 
		if(pt>=300 && pt<600) return 0.0295968 ; 
		if(pt>=600 && pt<10000) return 0.0527284 ; 
	}
}

return 1.0;

}

double BTag_MCEfficiency_bQCD(int flavor, double pt, double eta){

if(flavor==5) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.811052 ; 
		if(pt>=50 && pt<70) return 0.835717 ; 
		if(pt>=70 && pt<100) return 0.849549 ; 
		if(pt>=100 && pt<140) return 0.856016 ; 
		if(pt>=140 && pt<200) return 0.857548 ; 
		if(pt>=200 && pt<300) return 0.865497 ; 
		if(pt>=300 && pt<600) return 0.875698 ; 
		if(pt>=600 && pt<10000) return 0.837915 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.793602 ; 
		if(pt>=50 && pt<70) return 0.822111 ; 
		if(pt>=70 && pt<100) return 0.837186 ; 
		if(pt>=100 && pt<140) return 0.843131 ; 
		if(pt>=140 && pt<200) return 0.847002 ; 
		if(pt>=200 && pt<300) return 0.853699 ; 
		if(pt>=300 && pt<600) return 0.859497 ; 
		if(pt>=600 && pt<10000) return 0.812472 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.738744 ; 
		if(pt>=50 && pt<70) return 0.773496 ; 
		if(pt>=70 && pt<100) return 0.780866 ; 
		if(pt>=100 && pt<140) return 0.788358 ; 
		if(pt>=140 && pt<200) return 0.793432 ; 
		if(pt>=200 && pt<300) return 0.795082 ; 
		if(pt>=300 && pt<600) return 0.79771 ; 
		if(pt>=600 && pt<10000) return 0.741586 ; 
	}
}

if(flavor==4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.198887 ; 
		if(pt>=50 && pt<70) return 0.182735 ; 
		if(pt>=70 && pt<100) return 0.181083 ; 
		if(pt>=100 && pt<140) return 0.183282 ; 
		if(pt>=140 && pt<200) return 0.197171 ; 
		if(pt>=200 && pt<300) return 0.206793 ; 
		if(pt>=300 && pt<600) return 0.228625 ; 
		if(pt>=600 && pt<10000) return 0.238966 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.189248 ; 
		if(pt>=50 && pt<70) return 0.177034 ; 
		if(pt>=70 && pt<100) return 0.176295 ; 
		if(pt>=100 && pt<140) return 0.179956 ; 
		if(pt>=140 && pt<200) return 0.185377 ; 
		if(pt>=200 && pt<300) return 0.196279 ; 
		if(pt>=300 && pt<600) return 0.209888 ; 
		if(pt>=600 && pt<10000) return 0.22342 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.170749 ; 
		if(pt>=50 && pt<70) return 0.158496 ; 
		if(pt>=70 && pt<100) return 0.154094 ; 
		if(pt>=100 && pt<140) return 0.152209 ; 
		if(pt>=140 && pt<200) return 0.165958 ; 
		if(pt>=200 && pt<300) return 0.174828 ; 
		if(pt>=300 && pt<600) return 0.183736 ; 
		if(pt>=600 && pt<10000) return 0.209696 ; 
	}
}

if(flavor!=5 && flavor!=4) {
	if(fabs(eta)>=0 && fabs(eta)<0.6) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0142848 ; 
		if(pt>=50 && pt<70) return 0.00827488 ; 
		if(pt>=70 && pt<100) return 0.0067742 ; 
		if(pt>=100 && pt<140) return 0.0064532 ; 
		if(pt>=140 && pt<200) return 0.00689059 ; 
		if(pt>=200 && pt<300) return 0.00809407 ; 
		if(pt>=300 && pt<600) return 0.0121235 ; 
		if(pt>=600 && pt<10000) return 0.0261152 ; 
	}
	if(fabs(eta)>=0.6 && fabs(eta)<1.2) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0152294 ; 
		if(pt>=50 && pt<70) return 0.00909953 ; 
		if(pt>=70 && pt<100) return 0.00769217 ; 
		if(pt>=100 && pt<140) return 0.00733744 ; 
		if(pt>=140 && pt<200) return 0.00738141 ; 
		if(pt>=200 && pt<300) return 0.00923737 ; 
		if(pt>=300 && pt<600) return 0.0130379 ; 
		if(pt>=600 && pt<10000) return 0.0325017 ; 
	}
	if(fabs(eta)>=1.2 && fabs(eta)<2.5) {
		if(pt>=20 && pt<30) return 0 ; 
		if(pt>=30 && pt<50) return 0.0217575 ; 
		if(pt>=50 && pt<70) return 0.0145133 ; 
		if(pt>=70 && pt<100) return 0.0123139 ; 
		if(pt>=100 && pt<140) return 0.0126898 ; 
		if(pt>=140 && pt<200) return 0.0145814 ; 
		if(pt>=200 && pt<300) return 0.0208179 ; 
		if(pt>=300 && pt<600) return 0.03009 ; 
		if(pt>=600 && pt<10000) return 0.0602379 ; 
	}
}

return 1.0;

}

double TopTag_SF(float pt){

// scale factors taken from DAKX twiki

 if(pt>=400 & pt<480) return 1.00;
 if(pt>=480 & pt<600) return 0.98;	
 if(pt>=600 & pt<12000) return 0.99;	
 
 return 1;
}

double TopTag_Efficiency_TT(float pt){

 if(pt>=400 && pt<480) return 0.64376 ; 
 if(pt>=480 && pt<600) return 0.688251 ; 
 if(pt>=600 && pt<12000) return 0.701325 ;

 return 1;
}

double TopTag_Efficiency_ST(float pt){

 if(pt>=400 && pt<480) return 0.597507 ; 
 if(pt>=480 && pt<600) return 0.611694 ; 
 if(pt>=600 && pt<12000) return 0.600842 ; 
 
 return 1;
}

double TopTag_Efficiency_DIB(float pt){

 if(pt>=400 && pt<480) return 0.225214 ; 
 if(pt>=480 && pt<600) return 0.24963 ;  
 if(pt>=600 && pt<12000) return 0.208174 ; 

 return 1;
}

double TopTag_Efficiency_WJ(float pt){

 if(pt>=400 && pt<480) return 0.101516 ; 
 if(pt>=480 && pt<600) return 0.0822782 ; 
 if(pt>=600 && pt<12000) return 0.0842843 ; 

 return 1;
}

double TopTag_Efficiency_DY(float pt){
 
 if(pt>=400 && pt<480) return 0.168126 ; 
 if(pt>=480 && pt<600) return 0.14405 ; 
 if(pt>=600 && pt<12000) return 0.141325 ;

 return 1;
}

double TopTag_Efficiency_bQCD(float pt){

 if(pt>=400 && pt<480) return 0.3192 ; 
 if(pt>=480 && pt<600) return 0.274223 ; 
 if(pt>=600 && pt<12000) return 0.265134 ;

 return 1;
}

double SF_TOP(double alpha, double beta, double pt0, double pt1)
{
        double sfwt = sqrt(exp(alpha-beta*pt0) * exp(alpha-beta*pt1));
        return sfwt;
}

void reOrder(std::vector<double>& pt, std::vector<double>& eta, std::vector<double>& phi, std::vector<float>& ch) { //std::vector<TLorentzVector>& jetcorr) {
  for (unsigned int i=0; i<pt.size(); i++) {
    for(unsigned int j=i+1; j<pt.size(); j++) {
      if (pt[i]<pt[j]) {
        swap<double>(pt[i], pt[j]);
        swap<double>(eta[i], eta[j]);
        swap<double>(phi[i], phi[j]);
	swap<float>(ch[i], ch[j]);
        //swap<math::XYZTLorentzVector>(jetcorr[i], jetcorr[j]);
      }
    }
  }
}

void reOrderCollection(std::vector<float>& ch, std::vector<TLorentzVector>& kin){
  for (unsigned int k=0; k<kin.size(); k++) {
    for (unsigned int l=k+1; l<kin.size(); l++) {
      if (kin[k].Pt()<kin[l].Pt()) {
	swap<float>(ch[k], ch[l]);
	swap<TLorentzVector>(kin[k], kin[l]); 
      }
    }
  }
}

template <class Type>
void swap(Type &a, Type &b){
  Type tmp = a;
  a=b;
  b=tmp;
}

class Anal_Leptop_PROOF : public TSelector {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  static const int njetmx = 20;
  static const int njetmxAK8 = 10;
  static const int npartmx = 25;
  
  float weight;
  float weight_puwup;
  float weight_puwdown;
  float weight_btagwup;
  float weight_btagwdown;
  float event_pt_weight;
   
  bool isMC;
  bool SemiLeptt;
  bool DiLeptt;
  bool Hadtt;
  bool EE, EMU, MUMU, EJets, MUJets, TAUJets, TauTau, ETau, MuTau;
  
  bool isTT;
  bool isST;
  bool isDIB;
  bool isWJ;
  bool isDY;
  bool isbQCD;
  
  // Declaration of leaf types
  Int_t           irun;
  Int_t           ilumi;
  UInt_t          ievt;
  Int_t           nprim;
  Int_t           nprimi;
  Double_t        Rho;
  Double_t        event_weight;
  Float_t         qscale;
  Int_t           npu_vert;
  Int_t           trig_value;
  
  Bool_t          hlt_IsoMu24;
  Bool_t          hlt_Mu50;
  Bool_t          hlt_Ele50_PFJet165;
  Bool_t          hlt_Ele115;
  Bool_t          hlt_AK8PFJet500;
  Bool_t          hlt_Photon200;
  Bool_t          hlt_Mu37Ele27;
  Bool_t          hlt_Mu27Ele37;
  Bool_t          hlt_Mu37TkMu27;
  Bool_t          hlt_OldMu100;
  Bool_t          hlt_TkMu100;
  Bool_t          hlt_DoubleEle25;

  Int_t           ntrigobjs;
  Float_t         trigobjpt[20];   //[ntrigobjs]
  Float_t         trigobjeta[20];   //[ntrigobjs]
  Float_t         trigobjphi[20];   //[ntrigobjs]
  Float_t         trigobjmass[20];   //[ntrigobjs]                                                              
  Int_t           trigobjpdgId[20];   //[ntrigobjs]
  Bool_t          trigobjHLT[20];   //[ntrigobjs]
  Bool_t          trigobjL1[20];   //[ntrigobjs]
  Bool_t          trigobjBoth[20];   //[ntrigobjs]
  Int_t           trigobjIhlt[20];   //[ntrigobjs]
  Float_t         PFMET;
  Float_t         PFMETPhi;
  Float_t         MisEtSig;
  Float_t         sumEt;
  Int_t           npfjetAK8;
  Bool_t          pfjetAK8jetID[njetmxAK8];
  Bool_t          pfjetAK8jetID_tightlepveto [njetmxAK8];
  Float_t         pfjetAK8pt[njetmxAK8];
  Float_t         pfjetAK8y[njetmxAK8];
  Float_t         pfjetAK8phi[njetmxAK8];
  Float_t         pfjetAK8mass[njetmxAK8];
  Float_t         pfjetAK8JEC[njetmxAK8];
  Float_t         pfjetAK8btag_DeepCSV[njetmxAK8];
  Float_t 	  pfjetAK8matchAK4deepb[njetmxAK8];
  Float_t         pfjetAK8DeepTag_TvsQCD[njetmxAK8];
  Float_t         pfjetAK8DeepTag_WvsQCD[njetmxAK8];
  Float_t         pfjetAK8DeepTag_ZvsQCD[njetmxAK8];
  Float_t         pfjetAK8CHF[njetmxAK8];
  Float_t         pfjetAK8NHF[njetmxAK8];
  Float_t         pfjetAK8CEMF[njetmxAK8];
  Float_t         pfjetAK8NEMF[njetmxAK8];
  Float_t         pfjetAK8MUF[njetmxAK8];
  Float_t         pfjetAK8PHF[njetmxAK8];
  Float_t         pfjetAK8EEF[njetmxAK8];
  Float_t         pfjetAK8HFHF[njetmxAK8];
  Float_t	  pfjetAK8neuemfrac[njetmxAK8];
  Float_t	  pfjetAK8neunhadfrac[njetmxAK8];
  Float_t	  pfjetAK8HadF[njetmxAK8];
  Float_t	  pfjetAK8NHadF[njetmxAK8];
  Float_t	  pfjetAK8EmF[njetmxAK8];
  Int_t           pfjetAK8CHM[njetmxAK8];
  Int_t           pfjetAK8NHM[njetmxAK8];
  Int_t           pfjetAK8MUM[njetmxAK8];
  Int_t           pfjetAK8PHM[njetmxAK8];
  Int_t           pfjetAK8EEM[njetmxAK8];
  Int_t           pfjetAK8HFHM[njetmxAK8];
  Int_t           pfjetAK8Neucons[njetmxAK8];
  Int_t           pfjetAK8Chcons[njetmxAK8];
  Int_t		  pfjetAK8ncons[njetmxAK8];
  Float_t         pfjetAK8reso[njetmxAK8];
  Float_t         pfjetAK8resoup[njetmxAK8];
  Float_t         pfjetAK8resodn[njetmxAK8];
  Float_t         pfjetAK8jesup_pu[njetmxAK8];
  Float_t         pfjetAK8jesup_rel[njetmxAK8];
  Float_t         pfjetAK8jesup_scale[njetmxAK8];
  Float_t         pfjetAK8jesup_total[njetmxAK8];
  Float_t         pfjetAK8jesdn_pu[njetmxAK8];
  Float_t         pfjetAK8jesdn_rel[njetmxAK8];
  Float_t         pfjetAK8jesdn_scale[njetmxAK8];
  Float_t         pfjetAK8jesdn_total[njetmxAK8];
  Float_t         pfjetAK8chrad[njetmxAK8];
  Float_t         pfjetAK8pTD[njetmxAK8];
  Float_t         pfjetAK8sdmass[njetmxAK8];
  bool		  pfjetAK8haspfelectron[njetmxAK8];
  bool		  pfjetAK8haspfmuon[njetmxAK8];
  Int_t		  pfjetAK8muon[njetmxAK8];
  Int_t           pfjetAK8ele[njetmxAK8];

  Float_t         pfjetAK8tau1[njetmxAK8];
  Float_t         pfjetAK8tau2[njetmxAK8];
  Float_t         pfjetAK8tau3[njetmxAK8];
  Float_t         pfjetAK8tau21[njetmxAK8];
  Float_t         pfjetAK8tau32[njetmxAK8];
  Float_t         pfjetAK8sub1pt[njetmxAK8];
  Float_t         pfjetAK8sub1eta[njetmxAK8];
  Float_t         pfjetAK8sub1phi[njetmxAK8];
  Float_t         pfjetAK8sub1mass[njetmxAK8];
  Float_t         pfjetAK8sub1btag[njetmxAK8];
  Float_t         pfjetAK8sub1chhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub1neuhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub1emfrac[njetmxAK8];
  Float_t         pfjetAK8sub1phofrac[njetmxAK8];
  Float_t         pfjetAK8sub1mufrac[njetmxAK8];
  Float_t	  pfjetAK8sub1hadfrac[njetmxAK8];
  Float_t         pfjetAK8sub2pt[njetmxAK8];
  Float_t         pfjetAK8sub2eta[njetmxAK8];
  Float_t         pfjetAK8sub2phi[njetmxAK8];
  Float_t         pfjetAK8sub2mass[njetmxAK8];
  Float_t         pfjetAK8sub2btag[njetmxAK8];
  Float_t         pfjetAK8sub2chhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub2neuhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub2emfrac[njetmxAK8];
  Float_t         pfjetAK8sub2phofrac[njetmxAK8];
  Float_t         pfjetAK8sub2mufrac[njetmxAK8];
  Float_t	  pfjetAK8sub2hadfrac[njetmxAK8];
  Float_t	  pfjetAK8subhaddiff[njetmxAK8];
  Float_t	  pfjetAK8subemdiff[njetmxAK8];
  
  Bool_t 	  pfjetAK8_hasmatche[njetmxAK8];
  Bool_t          pfjetAK8_hasmatchmu[njetmxAK8];
  Float_t	  pfjetAK8stau21[njetmxAK8];
  Float_t	  pfjetAK8stau32[njetmxAK8];
  Float_t	  pfjetAK8subptdiff[njetmxAK8];
  Float_t	  pfjetAK8subpldiff[njetmxAK8];
  Float_t	  pfjetAK8subbtag[njetmxAK8];
  Bool_t 	  pfjetAK8haselectron[njetmxAK8];
  Bool_t 		   pfjetAK8hasmuon[njetmxAK8];
  Bool_t                   pfjetAK8hastau[njetmxAK8];
  Bool_t		   pfjetAK8hasb[njetmxAK8];
  Bool_t		   pfjetAK8hasqg[njetmxAK8];
  Bool_t 		   pfjetAK8hashadtop[njetmxAK8];
  Bool_t 		   pfjetAK8hasleptop[njetmxAK8];
  Bool_t		   pfjetAK8hastop[njetmxAK8];
  Float_t		   pfjetAK8re_tvsb[njetmxAK8];
  Float_t		   pfjetAK8rnu_tvsb[njetmxAK8];
  Float_t		   pfjetAK8rt[njetmxAK8];
  Float_t                  pfjetAK8rmu_tvsb[njetmxAK8];

  Float_t                  pfjetAK8elinsubpt[njetmxAK8];
  Float_t                  pfjetAK8elinsubeta[njetmxAK8];
  Float_t                  pfjetAK8elinsubphi[njetmxAK8];
  
  Float_t                  pfjetAK8elinsubmass[njetmxAK8];
  Float_t                  pfjetAK8elinsubjpt[njetmxAK8];
  Float_t                  pfjetAK8elinsubjeta[njetmxAK8];
  Float_t                  pfjetAK8elinsubjphi[njetmxAK8];
  
  Float_t                  pfjetAK8elinsubjmass[njetmxAK8];

  Float_t                  pfjetAK8muinsubpt[njetmxAK8];
  Float_t                  pfjetAK8muinsubeta[njetmxAK8];
  Float_t                  pfjetAK8muinsubphi[njetmxAK8];
  Float_t                  pfjetAK8muinsubmass[njetmxAK8];
  Float_t                  pfjetAK8muinsubjpt[njetmxAK8];
  Float_t                  pfjetAK8muinsubjeta[njetmxAK8];
  Float_t                  pfjetAK8muinsubjphi[njetmxAK8];
  Float_t                  pfjetAK8muinsubje[njetmxAK8];
  Float_t                  pfjetAK8muinsubjmass[njetmxAK8];
  Float_t                  pfjetAK8muinsubI0[njetmxAK8];
  Float_t                  pfjetAK8muinsubInear[njetmxAK8];
  Float_t                  pfjetAK8muinsubIfar[njetmxAK8];

  Int_t           npfjetAK4;
  Bool_t          pfjetAK4jetID[njetmx];
  Bool_t          pfjetAK4jetID_tightlepveto[njetmx];
  Float_t         pfjetAK4pt[njetmx];
  Float_t         pfjetAK4eta[njetmx];
  Float_t         pfjetAK4y[njetmx];
  Float_t         pfjetAK4phi[njetmx];
  Float_t         pfjetAK4mass[njetmx];
  Float_t         pfjetAK4JEC[njetmx];
  
  Float_t         pfjetAK4btag_DeepCSV[njetmx];
  Float_t         pfjetAK4btag_DeepFlav[njetmx];
  
  Float_t         pfjetAK4reso[njetmx];
  Float_t         pfjetAK4resoup[njetmx];
  Float_t         pfjetAK4resodn[njetmx];
  Float_t         pfjetAK4jesup_pu[njetmx];
  Float_t         pfjetAK4jesup_rel[njetmx];
  Float_t         pfjetAK4jesup_scale[njetmx];
  Float_t         pfjetAK4jesup_total[njetmx];
  Float_t         pfjetAK4jesdn_pu[njetmx];
  Float_t         pfjetAK4jesdn_rel[njetmx];
  Float_t         pfjetAK4jesdn_scale[njetmx];
  Float_t         pfjetAK4jesdn_total[njetmx];
  Int_t           pfjetAK4hadronflav[njetmx];
  Int_t           pfjetAK4partonflav[njetmx];
  Float_t         pfjetAK4qgl[njetmx];
   Float_t         pfjetAK4PUID[njetmx];
   Int_t           pfjetAK4GenMatch;
   Float_t         GENMET;
   Float_t         GENMETPhi;
   
   Int_t           ngenjetAK8;
   Float_t         genjetAK8pt[njetmx];
   
   Float_t         genjetAK8phi[njetmx];
   Float_t         genjetAK8mass[njetmx];
   Float_t         genjetAK8sdmass[njetmx];
   
   
   
   Int_t           ngenjetAK4;
   Float_t         genjetAK4pt[njetmx];
   
   Float_t         genjetAK4phi[njetmx];
   Float_t         genjetAK4mass[njetmx];
   
   Int_t           ngenparticles;
   Int_t           genpartstatus[npartmx];
   Int_t           genpartpdg[npartmx];
   Int_t           genpartmompdg[npartmx];
   Int_t           genpartdaugno[npartmx];
   Bool_t          genpartfromhard[npartmx];
   Bool_t          genpartfromhardbFSR[npartmx];
   Float_t         genpartpt[npartmx];
   Float_t         genparteta[npartmx];
   Float_t         genpartphi[npartmx];
   Float_t         genpartm[npartmx];
   
   Int_t           genpartpair[npartmx];
   Int_t           ngenelectrons;
   Float_t         genelectronpt[npartmx];
   Float_t         genelectroneta[npartmx];
   Float_t         genelectronphi[npartmx];
   Float_t         genelectronm[npartmx];
   Int_t	   nleptop;
   Int_t	   nhadtop;
   Int_t           ngenmuons;
   Float_t         genmuonpt[npartmx];
   Float_t         genmuoneta[npartmx];
   Float_t         genmuonphi[npartmx];
   Float_t         genmuonm[npartmx];
   Int_t           ngentaus;
   Float_t         gentaupt[npartmx];
   Float_t         gentaueta[npartmx];
   Float_t         gentauphi[npartmx];
   Float_t         gentaum[npartmx];
   Int_t           ngenqgs;
   Float_t         genqgpt[npartmx];
   Float_t         genqgeta[npartmx];
   Float_t         genqgphi[npartmx];
   Float_t         genqgm[npartmx];
   Int_t           ngenbs;
   Float_t         genbpt[npartmx];
   Float_t         genbeta[npartmx];
   Float_t         genbphi[npartmx];
   Float_t         genbm[npartmx];
   Int_t	   ngentops;
   Float_t	   gentoppt[npartmx];
   Float_t	   gentopeta[npartmx];
   Float_t	   gentopphi[npartmx];
   Float_t	   gentopm[npartmx];
   Int_t	   gentopid[npartmx];
   Int_t           nmuons;
   Bool_t          muonisPF[njetmx];
   Bool_t          muonisGL[njetmx];
   Bool_t          muonisTRK[njetmx];
   Bool_t          muonisLoose[njetmx];
   Bool_t          muonisGoodGL[njetmx];
   Bool_t          muonisMed[njetmx];
   Bool_t          muonisTight[njetmx];
   Bool_t          muonisHighPt[njetmx];
   Bool_t          muonisHighPttrk[njetmx];
   Bool_t          muonisMedPr[njetmx];
   Float_t         muonpt[njetmx];
   Float_t         muonp[njetmx];
   
   Float_t         muoneta[njetmx];
   Float_t         muonphi[njetmx];
   Float_t         muoncharge[njetmx];
   
   Float_t         muontrkvtx[njetmx];
   Float_t         muondz[njetmx];
   Float_t         muonpter[njetmx];
   Float_t         muonchi[njetmx];
   Int_t           muonndf[njetmx];
   Float_t         muonecal[njetmx];
   Float_t         muonhcal[njetmx];
   
   Float_t         muonpfiso[njetmx];
   
   Float_t         muonposmatch[njetmx];
   Float_t         muontrkink[njetmx];
   Float_t         muonsegcom[njetmx];

   Float_t         muonpixhit[njetmx];
   Float_t         muonmst[njetmx];
   Float_t         muontrklay[njetmx];
   Float_t         muonvalfrac[njetmx];
   Float_t         muonchiso[njetmx];
   Float_t         muonnhiso[njetmx];
   Float_t         muonphiso[njetmx];
   Float_t         muonminisoall[njetmx];
   Float_t         mudxy_sv[njetmx];
   Float_t         muonhit[njetmx];

   Int_t           nelecs;
   Float_t         elpt[njetmx];
   Float_t         eldxy_sv[njetmx];
   Float_t         eleta[njetmx];
   Float_t         elphi[njetmx];
   Float_t         elp[njetmx];
   Float_t         ele[njetmx];
   Float_t         elcharge[njetmx];
   Bool_t          elmvaid[njetmx];
   Bool_t          elmvaid_noIso[njetmx];

   Float_t         elhovere[njetmx];
   Float_t         elchi[njetmx];
   Int_t           elndf[njetmx];
   
   Float_t         eletain[njetmx];
   Float_t         elphiin[njetmx];

   Float_t         elsupcl_eta[njetmx];
   Float_t         elsupcl_phi[njetmx];
   Float_t         elsupcl_rawE[njetmx];
   Float_t         elfbrem[njetmx];
   
   Float_t         eleoverp[njetmx];
   Float_t         elietaieta[njetmx];
   Float_t         elmisshits[njetmx];
   
   Float_t         elpfiso[njetmx];
   
   
   //***al el ID variables****//
   Float_t         elsigmaieta[njetmx];
   Float_t         elsigmaiphi[njetmx];
   Float_t         elr9full[njetmx];
   Float_t         elsupcl_etaw[njetmx];
   Float_t         elsupcl_phiw[njetmx];
   Float_t         elhcaloverecal[njetmx];
   Float_t         elcloctftrkn[njetmx];
   Float_t         elcloctftrkchi2[njetmx];
   Float_t         ele1x5bye5x5[njetmx];
   Float_t         elnormchi2[njetmx];
   Float_t         elhitsmiss[njetmx];
   Float_t         eltrkmeasure[njetmx];
   Float_t         elconvtxprob[njetmx];
   Float_t         elecloverpout[njetmx];
   Float_t         elecaletrkmomentum[njetmx];
   Float_t         eldeltaetacltrkcalo[njetmx];
   Float_t         elsupcl_preshvsrawe[njetmx];
   Float_t         elpfisolsumphet[njetmx];
   Float_t         elpfisolsumchhadpt[njetmx];
   Float_t         elpfsiolsumneuhadet[njetmx];

   /*
     Int_t           nphotons;
     Float_t         phoe[njetmx];
     Float_t         phoeta[njetmx];
     Float_t         phophi[njetmx];
     Bool_t          phomvaid[njetmx];
     Float_t         phoe1by9[njetmx];
     Float_t         phoe9by25[njetmx];
     Float_t         photrkiso[njetmx];
     Float_t         phoemiso[njetmx];
     Float_t         phohadiso[njetmx];
     Float_t         phochhadiso[njetmx];
     Float_t         phoneuhadiso[njetmx];
     Float_t         phophoiso[njetmx];
     Float_t         phoPUiso[njetmx];
     Float_t         phohadbyem[njetmx];
     Float_t         phoietaieta[njetmx];
   */
   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ilumi;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_nprim;   //!
   TBranch        *b_nprimi;
   TBranch        *b_Rho;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_qscale;   //!
   TBranch        *b_npu_vert;   //!
   TBranch        *b_trig_value;   //!

   TBranch        *b_hlt_IsoMu24;
   TBranch        *b_hlt_Mu50;
   TBranch        *b_hlt_Ele50_PFJet165;
   TBranch        *b_hlt_Ele115;
   TBranch        *b_hlt_AK8PFJet500;
   TBranch        *b_hlt_Photon200;
   TBranch        *b_hlt_Mu37Ele27;
   TBranch        *b_hlt_Mu27Ele37;
   TBranch        *b_hlt_Mu37TkMu27;
   TBranch        *b_hlt_OldMu100;
   TBranch        *b_hlt_TkMu100;
   TBranch        *b_hlt_DoubleEle25;

   TBranch        *b_ntrigobjs;   //!
   TBranch        *b_trigobjpt;   //!
   TBranch        *b_trigobjeta;   //!
   TBranch        *b_trigobjphi;   //!
   TBranch        *b_trigobjmass;   //!                                                                       
   TBranch        *b_trigobjpdgId;   //!
   TBranch        *b_trigobjHLT;   //!
   TBranch        *b_trigobjL1;   //!
   TBranch        *b_trigobjBoth;   //!
   TBranch        *b_trigobjIhlt;   //!
   TBranch        *b_miset;   //!
   TBranch        *b_misphi;   //!
   TBranch        *b_misetsig;   //!
   TBranch        *b_sumEt;   //!
   TBranch        *b_npfjetAK8;   //!
   TBranch        *b_pfjetAK8jetID;   //!                                                    
   TBranch        *b_pfjetAK8jetID_tightlepveto;   //!
   TBranch        *b_pfjetAK8pt;   //!
   TBranch        *b_pfjetAK8y;   //!
   TBranch        *b_pfjetAK8phi;   //!
   TBranch        *b_pfjetAK8mass;   //!
   TBranch        *b_pfjetAK8JEC;   //!
   TBranch        *b_pfjetAK8btag_DeepCSV;   //!
   TBranch        *b_pfjetAK8DeepTag_TvsQCD;   //!
   TBranch        *b_pfjetAK8DeepTag_WvsQCD;   //!
   TBranch        *b_pfjetAK8DeepTag_ZvsQCD;   //!
   TBranch        *b_pfjetAK8CHF;   //!
   TBranch        *b_pfjetAK8NHF;   //!
   TBranch        *b_pfjetAK8CEMF;   //!
   TBranch        *b_pfjetAK8NEMF;   //!
   TBranch        *b_pfjetAK8MUF;   //!
   TBranch        *b_pfjetAK8PHF;   //!
   TBranch        *b_pfjetAK8EEF;   //!
   TBranch        *b_pfjetAK8HFHF;   //!
   TBranch        *b_pfjetAK8CHM;   //!
   TBranch        *b_pfjetAK8NHM;   //!
   TBranch        *b_pfjetAK8MUM;   //!
   TBranch        *b_pfjetAK8PHM;   //!
   TBranch        *b_pfjetAK8EEM;   //!
   TBranch        *b_pfjetAK8HFHM;   //!
   TBranch        *b_pfjetAK8Neucons;   //!
   TBranch        *b_pfjetAK8Chcons;   //!
   TBranch        *b_pfjetAK8reso;   //!
   TBranch        *b_pfjetAK8resoup;   //!
   TBranch        *b_pfjetAK8resodn;   //!
   TBranch        *b_pfjetAK8jesup_pu;   //!
   TBranch        *b_pfjetAK8jesup_rel;   //!
   TBranch        *b_pfjetAK8jesup_scale;   //!
   TBranch        *b_pfjetAK8jesup_total;   //!
   TBranch        *b_pfjetAK8jesdn_pu;   //!
   TBranch        *b_pfjetAK8jesdn_rel;   //!
   TBranch        *b_pfjetAK8jesdn_scale;   //!
   TBranch        *b_pfjetAK8jesdn_total;   //!
   TBranch        *b_pfjetAK8chrad;   //!
   TBranch        *b_pfjetAK8pTD;   //!
   TBranch        *b_pfjetAK8sdmass;   //!
   TBranch        *b_pfjetAK8tau1;   //!
   TBranch        *b_pfjetAK8tau2;   //!
   TBranch        *b_pfjetAK8tau3;   //!
   TBranch        *b_pfjetAK8sub1pt;   //!
   TBranch        *b_pfjetAK8sub1eta;   //!
   TBranch        *b_pfjetAK8sub1phi;   //!
   TBranch        *b_pfjetAK8sub1mass;   //!
   TBranch        *b_pfjetAK8sub1btag;   //!
   TBranch        *b_pfjetAK8sub1chhadfrac;   //!
   TBranch        *b_pfjetAK8sub1neuhadfrac;   //!
   TBranch        *b_pfjetAK8sub1emfrac;   //!
   TBranch        *b_pfjetAK8sub1phofrac;   //!
   TBranch        *b_pfjetAK8sub1mufrac;   //!
   TBranch        *b_pfjetAK8sub2pt;   //!
   TBranch        *b_pfjetAK8sub2eta;   //!
   TBranch        *b_pfjetAK8sub2phi;   //!
   TBranch        *b_pfjetAK8sub2mass;   //!
   TBranch        *b_pfjetAK8sub2btag;   //!
   TBranch        *b_pfjetAK8sub2chhadfrac;   //!
   TBranch        *b_pfjetAK8sub2neuhadfrac;   //!
   TBranch        *b_pfjetAK8sub2emfrac;   //!
   TBranch        *b_pfjetAK8sub2phofrac;   //!
   TBranch        *b_pfjetAK8sub2mufrac;   //!
   
   TBranch        *b_pfjetAK8elinsubpt;
   TBranch        *b_pfjetAK8elinsubeta;
   TBranch        *b_pfjetAK8elinsubphi;
   
   TBranch        *b_pfjetAK8elinsubmass;
   TBranch        *b_pfjetAK8elinsubjpt;
   TBranch        *b_pfjetAK8elinsubjeta;
   TBranch        *b_pfjetAK8elinsubjphi;
   
   TBranch        *b_pfjetAK8elinsubjmass;

   TBranch        *b_pfjetAK8muinsubpt;
   TBranch        *b_pfjetAK8muinsubeta;
   TBranch        *b_pfjetAK8muinsubphi;
   TBranch        *b_pfjetAK8muinsubmass;
   TBranch        *b_pfjetAK8muinsubjpt;
   TBranch        *b_pfjetAK8muinsubjeta;
   TBranch        *b_pfjetAK8muinsubjphi;
   TBranch        *b_pfjetAK8muinsubje;
   TBranch        *b_pfjetAK8muinsubjmass;
   TBranch        *b_pfjetAK8muinsubInear;
   TBranch        *b_pfjetAK8muinsubIfar;
   TBranch        *b_pfjetAK8muinsubI0;

   TBranch        *b_npfjetAK4;   //!
   TBranch        *b_pfjetAK4jetID;   //!                                                           
   TBranch        *b_pfjetAK4jetID_tightlepveto;   //!

   TBranch        *b_pfjetAK4pt;   //!
   TBranch        *b_pfjetAK4eta;   //!
   TBranch        *b_pfjetAK4y;   //!
   TBranch        *b_pfjetAK4phi;   //!
   TBranch        *b_pfjetAK4mass;   //!
   TBranch        *b_pfjetAK4JEC;   //!
   
   TBranch        *b_pfjetAK4btag_DeepCSV;   //!
   TBranch        *b_pfjetAK4btag_DeepFlav;   //!
   
   TBranch        *b_pfjetAK4reso;   //!
   TBranch        *b_pfjetAK4resoup;   //!
   TBranch        *b_pfjetAK4resodn;   //!
   TBranch        *b_pfjetAK4jesup_pu;   //!
   TBranch        *b_pfjetAK4jesup_rel;   //!
   TBranch        *b_pfjetAK4jesup_scale;   //!
   TBranch        *b_pfjetAK4jesup_total;   //!
   TBranch        *b_pfjetAK4jesdn_pu;   //!
   TBranch        *b_pfjetAK4jesdn_rel;   //!
   TBranch        *b_pfjetAK4jesdn_scale;   //!
   TBranch        *b_pfjetAK4jesdn_total;   //!
   TBranch        *b_pfjetAK4hadronflav;   //!
   TBranch        *b_pfjetAK4partonflav;   //!
   TBranch        *b_pfjetAK4qgl;   //!
   TBranch        *b_pfjetAK4PUID;   //!
   TBranch        *b_pfjetAK4GenMatch;   //!
   TBranch        *b_genmiset;   //!
   TBranch        *b_genmisphi;   //!
   TBranch        *b_genmisetsig;   //!
   TBranch        *b_ngenjetAK8;   //!
   TBranch        *b_genjetAK8pt;   //!
   
   TBranch        *b_genjetAK8phi;   //!
   TBranch        *b_genjetAK8mass;   //!
   TBranch        *b_genjetAK8sdmass;   //!
   
   TBranch        *b_ngenjetAK4;   //!
   TBranch        *b_genjetAK4pt;   //!
   TBranch        *b_genjetAK4phi;   //!
   TBranch        *b_genjetAK4mass;   //!
   
   TBranch        *b_ngenparticles;   //!
   TBranch        *b_genpartstatus;   //!
   TBranch        *b_genpartpdg;   //!
   TBranch        *b_genpartmompdg;   //!
   TBranch        *b_genpartdaugno;   //!
   TBranch        *b_genpartfromhard;   //!
   TBranch        *b_genpartfromhardbFSR;   //!
   TBranch        *b_genpartpt;   //!
   TBranch        *b_genparteta;   //!
   TBranch        *b_genpartphi;   //!
   TBranch        *b_genpartm;   //!
   
   TBranch        *b_nmuons;   //!
   TBranch        *b_muonisPF;   //!
   TBranch        *b_muonisGL;   //!
   TBranch        *b_muonisTRK;   //!
   TBranch        *b_muonisLoose;   //!
   TBranch        *b_muonisGoodGL;   //!
   TBranch        *b_muonisMed;   //!
   TBranch        *b_muonisTight;
   TBranch        *b_muonisHighPt;
   TBranch        *b_muonisHighPttrk;
   TBranch        *b_muonisMedPr;
   TBranch        *b_muonpt;   //!
   TBranch        *b_muonp;   //!
   
   TBranch        *b_muoneta;   //!
   TBranch        *b_muonphi;   //!
   TBranch        *b_muoncharge;
   
   TBranch        *b_muontrkvtx;   //!
   TBranch        *b_muondz;   //!
   TBranch        *b_muonpter;   //!
   TBranch        *b_muonchi;   //!
   TBranch        *b_muonndf;   //!
   TBranch        *b_muonecal;   //!
   TBranch        *b_muonhcal;   //!
   
   TBranch        *b_muonpfiso;   //!
   
   TBranch        *b_muonposmatch;   //!
   TBranch        *b_muontrkink;   //!
   TBranch        *b_muonsegcom;   //!

   TBranch        *b_muonpixhit;   //!
   TBranch        *b_muonmst;   //!
   TBranch        *b_muontrklay;   //!
   TBranch        *b_muonvalfrac;   //!
   TBranch        *b_muonchiso;
   TBranch        *b_muonnhiso;
   TBranch        *b_muonphiso;
   TBranch        *b_muonminisoall;
   TBranch        *b_mudxy_sv;
   TBranch        *b_muonhit;


   TBranch        *b_nelecs;   //!
   TBranch        *b_elpt;   //!
   TBranch        *b_eldxy_sv;
   TBranch        *b_eleta;   //!
   TBranch        *b_elphi;   //!
   TBranch        *b_elp;   //!
   TBranch        *b_ele;   //!
   TBranch        *b_elcharge;
   TBranch        *b_elmvaid;   //!
   TBranch        *b_elmvaid_noIso;   //!
   
   TBranch        *b_elhovere;   //!
   TBranch        *b_elchi;   //!
   TBranch        *b_elndf;   //!

   TBranch        *b_eletain;   //!
   TBranch        *b_elphiin;   //!

   TBranch        *b_elsupcl_eta;
   TBranch        *b_elsupcl_phi;
   TBranch        *b_elsupcl_rawE;
   TBranch        *b_elfbrem;   //!
   
   TBranch        *b_eleoverp;   //!
   TBranch        *b_elietaieta;   //!
   TBranch        *b_elmisshits;   //!
   TBranch        *b_elpfiso;   //!
   
   TBranch        *b_elsigmaieta;
   TBranch        *b_elsigmaiphi;
   TBranch        *b_elr9full;
   TBranch        *b_elsupcl_etaw;
   TBranch        *b_elsupcl_phiw;
   TBranch        *b_elhcaloverecal;
   TBranch        *b_elcloctftrkn;
   TBranch        *b_elcloctftrkchi2;
   TBranch        *b_ele1x5bye5x5;
   TBranch        *b_elnormchi2;
   TBranch        *b_elhitsmiss;
   TBranch        *b_eltrkmeasure;
   TBranch        *b_elconvtxprob;
   TBranch        *b_elecloverpout;
   TBranch        *b_elecaletrkmomentum;
   TBranch        *b_eldeltaetacltrkcalo;
   TBranch        *b_elsupcl_preshvsrawe;
   TBranch        *b_elpfisolsumphet;
   TBranch        *b_elpfisolsumchhadpt;
   TBranch        *b_elpfsiolsumneuhadet;

   /*
     TBranch        *b_nphotons;   //!
     TBranch        *b_phoe;   //!
     TBranch        *b_phoeta;   //!
     TBranch        *b_phophi;   //!
     TBranch        *b_phomvaid;   //!
     TBranch        *b_phoe1by9;   //!
     TBranch        *b_phoe9by25;   //!
     TBranch        *b_photrkiso;   //!
     TBranch        *b_phoemiso;   //!
     TBranch        *b_phohadiso;   //!
     TBranch        *b_phochhadiso;   //!
     TBranch        *b_phoneuhadiso;   //!
     TBranch        *b_phophoiso;   //!
     TBranch        *b_phoPUiso;   //!
     TBranch        *b_phohadbyem;   //!
     TBranch        *b_phoietaieta;   //!
   */
   float pfjetAK8pt_resoup[njetmxAK8], pfjetAK8mass_resoup[njetmxAK8], pfjetAK8pt_resodown[njetmxAK8], pfjetAK8mass_resodown[njetmxAK8]; 

 float selpfjetAK8pt[njetmxAK8], selpfjetAK8y[njetmxAK8], selpfjetAK8phi[njetmxAK8], selpfjetAK8mass[njetmxAK8];
 float selpfjetAK8matchAK4deepb[njetmxAK8];
 float selpfjetAK8DeepTag_TvsQCD[njetmxAK8], selpfjetAK8DeepTag_WvsQCD[njetmxAK8], selpfjetAK8DeepTag_ZvsQCD[njetmxAK8];
 float selpfjetAK8CHF[njetmxAK8], selpfjetAK8NHF[njetmxAK8], selpfjetAK8CEMF[njetmxAK8], selpfjetAK8NEMF[njetmxAK8], selpfjetAK8MUF[njetmxAK8], selpfjetAK8HOF[njetmxAK8], selpfjetAK8HadF[njetmxAK8], selpfjetAK8NHadF[njetmxAK8], selpfjetAK8EmF[njetmxAK8], selpfjetAK8neuemfrac[njetmxAK8], selpfjetAK8neunhadfrac[njetmxAK8];
 int selpfjetAK8CHM[njetmxAK8], selpfjetAK8NHM[njetmxAK8], selpfjetAK8EEM[njetmxAK8], selpfjetAK8MUM[njetmxAK8], selpfjetAK8Neucons[njetmxAK8], selpfjetAK8Chcons[njetmxAK8];
 float selpfjetAK8chrad[njetmxAK8], selpfjetAK8pTD[njetmxAK8], selpfjetAK8sdmass[njetmxAK8], selpfjetAK8tau21[njetmxAK8], selpfjetAK8tau32[njetmxAK8];

 float selpfjetAK8elinsubpt[njetmxAK8], selpfjetAK8elinsubeta[njetmxAK8], selpfjetAK8elinsubphi[njetmxAK8];
 float selpfjetAK8elinsubjpt[njetmxAK8], selpfjetAK8elinsubjeta[njetmxAK8], selpfjetAK8elinsubjphi[njetmxAK8];
 float selpfjetAK8matchedelID[njetmxAK8];
 float selpfjetAK8matchedelpt[njetmxAK8], selpfjetAK8matchedeleta[njetmxAK8], selpfjetAK8matchedelphi[njetmxAK8], selpfjetAK8matchedelE[njetmxAK8], selpfjetAK8matchedeldx[njetmxAK8], selpfjetAK8matchedeldxy_sv[njetmxAK8], selpfjetAK8matchedelcleta[njetmxAK8], selpfjetAK8matchedelclphi[njetmxAK8], selpfjetAK8matchedelclrawE[njetmxAK8]; 

 float selpfjetAK8matchedelsigmaieta[njetmxAK8];
 float selpfjetAK8matchedelsigmaiphi[njetmxAK8];
 float selpfjetAK8matchedelr9full[njetmxAK8];
 float selpfjetAK8matchedelsupcl_etaw[njetmxAK8];
 float selpfjetAK8matchedelsupcl_phiw[njetmxAK8];
 float selpfjetAK8matchedelhcaloverecal[njetmxAK8];
 float selpfjetAK8matchedelcloctftrkn[njetmxAK8];
 float selpfjetAK8matchedelcloctftrkchi2[njetmxAK8];
 float selpfjetAK8matchedele1x5bye5x5[njetmxAK8];
 float selpfjetAK8matchedelnormchi2[njetmxAK8];
 float selpfjetAK8matchedelhitsmiss[njetmxAK8];
 float selpfjetAK8matchedeltrkmeasure[njetmxAK8];
 float selpfjetAK8matchedelconvtxprob[njetmxAK8];
 float selpfjetAK8matchedelecloverpout[njetmxAK8];
 float selpfjetAK8matchedelecaletrkmomentum[njetmxAK8];
 float selpfjetAK8matchedeldeltaetacltrkcalo[njetmxAK8];
 float selpfjetAK8matchedelsupcl_preshvsrawe[njetmxAK8];
 float selpfjetAK8matchedelpfisolsumphet[njetmxAK8];
 float selpfjetAK8matchedelpfisolsumchhadpt[njetmxAK8];
 float selpfjetAK8matchedelpfisolsumneuhadet[njetmxAK8];
 float selpfjetAK8matchedeletain[njetmxAK8];
 float selpfjetAK8matchedelphiin[njetmxAK8];
 float selpfjetAK8matchedelfbrem[njetmxAK8];
 float selpfjetAK8matchedeleoverp[njetmxAK8];
 float selpfjetAK8matchedelhovere[njetmxAK8];
 float selpfjetAK8matchedelRho[njetmxAK8];

 float selpfjetAK8sub1mass[njetmxAK8], selpfjetAK8sub1btag[njetmxAK8], selpfjetAK8sub1hadfrac[njetmxAK8], selpfjetAK8sub1emfrac[njetmxAK8];
 float selpfjetAK8sub2mass[njetmxAK8], selpfjetAK8sub2btag[njetmxAK8], selpfjetAK8sub2hadfrac[njetmxAK8], selpfjetAK8sub2emfrac[njetmxAK8];
 float selpfjetAK8subbtag[njetmxAK8], selpfjetAK8subhaddiff[njetmxAK8], selpfjetAK8subemdiff[njetmxAK8], selpfjetAK8subptdiff[njetmxAK8];
 
 bool selpfjetAK8haselectron[njetmxAK8], selpfjetAK8hasmuon[njetmxAK8], selpfjetAK8hasqg[njetmxAK8], selpfjetAK8hasb[njetmxAK8], selpfjetAK8hashadtop[njetmxAK8], selpfjetAK8hasleptop[njetmxAK8], selpfjetAK8hastop[njetmxAK8];
 bool pfjetAK8hashadtop_alldecay[njetmxAK8],pfjetAK8hasleptop_alldecay[njetmxAK8];
 bool selpfjetAK8hashadtop_alldecay[njetmxAK8],selpfjetAK8hasleptop_alldecay[njetmxAK8];
 float selpfjetAK8re_tvsb[njetmxAK8], selpfjetAK8rnu_tvsb[njetmxAK8];
 
 int npfjetAK8_thad;
 int npfjetAK8_te;
 int npfjetAK8_tmu;
 int npfjetAK8_qg;
 int npfjetAK8_b;
 int npfjetAK8_all;
 
 float puWeight, puWeightUp, puWeightDown;
 
 static const int nmaxjet = 10;
 /*
 static const int nhist = 75;
   const char *branchnames[nhist] = {
	  "pfjetAK8pt","pfjetAK8y","pfjetAK8phi","pfjetAK8mass",
	  "pfjetAK8btag_CMVA","pfjetAK8btag_CSV","pfjetAK8btag_DeepCSV","pfjetAK8matchAK4deepb",
	  "pfjetAK8DeepTag_TvsQCD","pfjetAK8DeepTag_WvsQCD","pfjetAK8DeepTag_ZvsQCD",
	  "pfjetAK8CHF","pfjetAK8NHF","pfjetAK8CEMF","pfjetAK8NEMF","pfjetAK8MUF","pfjetAK8HOF","pfjetAK8HadF","pfjetAK8NHadF","pfjetAK8EmF","pfjetAK8neuemfrac","pfjetAK8neunhadfrac",
	  "pfjetAK8chrad","pfjetAK8pTD","pfjetAK8sdmass","pfjetAK8tau21","pfjetAK8tau32",
	  "pfjetAK8sub1mass","pfjetAK8sub1btag","pfjetAK8sub1hadfrac","pfjetAK8sub1emfrac",
	  "pfjetAK8sub2mass","pfjetAK8sub2btag","pfjetAK8sub2hadfrac","pfjetAK8sub2emfrac",
	  "pfjetAK8subhaddiff","pfjetAK8subemdiff","pfjetAK8subptdiff","pfjetAK8subbtag",
	  "pfjetAK8_leppt","pfjetAK8_bpt","pfjetAK8_nupt","pfjetAK8_Rnew","pfjetAK8_bbyW_E",
	    "pfjetAK8_Kfactor","pfjetAK8_re_tvsb","pfjetAK8_rnu_tvsb","pfjetAK8_eldxy_sv",
          "pfjetAK8_sigmaieta", "pfjetAK8_sigmaiphi","pfjetAK8_r9full","pfjetAK8_etaw",
          "pfjetAK8_phiw","pfjetAK8_hcaloverecal","pfjetAK8_cloctftrkn","pfjetAK8_cloctftrkchi2",
          "pfjetAK8_e1x5bye5x5","pfjetAK8_normchi2","pfjetAK8_hitsmiss","pfjetAK8_trkmeasure",
          "pfjetAK8_ecloverpout","pfjetAK8_ecaletrkmomentum","pfjetAK8_deltaetacltrkcalo",
          "pfjetAK8_preshvsrawe","pfjetAK8_pfisolsumphet","pfjetAK8_pfisolsumchhadpt",
          "pfjetAK8_pfisolsumneuhadet","pfjetAK8_etain","pfjetAK8_phiin","pfjetAK8_fbrem",
          "pfjetAK8_eoverp","pfjetAK8_hovere","pfjetAK8_Rho","pfjetAK8_elpt","pfjetAK8_supeleta" 
	  };
  	  
   static const int nhist2D = 24;

   const char *branch2Dnames[nhist2D] = {"h2d_el_subjet", "h2d_el_subjet_gRe", "h2d_el_subjet_lRe", "h2d_el_subjeteta", "h2d_el_subjeteta_gRe", "h2d_el_subjeteta_lRe", "h2d_elcl_subjeteta", "h2d_elcl_subjeteta_gRe", "h2d_elcl_subjeteta_lRe", "h2d_elcl_subjete", "h2d_elcl_subjete_gRe", "h2d_elcl_subjete_lRe", "h2d_el_subjetM", "h2d_el_subjetM_gRe", "h2d_el_subjetM_lRe", "h2d_el_subjetetaM", "h2d_el_subjetetaM_gRe", "h2d_el_subjetetaM_lRe", "h2d_elcl_subjetetaM", "h2d_elcl_subjetetaM_gRe", "h2d_elcl_subjetetaM_lRe", "h2d_elcl_subjeteM", "h2d_elcl_subjeteM_gRe", "h2d_elcl_subjeteM_lRe"};
   
   static const int nhist2Dnew = 8;
   const char *branch2Dnamesnew[nhist2Dnew] = {"h2d_eldxy_subelpt", "h2d_eldxy_sv_subelpt", "h2d_eldxy_eldxy_sv", "h2d_eldxy_elID", "h2d_eldxy_sv_elID", "h2d_eldxy_sv_elcleta", "h2d_eldxy_sv_elclrawE", "h2d_eldxy_sv_elE"};

  double hist_low[nhist] = {
	  400,-2.5,-3.14,0,
	  0,0,0,0,
	  0,0,0,
	  0,0,0,0,0,0,0,0,0,0,0,
	  -0.25,0,0,0,0,
	  0,0,0,0,
	  0,0,0,0,
	  0,0,0,0,
	  100,100,100,
	  0,0,-2,
	  -1,-1,-0.2,0,0,0,0,0,0,-5,0,-1,0,0,0,0,-2.5,-2.5,0,0,0,0,-4,-2.5,-100,0,0,0,0,-3
	  };
  double hist_lowx2D[nhist2Dnew] = {-0.01,-0.01,-0.01,-0.01,-0.01,-0.01,-0.01,-0.01};
  double hist_highx2D[nhist2Dnew] = {0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02};
  double hist_lowy2D[nhist2Dnew] = {0,0,-0.01,-0.01,-0.01,-5.0,0,0};
  double hist_highy2D[nhist2Dnew] = {2000,2000,0.02,1.05,1.05,5.0,2000,2000};
  int hist_nbins2Dnew[nhist2Dnew] = {200,200,200,200,200,200,200,200};

  double hist_low2D[nhist2D] = {
    0,0,0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,0,0,0,0,0,0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,0,0,0};

  double hist_up[nhist] = {
	  3100,2.5,3.14,300,
	  1,1,1,1,
	  1,1,1,
	  1,1,1,1,1,1,1,1,1,1,1,
	  0.25,1,300,1,1,
	  300,1,1,1,
	  300,1,1,1,
	  1,1,1,1,
	  3100,3100,3100,
	  3.14,2.5,5,
	  1,1,0.2,0.1,0.1,50,1,1,500,20,200,1,2000,10,20,10000,2.5,2.5,2,1500,1500,1000,4,2.5,100,1000,500,100,1500,3
	  };

  double hist_up2D[nhist2D] = {
    2000,2000,2000,5.0,5.0,5.0,5.0,5.0,5.0,2000,2000,2000,2000,2000,2000,5.0,5.0,5.0,5.0,5.0,5.0,2000,2000,2000};
  


  int hist_nbins[nhist] = {
	  25,25,25,25,
	  20,20,20,20,
	  20,20,20,
	  25,25,25,25,25,25,25,25,25,25,25,
	  50,25,25,25,25,
	  25,25,25,25,
	  25,25,25,25,
	  25,25,25,25,
	  30,30,30,
	  30,25,35,
	  1000,16,2000,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200
  };


  int hist_nbins2D[nhist2D] = {200,200,200,100,100,100,100,100,100,200,200,200,200,200,200,100,100,100,100,100,100,200,200,200};
*/
  double pu_rat18[100] =    {15.0025,9.25619,7.25497,5.06682,3.81133,3.00109,2.47446,2.15125,1.91921,1.76245,1.65799,1.5848,1.53433,1.49373,1.46908,1.44313,1.42046,1.40154,1.37988,1.35679,1.3361,1.31139,1.28604,1.26065,1.23868,1.21481,1.19261,1.17143,1.15268,1.13366,1.11664,1.10102,1.08379,1.06829,1.05251,1.03584,1.01745,1.00004,0.980328,0.958154,0.935215,0.910953,0.885202,0.856835,0.827565,0.797445,0.765605,0.733814,0.701484,0.668361,0.634682,0.60224,0.570207,0.537462,0.505992,0.475353,0.445231,0.417069,0.390368,0.363391,0.339587,0.316478,0.293616,0.272703,0.253076,0.23398,0.216635,0.200789,0.185823,0.171907,0.159236,0.148382,0.13732,0.127015,0.11736,0.108589,0.0994979,0.0923745,0.0861297,0.0803695,0.0741731,0.0689201,0.0636846,0.0574831,0.053418,0.0491854,0.0468157,0.0431244,0.0407697,0.0352331,0.0327207,0.0311201,0.0284833,0.0257933,0.0234082,0.0200422,0.0199332,0.0189272,0.020733,0.0166171};
  double pu_rat18_up[100] = {0,11.3701,49.1593,16.3978,10.4484,7.79227,5.70396,4.15872,3.02768,2.28549,1.82582,1.52983,1.3595,1.2554,1.19605,1.1684,1.16115,1.17185,1.18964,1.20936,1.22873,1.23491,1.23159,1.21107,1.18259,1.14644,1.11133,1.08136,1.05384,1.03331,1.01987,1.01367,1.01107,1.01298,1.01865,1.02593,1.03512,1.0447,1.05099,1.0554,1.05447,1.04466,1.02824,1.00332,0.965566,0.923431,0.871249,0.814665,0.752156,0.689408,0.624858,0.564,0.505617,0.452167,0.402,0.359344,0.321227,0.285921,0.258403,0.233682,0.210464,0.192413,0.174424,0.159861,0.146181,0.131623,0.119227,0.10899,0.0963316,0.086803,0.0773651,0.0712667,0.0629173,0.0552031,0.0481823,0.0455058,0.0376989,0.0339163,0.0298286,0.0264131,0.0255965,0.0179475,0.0169746,0.0136435,0.0117583,0.00988318,0.00674005,0.00661599,0.00316237,0.00149674,0.0010104,0.00106782,0.000384941,0.000591271,0.000423128,0.000165822,7.60044e-05,4.96232e-05,7.51979e-05,1.05862e-05};
  double pu_rat18_dn[100] = {0,15.0557,67.8751,22.3278,14.1211,10.4821,7.88069,5.86513,4.31762,3.35551,2.78627,2.40097,2.16428,2.00485,1.9056,1.85092,1.82051,1.80608,1.78719,1.75544,1.71117,1.64481,1.57234,1.49261,1.42092,1.35612,1.3043,1.26517,1.23118,1.20443,1.18302,1.16596,1.14834,1.13047,1.11055,1.08517,1.05388,1.01479,0.96502,0.907499,0.841466,0.767187,0.68971,0.610695,0.530471,0.45611,0.385995,0.32355,0.268127,0.221267,0.181416,0.149012,0.122387,0.100955,0.0832931,0.0694147,0.0579993,0.0482614,0.0406839,0.0341693,0.0284128,0.0238208,0.0196651,0.0163071,0.0134164,0.0108213,0.00875349,0.00713274,0.00561523,0.00450669,0.00357902,0.00293888,0.00231295,0.00180802,0.00140385,0.00117654,0.000861839,0.000682485,0.000525487,0.000404909,0.00033922,0.000204219,0.000164688,0.000112084,8.12391e-05,5.70485e-05,3.2298e-05,2.61592e-05,1.02574e-05,3.96059e-06,2.16985e-06,1.85204e-06,5.36884e-07,6.60936e-07,3.78607e-07,1.19189e-07,4.4536e-08,2.4673e-08,3.47283e-08,5.35281e-09};
  
  static const int nobshist = 23;
  
  const char *obsnames[nobshist] = {
	  "pt","y","mass",
	  "NHad","neuhad","sdmass","chrad","subhaddiff","tau21",
	  "DAK8_topvsQCD","bbyW_E","Kfactor",
	  "re","rnu",
	  "hadsdmass",
	  "haspfelectron",
	  "haspfmuon",
	  "deltaR_muon",
	  "PDGId","Numdaught","Isleptop","Bidaughtid",
	  "rt"
	  };
	  
 double obs_low[nobshist] = {400,-2.5,0,0,0,0,-2.5,0,0,0,0,-2,-1,-1,0,-0.5,-0.5,0,-0.5,-0.5,-0.5,-0.5,-1}; 
 double obs_up[nobshist] = {3100,2.5,300,1,1,300,2.5,1,1,1,2.5,4.5,1,1,300,1.5,1.5,3,19.5,5.5,1.5,12.5,1}; 
 int obs_nbins[nobshist] = {25,25,25,25,25,25,25,25,25,20,25,25,40,40,25,2,2,10,20,6,2,13,25};

 TH1D *hist_obs[nobshist] ;
 TH1D *hist_obspuwup[nobshist] ;
 TH1D *hist_obspuwdown[nobshist] ;
 TH1D *hist_obsbtagwup[nobshist] ;
 TH1D *hist_obsbtagwdown[nobshist] ;

 TH1D *hist_obswel[nobshist] ;
 TH1D *hist_obswelpuwup[nobshist] ;
 TH1D *hist_obswelpuwdown[nobshist] ;
 TH1D *hist_obswelbtagwup[nobshist] ;
 TH1D *hist_obswelbtagwdown[nobshist] ;

 TH1D *hist_obsnoel[nobshist] ;
 TH1D *hist_obsnoelpuwup[nobshist] ;
 TH1D *hist_obsnoelpuwdown[nobshist] ;
 TH1D *hist_obsnoelbtagwup[nobshist] ;
 TH1D *hist_obsnoelbtagwdown[nobshist] ;

 TH1D *hist_obsljesup[3]; 
 TH1D *hist_obsljesdown[3];
 TH1D *hist_obsljerup[3]; 
 TH1D *hist_obsljerdown[3];

 TH1D *hist_obswelljesup[3];
 TH1D *hist_obswelljesdown[3];
 TH1D *hist_obswelljerup[3];
 TH1D *hist_obswelljerdown[3];

 TH1D *hist_obsnoelljesup[3];
 TH1D *hist_obsnoelljesdown[3];
 TH1D *hist_obsnoelljerup[3];
 TH1D *hist_obsnoelljerdown[3];

 TH1D *hist_obs_mu[nobshist];
 TH1D *hist_obs_mupuwup[nobshist];
 TH1D *hist_obs_mupuwdown[nobshist];
 TH1D *hist_obs_mubtagwup[nobshist] ;
 TH1D *hist_obs_mubtagwdown[nobshist] ;

 TH1D *hist_obs_wmumu[nobshist];
 TH1D *hist_obs_wmumupuwup[nobshist];
 TH1D *hist_obs_wmumupuwdown[nobshist];
 TH1D *hist_obs_wmumubtagwup[nobshist] ;
 TH1D *hist_obs_wmumubtagwdown[nobshist] ;

 TH1D *hist_obs_muljesup[3];
 TH1D *hist_obs_muljesdown[3];
 TH1D *hist_obs_muljerup[3];
 TH1D *hist_obs_muljerdown[3];

 TH1D *hist_obs_wmumuljesup[3];
 TH1D *hist_obs_wmumuljesdown[3];
 TH1D *hist_obs_wmumuljerup[3];
 TH1D *hist_obs_wmumuljerdown[3];

 const char *obsnames_mu[nobshist] = {
   "pt","y","mass",
   "NHad","neuhad","sdmass","chrad","subhaddiff","tau21",
   "DAK8_topvsQCD","bbyW_E","Kfactor",
   "re","rmu",
   "hadsdmass",
   "haspfelectron",
   "haspfmuon",
   "deltaR_electron",
   "PDGId","Numdaught","Isleptop","Bidaughtid",
   "rt"
 };

 
 TH2D *hist_2D_msd_deepak8;
 TH2D *hist_2D_bpass_flavb;
 TH2D *hist_2D_bpass_flavc;
 TH2D *hist_2D_bpass_flavq;
 TH2D *hist_2D_ball_flavb;
 TH2D *hist_2D_ball_flavc;
 TH2D *hist_2D_ball_flavq;
 /*
   TH1D *hist_top_deepak8_pass;
   TH1D *hist_top_deepak8;
 
   TH1D *hist_th_pt;
   TH1D *hist_th_y;
   TH1D *hist_th_sdmass;
   TH1D *hist_th_tau32;
   TH1D *hist_th_deepak8;
   
   TH1D *hist_dilep_tag;
   TH1D *hist_semilep_tag;
 */
 TH1D *hist_count;
 
 static const int noptbins = 32;
 double ptbins[noptbins+1] = {395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103};
  
 
  static const int nobptbins = 9;
  double bptbins[nobptbins+1] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
  
  static const int nobetabins = 3;
  double betabins[nobetabins+1] = {0, 0.6, 1.2, 2.5};
 
  static const int notptbins = 3;
  double tptbins[notptbins+1] = {400, 480, 600, 1200};
 
  float DAK8_topcut = 0.470; //0.470 ;//1% mistag rate (0.920 // 0.1% mistag) 
  float deep_btag_cut = 0.2770; //0.0494//loose //(0.2770//medium) (0.7264 //tight) => These are for Autumn18
  //for UL18 => 0.0490: loose, 0.2783: medium, 0.7100: tight 
  float re_cut = 0.3;
  float rt_cut = 0.7;//0.725;
  
  static const int nre = 100;
  
  TProofOutputFile *OutFile;
  TFile *fileOut;
   
  TTree *Tout ;
  TTree *Tout1 ;
  TTree *Tout2 ;
  TTree *Tout3 ;
  TTree *Tout4 ;
  TTree *Tout5 ;
  TTree *Tout6 ;
  

  TH1D *hist_dReb;
  TH1D *hist_dRnue;
  TH1D *hist_dRnub;
  TH1D *hist_dRet;
  TH1D *hist_dRnut;
  TH1D *hist_dRbt;

  TH2D *hist2D_dReb;
  TH2D *hist2D_dRnue;
  TH2D *hist2D_dRnub;

  TH2D *hist2D_dRet;
  TH2D *hist2D_dRnut;
  TH2D *hist2D_dRbt;
  
  TH1D *hist_init[18];
  TH1D *hist_binit[6];
  TH1D *hist_dRdPhi[14];
  TH1D *hist_initpuwup[18];
  TH1D *hist_initpuwdown[18];
  TH1D *hist_initbtagwup[18];
  TH1D *hist_initbtagwdown[18];

  const char *initnames[18] = {
    "nmu","nel","PFMET","N_PV_sel","NJets_AK4","NBJets_AK4","NJets_AK8","mll","l1pt","l1eta","l1phi","l2pt","l2eta","l2phi","bjetpt","bjeteta","bjetphi","NJets_balAK8"};//,"dRmuelcand","dPhimuelcand","Mmuelcand","Ptmuelcand","njets_AK4_out"};

  const char *titlenames[18] = {"nmu","nel","PFMET","# of Primary Vertices","# of AK4 jets","# of b tagged AK4 jets","# of AK8 jets","mll","l1pt","l1eta","l1phi","l2pt","l2eta","l2phi","bjet pT","bjet eta","bjet phi","# of bal.AK8 jets wrt mu"};//"dR between top candidates","dPhi between top candidates","Invariant mass of top candidates","pT of top candidates","# of AK4 jets outside top candidates"}; 

  double ini_low[18] = {0.5,0.5,0.0,-0.1,0.5,0.5,0.5,0.0,25.0,-2.5,-5.0,25.0,-2.5,-5.0,25.0,-2.5,-5.0,0.5};//0.0,-3.0,0.0,400.0,-0.5};
  double ini_up[18] = {10.5,10.5,500,99.9,10.5,10.5,10.5,500.0,1000.0,2.5,5.0,1000.0,2.5,5.0,1000.0,2.5,5.0,10.5};//5.0,3.0,1500.0,3100.0,10.5};
  int ini_nbins[18] = {10,10,25,100,10,10,10,25,50,25,50,50,25,50,50,25,50,10};//,50,50,50,25,10};
  const char *binitnames[6] = {"bNJets_AK4","bNBJets_AK4","jet1pt","jet2pt","bjet1pt","bjet2pt"}; 
  const char *btitlenames[6] = {"# of AK4 jets with boost", "# of b tagged AK4 jets with boost", "Leading AK4 jet pT with boost", "Subleading AK4 jet pT with boost", "Leading bjet pT with boost", "Subleading bjet pT with boost"};

  double ini_blow[6] = {-0.1,-0.1,25.0,25.0,25.0,25.0};
  double ini_bup[6] = {10.5,10.5,1000.0,1000.0,1000.0,1000.0};
  int ini_bnbins[6] = {10,10,50,50,50,50};


  const char *drnames[14] = {"dR_ljmu", "dPhi_ljmu", "dR_ljel", "dPhi_ljel", "dR_ljmu_mu", "dPhi_ljmu_mu", "dR_ljmu_el", "dPhi_ljmu_el", "dR_ljel_mu", "dPhi_ljel_mu", "dR_ljel_el", "dPhi_ljel_el", "dR_muel", "dPhi_muel"}; 

  const char *drtitlenames[14] = {"dR lj,#mu", "dPhi lj,#mu", "dR lj,el", "dPhi lj,el", "dR lj#mu,#mu", "dPhi lj#mu,#mu", "dR lj#mu,el", "dPhi_lj#mu,el", "dR ljel,#mu", "dPhi ljel,#mu", "dR ljel,el", "dPhi ljel,el", "dR #mu el", "dPhi #mu el"};

  double dr_low[14] = {0.0,-3.0,0.0,-3.0,0.0,-3.0,0.0,-3.0,0.0,-3.0,0.0,-3.0,0.0,-3.0};
  double dr_high[14] = {5.0,3.0,5.0,3.0,5.0,3.0,5.0,3.0,5.0,3.0,5.0,3.0,5.0,3.0};


  TH1D *hist_npv;
  TH1D *hist_npv_nopuwt;


  TH1D *hist_event_count;
  
  float in_pfjetAK8NHadF;
  float in_pfjetAK8neunhadfrac;
  float in_pfjetAK8subhaddiff;
  float in_pfjetAK8tau21;
  float in_pfjetAK8chrad;
  float in_pfjetAK8sdmass;
  float in_pfjetAK8eldxy_sv;
  float in_pfjetAK8matchedelcleta;
  float in_pfjetAK8matchedelpt;
  float in_pfjetAK8matchedelsigmaieta;
  float in_pfjetAK8matchedelsigmaiphi;
  float in_pfjetAK8matchedelr9full;
  float in_pfjetAK8matchedelsupcl_etaw;
  float in_pfjetAK8matchedelsupcl_phiw;
  float in_pfjetAK8matchedelhcaloverecal;
  float in_pfjetAK8matchedelcloctftrkn;
  float in_pfjetAK8matchedelcloctftrkchi2;
  float in_pfjetAK8matchedele1x5bye5x5;
  float in_pfjetAK8matchedelnormchi2;
  float in_pfjetAK8matchedelhitsmiss;
  float in_pfjetAK8matchedeltrkmeasure;
  float in_pfjetAK8matchedelecloverpout;
  float in_pfjetAK8matchedelecaletrkmomentum;
  float in_pfjetAK8matchedeldeltaetacltrkcalo;
  float in_pfjetAK8matchedelsupcl_preshvsrawe;
  float in_pfjetAK8matchedelpfisolsumphet;
  float in_pfjetAK8matchedelpfisolsumchhadpt;
  float in_pfjetAK8matchedelpfisolsumneuhadet;
  float in_pfjetAK8matchedeletain;
  float in_pfjetAK8matchedelphiin;
  float in_pfjetAK8matchedelfbrem;
  float in_pfjetAK8matchedeleoverp;
  float in_pfjetAK8matchedelhovere;
  float in_pfjetAK8matchedelRho;
  
  float in_pfjetAK8matchedmuonchi;
  float in_pfjetAK8matchedmuonposmatch;
  float in_pfjetAK8matchedmuontrkink;
  float in_pfjetAK8matchedmuonsegcom;
  float in_pfjetAK8matchedmuonhit;
  float in_pfjetAK8matchedmuonmst;
  float in_pfjetAK8matchedmuontrkvtx;
  float in_pfjetAK8matchedmuondz;
  float in_pfjetAK8matchedmuonpixhit;
  float in_pfjetAK8matchedmuontrklay;
  float in_pfjetAK8matchedmuonvalfrac;
  float in_pfjetAK8muinsubptrat;
  float in_pfjetAK8muinsubmassrat;
  float in_pfjetAK8muinsubinvmass;
  float in_pfjetAK8muinsubIfarbyI0;
  float in_pfjetAK8muinsubInearbyI0;
  
  TMVA::Reader *reader1;
  TMVA::Reader *reader2;
  TMVA::Reader *reader3;
  TMVA::Reader *reader4;

  //  TString dir = "/home/deroy/t3store3/CMSSW_10_5_0/src/BDTResponse_validator/Analysis/newvar_sv/Signal/";
  TString dir = "/home/deroy/t3store3/Muon_MuEl/";
  //TString weightfile1 = dir + TString("TMVAClassification_BDTG_elIDvarv3.weights.xml");
  TString weightfile1 = dir + TString("TMVAClassification_BDTG_elIDvar_Jan2021Corr_TTbarUL18.weights.xml");  
  TString weightfile2 = dir + TString("TMVAClassification_BDTG_rnu.weights.xml");
  //TString weightfile2 = dir + TString("TMVAClassification_BDTG_muIDvar_Jan2021Corr_TTbarUL18.weights.xml");
  TString weightfile3 = dir + TString("TMVAClassification_only7varsnomatchel_BDTG.weights.xml");
  TString weightfile4 = dir + TString("TMVAClassification_BDTG_muIDvar_Jan2021Corr_TTbarUL18.weights.xml");

  float ptcut = 400;

   Anal_Leptop_PROOF(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Anal_Leptop_PROOF() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Anal_Leptop_PROOF,0);
};

#endif

#ifdef Anal_Leptop_PROOF_cxx
void Anal_Leptop_PROOF::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   
   fChain->SetBranchAddress("irun", &irun, &b_irun);
   fChain->SetBranchAddress("ilumi", &ilumi, &b_ilumi);
   fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
   fChain->SetBranchAddress("nprim", &nprim, &b_nprim);
   fChain->SetBranchAddress("nprimi", &nprimi, &b_nprimi);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
   fChain->SetBranchAddress("npu_vert", &npu_vert, &b_npu_vert);
   fChain->SetBranchAddress("trig_value", &trig_value, &b_trig_value);
   fChain->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24,&b_hlt_IsoMu24);
   fChain->SetBranchAddress("hlt_Mu50",&hlt_Mu50,&b_hlt_Mu50);
   fChain->SetBranchAddress("hlt_Ele50_PFJet165",&hlt_Ele50_PFJet165,&b_hlt_Ele50_PFJet165);
   fChain->SetBranchAddress("hlt_Ele115",&hlt_Ele115,&b_hlt_Ele115);
   fChain->SetBranchAddress("hlt_AK8PFJet500", &hlt_AK8PFJet500, &b_hlt_AK8PFJet500);
   fChain->SetBranchAddress("hlt_Photon200", &hlt_Photon200, &b_hlt_Photon200);
   fChain->SetBranchAddress("hlt_Mu37Ele27", &hlt_Mu37Ele27, &b_hlt_Mu37Ele27);
   fChain->SetBranchAddress("hlt_Mu27Ele37", &hlt_Mu27Ele37, &b_hlt_Mu27Ele37);
   fChain->SetBranchAddress("hlt_Mu37TkMu27", &hlt_Mu37TkMu27, &b_hlt_Mu37TkMu27);
   fChain->SetBranchAddress("hlt_OldMu100", &hlt_OldMu100, &b_hlt_OldMu100);
   fChain->SetBranchAddress("hlt_TkMu100", &hlt_TkMu100, &b_hlt_TkMu100);
   fChain->SetBranchAddress("hlt_DoubleEle25", &hlt_DoubleEle25, &b_hlt_DoubleEle25);
   fChain->SetBranchAddress("ntrigobjs", &ntrigobjs, &b_ntrigobjs);
   fChain->SetBranchAddress("trigobjpt", trigobjpt, &b_trigobjpt);
   fChain->SetBranchAddress("trigobjeta", trigobjeta, &b_trigobjeta);
   fChain->SetBranchAddress("trigobjphi", trigobjphi, &b_trigobjphi);
   fChain->SetBranchAddress("trigobjmass", trigobjmass, &b_trigobjmass);
   fChain->SetBranchAddress("trigobjpdgId", trigobjpdgId, &b_trigobjpdgId);
   fChain->SetBranchAddress("trigobjHLT", trigobjHLT, &b_trigobjHLT);
   fChain->SetBranchAddress("trigobjL1", trigobjL1, &b_trigobjL1);
   fChain->SetBranchAddress("trigobjBoth", trigobjBoth, &b_trigobjBoth);
   fChain->SetBranchAddress("trigobjIhlt", trigobjIhlt, &b_trigobjIhlt);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_miset);
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_misphi);
   fChain->SetBranchAddress("MisEtSig", &MisEtSig, &b_misetsig);
   fChain->SetBranchAddress("sumEt", &sumEt, &b_sumEt);
   fChain->SetBranchAddress("npfjetAK8", &npfjetAK8, &b_npfjetAK8);
   fChain->SetBranchAddress("pfjetAK8jetID_tightlepveto", pfjetAK8jetID_tightlepveto, &b_pfjetAK8jetID_tightlepveto);
   fChain->SetBranchAddress("pfjetAK8jetID", pfjetAK8jetID, &b_pfjetAK8jetID);
   fChain->SetBranchAddress("pfjetAK8pt", pfjetAK8pt, &b_pfjetAK8pt);
   fChain->SetBranchAddress("pfjetAK8y", pfjetAK8y, &b_pfjetAK8y);
   fChain->SetBranchAddress("pfjetAK8phi", pfjetAK8phi, &b_pfjetAK8phi);
   fChain->SetBranchAddress("pfjetAK8mass", pfjetAK8mass, &b_pfjetAK8mass);
   fChain->SetBranchAddress("pfjetAK8JEC", pfjetAK8JEC, &b_pfjetAK8JEC);
   fChain->SetBranchAddress("pfjetAK8btag_DeepCSV", pfjetAK8btag_DeepCSV, &b_pfjetAK8btag_DeepCSV);

   fChain->SetBranchAddress("pfjetAK8DeepTag_TvsQCD", pfjetAK8DeepTag_TvsQCD, &b_pfjetAK8DeepTag_TvsQCD);
   fChain->SetBranchAddress("pfjetAK8DeepTag_WvsQCD", pfjetAK8DeepTag_WvsQCD, &b_pfjetAK8DeepTag_WvsQCD);
   fChain->SetBranchAddress("pfjetAK8DeepTag_ZvsQCD", pfjetAK8DeepTag_ZvsQCD, &b_pfjetAK8DeepTag_ZvsQCD);
   fChain->SetBranchAddress("pfjetAK8CHF", pfjetAK8CHF, &b_pfjetAK8CHF);
   fChain->SetBranchAddress("pfjetAK8NHF", pfjetAK8NHF, &b_pfjetAK8NHF);
   fChain->SetBranchAddress("pfjetAK8CEMF", pfjetAK8CEMF, &b_pfjetAK8CEMF);
   fChain->SetBranchAddress("pfjetAK8NEMF", pfjetAK8NEMF, &b_pfjetAK8NEMF);
   fChain->SetBranchAddress("pfjetAK8MUF", pfjetAK8MUF, &b_pfjetAK8MUF);
   fChain->SetBranchAddress("pfjetAK8PHF", pfjetAK8PHF, &b_pfjetAK8PHF);
   fChain->SetBranchAddress("pfjetAK8EEF", pfjetAK8EEF, &b_pfjetAK8EEF);
   fChain->SetBranchAddress("pfjetAK8HFHF", pfjetAK8HFHF, &b_pfjetAK8HFHF);
   fChain->SetBranchAddress("pfjetAK8CHM", pfjetAK8CHM, &b_pfjetAK8CHM);
   fChain->SetBranchAddress("pfjetAK8NHM", pfjetAK8NHM, &b_pfjetAK8NHM);
   fChain->SetBranchAddress("pfjetAK8MUM", pfjetAK8MUM, &b_pfjetAK8MUM);
   fChain->SetBranchAddress("pfjetAK8PHM", pfjetAK8PHM, &b_pfjetAK8PHM);
   fChain->SetBranchAddress("pfjetAK8EEM", pfjetAK8EEM, &b_pfjetAK8EEM);
   fChain->SetBranchAddress("pfjetAK8HFHM", pfjetAK8HFHM, &b_pfjetAK8HFHM);
   fChain->SetBranchAddress("pfjetAK8Neucons", pfjetAK8Neucons, &b_pfjetAK8Neucons);
   fChain->SetBranchAddress("pfjetAK8Chcons", pfjetAK8Chcons, &b_pfjetAK8Chcons);
   fChain->SetBranchAddress("pfjetAK8reso", pfjetAK8reso, &b_pfjetAK8reso);
   fChain->SetBranchAddress("pfjetAK8resoup", pfjetAK8resoup, &b_pfjetAK8resoup);
   fChain->SetBranchAddress("pfjetAK8resodn", pfjetAK8resodn, &b_pfjetAK8resodn);
   fChain->SetBranchAddress("pfjetAK8jesup_pu", pfjetAK8jesup_pu, &b_pfjetAK8jesup_pu);
   fChain->SetBranchAddress("pfjetAK8jesup_rel", pfjetAK8jesup_rel, &b_pfjetAK8jesup_rel);
   fChain->SetBranchAddress("pfjetAK8jesup_scale", pfjetAK8jesup_scale, &b_pfjetAK8jesup_scale);
   fChain->SetBranchAddress("pfjetAK8jesup_total", pfjetAK8jesup_total, &b_pfjetAK8jesup_total);
   fChain->SetBranchAddress("pfjetAK8jesdn_pu", pfjetAK8jesdn_pu, &b_pfjetAK8jesdn_pu);
   fChain->SetBranchAddress("pfjetAK8jesdn_rel", pfjetAK8jesdn_rel, &b_pfjetAK8jesdn_rel);
   fChain->SetBranchAddress("pfjetAK8jesdn_scale", pfjetAK8jesdn_scale, &b_pfjetAK8jesdn_scale);
   fChain->SetBranchAddress("pfjetAK8jesdn_total", pfjetAK8jesdn_total, &b_pfjetAK8jesdn_total);
   fChain->SetBranchAddress("pfjetAK8chrad", pfjetAK8chrad, &b_pfjetAK8chrad);
   fChain->SetBranchAddress("pfjetAK8pTD", pfjetAK8pTD, &b_pfjetAK8pTD);
   fChain->SetBranchAddress("pfjetAK8sdmass", pfjetAK8sdmass, &b_pfjetAK8sdmass);
   fChain->SetBranchAddress("pfjetAK8tau1", pfjetAK8tau1, &b_pfjetAK8tau1);
   fChain->SetBranchAddress("pfjetAK8tau2", pfjetAK8tau2, &b_pfjetAK8tau2);
   fChain->SetBranchAddress("pfjetAK8tau3", pfjetAK8tau3, &b_pfjetAK8tau3);
   fChain->SetBranchAddress("pfjetAK8sub1pt", pfjetAK8sub1pt, &b_pfjetAK8sub1pt);
   fChain->SetBranchAddress("pfjetAK8sub1eta", pfjetAK8sub1eta, &b_pfjetAK8sub1eta);
   fChain->SetBranchAddress("pfjetAK8sub1phi", pfjetAK8sub1phi, &b_pfjetAK8sub1phi);
   fChain->SetBranchAddress("pfjetAK8sub1mass", pfjetAK8sub1mass, &b_pfjetAK8sub1mass);
   fChain->SetBranchAddress("pfjetAK8sub1btag", pfjetAK8sub1btag, &b_pfjetAK8sub1btag);
   fChain->SetBranchAddress("pfjetAK8sub1chhadfrac", pfjetAK8sub1chhadfrac, &b_pfjetAK8sub1chhadfrac);
   fChain->SetBranchAddress("pfjetAK8sub1neuhadfrac", pfjetAK8sub1neuhadfrac, &b_pfjetAK8sub1neuhadfrac);
   fChain->SetBranchAddress("pfjetAK8sub1emfrac", pfjetAK8sub1emfrac, &b_pfjetAK8sub1emfrac);
   fChain->SetBranchAddress("pfjetAK8sub1phofrac", pfjetAK8sub1phofrac, &b_pfjetAK8sub1phofrac);
   fChain->SetBranchAddress("pfjetAK8sub1mufrac", pfjetAK8sub1mufrac, &b_pfjetAK8sub1mufrac);
   fChain->SetBranchAddress("pfjetAK8sub2pt", pfjetAK8sub2pt, &b_pfjetAK8sub2pt);
   fChain->SetBranchAddress("pfjetAK8sub2eta", pfjetAK8sub2eta, &b_pfjetAK8sub2eta);
   fChain->SetBranchAddress("pfjetAK8sub2phi", pfjetAK8sub2phi, &b_pfjetAK8sub2phi);
   fChain->SetBranchAddress("pfjetAK8sub2mass", pfjetAK8sub2mass, &b_pfjetAK8sub2mass);
   fChain->SetBranchAddress("pfjetAK8sub2btag", pfjetAK8sub2btag, &b_pfjetAK8sub2btag);
   fChain->SetBranchAddress("pfjetAK8sub2chhadfrac", pfjetAK8sub2chhadfrac, &b_pfjetAK8sub2chhadfrac);
   fChain->SetBranchAddress("pfjetAK8sub2neuhadfrac", pfjetAK8sub2neuhadfrac, &b_pfjetAK8sub2neuhadfrac);
   fChain->SetBranchAddress("pfjetAK8sub2emfrac", pfjetAK8sub2emfrac, &b_pfjetAK8sub2emfrac);
   fChain->SetBranchAddress("pfjetAK8sub2phofrac", pfjetAK8sub2phofrac, &b_pfjetAK8sub2phofrac);
   fChain->SetBranchAddress("pfjetAK8sub2mufrac", pfjetAK8sub2mufrac, &b_pfjetAK8sub2mufrac);
   
   fChain->SetBranchAddress("pfjetAK8elinsubpt", pfjetAK8elinsubpt, &b_pfjetAK8elinsubpt);
   fChain->SetBranchAddress("pfjetAK8elinsubeta", pfjetAK8elinsubeta, &b_pfjetAK8elinsubeta);
   fChain->SetBranchAddress("pfjetAK8elinsubphi", pfjetAK8elinsubphi, &b_pfjetAK8elinsubphi);
   fChain->SetBranchAddress("pfjetAK8elinsubmass", pfjetAK8elinsubmass, &b_pfjetAK8elinsubmass);

   fChain->SetBranchAddress("pfjetAK8elinsubjpt", pfjetAK8elinsubjpt, &b_pfjetAK8elinsubjpt);
   fChain->SetBranchAddress("pfjetAK8elinsubjeta", pfjetAK8elinsubjeta, &b_pfjetAK8elinsubjeta);
   fChain->SetBranchAddress("pfjetAK8elinsubjphi", pfjetAK8elinsubjphi, &b_pfjetAK8elinsubjphi);
   fChain->SetBranchAddress("pfjetAK8elinsubjmass", pfjetAK8elinsubjmass, &b_pfjetAK8elinsubjmass);

   fChain->SetBranchAddress("pfjetAK8muinsubpt", pfjetAK8muinsubpt, &b_pfjetAK8muinsubpt);
   fChain->SetBranchAddress("pfjetAK8muinsubeta", pfjetAK8muinsubeta, &b_pfjetAK8muinsubeta);
   fChain->SetBranchAddress("pfjetAK8muinsubphi", pfjetAK8muinsubphi, &b_pfjetAK8muinsubphi);
   fChain->SetBranchAddress("pfjetAK8muinsubmass", pfjetAK8muinsubmass, &b_pfjetAK8muinsubmass);
   fChain->SetBranchAddress("pfjetAK8muinsubIfar", pfjetAK8muinsubIfar, &b_pfjetAK8muinsubIfar);
   fChain->SetBranchAddress("pfjetAK8muinsubInear", pfjetAK8muinsubInear, &b_pfjetAK8muinsubInear);
   fChain->SetBranchAddress("pfjetAK8muinsubI0", pfjetAK8muinsubI0, &b_pfjetAK8muinsubI0);


   fChain->SetBranchAddress("npfjetAK4", &npfjetAK4, &b_npfjetAK4);
   fChain->SetBranchAddress("pfjetAK4jetID", pfjetAK4jetID, &b_pfjetAK4jetID);
   fChain->SetBranchAddress("pfjetAK4jetID_tightlepveto", pfjetAK4jetID_tightlepveto, &b_pfjetAK4jetID_tightlepveto);
   fChain->SetBranchAddress("pfjetAK4pt", pfjetAK4pt, &b_pfjetAK4pt);
   fChain->SetBranchAddress("pfjetAK4eta", pfjetAK4eta, &b_pfjetAK4eta);
   fChain->SetBranchAddress("pfjetAK4y", pfjetAK4y, &b_pfjetAK4y);
   fChain->SetBranchAddress("pfjetAK4phi", pfjetAK4phi, &b_pfjetAK4phi);
   fChain->SetBranchAddress("pfjetAK4mass", pfjetAK4mass, &b_pfjetAK4mass);
   fChain->SetBranchAddress("pfjetAK4JEC", pfjetAK4JEC, &b_pfjetAK4JEC);
   
   fChain->SetBranchAddress("pfjetAK4btag_DeepCSV", pfjetAK4btag_DeepCSV, &b_pfjetAK4btag_DeepCSV);
   fChain->SetBranchAddress("pfjetAK4btag_DeepFlav", pfjetAK4btag_DeepFlav, &b_pfjetAK4btag_DeepFlav);
   
   fChain->SetBranchAddress("pfjetAK4reso", pfjetAK4reso, &b_pfjetAK4reso);
   fChain->SetBranchAddress("pfjetAK4resoup", pfjetAK4resoup, &b_pfjetAK4resoup);
   fChain->SetBranchAddress("pfjetAK4resodn", pfjetAK4resodn, &b_pfjetAK4resodn);
   fChain->SetBranchAddress("pfjetAK4jesup_pu", pfjetAK4jesup_pu, &b_pfjetAK4jesup_pu);
   fChain->SetBranchAddress("pfjetAK4jesup_rel", pfjetAK4jesup_rel, &b_pfjetAK4jesup_rel);
   fChain->SetBranchAddress("pfjetAK4jesup_scale", pfjetAK4jesup_scale, &b_pfjetAK4jesup_scale);
   fChain->SetBranchAddress("pfjetAK4jesup_total", pfjetAK4jesup_total, &b_pfjetAK4jesup_total);
   fChain->SetBranchAddress("pfjetAK4jesdn_pu", pfjetAK4jesdn_pu, &b_pfjetAK4jesdn_pu);
   fChain->SetBranchAddress("pfjetAK4jesdn_rel", pfjetAK4jesdn_rel, &b_pfjetAK4jesdn_rel);
   fChain->SetBranchAddress("pfjetAK4jesdn_scale", pfjetAK4jesdn_scale, &b_pfjetAK4jesdn_scale);
   fChain->SetBranchAddress("pfjetAK4jesdn_total", pfjetAK4jesdn_total, &b_pfjetAK4jesdn_total);
   fChain->SetBranchAddress("pfjetAK4hadronflav", pfjetAK4hadronflav, &b_pfjetAK4hadronflav);
   fChain->SetBranchAddress("pfjetAK4partonflav", pfjetAK4partonflav, &b_pfjetAK4partonflav);
   fChain->SetBranchAddress("pfjetAK4qgl", pfjetAK4qgl, &b_pfjetAK4qgl);
   fChain->SetBranchAddress("pfjetAK4PUID", pfjetAK4PUID, &b_pfjetAK4PUID);
   fChain->SetBranchAddress("pfjetAK4GenMatch", &pfjetAK4GenMatch, &b_pfjetAK4GenMatch);
   fChain->SetBranchAddress("GENMET", &GENMET, &b_genmiset);
   fChain->SetBranchAddress("GENMETPhi", &GENMETPhi, &b_genmisphi);
   
   fChain->SetBranchAddress("ngenjetAK8", &ngenjetAK8, &b_ngenjetAK8);
   fChain->SetBranchAddress("genjetAK8pt", genjetAK8pt, &b_genjetAK8pt);
   
   fChain->SetBranchAddress("genjetAK8phi", genjetAK8phi, &b_genjetAK8phi);
   fChain->SetBranchAddress("genjetAK8mass", genjetAK8mass, &b_genjetAK8mass);
   fChain->SetBranchAddress("genjetAK8sdmass", genjetAK8sdmass, &b_genjetAK8sdmass);
   
   fChain->SetBranchAddress("ngenjetAK4", &ngenjetAK4, &b_ngenjetAK4);
   fChain->SetBranchAddress("genjetAK4pt", genjetAK4pt, &b_genjetAK4pt);
   fChain->SetBranchAddress("genjetAK4phi", genjetAK4phi, &b_genjetAK4phi);
   fChain->SetBranchAddress("genjetAK4mass", genjetAK4mass, &b_genjetAK4mass);
   
   fChain->SetBranchAddress("ngenparticles", &ngenparticles, &b_ngenparticles);
   fChain->SetBranchAddress("genpartstatus", genpartstatus, &b_genpartstatus);
   fChain->SetBranchAddress("genpartpdg", genpartpdg, &b_genpartpdg);
   fChain->SetBranchAddress("genpartmompdg", genpartmompdg, &b_genpartmompdg);
   fChain->SetBranchAddress("genpartdaugno", genpartdaugno, &b_genpartdaugno);
   fChain->SetBranchAddress("genpartfromhard", genpartfromhard, &b_genpartfromhard);
   fChain->SetBranchAddress("genpartfromhardbFSR", genpartfromhardbFSR, &b_genpartfromhardbFSR);
   fChain->SetBranchAddress("genpartpt", genpartpt, &b_genpartpt);
   fChain->SetBranchAddress("genparteta", genparteta, &b_genparteta);
   fChain->SetBranchAddress("genpartphi", genpartphi, &b_genpartphi);
   fChain->SetBranchAddress("genpartm", genpartm, &b_genpartm);
   
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("muonisPF", muonisPF, &b_muonisPF);
   fChain->SetBranchAddress("muonisGL", muonisGL, &b_muonisGL);
   fChain->SetBranchAddress("muonisTRK", muonisTRK, &b_muonisTRK);
   fChain->SetBranchAddress("muonisLoose", muonisLoose, &b_muonisLoose);
   fChain->SetBranchAddress("muonisGoodGL", muonisGoodGL, &b_muonisGoodGL);
   fChain->SetBranchAddress("muonisMed", muonisMed, &b_muonisMed);
   fChain->SetBranchAddress("muonisTight", muonisTight, &b_muonisTight);
   fChain->SetBranchAddress("muonisHighPt", muonisHighPt, &b_muonisHighPt);
   fChain->SetBranchAddress("muonisHighPttrk", muonisHighPttrk, &b_muonisHighPttrk);
   fChain->SetBranchAddress("muonisMedPr", muonisMedPr, &b_muonisMedPr);
   fChain->SetBranchAddress("muonpt", muonpt, &b_muonpt);
   fChain->SetBranchAddress("muonp", muonp, &b_muonp);
   fChain->SetBranchAddress("muoneta", muoneta, &b_muoneta);
   fChain->SetBranchAddress("muonphi", muonphi, &b_muonphi);
   fChain->SetBranchAddress("muoncharge", muoncharge, &b_muoncharge);
   fChain->SetBranchAddress("muontrkvtx", muontrkvtx, &b_muontrkvtx);
   fChain->SetBranchAddress("muondz", muondz, &b_muondz);
   fChain->SetBranchAddress("muonpter", muonpter, &b_muonpter);
   fChain->SetBranchAddress("muonchi", muonchi, &b_muonchi);
   fChain->SetBranchAddress("muonndf", muonndf, &b_muonndf);
   fChain->SetBranchAddress("muonecal", muonecal, &b_muonecal);
   fChain->SetBranchAddress("muonhcal", muonhcal, &b_muonhcal);
   fChain->SetBranchAddress("muonpfiso", muonpfiso, &b_muonpfiso);
   fChain->SetBranchAddress("muonposmatch", muonposmatch, &b_muonposmatch);
   fChain->SetBranchAddress("muontrkink", muontrkink, &b_muontrkink);
   fChain->SetBranchAddress("muonsegcom", muonsegcom, &b_muonsegcom);
   fChain->SetBranchAddress("muonpixhit", muonpixhit, &b_muonpixhit);
   fChain->SetBranchAddress("muonmst", muonmst, &b_muonmst);
   fChain->SetBranchAddress("muontrklay", muontrklay, &b_muontrklay);
   fChain->SetBranchAddress("muonvalfrac", muonvalfrac, &b_muonvalfrac);
   fChain->SetBranchAddress("muonchiso", muonchiso, &b_muonchiso);
   fChain->SetBranchAddress("muonnhiso", muonnhiso, &b_muonnhiso);
   fChain->SetBranchAddress("muonphiso", muonphiso, &b_muonphiso);
   fChain->SetBranchAddress("muonminisoall", muonminisoall, &b_muonminisoall);
   fChain->SetBranchAddress("mudxy_sv", mudxy_sv, &b_mudxy_sv);
   fChain->SetBranchAddress("muonhit", muonhit, &b_muonhit);

   fChain->SetBranchAddress("nelecs", &nelecs, &b_nelecs);
   fChain->SetBranchAddress("elpt", elpt, &b_elpt);
   fChain->SetBranchAddress("eldxy_sv", eldxy_sv, &b_eldxy_sv);
   fChain->SetBranchAddress("eleta", eleta, &b_eleta);
   fChain->SetBranchAddress("elphi", elphi, &b_elphi);
   fChain->SetBranchAddress("elp", elp, &b_elp);
   fChain->SetBranchAddress("ele", ele, &b_ele);
   fChain->SetBranchAddress("elcharge", elcharge, &b_elcharge);
   fChain->SetBranchAddress("elmvaid", elmvaid, &b_elmvaid);
   fChain->SetBranchAddress("elmvaid_noIso", elmvaid_noIso, &b_elmvaid_noIso);

   fChain->SetBranchAddress("elhovere", elhovere, &b_elhovere);
   fChain->SetBranchAddress("elchi", elchi, &b_elchi);
   fChain->SetBranchAddress("elndf", elndf, &b_elndf);
   
   fChain->SetBranchAddress("eletain", eletain, &b_eletain);
   fChain->SetBranchAddress("elphiin", elphiin, &b_elphiin);
   
   fChain->SetBranchAddress("elsupcl_eta", elsupcl_eta, &b_elsupcl_eta);
   fChain->SetBranchAddress("elsupcl_phi", elsupcl_phi, &b_elsupcl_phi);
   fChain->SetBranchAddress("elsupcl_rawE", elsupcl_rawE, &b_elsupcl_rawE);
   fChain->SetBranchAddress("elfbrem", elfbrem, &b_elfbrem);
   
   fChain->SetBranchAddress("eleoverp", eleoverp, &b_eleoverp);
   fChain->SetBranchAddress("elietaieta", elietaieta, &b_elietaieta);
   fChain->SetBranchAddress("elmisshits", elmisshits, &b_elmisshits);
   
   fChain->SetBranchAddress("elpfiso", elpfiso, &b_elpfiso);
   
   
   
   fChain->SetBranchAddress("elsigmaieta", elsigmaieta, &b_elsigmaieta);
   fChain->SetBranchAddress("elsigmaiphi", elsigmaiphi, &b_elsigmaiphi);
   fChain->SetBranchAddress("elr9full", elr9full, &b_elr9full);
   fChain->SetBranchAddress("elsupcl_etaw", elsupcl_etaw, &b_elsupcl_etaw); 
   fChain->SetBranchAddress("elsupcl_phiw", elsupcl_phiw, &b_elsupcl_phiw);
   fChain->SetBranchAddress("elhcaloverecal", elhcaloverecal, &b_elhcaloverecal);
   fChain->SetBranchAddress("elcloctftrkn",  elcloctftrkn, &b_elcloctftrkn);
   fChain->SetBranchAddress("elcloctftrkchi2", elcloctftrkchi2, &b_elcloctftrkchi2);
   fChain->SetBranchAddress("ele1x5bye5x5", ele1x5bye5x5, &b_ele1x5bye5x5);
   fChain->SetBranchAddress("elnormchi2", elnormchi2, &b_elnormchi2);
   fChain->SetBranchAddress("elhitsmiss", elhitsmiss, &b_elhitsmiss);
   fChain->SetBranchAddress("eltrkmeasure", eltrkmeasure, &b_eltrkmeasure);
   fChain->SetBranchAddress("elconvtxprob", elconvtxprob, &b_elconvtxprob);
   fChain->SetBranchAddress("elecloverpout", elecloverpout, &b_elecloverpout);
   fChain->SetBranchAddress("elecaletrkmomentum", elecaletrkmomentum, &b_elecaletrkmomentum);
   fChain->SetBranchAddress("eldeltaetacltrkcalo", eldeltaetacltrkcalo, &b_eldeltaetacltrkcalo);
   fChain->SetBranchAddress("elsupcl_preshvsrawe", elsupcl_preshvsrawe, &b_elsupcl_preshvsrawe);
   fChain->SetBranchAddress("elpfisolsumphet", elpfisolsumphet, &b_elpfisolsumphet);
   fChain->SetBranchAddress("elpfisolsumchhadpt", elpfisolsumchhadpt, &b_elpfisolsumchhadpt);
   fChain->SetBranchAddress("elpfsiolsumneuhadet", elpfsiolsumneuhadet, &b_elpfsiolsumneuhadet);
   /*
   fChain->SetBranchAddress("nphotons", &nphotons, &b_nphotons);
   fChain->SetBranchAddress("phoe", phoe, &b_phoe);
   fChain->SetBranchAddress("phoeta", phoeta, &b_phoeta);
   fChain->SetBranchAddress("phophi", phophi, &b_phophi);
   fChain->SetBranchAddress("phomvaid", phomvaid, &b_phomvaid);
   fChain->SetBranchAddress("phoe1by9", phoe1by9, &b_phoe1by9);
   fChain->SetBranchAddress("phoe9by25", phoe9by25, &b_phoe9by25);
   fChain->SetBranchAddress("photrkiso", photrkiso, &b_photrkiso);
   fChain->SetBranchAddress("phoemiso", phoemiso, &b_phoemiso);
   fChain->SetBranchAddress("phohadiso", phohadiso, &b_phohadiso);
   fChain->SetBranchAddress("phochhadiso", phochhadiso, &b_phochhadiso);
   fChain->SetBranchAddress("phoneuhadiso", phoneuhadiso, &b_phoneuhadiso);
   fChain->SetBranchAddress("phophoiso", phophoiso, &b_phophoiso);
   fChain->SetBranchAddress("phoPUiso", phoPUiso, &b_phoPUiso);
   fChain->SetBranchAddress("phohadbyem", phohadbyem, &b_phohadbyem);
   fChain->SetBranchAddress("phoietaieta", phoietaieta, &b_phoietaieta);
   */
}

Bool_t Anal_Leptop_PROOF::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Anal_Leptop_PROOF_cxx
