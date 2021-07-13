#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>

#include <iostream>
#include <fstream>

using namespace std;

class AK4Jet {

 public:
  
  int  jetid;
  int  jetid_tightlepveto;
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  int    hadronFlavour;                                                                        
  int    partonFlavour;
  float  btag_DeepFlav;
  float  btag_DeepCSV;
  float  puid;
  float  qgl;
  
  /*
    float  reso;
    float  reso_up;
    float  reso_dn;
    float jesup_pu;
    float jesup_rel;
    float jesup_scale;
    float jesup_total;
    float jesdn_pu;
    float jesdn_rel;
    float jesdn_scale;
    float jesdn_total;
    int   genmatch;
  */
};

class AK8Jet {

 public:
  
  int  jetid;
  int  jetid_tightlepveto;
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  float pt_resoup;
  float mass_resoup;
  float pt_resodn;
  float mass_resodn;
  float jesup_total;
  float jesdn_total;
  float chrad;
  float tau21;
  float tau32;
  float DeepTag_TvsQCD;
  float DeepTag_WvsQCD;
  float DeepTag_ZvsQCD;
  float btag_DeepCSV;
  float CHF;
  float NHF;
  float CEMF;
  float NEMF;
  float MUF;
  float PHF;
  float HadF;
  float NHadF;
  float EmF;
  float EEM;
  float ncons;
  float Chcons;
  float neuemfrac;
  float neunhadfrac;
  float sdmass;
  float sub1pt;
  float sub1eta;
  float sub1phi;
  float sub1mass;
  float sub1btag;
  float sub1hadfrac;
  float sub1emfrac;
  float sub2pt;
  float sub2eta;
  float sub2phi;
  float sub2mass;
  float sub2btag;
  float sub2hadfrac;
  float sub2emfrac;
  float subbtag;
  float subhaddiff;
  float subemdiff;
  float subptdiff;
  float elinsubpt;
  float elinsubeta;
  float elinsubphi;
  float elinsubjpt;
  float elinsubjeta;
  float elinsubjphi;
  float elinsubjmass;
  float muinsubpt;
  float muinsubeta;
  float muinsubphi;
  float muinsubjpt;
  float muinsubjeta;
  float muinsubjphi;
  float muinsubjmass;
  float muinsubI0;
  float muinsubInear;
  float muinsubIfar;
  float matchAK4deepb;
  float re_tvsb;
  float rmu_tvsb;

  bool haselectron;
  bool hasmuon;
  bool hastau;
  bool hasqg; 
  bool hasb;
  bool hasleptop; 
  bool hashadtop;
  bool hastop;
  bool hasleptop_alldecay;
  bool hashadtop_alldecay;
  bool haspfelectron;
  bool haspfmuon;
  bool hasmatchmu;
  bool hasmatche;
};


class Muon {

 public:

  float pt;
  float eta;
  float phi;
  float charge;
  float trkvtx;
  float dz;
  bool ip;
  bool isTRK;
  bool isGL;
  bool isPF;
  bool isLoose;
  bool isGoodGL;
  bool isMed;
  bool isMedPr;
  bool isTight;
  bool isHighPt;
  bool isHighPttrk;
  float minisoall;
  float chi;
  float posmatch;
  float trkink;
  float segcom;
  float hit;
  float mst;
  float pixhit;
  float trklay;
  float valfrac;
  float pfiso;
  float p;
  float mudxy_sv;
  

};

class Electron {

 public:

  float pt;
  float eta;
  float phi;
  float charge;
  bool id;
  bool Fallv2WP80;
  bool id_noIso;
  bool Fallv2WP80_noIso;
  float p;
  float dxy;
  float dz;
  bool ip;
  float pfiso;
  float eldxy_sv;
  float supcl_eta;
  float supcl_phi;
  float supcl_rawE;
  float sigmaieta;
  float sigmaiphi;
  float r9full;
  float supcl_etaw;
  float supcl_phiw;
  float hcaloverecal;
  float cloctftrkn;
  float cloctftrkchi2;
  float e1x5bye5x5;
  float normchi2;
  float hitsmiss;
  float trkmeasure;
  float ecloverpout;
  float ecaletrkmomentum;
  float deltaetacltrkcalo;
  float supcl_preshvsrawe;
  float pfisolsumphet;
  float pfisolsumchhadpt;
  float pfsiolsumneuhadet;
  float etain;
  float phiin;
  float fbrem;
  float eoverp;
  float hovere;

};

class Lepton {
 
 public:

  float pt;
  float eta;
  float phi;
  float charge;
  int lepton_id;
};

class AK4GenJet {

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

};

class AK8GenJet {

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

} ;

class GenParton{

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;
  int genPartIdxMother;
  int pdgId;
  int status;
  int statusFlags;

  TLorentzVector p4;

} ;
bool AK4Jet_sort_by_pt(AK4Jet i1, AK4Jet i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<AK4Jet> & objs) {
  sort(objs.begin(), objs.end(), AK4Jet_sort_by_pt);
}
bool AK8Jet_sort_by_pt(AK8Jet i1, AK8Jet i2)                                                        
{                                                                                                   
  return (i1.pt > i2.pt);                                                                           
}                                                                                                   
void sorted_by_pt(vector<AK8Jet> & objs) {                                                          
  sort(objs.begin(), objs.end(), AK8Jet_sort_by_pt);                                                
}
bool Muon_sort_by_pt(Muon i1, Muon i2)                                                                           
{                                                                                                                
  return (i1.pt > i2.pt);                                                                                        
}                                                                                                                
void sorted_by_pt(vector<Muon> & objs) {                                                                         
  sort(objs.begin(), objs.end(), Muon_sort_by_pt);                                                               
}
bool Electron_sort_by_pt(Electron i1, Electron i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Electron> & objs) {
  sort(objs.begin(), objs.end(), Electron_sort_by_pt);
}
bool Lepton_sort_by_pt(Lepton i1, Lepton i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Lepton> & objs) {
  sort(objs.begin(), objs.end(), Lepton_sort_by_pt);
}
/*
  float compute_HT(vector<AK4Jet>  & objs, float ptcut, float etacut){
  
  float HT = 0;
  
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
  
    if(objs[iobs].pt > ptcut && abs(objs[iobs].eta)<=etacut){
    
    HT += objs[iobs].pt;
      
    }
    }
    
    return HT;
    }
    
      
    bool AK4Jet_sort_by_DeepFlav(AK4Jet i1, AK4Jet i2)
    {
    return (i1.btagDeepFlavB > i2.btagDeepFlavB);
    }
    void sorted_by_DeepFlav(vector<AK4Jet> & objs) {
    sort(objs.begin(), objs.end(), AK4Jet_sort_by_DeepFlav);
    }
    
    bool AK8Jet_sort_by_DeepAK8_Htag(AK8Jet i1, AK8Jet i2)
    {
    return (i1.deepTagMD_bbvsLight > i2.deepTagMD_bbvsLight);
    }
    void sorted_by_DeepAK8_Htag(vector<AK8Jet> & objs) {
    sort(objs.begin(), objs.end(), AK8Jet_sort_by_DeepAK8_Htag);
    }
    
    bool GenAK4Jet_sort_by_pt(AK4GenJet i1, AK4GenJet i2)
    {
    return (i1.pt > i2.pt);
    }
    void sorted_by_pt(vector<AK4GenJet> & objs) {
    sort(objs.begin(), objs.end(), GenAK4Jet_sort_by_pt);
    }
    
    bool GenAK8Jet_sort_by_pt(AK8GenJet i1, AK8GenJet i2)
    {
    return (i1.pt > i2.pt);
    }
    void sorted_by_pt(vector<AK8GenJet> & objs) {
    sort(objs.begin(), objs.end(), GenAK8Jet_sort_by_pt);
    }
*/


