#define Anal_Leptop_PROOF_cxx
#include "Anal_Leptop_PROOF.h"

void Anal_Leptop_PROOF::getmuons(std::vector<Muon> &vmuons, float ptcut=25, float etacut=2.5, int maxsize=njetmx, float dxy_cut=0.2, float dz_cut=0.5)
{
  
  for(int mu=0; mu<(nmuons); mu++){
    
    if(muonpt[mu]<ptcut) continue; 
    if(fabs(muoneta[mu])>etacut)  continue; 
    
    bool mu_id = Muon_TightID(muonisGL[mu],muonisPF[mu],
			      muonchi[mu],muonhit[mu],muonmst[mu],
			      muontrkvtx[mu],muondz[mu],muonpixhit[mu],muontrklay[mu]);
    if(!mu_id) continue;
    //bool mu_iso = Muon_Iso_ID(muonpfiso[mu]);
    //if(!mu_iso) continue;
    
    if(fabs(muontrkvtx[mu])>dxy_cut || fabs(muondz[mu])>dz_cut) continue;
    
    Muon vmuon;
    vmuon.pt = muonpt[mu];
    vmuon.eta = muoneta[mu];
    vmuon.phi = muonphi[mu];
    vmuon.mass = 0.105658;
    vmuon.charge = muoncharge[mu];
    vmuon.trkvtx = muontrkvtx[mu];
    vmuon.dz = muondz[mu];
    vmuon.isTRK = muonisTRK[mu];
    vmuon.isGL = muonisGL[mu];
    vmuon.isPF = muonisPF[mu];
    vmuon.isLoose = muonisLoose[mu];
    vmuon.isGoodGL = muonisGoodGL[mu];
    vmuon.isMed = muonisMed[mu];
    vmuon.isMedPr = muonisMedPr[mu];
    vmuon.isTight = muonisTight[mu];
    vmuon.isHighPt = muonisHighPt[mu];
    vmuon.isHighPttrk = muonisHighPttrk[mu];
    vmuon.minisoall = muonminisoall[mu];
    vmuon.chi = muonchi[mu];
    vmuon.posmatch = muonposmatch[mu];
    vmuon.trkink = muontrkink[mu];
    vmuon.segcom = muonsegcom[mu];
    vmuon.hit = muonhit[mu];
    vmuon.mst = muonmst[mu];
    vmuon.pixhit = muonpixhit[mu];
    vmuon.trklay = muontrklay[mu];
    vmuon.valfrac = muonvalfrac[mu];
    vmuon.pfiso = muonpfiso[mu];   
    vmuon.p = muonp[mu];
    vmuon.mudxy_sv = mudxy_sv[mu];
    
    vmuon.p4.SetPtEtaPhiM(vmuon.pt,vmuon.eta,vmuon.phi,vmuon.mass);
    
    vmuons.push_back(vmuon);
    
    if(int(vmuons.size())>=maxsize) break;
  }
  
  sorted_by_pt(vmuons);
	
}

void Anal_Leptop_PROOF::getelectrons(std::vector<Electron> &velectrons, float ptcut=25, float etacut=2.5, int maxsize=njetmx, float dxy_cut=0.05, float dz_cut=0.1)
{
  
 for(int ie=0; ie<(nelecs); ie++) {
   
   if (elpt[ie]<ptcut) continue;
   if(fabs(eleta[ie])>etacut)  continue; 
   
   // if(!elmvaid[ie]) continue;
   if(!elmvaid_noIso[ie]) continue;
    
   if(!((fabs(elsupcl_eta[ie])<1.4442 && fabs(eldxytrk[ie])<dxy_cut && fabs(eldztrk[ie])<dz_cut)||(fabs(elsupcl_eta[ie])>1.5660 && fabs(elsupcl_eta[ie])<2.5 && fabs(eldxytrk[ie])<(2*dxy_cut) && fabs(eldztrk[ie])<(2*dz_cut)))) continue;
   
   Electron velectron;
   
    velectron.pt = elpt[ie];
    velectron.eta = eleta[ie];
    velectron.phi = elphi[ie];
    velectron.mass = 0.000511;
    velectron.charge = elcharge[ie];
    velectron.id = elmvaid[ie];
    velectron.Fallv2WP80 = elmvaid_Fallv2WP80[ie];
    velectron.id_noIso = elmvaid_noIso[ie];
    velectron.Fallv2WP80_noIso = elmvaid_Fallv2WP80_noIso[ie];
    velectron.p = elp[ie];
    velectron.dxy = eldxytrk[ie];
    velectron.dz = eldztrk[ie];
    velectron.pfiso = elpfiso[ie];
    velectron.eldxy_sv = eldxy_sv[ie];
    velectron.supcl_eta = elsupcl_eta[ie];
    velectron.supcl_phi = elsupcl_phi[ie];
    velectron.supcl_rawE = elsupcl_rawE[ie];
    velectron.sigmaieta = elsigmaieta[ie];
    velectron.sigmaiphi = elsigmaiphi[ie];
    velectron.r9full = elr9full[ie];
    velectron.supcl_etaw = elsupcl_etaw[ie];
    velectron.supcl_phiw = elsupcl_phiw[ie];
    velectron.hcaloverecal = elhcaloverecal[ie];
    velectron.cloctftrkn = elcloctftrkn[ie];
    velectron.cloctftrkchi2 = elcloctftrkchi2[ie];
    velectron.e1x5bye5x5 = ele1x5bye5x5[ie];
    velectron.normchi2 = elnormchi2[ie];
    velectron.hitsmiss = elhitsmiss[ie];
    velectron.trkmeasure = eltrkmeasure[ie];
    velectron.ecloverpout = elecloverpout[ie];
    velectron.ecaletrkmomentum = elecaletrkmomentum[ie];
    velectron.deltaetacltrkcalo = eldeltaetacltrkcalo[ie];
    velectron.supcl_preshvsrawe = elsupcl_preshvsrawe[ie];
    velectron.pfisolsumphet = elpfisolsumphet[ie];
    velectron.pfisolsumchhadpt = elpfisolsumchhadpt[ie];
    velectron.pfsiolsumneuhadet = elpfsiolsumneuhadet[ie];
    velectron.etain = eletain[ie];
    velectron.phiin = elphiin[ie];
    velectron.fbrem = elfbrem[ie];
    velectron.eoverp = eleoverp[ie];
    velectron.hovere = elhovere[ie];
    
    velectron.p4.SetPtEtaPhiM(velectron.pt,velectron.eta,velectron.phi,velectron.mass);
    
    velectrons.push_back(velectron);
    
    if(int(velectrons.size()) >= maxsize) break;
    
 }
 
 sorted_by_pt(velectrons);	
	
}

void Anal_Leptop_PROOF::getLeptons(std::vector<Lepton> &vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, float pt_cut=30)
{ 
  for(unsigned imu=0; imu<vmuons.size(); imu++){
    if(vmuons[imu].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = vmuons[imu].pt;
    vlepton.eta = vmuons[imu].eta;
    vlepton.phi = vmuons[imu].phi;
    vlepton.mass = vmuons[imu].mass;
    vlepton.charge = vmuons[imu].charge;
    vlepton.lepton_id = 1;
    vlepton.pdgId = 13;
    vlepton.p4 = vmuons[imu].p4;
    vlepton.indexemu = imu; 
    vleptons.push_back(vlepton);
  }
  for(unsigned ie=0; ie<velectrons.size(); ie++){
    if(velectrons[ie].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = velectrons[ie].pt;
    vlepton.eta = velectrons[ie].eta;
    vlepton.phi = velectrons[ie].phi;
    vlepton.mass = velectrons[ie].mass;
    vlepton.charge = velectrons[ie].charge;
    vlepton.lepton_id = 2;
    vlepton.pdgId = 11;
    vlepton.p4 = velectrons[ie].p4;
    vlepton.indexemu=ie;
    vleptons.push_back(vlepton);
  }
  sorted_by_pt(vleptons);
}

void Anal_Leptop_PROOF::getAK4jets(std::vector<AK4Jet> &Jets, float ptcut=30, float etacut=2.5, bool isMC=false, int maxsize=njetmx)
{
  
  for(int ijet=0; ijet<(npfjetAK4); ijet++){
    
    AK4Jet sJet; 
    
    if(!pfjetAK4jetID[ijet]) continue;
    
    pfjetAK4pt[ijet] *= pfjetAK4JEC[ijet] ;
    pfjetAK4mass[ijet] *= pfjetAK4JEC[ijet];
    
    if(isMC){
      pfjetAK4pt[ijet] *= (1+pfjetAK4reso[ijet]) ;
      pfjetAK4mass[ijet] *= (1+pfjetAK4reso[ijet]) ;
    }
    
    if(fabs(pfjetAK4eta[ijet])>etacut) continue;
    if(pfjetAK4pt[ijet]<ptcut) continue;
    
    sJet.pt = pfjetAK4pt[ijet];
    sJet.mass = pfjetAK4mass[ijet];
    sJet.eta = pfjetAK4eta[ijet];
    sJet.y = pfjetAK4y[ijet];
    sJet.phi = pfjetAK4phi[ijet];
    sJet.p4.SetPtEtaPhiM(sJet.pt,sJet.eta,sJet.phi,sJet.mass);
  
    sJet.jetid = pfjetAK4jetID[ijet];
    sJet.jetid_tightlepveto = pfjetAK4jetID_tightlepveto[ijet];
    sJet.hadronFlavour = pfjetAK4hadronflav[ijet];
    sJet.partonFlavour = pfjetAK4partonflav[ijet];
    sJet.btag_DeepFlav = pfjetAK4btag_DeepFlav[ijet];
    sJet.btag_DeepCSV = pfjetAK4btag_DeepCSV[ijet];
    sJet.puid = pfjetAK4PUID[ijet];
    sJet.qgl = pfjetAK4qgl[ijet];
		
    Jets.push_back(sJet);
    
    if(int(Jets.size())>=maxsize) break;
    
  }

  sorted_by_pt(Jets);
			
}

void Anal_Leptop_PROOF::getAK8jets(std::vector<AK8Jet> &LJets, float ptcut=200, float etacut=2.5, bool isMC=false, int maxsize=njetmxAK8)
{

  for(int ijet=0; ijet<(npfjetAK8); ijet++){
  
    if(!pfjetAK8jetID[ijet]) continue;
    
    pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet] ;
    pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
    
    if(isMC){
      
      pfjetAK8pt[ijet] *= (1+pfjetAK8reso[ijet]) ;
      pfjetAK8mass[ijet] *= (1+pfjetAK8reso[ijet]) ;
      
      pfjetAK8pt_resoup[ijet] = pfjetAK8pt[ijet]*(1+pfjetAK8resoup[ijet]);
      pfjetAK8mass_resoup[ijet] = pfjetAK8mass[ijet]*(1+pfjetAK8resoup[ijet]);
      pfjetAK8pt_resodown[ijet] = pfjetAK8pt[ijet]*(1+pfjetAK8resodn[ijet]);
      pfjetAK8mass_resodown[ijet] = pfjetAK8mass[ijet]*(1+pfjetAK8resodn[ijet]);
    }
    
    if(fabs(pfjetAK8eta[ijet])>etacut) continue;
    if(pfjetAK8pt[ijet] < ptcut) continue;
    
    AK8Jet LJet;
    
    LJet.jetid = pfjetAK8jetID[ijet];
    LJet.jetid_tightlepveto = pfjetAK8jetID_tightlepveto[ijet];
    LJet.pt = pfjetAK8pt[ijet];
    LJet.eta = pfjetAK8eta[ijet];
    LJet.mass = pfjetAK8mass[ijet];
    LJet.phi = pfjetAK8phi[ijet];
    LJet.y = pfjetAK8y[ijet];
    LJet.p4.SetPtEtaPhiM(LJet.pt,LJet.eta,LJet.phi,LJet.mass);
    LJet.pt_resoup = pfjetAK8pt_resoup[ijet];
    LJet.mass_resoup = pfjetAK8mass_resoup[ijet];
    LJet.pt_resodn = pfjetAK8pt_resodown[ijet];
    LJet.mass_resodn = pfjetAK8mass_resodown[ijet];
    LJet.jesup_total = pfjetAK8jesup_total[ijet];
    LJet.jesdn_total = pfjetAK8jesdn_total[ijet];
    
    LJet.chrad = pfjetAK8chrad[ijet];
    LJet.tau21 = pfjetAK8tau2[ijet]*1./max(float(1.e-6),pfjetAK8tau1[ijet]);
    LJet.tau32 = pfjetAK8tau3[ijet]*1./max(float(1.e-6),pfjetAK8tau2[ijet]);
    LJet.DeepTag_TvsQCD = pfjetAK8DeepTag_TvsQCD[ijet];
    LJet.DeepTag_WvsQCD = pfjetAK8DeepTag_WvsQCD[ijet];
    LJet.DeepTag_ZvsQCD = pfjetAK8DeepTag_ZvsQCD[ijet];
    LJet.btag_DeepCSV = pfjetAK8btag_DeepCSV[ijet];
    
    LJet.CHF = pfjetAK8CHF[ijet];
    LJet.NHF = pfjetAK8NHF[ijet];
    LJet.CEMF = pfjetAK8CEMF[ijet];
    LJet.NEMF = pfjetAK8NEMF[ijet];
    LJet.MUF = pfjetAK8MUF[ijet];
    LJet.PHF = pfjetAK8PHF[ijet];
    pfjetAK8HadF[ijet] = (pfjetAK8NHF[ijet]+pfjetAK8CHF[ijet]);
    LJet.NHadF = (1.- pfjetAK8HadF[ijet]);
    
    LJet.EMF = (pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]);
    LJet.EEM = pfjetAK8EEM[ijet];
    
    LJet.neucons = pfjetAK8Neucons[ijet];
    LJet.chcons = pfjetAK8Chcons[ijet];
    
    LJet.neuemfrac = (pfjetAK8NEMF[ijet]*1.)/max(float(1.e-6),(pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]));
    LJet.neunhadfrac = (pfjetAK8NEMF[ijet]*1.)/max(1.e-6,(1.- pfjetAK8HadF[ijet]));
    
    LJet.sdmass = pfjetAK8sdmass[ijet];
    LJet.sub1pt = pfjetAK8sub1pt[ijet];
    LJet.sub1eta = pfjetAK8sub1eta[ijet];
    LJet.sub1phi = pfjetAK8sub1phi[ijet];
    LJet.sub1mass = pfjetAK8sub1mass[ijet];
    LJet.sub1btag = pfjetAK8sub1btag[ijet];
    LJet.sub1hadfrac = pfjetAK8sub1chhadfrac[ijet]+pfjetAK8sub1neuhadfrac[ijet];
    LJet.sub1emfrac = pfjetAK8sub1emfrac[ijet];
    
    LJet.sub2pt = pfjetAK8sub2pt[ijet];
    LJet.sub2eta = pfjetAK8sub2eta[ijet];
    LJet.sub2phi = pfjetAK8sub2phi[ijet];
    LJet.sub2mass = pfjetAK8sub2mass[ijet];
    LJet.sub2btag = pfjetAK8sub2btag[ijet];
    LJet.sub2hadfrac = pfjetAK8sub2chhadfrac[ijet]+pfjetAK8sub2neuhadfrac[ijet];
    LJet.sub2emfrac = pfjetAK8sub2emfrac[ijet];
    
    LJet.subbtag = max(pfjetAK8sub1btag[ijet],pfjetAK8sub2btag[ijet]);
    
    LJet.subhaddiff = pfjetAK8subhaddiff[ijet];
    LJet.subemdiff = pfjetAK8subemdiff[ijet];
    LJet.subptdiff = pfjetAK8subptdiff[ijet];

    LJet.elinsubpt = pfjetAK8elinsubpt[ijet];
    LJet.elinsubeta = pfjetAK8elinsubeta[ijet];
    LJet.elinsubphi = pfjetAK8elinsubphi[ijet];
        
    LJet.elinsubjpt = pfjetAK8elinsubjpt[ijet];
    LJet.elinsubjeta = pfjetAK8elinsubjeta[ijet];
    LJet.elinsubjphi = pfjetAK8elinsubjphi[ijet];
    LJet.elinsubjmass = pfjetAK8elinsubjmass[ijet];
    
    LJet.muinsubpt = pfjetAK8muinsubpt[ijet];
    LJet.muinsubeta = pfjetAK8muinsubeta[ijet];
    LJet.muinsubphi = pfjetAK8muinsubphi[ijet];
    LJet.muinsubjpt = pfjetAK8muinsubjpt[ijet];
    LJet.muinsubjeta = pfjetAK8muinsubjeta[ijet];
    LJet.muinsubjphi = pfjetAK8muinsubjphi[ijet];
    LJet.muinsubjmass = pfjetAK8muinsubjmass[ijet];
    
    LJet.muinsubI0 = pfjetAK8muinsubI0[ijet];
    LJet.muinsubInear = pfjetAK8muinsubInear[ijet];
    LJet.muinsubIfar = pfjetAK8muinsubIfar[ijet];
    
    LJet.haselectron = LJet.hasmuon = LJet.hastau = LJet.hasqg = LJet.hasb = LJet.hasleptop = LJet.hashadtop = LJet.hastop = LJet.hasmatchmu = LJet.hasmatche = false;
    LJet.hasleptop_alldecay = LJet.hashadtop_alldecay = false;
    
    LJet.haspfelectron = LJet.haspfmuon = false;
    LJet.matchAK4deepb = LJet.re_tvsb = LJet.rmu_tvsb = -99;
    
    if(LJet.elinsubpt > 0. && LJet.elinsubjpt > 0.){
      if (delta2R(LJet.elinsubeta,LJet.elinsubphi,LJet.y,LJet.phi) < 0.8) {
	LJet.haspfelectron = true;
      }
    }
    
    if (LJet.muinsubpt > 0. && LJet.muinsubjpt > 0.) {
      if (delta2R(LJet.muinsubeta,LJet.muinsubphi,LJet.y,LJet.phi) < 0.8) {
	LJet.haspfmuon = true;
      }
    }
    
    LJets.push_back(LJet);
    
    if(int(LJets.size())>=maxsize) break;
  }
  
  sorted_by_pt(LJets);
  
}

void Anal_Leptop_PROOF::LeptonJet_cleaning(std::vector<AK4Jet> &Jets, std::vector<Lepton> Leptons, float dR_cut=0.4, float ptcut=30, float etacut=2.5)
{
  
   auto jet = Jets.begin();
   while (jet != Jets.end()){
	
	for(auto & lepton: Leptons){
		
		if(delta2R(jet->eta, jet->phi, lepton.eta, lepton.phi)<dR_cut){
			jet->p4 -= lepton.p4;
			jet->pt = (jet->p4).Pt();
			jet->eta = (jet->p4).Eta();
			jet->y = (jet->p4).Rapidity();
			jet->phi = (jet->p4).Phi();
			jet->mass = (jet->p4).M();
		}
	}
	
	if(jet->pt<ptcut || fabs(jet->eta)>etacut){
		Jets.erase(jet);
		}
	else{
		++jet;
		}
  }	

  sorted_by_pt(Jets);

}

void Anal_Leptop_PROOF::getPartons(std::vector<GenParton> &GenPartons, int maxsize=npartmx)
{
  
  for(int igen=0; igen<(ngenparticles); igen++){
	 
	GenParton parton; 
	
	parton.pt = genpartpt[igen];
	parton.eta = genparteta[igen];
	parton.phi = genpartphi[igen];
	parton.mass = genpartm[igen];  
	parton.p4.SetPtEtaPhiM(parton.pt,parton.eta,parton.phi,parton.mass);
    
    parton.status = genpartstatus[igen];
	parton.pdgId = genpartpdg[igen];
	parton.mompdgId = genpartmompdg[igen];
	parton.grmompdgId = genpartgrmompdg[igen]; 
	
	parton.fromhard = genpartfromhard[igen];
	parton.fromhardbFSR = genpartfromhardbFSR[igen];
	parton.isPromptFinalState = genpartisPromptFinalState[igen];
	parton.isLastCopyBeforeFSR = genpartisLastCopyBeforeFSR[igen]; 
		
    GenPartons.push_back(parton);
    
    if(int(GenPartons.size())>=maxsize) break;
    
  }

//	GenPartons.sorted_by_pt();	
	
}

void Anal_Leptop_PROOF::getLHETops(std::vector<GenParton> &LHETops, std::vector<GenParton> GenPartons)
{
  
  for(auto & part: GenPartons){
	  
	  if(abs(part.status)!=22) continue;
      if(!(part.fromhard)) continue;
      if(abs(part.pdgId)!=6) continue;
	  
	  LHETops.push_back(part);
	  
	  }
}

void Anal_Leptop_PROOF::getGENTops(vector<TopQuark> &gentops, vector<GenParton> genpartons)  // with daughters after shower
{     
    vector<GenParton> W_dau;
    vector<GenParton> t_bp;
  
    for(unsigned igen=0; igen<genpartons.size(); igen++){
      
      if(!(genpartons[igen].status==23 || genpartons[igen].status==1)) continue;
      if(!(genpartons[igen].fromhard)) continue;
      
      if(!((abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<=5)||(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16))) continue;
      if(!(abs(genpartons[igen].mompdgId)==6 || abs(genpartons[igen].mompdgId)==24)) continue;
      if(abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)!=6) continue;
            
      if(abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<=5 && abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==6)   {  W_dau.push_back(genpartons[igen]); }
      if(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16 && abs(genpartons[igen].mompdgId)==24 && abs(genpartons[igen].grmompdgId)==6) {  W_dau.push_back(genpartons[igen]); }
      if(abs(genpartons[igen].pdgId)==5 && abs(genpartons[igen].mompdgId)==6) {  t_bp.push_back(genpartons[igen]); }
    }
    
    for(unsigned ipart=0; ipart<W_dau.size(); ipart++){
      
      if(!((abs(W_dau[ipart].pdgId)==11 || abs(W_dau[ipart].pdgId)==13 || abs(W_dau[ipart].pdgId)==15 || abs(W_dau[ipart].pdgId)==1 || abs(W_dau[ipart].pdgId)==3))) continue;
      
      unsigned partner = -1;
      unsigned match_b = -1;
      
      for(unsigned jpart=(ipart+1); jpart<W_dau.size(); jpart++){
	if((W_dau[ipart].mompdgId==W_dau[jpart].mompdgId) && (W_dau[ipart].grmompdgId==W_dau[jpart].grmompdgId) && (W_dau[ipart].pdgId*W_dau[jpart].pdgId<0)){
	  partner = jpart;
	  break;
	}
      }
      
      for(unsigned ib=0; ib<t_bp.size(); ib++){
	if(t_bp[ib].mompdgId==W_dau[ipart].grmompdgId){
	  match_b = ib;
	  break;
	}
      }
      
      GenParton q1, q2, b;
      
      q1 = W_dau[ipart];
      
      if(int(partner)>=0 && partner<W_dau.size()){
	q2 = W_dau[partner];
      }
      if(int(match_b)>=0 && match_b<t_bp.size()){
	b = t_bp[match_b];
      }
      
      if(int(partner)>=0 && partner<W_dau.size() && int(match_b)>=0 && match_b<t_bp.size()){
	
	TopQuark topQ;
	
	topQ.p4 = (b.p4 + q1.p4 + q2.p4);
	topQ.daughter.push_back(q1);
	topQ.daughter.push_back(q2);
	topQ.daughter.push_back(b);
	
	gentops.push_back(topQ);	
      }
    }
}



void Anal_Leptop_PROOF::TopAssignment_toJet(std::vector<AK8Jet> &LJets, std::vector<GenParton> lhetops, std::vector<TopQuark> gentops)
{
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
    for(unsigned itop=0; itop<lhetops.size(); itop++)
      {
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,lhetops[itop].eta,lhetops[itop].phi)<0.6){
	  LJets[ijet].hastop = true;
	  break;
	}
      }
    
    for(unsigned itop=0; itop<gentops.size(); itop++){
      
      if(abs(gentops[itop].daughter[0].pdgId)==11 || abs(gentops[itop].daughter[0].pdgId)==13 || abs(gentops[itop].daughter[0].pdgId)==15){
	
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].p4.Rapidity(),gentops[itop].p4.Phi())<0.8){
	  LJets[ijet].hasleptop = true;	
	}
	
	if((delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[0].eta,gentops[itop].daughter[0].phi)<0.8)
	   && (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[1].eta,gentops[itop].daughter[1].phi)<0.8)
	   && (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[2].eta,gentops[itop].daughter[2].phi)<0.8))
	  {
	    LJets[ijet].hasleptop_alldecay = true;
	  }
	
      }
      
      else if(abs(gentops[itop].daughter[0].pdgId)==1 || abs(gentops[itop].daughter[0].pdgId)==3){
	
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].p4.Rapidity(),gentops[itop].p4.Phi())<0.8){
	  LJets[ijet].hashadtop = true;	
	}
	
	if((delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[0].eta,gentops[itop].daughter[0].phi)<0.8)
	   && (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[1].eta,gentops[itop].daughter[1].phi)<0.8)
	   && (delta2R(LJets[ijet].y,LJets[ijet].phi,gentops[itop].daughter[2].eta,gentops[itop].daughter[2].phi)<0.8))
	  {
	    LJets[ijet].hashadtop_alldecay = true;
	  }
	
      }
    }
      
  }//ijet
	
}

void Anal_Leptop_PROOF::AssignGen(std::vector<AK8Jet> &LJets, std::vector<GenParton> GenPartons){
  
  for(auto & ljet: LJets){
    
    for(auto & part: GenPartons){
      
      if(abs(part.status)!=23 && part.status!=1) continue;
      if(!(part.fromhard)) continue;
      
      if(abs(part.pdgId)==11 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.7){
	ljet.haselectron  = true;
	break;
      }
      
      if(abs(part.pdgId)==13 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.7){
	ljet.hasmuon  = true;
	break;
      }
      
      if(abs(part.pdgId)==15 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.7){
	ljet.hastau  = true;
	break;
      }
      
      if(((abs(part.pdgId)>=1 && abs(part.pdgId)<5) || abs(part.pdgId)==21) && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.7){
	ljet.hasqg  = true;
	break;
      }
      
      if(abs(part.pdgId)==5 && delta2R(ljet.y,ljet.phi,part.eta,part.phi)<0.8){
	ljet.hasb  = true;
	break;
      } 
    }
  }
}

bool Anal_Leptop_PROOF::isBJet(AK4Jet jet, float btag_cut)
{
if(jet.btag_DeepFlav >= btag_cut){
	return true;
}
else{
	return false;
	}	
}


void Anal_Leptop_PROOF::ReadTagger(std::vector<AK8Jet> &LJets, std::vector<Lepton> vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, TMVA::Reader *reader_electron, TMVA::Reader *reader_muon){
  
  //for (int ij=0; ij<min(int(LJets.size()),2); ij++) {
  for (int ij=0; ij<int(LJets.size()); ij++) {
    //Muon
    in_mupfjetAK8NHadF = -999;
    in_mupfjetAK8neunhadfrac = -999;
    in_mupfjetAK8subhaddiff = -999;
    in_mupfjetAK8tau21 = -999;
    in_mupfjetAK8chrad = -999;
    in_mupfjetAK8sdmass = -999;
    in_pfjetAK8matchedmuonchi = -999;
    in_pfjetAK8matchedmuonposmatch = -999;
    in_pfjetAK8matchedmuontrkink = -999;
    in_pfjetAK8matchedmuonsegcom = -999;
    in_pfjetAK8matchedmuonhit = -999;
    in_pfjetAK8matchedmuonmst = -999;
    in_pfjetAK8matchedmuontrkvtx = -999;
    in_pfjetAK8matchedmuondz = -999;
    in_pfjetAK8matchedmuonpixhit = -999;
    in_pfjetAK8matchedmuontrklay = -999;
    in_pfjetAK8matchedmuonvalfrac = -999;
    in_pfjetAK8muinsubptrat = -999;
    in_pfjetAK8muinsubmassrat = -999;
    in_pfjetAK8muinsubinvmass = -999;
    in_pfjetAK8muinsubIfarbyI0 = -999;
    in_pfjetAK8muinsubInearbyI0 = -999;
    //electron
    in_pfjetAK8NHadF = -999;
    in_pfjetAK8neunhadfrac = -999;
    in_pfjetAK8subhaddiff = -999;
    in_pfjetAK8tau21 = -999;
    in_pfjetAK8chrad = -999;
    in_pfjetAK8sdmass = -999;
    in_pfjetAK8matchedelcleta = -999;
    in_pfjetAK8matchedelpt = -999;
    in_pfjetAK8matchedelsigmaieta = -999;
    in_pfjetAK8matchedelsigmaiphi = -999;
    in_pfjetAK8matchedelr9full = -999;
    in_pfjetAK8matchedelsupcl_etaw = -999;
    in_pfjetAK8matchedelsupcl_phiw = -999;
    in_pfjetAK8matchedelhcaloverecal = -999;
    in_pfjetAK8matchedelcloctftrkn = -999;
    in_pfjetAK8matchedelcloctftrkchi2 = -999;
    in_pfjetAK8matchedele1x5bye5x5 = -999;
    in_pfjetAK8matchedelnormchi2 = -999;
    in_pfjetAK8matchedelhitsmiss = -999;
    in_pfjetAK8matchedeltrkmeasure = -999;
    in_pfjetAK8matchedelecloverpout = -999;
    in_pfjetAK8matchedelecaletrkmomentum = -999;
    in_pfjetAK8matchedeldeltaetacltrkcalo = -999;
    in_pfjetAK8matchedelsupcl_preshvsrawe = -999;
    in_pfjetAK8matchedelpfisolsumphet = -999;
    in_pfjetAK8matchedelpfisolsumchhadpt = -999;
    in_pfjetAK8matchedelpfisolsumneuhadet = -999;
    in_pfjetAK8matchedeletain = -999;
    in_pfjetAK8matchedelphiin = -999;
    in_pfjetAK8matchedelfbrem = -999;
    in_pfjetAK8matchedeleoverp = -999;
    in_pfjetAK8matchedelhovere = -999;
    in_pfjetAK8matchedelRho = (Rho);
    
    LJets[ij].hasmatche = false;
    LJets[ij].hasmatchmu = false;
    
    if (int(vleptons.size()) >ij) {
      
      if (vleptons[ij].pdgId==11) { //Electron
	
	in_pfjetAK8NHadF = LJets[ij].NHadF;
	in_pfjetAK8neunhadfrac = LJets[ij].neunhadfrac;
	in_pfjetAK8subhaddiff = LJets[ij].subhaddiff;
	in_pfjetAK8tau21 = LJets[ij].tau21;
	in_pfjetAK8chrad = LJets[ij].chrad;
	in_pfjetAK8sdmass = LJets[ij].sdmass;
	
	LJets[ij].hasmatche = true;
	int iel = vleptons[ij].indexemu;
	
	in_pfjetAK8matchedelcleta = velectrons[iel].supcl_eta; //elsupcl_eta[nearest];
	in_pfjetAK8matchedelpt = fabs(velectrons[iel].pt); //elpt[nearest]);
	in_pfjetAK8matchedelsigmaieta = velectrons[iel].eta; // elsigmaieta[nearest];
	in_pfjetAK8matchedelsigmaiphi = velectrons[iel].phi; // elsigmaiphi[nearest];
	in_pfjetAK8matchedelr9full = velectrons[iel].r9full; // elr9full[nearest];
	in_pfjetAK8matchedelsupcl_etaw = velectrons[iel].supcl_etaw; // elsupcl_etaw[nearest];
	in_pfjetAK8matchedelsupcl_phiw = velectrons[iel].supcl_phiw; // elsupcl_phiw[nearest];
	in_pfjetAK8matchedelhcaloverecal = velectrons[iel].hcaloverecal; // elhcaloverecal[nearest];
	in_pfjetAK8matchedelcloctftrkn = velectrons[iel].cloctftrkn; // elcloctftrkn[nearest];
	in_pfjetAK8matchedelcloctftrkchi2 = velectrons[iel].cloctftrkchi2; // elcloctftrkchi2[nearest];
	in_pfjetAK8matchedele1x5bye5x5 = velectrons[iel].e1x5bye5x5; // ele1x5bye5x5[nearest];
	in_pfjetAK8matchedelnormchi2 = velectrons[iel].normchi2; // elnormchi2[nearest];
	in_pfjetAK8matchedelhitsmiss = velectrons[iel].hitsmiss; // elhitsmiss[iel];
	in_pfjetAK8matchedeltrkmeasure = velectrons[iel].trkmeasure; // eltrkmeasure[nearest];
	in_pfjetAK8matchedelecloverpout = velectrons[iel].ecloverpout; // elecloverpout[nearest];
	in_pfjetAK8matchedelecaletrkmomentum = velectrons[iel].ecaletrkmomentum; //elecaletrkmomentum[nearest];
	in_pfjetAK8matchedeldeltaetacltrkcalo = velectrons[iel].deltaetacltrkcalo; // eldeltaetacltrkcalo[nearest];
	in_pfjetAK8matchedelsupcl_preshvsrawe = velectrons[iel].supcl_preshvsrawe; // elsupcl_preshvsrawe[nearest];
	in_pfjetAK8matchedelpfisolsumphet = velectrons[iel].pfisolsumphet; // elpfisolsumphet[nearest];
	in_pfjetAK8matchedelpfisolsumchhadpt = velectrons[iel].pfisolsumchhadpt; // elpfisolsumchhadpt[nearest];
	in_pfjetAK8matchedelpfisolsumneuhadet = velectrons[iel].pfsiolsumneuhadet; //elpfsiolsumneuhadet[nearest];
	in_pfjetAK8matchedeletain = velectrons[iel].etain; // eletain[nearest];
	in_pfjetAK8matchedelphiin = velectrons[iel].phiin; // elphiin[nearest];
	in_pfjetAK8matchedelfbrem = velectrons[iel].fbrem; //elfbrem[nearest];
	in_pfjetAK8matchedeleoverp = velectrons[iel].eoverp; // eleoverp[nearest];
	in_pfjetAK8matchedelhovere = velectrons[iel].hovere; //elhovere[nearest];
	
	LJets[ij].re_tvsb = reader_electron->EvaluateMVA("BDTG method");
	
      } else { // Muon
	
	in_mupfjetAK8NHadF = LJets[ij].NHadF;
	in_mupfjetAK8neunhadfrac = LJets[ij].neunhadfrac;
	in_mupfjetAK8subhaddiff = LJets[ij].subhaddiff;
	in_mupfjetAK8tau21 = LJets[ij].tau21;
	in_mupfjetAK8chrad = LJets[ij].chrad;
	in_mupfjetAK8sdmass = LJets[ij].sdmass;
	
	LJets[ij].hasmatchmu = true;
	int imu = vleptons[ij].indexemu;
	in_pfjetAK8matchedmuonchi = vmuons[imu].chi; //GMA //muonchi[nearestmu];
	in_pfjetAK8matchedmuonposmatch = vmuons[imu].posmatch; //muonposmatch[nearestmu];
	in_pfjetAK8matchedmuontrkink = vmuons[imu].trkink; // muontrkink[nearestmu];
	in_pfjetAK8matchedmuonsegcom = vmuons[imu].segcom; // muonsegcom[nearestmu];
	in_pfjetAK8matchedmuonhit = vmuons[imu].hit; // muonhit[nearestmu];
	in_pfjetAK8matchedmuonmst = vmuons[imu].mst; //muonmst[nearestmu];
	in_pfjetAK8matchedmuontrkvtx = vmuons[imu].trkvtx; //muontrkvtx[nearestmu];
	in_pfjetAK8matchedmuondz = vmuons[imu].dz; //muondz[nearestmu];
	in_pfjetAK8matchedmuonpixhit = vmuons[imu].pixhit; // muonpixhit[nearestmu];
	in_pfjetAK8matchedmuontrklay = vmuons[imu].trklay; // muontrklay[nearestmu];
	in_pfjetAK8matchedmuonvalfrac = vmuons[imu].valfrac; // muonvalfrac[nearestmu];
	TLorentzVector lep, lepbj, bj;
	lep.SetPtEtaPhiE(LJets[ij].muinsubpt, LJets[ij].muinsubeta, LJets[ij].muinsubphi, 0.105658);
	bj.SetPtEtaPhiM(LJets[ij].muinsubjpt, LJets[ij].muinsubjeta, LJets[ij].muinsubjphi, LJets[ij].muinsubjmass);
	lepbj = lep + bj;
	
	in_pfjetAK8muinsubptrat = (lep.Pt())/(max(1.e-6,lepbj.Pt()));
	in_pfjetAK8muinsubmassrat = (bj.M())/(max(1.e-6,lepbj.M()));
	in_pfjetAK8muinsubinvmass = lepbj.M();
	in_pfjetAK8muinsubIfarbyI0 = LJets[ij].muinsubIfar/max(float(1.e-6),LJets[ij].muinsubI0);
	in_pfjetAK8muinsubInearbyI0 = LJets[ij].muinsubInear/max(float(1.e-6),LJets[ij].muinsubI0);
	
	LJets[ij].rmu_tvsb = reader_muon->EvaluateMVA("BDTG method");
      } //else of if (vleptons[ij].lepton_id==1)
    } // if (vleptons.size() >=ij)
  }//for (int ij=0; ij<min(LJets.size(),2); ij++)	
}

void Anal_Leptop_PROOF::Match_trigger(vector<bool> double_hlts, vector<bool> single_hlts,vector<vector<float>> double_pt_cuts, vector<float> single_pt_cuts, vector<vector<int>> double_pids, vector<int> single_pids, vector<float> single_other_pt_cuts, vector<int> single_other_pids, vector<std::pair<int,TLorentzVector>> TrigRefObj,Lepton lepcand_1, Lepton lepcand_2, vector<AK4Jet> Jets, bool &trig_threshold_pass, bool &trig_matching_pass, vector<TH1D*> &hist_init,vector<TH2D*> &hist2d)
{
  
  if(double_hlts.size()<1 && single_hlts.size()<1) { 
    return; 
  }
  
  else{
    
    // checking if offline objects passed trigger thresholds //
    
    bool double_trig_pass(false), single_trig_pass(false);
    bool any_double_hlt_pass = false;
    
    if(double_hlts.size()>0){
      for(unsigned ihlt=0; ihlt<double_hlts.size(); ihlt++){
	if (double_hlts[ihlt]){ any_double_hlt_pass = true; }
	if (double_hlts[ihlt] && ((lepcand_1.pt>double_pt_cuts[ihlt][0] && lepcand_1.pdgId==double_pids[ihlt][0] && lepcand_2.pt>double_pt_cuts[ihlt][1] && lepcand_2.pdgId==double_pids[ihlt][1]) || (lepcand_1.pt>double_pt_cuts[ihlt][1] && lepcand_1.pdgId==double_pids[ihlt][1] && lepcand_2.pt>double_pt_cuts[ihlt][0] && lepcand_2.pdgId==double_pids[ihlt][0]))) { 
	  trig_threshold_pass = true; 
	  double_trig_pass = true;
	  break; 
	}
      }
    }
    
    int fired_single_trig = -1;
    
    if(!any_double_hlt_pass){
      if(single_hlts.size()>0){
	for(unsigned ihlt=0; ihlt<single_hlts.size(); ihlt++){
	  if (single_hlts[ihlt] && ((lepcand_1.pt>single_pt_cuts[ihlt] && lepcand_1.pdgId==single_pids[ihlt]) || (lepcand_2.pt>single_pt_cuts[ihlt] && lepcand_2.pdgId==single_pids[ihlt])) && single_other_pt_cuts[ihlt]>0 && Jets.size()>0 && Jets[0].pt>single_other_pt_cuts[ihlt])  { 
	    trig_threshold_pass = true; 
	    single_trig_pass = true;
	    fired_single_trig = ihlt;
	    break; 
	  }
	  else if (single_hlts[ihlt] && ((lepcand_1.pt>single_pt_cuts[ihlt] && lepcand_1.pdgId==single_pids[ihlt]) || (lepcand_2.pt>single_pt_cuts[ihlt] && lepcand_2.pdgId==single_pids[ihlt])) && single_other_pt_cuts[ihlt]<0)
	    {
	      trig_threshold_pass = true; 
	      single_trig_pass = true;
	      fired_single_trig = ihlt;
	      break; 
	    }
	}
      }
    }
    
    // check if offline objects match to trigger objects //
    
    bool lep1_match = false;
    bool lep2_match = false;
    bool jet_match = false;
    
    float ptratmin_lep1(0.3), ptratmin_lep2(0.3), ptratmin_jet(0.3);
    
    for (uint trv=0; trv<TrigRefObj.size(); trv++) {
      bool mutrobj(false), eltrobj(false), jettrobj(false);
      if (abs(TrigRefObj[trv].first)==13) mutrobj=true;
      else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() < 1.e-3) eltrobj=true;
      else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.e-3) jettrobj=true;
      //else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.) jettrobj=true;
      //else if (abs(TrigRefObj[trv].first)==0) eltrobj=true;
      
      TVector3 Trv = TrigRefObj[trv].second.Vect();
      
      if (mutrobj||eltrobj) {
	
	TVector3 flep1v = lepcand_1.p4.Vect();
	TVector3 flep2v = lepcand_2.p4.Vect();
	
	double tmprat1 = fabs((flep1v-Trv).Mag()/max(1.e-6,flep1v.Mag()));
	double tmprat2 = fabs((flep2v-Trv).Mag()/max(1.e-6,flep2v.Mag()));
	
	hist_init[0]->Fill(tmprat1, weight);
	hist_init[1]->Fill(tmprat2, weight);
	
	if(lepcand_1.pdgId==11||lepcand_2.pdgId==11) {
	  TVector3 felv = (lepcand_1.pdgId==11) ? flep1v : flep2v;
	  double tmpratel = (lepcand_1.pdgId==11) ? tmprat1 : tmprat2;
	  hist2d[0]->Fill(felv.Mag(), tmpratel, weight);
	}
	if(lepcand_1.pdgId==13||lepcand_2.pdgId==13) {
          TVector3 fmuv = (lepcand_1.pdgId==13) ? flep1v : flep2v;
	  double tmpratmu = (lepcand_1.pdgId==13) ? tmprat1 : tmprat2;
          hist2d[0]->Fill(fmuv.Mag(), tmpratmu, weight);
        }
	if (tmprat1 < ptratmin_lep1) {
	  if((mutrobj && lepcand_1.pdgId==13)||(eltrobj && lepcand_1.pdgId==11)){
	    lep1_match = true;
	    ptratmin_lep1 = tmprat1;
	  }
	}
	
	if (tmprat2 < ptratmin_lep2) {
	  if((mutrobj && lepcand_2.pdgId==13)||(eltrobj && lepcand_2.pdgId==11)){
	    lep2_match = true;
	    ptratmin_lep2 = tmprat2;
	  }
	}
      }
		
      if(jettrobj){
	
	double tmprat = fabs((Jets[0].p4.Vect()-Trv).Mag()/max(1.e-6,(Jets[0].p4.Vect().Mag())));
	if (tmprat < ptratmin_jet) {
	  jet_match = true;
	  ptratmin_jet = tmprat;
	}
	
      }
    }
    
    //if (((double_hlts.size()>0 && double_trig_pass && lep1_match && lep2_match) || (single_hlts.size()>0 && single_trig_pass && lep1_match))) { trig_matching_pass = true; }
    if (double_hlts.size()>0 && double_trig_pass && lep1_match && lep2_match)  { trig_matching_pass = true; }
    else if (single_hlts.size()>0 && single_trig_pass && (lep1_match||lep2_match) && single_other_pt_cuts[fired_single_trig]>0 && jet_match) { trig_matching_pass = true; }
    else if (single_hlts.size()>0 && single_trig_pass && (lep1_match||lep2_match) && single_other_pt_cuts[fired_single_trig]<0) { trig_matching_pass = true; }
    
  }
  
}
