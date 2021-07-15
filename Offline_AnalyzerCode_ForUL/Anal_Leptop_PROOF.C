#define Anal_Leptop_PROOF_cxx
#include "Anal_Leptop_PROOF.h"
#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TProofServ.h>
#include "Objects.h"

//#define E_MU_TTBar
//#define E_E_TTBar
#define MU_MU_TTBar

void Anal_Leptop_PROOF::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void Anal_Leptop_PROOF::SlaveBegin(TTree * /*tree*/)
{
  //The SlaveBegin() function is called after the Begin() function.
  //When running with PROOF SlaveBegin() is called on each slave server.
  //The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  OutFile = new TProofOutputFile("Summer20UL18_TTBar_DiLeptonic_mumuoutput.root");
  
  fileOut = OutFile->OpenFile("RECREATE");
  if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",
	      OutFile->GetDir(), OutFile->GetFileName());
    }
  
  isMC = true;
  isTT = true;
  isST = false;
  isDIB = false;
  isWJ = false;
  isDY = false;
  isbQCD = false;
  
  Tout = new TTree("leptop","leptop");
  Tout->Branch("event_pt_weight",&event_pt_weight,"event_pt_weight/F");
  Tout->Branch("weight",&weight,"weight/F");

  Tnewvar = new TTree("newvars","newvars");
  Tnewvar->Branch("M_l1l2",&M_l1l2,"M_l1l2/F");
  Tnewvar->Branch("rat_l1pt_l2pt",&rat_l1pt_l2pt,"rat_l1pt_l2pt/F");
  Tnewvar->Branch("deltaPhi_l1l2",&deltaPhi_l1l2,"deltaPhi_l1l2/F");
  Tnewvar->Branch("l1pt_nearjet",&l1pt_nearjet,"l1pt_nearjet/F");
  Tnewvar->Branch("l2pt_nearjet",&l2pt_nearjet,"l2pt_nearjet/F");
  Tnewvar->Branch("met_pt",&met_pt,"met_pt/F");
  Tnewvar->Branch("met_eta",&met_eta,"met_eta/F");
  Tnewvar->Branch("M_bl1",&M_bl1,"M_bl1/F"); 
  Tnewvar->Branch("M_bl2",&M_bl2,"M_bl2/F");
  Tnewvar->Branch("M_jl1",&M_jl1,"M_jl1/F");
  Tnewvar->Branch("M_jl2",&M_jl2,"M_jl2/F");
  Tnewvar->Branch("delta_phil1_met",&delta_phil1_met,"delta_phil1_met/F");
  Tnewvar->Branch("delta_phil2_met",&delta_phil2_met,"delta_phil2_met/F");
  Tnewvar->Branch("delta_phibl1_met",&delta_phibl1_met,"delta_phibl1_met/F"); 
  Tnewvar->Branch("delta_phibl2_met",&delta_phibl2_met,"delta_phibl2_met/F");
  Tnewvar->Branch("rat_metpt_ak4pt",&rat_metpt_ak4pt,"rat_metpt_ak4pt/F");
  Tnewvar->Branch("rat_metpt_ak8pt",&rat_metpt_ak8pt,"rat_metpt_ak8pt/F");
  Tnewvar->Branch("rat_metpt_eventHT",&rat_metpt_eventHT,"rat_metpt_eventHT/F");
  Tnewvar->Branch("mt_of_l1met",&mt_of_l1met,"mt_of_l1met/F");
  Tnewvar->Branch("mt_of_l2met",&mt_of_l2met,"mt_of_l2met/F");
  Tnewvar->Branch("no_ak4jets",&no_ak4jets,"no_ak4jets/F");
  Tnewvar->Branch("no_ak4bjets",&no_ak4bjets,"no_ak4bjets/F");
  Tnewvar->Branch("no_ak8jets",&no_ak8jets,"no_ak8jets/F");
  Tnewvar->Branch("EventHT",&EventHT,"EventHT/F");
  Tnewvar->Branch("extra_ak4j",&extra_ak4j,"extra_ak4j/F");
  Tnewvar->Branch("ptsum_extra_ak4",&ptsum_extra_ak4,"ptsum_extra_ak4/F");
  Tnewvar->Branch("extra_ak4jqgl",&extra_ak4jqgl,"extra_ak4jqgl/F");
  Tnewvar->Branch("extra_ak4jdeepb",&extra_ak4jdeepb,"extra_ak4jdeepb/F");
  Tnewvar->Branch("rat_extra_ak4jpt_lpt",&rat_extra_ak4jpt_lpt,"rat_extra_ak4jpt_lpt/F");
  Tnewvar->Branch("ak81pt",&ak81pt,"ak81pt/F"); 
  Tnewvar->Branch("ak81y",&ak81y,"ak81y/F");
  Tnewvar->Branch("ak81mass",&ak81mass,"ak81mass/F");
  Tnewvar->Branch("ak81sdmass",&ak81sdmass,"ak81sdmass/F");
  Tnewvar->Branch("ak81deep_tvsqcd",&ak81deep_tvsqcd,"ak81deep_tvsqcd/F");
  Tnewvar->Branch("ak81deep_wvsqcd",&ak81deep_wvsqcd,"ak81deep_wvsqcd/F");
  Tnewvar->Branch("ak82pt",&ak82pt,"ak82pt/F");
  Tnewvar->Branch("ak82y",&ak82y,"ak82y/F");
  Tnewvar->Branch("ak82mass",&ak82mass,"ak82mass/F");
  Tnewvar->Branch("ak82sdmass",&ak82sdmass,"ak82sdmass/F");
  Tnewvar->Branch("ak82deep_tvsqcd",&ak82deep_tvsqcd,"ak82deep_tvsqcd/F");
  Tnewvar->Branch("ak82deep_wvsqcd",&ak82deep_wvsqcd,"ak82deep_wvsqcd/F");
  Tnewvar->Branch("delta_phibl1bl2",&delta_phibl1bl2,"delta_phibl1bl2/F"); 
  //Tnewvar->Branch("delta_phijl1jl2",&delta_phijl1jl2,"delta_phijl1jl2/F");
  Tnewvar->Branch("deltaR_l1l2",&deltaR_l1l2,"deltaR_l1l2/F");
  Tnewvar->Branch("deltaR_l1b1",&deltaR_l1b1,"deltaR_l1b1/F");
  Tnewvar->Branch("deltaR_l2b1",&deltaR_l2b1,"deltaR_l2b1/F");
  Tnewvar->Branch("deltaR_l1b2",&deltaR_l1b2,"deltaR_l1b2/F");
  Tnewvar->Branch("deltaR_l2b2",&deltaR_l2b2,"deltaR_l2b2/F");
  Tnewvar->Branch("deltaR_l1j1",&deltaR_l1j1,"deltaR_l1j1/F");
  Tnewvar->Branch("deltaR_l2j1",&deltaR_l2j1,"deltaR_l2j1/F");
  Tnewvar->Branch("deltaR_l1j2",&deltaR_l1j2,"deltaR_l1j2/F");
  Tnewvar->Branch("deltaR_l2j2",&deltaR_l2j2,"deltaR_l2j2/F");
    
  char name[1000];
  
  for(int binit=0; binit<6; binit++){
    char bnamein[1000]; 
    char btitlein[1000];
    sprintf(bnamein,"hist_%s",binitnames[binit]);
    sprintf(btitlein,"%s",btitlenames[binit]);
    hist_binit[binit] = new TH1D(bnamein,btitlein,ini_bnbins[binit],ini_blow[binit],ini_bup[binit]);
    hist_binit[binit]->Sumw2();
  }
  
  for(int init=0; init<17; init++){
    char namein[1000];// nameinup[1000], nameindown[1000];
    char titlein[1000];
    
    if (init <= 2 || init > 6) {
      
      sprintf(namein,"hist_%s",initnames[init]);
      sprintf(titlein,"%s",titlenames[init]);
      hist_init[init] = new TH1D(namein,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init[init]->Sumw2();
      /*
	sprintf(nameinup,"hist_%s_puwup",initnames[init]);
	hist_initpuwup[init] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initpuwup[init]->Sumw2();
	sprintf(nameindown,"hist_%s_puwdown",initnames[init]);
	hist_initpuwdown[init] = new TH1D(nameindown,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initpuwdown[init]->Sumw2();
	sprintf(nameinup,"hist_%s_btagwup",initnames[init]);
	hist_initbtagwup[init] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initbtagwup[init]->Sumw2();
	sprintf(nameindown,"hist_%s_btagwdown",initnames[init]);
	hist_initbtagwdown[init] = new TH1D(nameindown,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initbtagwdown[init]->Sumw2();
      */
    }
    else {
      sprintf(namein,"%s",initnames[init]);
      sprintf(titlein,"%s",titlenames[init]);
      hist_init[init] = new TH1D(namein,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init[init]->Sumw2();
      /*
	sprintf(nameinup,"hist_%s_puwup",initnames[init]);
	hist_initpuwup[init] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initpuwup[init]->Sumw2();
	sprintf(nameindown,"hist_%s_puwdown",initnames[init]);
	hist_initpuwdown[init] = new TH1D(nameindown,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initpuwdown[init]->Sumw2();
	sprintf(nameinup,"hist_%s_btagwup",initnames[init]);
	hist_initbtagwup[init] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initbtagwup[init]->Sumw2();
	sprintf(nameindown,"hist_%s_btagwdown",initnames[init]);
	hist_initbtagwdown[init] = new TH1D(nameindown,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
	hist_initbtagwdown[init]->Sumw2();
      */
    }
  }
  for(int lvar=0; lvar<45; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Obs_%s",obsnames[lvar]);
    hist_obs[lvar] = new TH1D(lnamein,lnamein,obs_nbins[lvar],obs_low[lvar],obs_up[lvar]);
    hist_obs[lvar]->Sumw2();

  }
  /*
    for(int nvar=0; nvar<15; nvar++){
    char namein_nvar[1000], 
    char titlein_nvar[1000];
    
    sprintf(namein_nvar,"hist_%s",new_var_names[nvar]);
    sprintf(titlein_nvar,"%s",new_var_title[nvar]);
    hist_new_var[nvar] = new TH1D(namein_nvar,titlein_nvar,new_var_nbins[nvar],new_var_low[nvar],new_var_up[nvar]);
    hist_new_var[nvar]->Sumw2();
    }
  */

  char title[1000];
  sprintf(name,"N_PV");
  sprintf(title,"# of Primary Vertices");
  hist_npv = new TH1D(name,title,100,-0.1,99.9);//80,-0.1,79.9);
  hist_npv->Sumw2();
    
  sprintf(name,"N_PV_nopuwt");
  sprintf(title,"# of Primary Vertices");
  hist_npv_nopuwt = new TH1D(name,title,100,-0.1,99.9);
  hist_npv_nopuwt->Sumw2();
  
  sprintf(name,"Counter_event");
  hist_event_count = new TH1D(name,title,2,-0.5,1.5);
  hist_event_count->Sumw2();
        
  hist_2D_msd_deepak8 = new TH2D("hist_2D_msd_deepak8","hist_2D_msd_deepak8",25,0,300,25,0,1);
  hist_2D_msd_deepak8->Sumw2();
  hist_2D_bpass_flavb = new TH2D("h2d_btagpass_flavb","h2d_btagpass_flavb",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavb->Sumw2();
  hist_2D_bpass_flavc = new TH2D("h2d_btagpass_flavc","h2d_btagpass_flavc",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavc->Sumw2();
  hist_2D_bpass_flavq = new TH2D("h2d_btagpass_flavq","h2d_btagpass_flavq",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_bpass_flavq->Sumw2();
  hist_2D_ball_flavb = new TH2D("h2d_flavb","h2d_flavb",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavb->Sumw2();
  hist_2D_ball_flavc = new TH2D("h2d_flavc","h2d_flavc",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavc->Sumw2();
  hist_2D_ball_flavq = new TH2D("h2d_flavq","h2d_flavq",nobptbins,bptbins,nobetabins,betabins);
  hist_2D_ball_flavq->Sumw2();  
  
  hist_count = new TH1D("Counter","Counter",21,0.5,21.5);
  hist_count->Sumw2();  
  
  reader1 = new TMVA::Reader( "BDTG_Re" );
  reader1->AddVariable( "selpfjetAK8NHadF", &in_pfjetAK8NHadF);
  reader1->AddVariable( "selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
  reader1->AddVariable( "selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
  reader1->AddVariable( "selpfjetAK8tau21", &in_pfjetAK8tau21);
  reader1->AddVariable( "selpfjetAK8chrad", &in_pfjetAK8chrad);
  reader1->AddVariable( "selpfjetAK8sdmass", &in_pfjetAK8sdmass);
  reader1->AddVariable( "selpfjetAK8matchedeldxy_sv", &in_pfjetAK8eldxy_sv);
  reader1->AddVariable( "selpfjetAK8matchedelcleta", &in_pfjetAK8matchedelcleta);
  reader1->AddVariable("selpfjetAK8matchedelpt", &in_pfjetAK8matchedelpt);
  reader1->AddVariable("selpfjetAK8matchedelsigmaieta", &in_pfjetAK8matchedelsigmaieta);
  reader1->AddVariable("selpfjetAK8matchedelsigmaiphi", &in_pfjetAK8matchedelsigmaiphi);
  reader1->AddVariable("selpfjetAK8matchedelr9full", &in_pfjetAK8matchedelr9full);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_etaw", &in_pfjetAK8matchedelsupcl_etaw);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_phiw", &in_pfjetAK8matchedelsupcl_phiw);
  reader1->AddVariable("selpfjetAK8matchedelhcaloverecal", &in_pfjetAK8matchedelhcaloverecal);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkn", &in_pfjetAK8matchedelcloctftrkn);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkchi2", &in_pfjetAK8matchedelcloctftrkchi2);
  reader1->AddVariable("selpfjetAK8matchedele1x5bye5x5", &in_pfjetAK8matchedele1x5bye5x5);
  reader1->AddVariable("selpfjetAK8matchedelnormchi2", &in_pfjetAK8matchedelnormchi2);
  reader1->AddVariable("selpfjetAK8matchedelhitsmiss", &in_pfjetAK8matchedelhitsmiss);
  reader1->AddVariable("selpfjetAK8matchedeltrkmeasure", &in_pfjetAK8matchedeltrkmeasure);
  reader1->AddVariable("selpfjetAK8matchedelecloverpout", &in_pfjetAK8matchedelecloverpout);
  reader1->AddVariable("selpfjetAK8matchedelecaletrkmomentum", &in_pfjetAK8matchedelecaletrkmomentum);
  reader1->AddVariable("selpfjetAK8matchedeldeltaetacltrkcalo", &in_pfjetAK8matchedeldeltaetacltrkcalo);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_preshvsrawe", &in_pfjetAK8matchedelsupcl_preshvsrawe);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumphet", &in_pfjetAK8matchedelpfisolsumphet);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumchhadpt", &in_pfjetAK8matchedelpfisolsumchhadpt);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumneuhadet", &in_pfjetAK8matchedelpfisolsumneuhadet);
  reader1->AddVariable("selpfjetAK8matchedeletain", &in_pfjetAK8matchedeletain);
  reader1->AddVariable("selpfjetAK8matchedelphiin", &in_pfjetAK8matchedelphiin);
  reader1->AddVariable("selpfjetAK8matchedelfbrem", &in_pfjetAK8matchedelfbrem);
  reader1->AddVariable("selpfjetAK8matchedeleoverp", &in_pfjetAK8matchedeleoverp);
  reader1->AddVariable("selpfjetAK8matchedelhovere", &in_pfjetAK8matchedelhovere);
  reader1->AddVariable("selpfjetAK8matchedelRho", &in_pfjetAK8matchedelRho);
  reader1->BookMVA("BDTG method", weightfile1);
  
  reader3 = new TMVA::Reader( "BDTG_Rt" );
  reader3->AddVariable( "selpfjetAK8NHadF", &in_pfjetAK8NHadF);
  reader3->AddVariable( "selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
  reader3->AddVariable( "selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
  reader3->AddVariable( "selpfjetAK8tau21", &in_pfjetAK8tau21);
  reader3->AddVariable( "selpfjetAK8chrad", &in_pfjetAK8chrad);
  reader3->AddVariable( "selpfjetAK8sdmass", &in_pfjetAK8sdmass);
  reader3->BookMVA("BDTG method", weightfile3);
  
  reader4 = new TMVA::Reader( "BDTG_Rmu" );
  reader4->AddVariable( "selpfjetAK8NHadF", &in_pfjetAK8NHadF);
  reader4->AddVariable( "selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
  reader4->AddVariable( "selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
  reader4->AddVariable( "selpfjetAK8tau21", &in_pfjetAK8tau21);
  reader4->AddVariable( "selpfjetAK8chrad", &in_pfjetAK8chrad);
  reader4->AddVariable( "selpfjetAK8sdmass", &in_pfjetAK8sdmass);
  reader4->AddVariable("selpfjetAK8matchedmuonchi", &in_pfjetAK8matchedmuonchi);
  reader4->AddVariable("selpfjetAK8matchedmuonposmatch", &in_pfjetAK8matchedmuonposmatch);
  reader4->AddVariable("selpfjetAK8matchedmuontrkink", &in_pfjetAK8matchedmuontrkink);
  reader4->AddVariable("selpfjetAK8matchedmuonsegcom", &in_pfjetAK8matchedmuonsegcom);
  reader4->AddVariable("selpfjetAK8matchedmuonhit", &in_pfjetAK8matchedmuonhit);
  reader4->AddVariable("selpfjetAK8matchedmuonmst", &in_pfjetAK8matchedmuonmst);
  reader4->AddVariable("selpfjetAK8matchedmuontrkvtx", &in_pfjetAK8matchedmuontrkvtx);
  reader4->AddVariable("selpfjetAK8matchedmuondz", &in_pfjetAK8matchedmuondz);
  reader4->AddVariable("selpfjetAK8matchedmuonpixhit", &in_pfjetAK8matchedmuonpixhit);
  reader4->AddVariable("selpfjetAK8matchedmuontrklay", &in_pfjetAK8matchedmuontrklay);
  reader4->AddVariable("selpfjetAK8matchedmuonvalfrac", &in_pfjetAK8matchedmuonvalfrac);
  reader4->AddVariable("selpfjetAK8muinsubptrat", &in_pfjetAK8muinsubptrat);
  reader4->AddVariable("selpfjetAK8muinsubmassrat", &in_pfjetAK8muinsubmassrat);
  reader4->AddVariable("selpfjetAK8muinsubinvmass", &in_pfjetAK8muinsubinvmass);
  reader4->AddVariable("selpfjetAK8muinsubIfarbyI0", &in_pfjetAK8muinsubIfarbyI0);
  reader4->AddVariable("selpfjetAK8muinsubInearbyI0", &in_pfjetAK8muinsubInearbyI0);
  reader4->BookMVA("BDTG method", weightfile4);
}

Bool_t Anal_Leptop_PROOF::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Anal_Leptop_PROOF::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  GetEntry(entry);
  if(isMC){
    weight = event_weight;
  }else{
    weight = 1;
  }
  Tout->Fill();
  
  int ngenelc = 0;
  int ngenmu = 0;
  int ngentau = 0;
  int ngenqg = 0;
  int ngenb = 0;
  
  TLorentzVector leptop4v[2];
  TLorentzVector leptop4v_daught[3][2];
  TLorentzVector hadtop4v[2];
  TLorentzVector hadtop4v_daught[3][2];
  int leptop_id_daught[2];
  
  if(isMC){
    
    for(int igen=0; igen<ngenparticles; igen++){
      if(abs(genpartstatus[igen])!=23 && genpartstatus[igen]!=1) continue;
      if(!(genpartfromhard[igen]) /*|| !(genpartfromhardbFSR[igen])*/) continue;
      
      if(abs(genpartpdg[igen])==11) {
        genelectronpt[ngenelc] = genpartpt[igen];
        genelectroneta[ngenelc] = genparteta[igen];
        genelectronphi[ngenelc] = genpartphi[igen];
        genelectronm[ngenelc] = genpartm[igen];
        ngenelc++;
      }
      if(abs(genpartpdg[igen])==13) {
        genmuonpt[ngenmu] = genpartpt[igen];
	genmuoneta[ngenmu] = genparteta[igen];
        genmuonphi[ngenmu] = genpartphi[igen];
        genmuonm[ngenmu] = genpartm[igen];
        ngenmu++;
      }
      if(abs(genpartpdg[igen])==15) {
	gentaupt[ngentau] = genpartpt[igen];
        gentaueta[ngentau] = genparteta[igen];
        gentauphi[ngentau] = genpartphi[igen];
        gentaum[ngentau] = genpartm[igen];
        ngentau++;
      }
      if((abs(genpartpdg[igen])>=1 && abs(genpartpdg[igen])<5)||(abs(genpartpdg[igen])==21)) {
        genqgpt[ngenqg] = genpartpt[igen];
        genqgeta[ngenqg] = genparteta[igen];
        genqgphi[ngenqg] = genpartphi[igen];
        genqgm[ngenqg] = genpartm[igen];
        ngenqg++;
      }
      if(abs(genpartpdg[igen])==5) {
        genbpt[ngenb] = genpartpt[igen];
        genbeta[ngenb] = genparteta[igen];
        genbphi[ngenb] = genpartphi[igen];
        genbm[ngenb] = genpartm[igen];
        ngenb++;
      }
    }
    ngenelectrons = ngenelc;
    ngenmuons = ngenmu;
    ngentaus = ngentau;
    ngenqgs = ngenqg;
    ngenbs = ngenb;
    int ngent = 0;
    
    for(int igen=0; igen<ngenparticles; igen++){
      
      if(abs(genpartstatus[igen])!=22) continue;
      if(!(genpartfromhard[igen])) continue;
      if(abs(genpartpdg[igen])!=6) continue;
      
      gentoppt[ngent] = genpartpt[igen];
      gentopeta[ngent] = genparteta[igen];
      gentopphi[ngent] = genpartphi[igen];
      gentopm[ngent] = genpartm[igen];
      gentopid[ngent] = igen;
      
      ngent++;
    }
    
    ngentops = ngent;
  }

  int top_dp[6];
  int idp = 0;
  int ndq = 0; int ndb = 0; int ndl = 0;
  int top_dqp[4] = {-1,-1,-1,-1};
  int top_dbp[2] = {-1,-1};
  int top_dlp[4] = {-1,-1,-1,-1};
  
  for(int igen=0; igen<ngenparticles; igen++){
    
    if(!(genpartstatus[igen]==23 || genpartstatus[igen]==1)) continue;
    if(!(genpartfromhard[igen])) continue;
    
    if(!((abs(genpartpdg[igen])>=1 && abs(genpartpdg[igen])<=5)||(abs(genpartpdg[igen])>=11 && abs(genpartpdg[igen])<=16))) continue;
    if(!(abs(genpartmompdg[igen])==6 || abs(genpartmompdg[igen])==24)) continue;
    top_dp[idp] = igen;
    if(abs(genpartpdg[igen])>=1 && abs(genpartpdg[igen])<5 && abs(genpartmompdg[igen])==24) {  top_dqp[ndq] = igen;  ndq++; }
    if(abs(genpartpdg[igen])>=11 && abs(genpartpdg[igen])<=16 && abs(genpartmompdg[igen])==24) {  top_dlp[ndl] = igen;  ndl++; }
    if(abs(genpartpdg[igen])==5 && abs(genpartmompdg[igen])==6) {  top_dbp[ndb] = igen;  ndb++; }
    idp++;
  }
  
  if(isMC && isTT){
    
    nleptop = nhadtop = -1;
    
    if(ndl==0) { nleptop = 0; }
    if(ndl==2) { nleptop = 1; }
    if(ndl==4) { nleptop = 2; }
    
    if(ndq==0) { nhadtop = 0; }
    if(ndq==2) { nhadtop = 1; }
    if(ndq==4) { nhadtop = 2; }
    
    for(int ilep=0; ilep<ndl; ilep++){
      for(int jlep=(ilep+1); jlep<ndl; jlep++){
	if((abs(abs(genpartpdg[top_dlp[ilep]])-abs(genpartpdg[top_dlp[jlep]]))==1) && (genpartpdg[top_dlp[ilep]]*genpartpdg[top_dlp[jlep]]<0) && (genpartmompdg[top_dlp[ilep]]==genpartmompdg[top_dlp[jlep]]))
	  {
	    genpartpair[top_dlp[ilep]] = top_dlp[jlep];
	    genpartpair[top_dlp[jlep]] = top_dlp[ilep];
	  }
      }
    }
    
    for(int iq=0; iq<ndq; iq++){
      for(int jq=(iq+1); jq<ndq; jq++){
	if( ( (abs(abs(genpartpdg[top_dqp[iq]])-abs(genpartpdg[top_dqp[jq]]))==1)||(abs(abs(genpartpdg[top_dqp[iq]])-abs(genpartpdg[top_dqp[jq]]))==3) ) && (genpartpdg[top_dqp[iq]]*genpartpdg[top_dqp[jq]]<0) && (genpartmompdg[top_dqp[iq]]==genpartmompdg[top_dqp[jq]]) )
	  {
	    genpartpair[top_dqp[iq]] = top_dqp[jq];
	    genpartpair[top_dqp[jq]] = top_dqp[iq];
	  }
      }
    }
    
    int ileptop = 0;
    int ihadtop = 0;
    
    if(nleptop>0){
      
      for(int ilep=0; ilep<ndl; ilep++){
	if(!((abs(genpartpdg[top_dlp[ilep]])==11 || abs(genpartpdg[top_dlp[ilep]])==13 || abs(genpartpdg[top_dlp[ilep]])==15) && genpartpair[top_dlp[ilep]]>=0)) continue;
	
	for(int ib=0; ib<ndb; ib++){
	  
	  if(genpartmompdg[top_dbp[ib]]*genpartmompdg[top_dlp[ilep]]>0){
	    TLorentzVector b4; b4.SetPtEtaPhiM(genpartpt[top_dbp[ib]],genparteta[top_dbp[ib]],genpartphi[top_dbp[ib]],genpartm[top_dbp[ib]]);
	    int ipar = top_dlp[ilep];
	    if(ipar>=0){
	      TLorentzVector q1; q1.SetPtEtaPhiM(genpartpt[ipar],genparteta[ipar],genpartphi[ipar],genpartm[ipar]);
	      int jpar = genpartpair[top_dlp[ilep]];
	      TLorentzVector q2; q2.SetPtEtaPhiM(genpartpt[jpar],genparteta[jpar],genpartphi[jpar],genpartm[jpar]);
	      leptop4v[ileptop] = (b4+q1+q2);
	      leptop4v_daught[0][ileptop] = q1;
	      leptop4v_daught[1][ileptop] = q2;
	      leptop4v_daught[2][ileptop] = b4;
	      leptop_id_daught[ileptop] = genpartpdg[top_dlp[ilep]];
	      
	      ileptop++;
	      break;
	    }
	  }
	  
	}
	if(ileptop>=nleptop) break;
      }
    }
    
    if(nhadtop>0){
      
      for(int iq=0; iq<ndq; iq++){
	
	if(!(( (abs(genpartpdg[top_dqp[iq]])==1) || (abs(genpartpdg[top_dqp[iq]])==3) ) && genpartpair[top_dqp[iq]]>=0)) continue;
	
	for(int ib=0; ib<ndb; ib++){
	  
	  if(genpartmompdg[top_dbp[ib]]*genpartmompdg[top_dqp[iq]]>0){
	    TLorentzVector b4; b4.SetPtEtaPhiM(genpartpt[top_dbp[ib]],genparteta[top_dbp[ib]],genpartphi[top_dbp[ib]],genpartm[top_dbp[ib]]);
	    int ipar = top_dqp[iq];
	    if(ipar>=0){
	      TLorentzVector q1; q1.SetPtEtaPhiM(genpartpt[ipar],genparteta[ipar],genpartphi[ipar],genpartm[ipar]);
	      int jpar = genpartpair[top_dqp[iq]];
	      TLorentzVector q2; q2.SetPtEtaPhiM(genpartpt[jpar],genparteta[jpar],genpartphi[jpar],genpartm[jpar]);
	      hadtop4v[ihadtop] = (b4+q1+q2);
	      hadtop4v_daught[0][ihadtop] = q1;
	      hadtop4v_daught[1][ihadtop] = q2;
	      hadtop4v_daught[2][ihadtop] = b4;
	      
	      ihadtop++;
	      break;
	    }
	  }
	  
	}
	if(ihadtop>=nhadtop) break;
      }
      
    }
    
    nleptop = ileptop;
    nhadtop = ihadtop;
    
  }//isMC
  
  DiLeptt = SemiLeptt = Hadtt = EE = MUMU = EMU = EJets = MUJets = TauTau = ETau = MuTau = false;
  
  if(nleptop==2 && nhadtop==0) { DiLeptt = true; }
  if(nleptop==1 && nhadtop==1) { SemiLeptt = true; }
  if(nleptop==0 && nhadtop==2) { Hadtt = true; }
  
  if(DiLeptt && abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==11) { EE = true; }
  if(DiLeptt && abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==13) { MUMU = true; }
  if(DiLeptt && abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==15) { TauTau = true; }
  if(DiLeptt && ((abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==13) || (abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==11)) ) { EMU = true; }
  if(DiLeptt && ((abs(leptop_id_daught[0])==11 && abs(leptop_id_daught[1])==15) || (abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==11)) ) { ETau = true; }
  if(DiLeptt && ((abs(leptop_id_daught[0])==13 && abs(leptop_id_daught[1])==15) || (abs(leptop_id_daught[0])==15 && abs(leptop_id_daught[1])==13)) ) { MuTau = true; }
  
  
  if(SemiLeptt && abs(leptop_id_daught[0])==11) { EJets = true; }
  if(SemiLeptt && abs(leptop_id_daught[0])==13) { MUJets = true; }
  if(SemiLeptt && abs(leptop_id_daught[0])==15) { TAUJets = true; }
  
  bool boosted = ((SemiLeptt && (leptop4v[0].Pt()>400)) || (DiLeptt && ((abs(leptop_id_daught[0])==11 && (leptop4v[0].Pt()>400)) || (abs(leptop_id_daught[1])==11 && (leptop4v[1].Pt()>400))))) ; 
  
  if (isMC && isTT) {
#ifdef E_MU_TTBar
    //if(!(DiLeptt && EMU)) return kFALSE; //for signal EMU
    if((DiLeptt && EMU)) return kFALSE; //for non-signal EMU TTbar
#elif defined(E_E_TTBar)
    //if(!(DiLeptt && EE)) return kFALSE; //for signal EE
    if((DiLeptt && EE)) return kFALSE; //for non-signal EE TTbar
#elif defined(MU_MU_TTBar) 
    if(!(DiLeptt && MUMU)) return kFALSE; //for signal MUMU
    //if((DiLeptt && MUMU)) return kFALSE; //for non-signal MUMU TTbar
#endif
  }
  
  hist_event_count->Fill(1,weight);
  int nbjet_cut = -1;
  
#ifdef E_MU_TTBar
  nbjet_cut = 1;
#elif defined(E_E_TTBar)
  nbjet_cut = 1;
#elif defined(MU_MU_TTBar)
  nbjet_cut = 1;
#elif defined(E_Jets_TTBar)
  nbjet_cut = 1;
#endif
    
  //******looping first on AK4 jets as used by AK8, el and muon*********//
  
  int nbjetAK4 = 0;
  int nbjetAK4_lead = 0;
  
  float btagwt = 1.;
  float btagwtup = 1.;
  float btagwtdown = 1.;
  float btag_eff = 1;
  float Event_HT = 0;
  
  vector <AK4Jet> Jets;
  std::vector<TLorentzVector> bjv;
  
  for(int ijet=0; ijet<npfjetAK4; ijet++){

    if(pfjetAK4jetID[ijet]==0) continue;
    
    if(fabs(pfjetAK4eta[ijet])>2.5) continue;
    pfjetAK4pt[ijet] *= pfjetAK4JEC[ijet] ;
    pfjetAK4mass[ijet] *= pfjetAK4JEC[ijet];
    
    if(isMC){
      pfjetAK4pt[ijet] *= (1+pfjetAK4reso[ijet]) ;
      pfjetAK4mass[ijet] *= (1+pfjetAK4reso[ijet]) ;
    }
    if(pfjetAK4pt[ijet]<30.) continue;

    Event_HT += pfjetAK4pt[ijet];
    AK4Jet sJet;
    sJet.jetid = pfjetAK4jetID[ijet];
    sJet.jetid_tightlepveto = pfjetAK4jetID_tightlepveto[ijet];
    sJet.pt = pfjetAK4pt[ijet];
    sJet.mass = pfjetAK4mass[ijet];
    sJet.eta = pfjetAK4eta[ijet];
    sJet.y = pfjetAK4y[ijet];
    sJet.phi = pfjetAK4phi[ijet];
    sJet.hadronFlavour = pfjetAK4hadronflav[ijet];
    sJet.partonFlavour = pfjetAK4partonflav[ijet];
    sJet.btag_DeepFlav = pfjetAK4btag_DeepFlav[ijet];
    sJet.btag_DeepCSV = pfjetAK4btag_DeepCSV[ijet];
    sJet.puid = pfjetAK4PUID[ijet];
    sJet.qgl = pfjetAK4qgl[ijet];
    
    Jets.push_back(sJet);
    
    if(Jets.size()>=njetmx) break;
    
  }
  npfjetAK4 = Jets.size();
  sorted_by_pt(Jets);
  
  for(unsigned ijet=0; ijet<Jets.size(); ijet++){
    
    if(Jets[ijet].btag_DeepFlav > deep_btag_cut) {
      if(abs(Jets[ijet].hadronFlavour)==5){  hist_2D_bpass_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
      if(abs(Jets[ijet].hadronFlavour)==4){  hist_2D_bpass_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
      if(abs(Jets[ijet].hadronFlavour)!=5 && abs(Jets[ijet].hadronFlavour)!=4){  
	hist_2D_bpass_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); 
      }
    }
	if(abs(Jets[ijet].hadronFlavour)==5){ hist_2D_ball_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
	if(abs(Jets[ijet].hadronFlavour)==4){ hist_2D_ball_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
	if(abs(Jets[ijet].hadronFlavour)!=5 && abs(Jets[ijet].hadronFlavour)!=4){  hist_2D_ball_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
	
	if(Jets[ijet].btag_DeepFlav > deep_btag_cut) { nbjetAK4++;
	TLorentzVector bv;
      bv.SetPtEtaPhiM(Jets[ijet].pt,Jets[ijet].eta,Jets[ijet].phi,Jets[ijet].mass);
      bjv.push_back(bv);
      }
      if((Jets[ijet].btag_DeepFlav > deep_btag_cut) && (ijet==0||ijet==1)) { nbjetAK4_lead++; }
      
      /*btag sf along with its uncertainty need to be updated for UL*

	if(isMC){
	
	if(Jets[ijet].btag_DeepFlav > deep_btag_cut) {
	
	btagwt *= max(float(0),float(BTag_SF(Jets[ijet].hadronFlavour,"noSyst",Jets[ijet].pt)));
	btagwtup *= max(float(0),float(BTag_SF(Jets[ijet].hadronFlavour,"up",Jets[ijet].pt)));
	btagwtdown *= max(float(0),float(BTag_SF(Jets[ijet].hadronFlavour,"down",Jets[ijet].pt)));
	} else{
	
	if(isTT){       btag_eff = BTag_MCEfficiency_TT(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);      }
	if(isST){       btag_eff = BTag_MCEfficiency_ST(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);      }
	if(isDIB){      btag_eff = BTag_MCEfficiency_DIB(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);     }
	if(isWJ){       btag_eff = BTag_MCEfficiency_WJ(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);      }
	if(isDY){       btag_eff = BTag_MCEfficiency_DY(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);      }
	if(isbQCD){     btag_eff = BTag_MCEfficiency_bQCD(abs(Jets[ijet].hadronFlavour),Jets[ijet].pt,pfjetAK4eta[fjet]);    }
	
	btagwt *= max(float(0),float((1. - btag_eff*BTag_SF(Jets[ijet].hadronFlavour,"noSyst",Jets[ijet].pt))*1./(1. - btag_eff)) );
	btagwtup *= max(float(0),float((1. - btag_eff*BTag_SF(Jets[ijet].hadronFlavour,"up",Jets[ijet].pt))*1./(1. - btag_eff)) );
	btagwtdown *= max(float(0),float((1. - btag_eff*BTag_SF(Jets[ijet].hadronFlavour,"down",Jets[ijet].pt))*1./(1. - btag_eff)) );
	}
      */
  }
  
  //*****//
  vector <Muon> vmuons;
  
  for(int mu=0; mu<nmuons; mu++){
    
    if(muonpt[mu]<25.) continue; 
    
    if(fabs(muoneta[mu])>2.5)  continue; 
    bool mu_id = Muon_TightID(muonisGL[mu],muonisPF[mu],
			      muonchi[mu],muonhit[mu],muonmst[mu],
			      muontrkvtx[mu],muondz[mu],muonpixhit[mu],muontrklay[mu]);
    bool mu_iso = Muon_Iso_ID(muonpfiso[mu]);
    
    if(!mu_id) continue;
    //if(!mu_iso) continue;
    
    Muon vmuon;
    vmuon.pt = muonpt[mu];
    vmuon.eta = muoneta[mu];
    vmuon.phi = muonphi[mu];
    vmuon.charge = muoncharge[mu];
    vmuon.trkvtx = muontrkvtx[mu];
    vmuon.dz = muondz[mu];
    vmuon.ip = (fabs(muontrkvtx[mu])<0.2 && fabs(muondz[mu])<0.5);
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

    vmuons.push_back(vmuon);
    
    if(vmuons.size()>=njetmx) break;
  }
  
  nmuons = vmuons.size();
  
  sorted_by_pt(vmuons);

  vector <AK8Jet> LJets;

  for(int ijet=0; ijet<npfjetAK8; ijet++){
    
    if(!pfjetAK8jetID[ijet]) continue;

    if(fabs(pfjetAK8y[ijet])>2.5) continue;
    
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
    
    if(pfjetAK8pt[ijet] < ptcut) continue;
    
    AK8Jet LJet;
    LJet.jetid = pfjetAK8jetID[ijet];
    LJet.jetid_tightlepveto = pfjetAK8jetID_tightlepveto[ijet];
    LJet.pt = pfjetAK8pt[ijet];
    LJet.eta = pfjetAK8eta[ijet];
    LJet.mass = pfjetAK8mass[ijet];
    LJet.phi = pfjetAK8phi[ijet];
    LJet.y = pfjetAK8y[ijet];
    LJet.pt_resoup = pfjetAK8pt_resoup[ijet];
    LJet.mass_resoup = pfjetAK8mass_resoup[ijet];
    LJet.pt_resodn = pfjetAK8pt_resodown[ijet];
    LJet.mass_resodn = pfjetAK8mass_resodown[ijet];
    LJet.jesup_total = pfjetAK8jesup_total[ijet];
    LJet.jesdn_total = pfjetAK8jesdn_total[ijet];

    LJet.chrad = pfjetAK8chrad[ijet];
    LJet.tau21 = pfjetAK8tau2[ijet]*1./pfjetAK8tau1[ijet];
    LJet.tau32 = pfjetAK8tau3[ijet]*1./pfjetAK8tau2[ijet];
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
    LJet.HadF = (pfjetAK8NHF[ijet]+pfjetAK8CHF[ijet]);
    pfjetAK8HadF[ijet] = (pfjetAK8NHF[ijet]+pfjetAK8CHF[ijet]);
    LJet.NHadF = (1.- pfjetAK8HadF[ijet]);

    LJet.EMF = (pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]);
    LJet.EEM = pfjetAK8EEM[ijet];

    LJet.neucons = pfjetAK8Neucons[ijet];
    LJet.chcons = pfjetAK8Chcons[ijet];

    LJet.neuemfrac = (pfjetAK8NEMF[ijet]*1.)/(pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]);
    LJet.neunhadfrac = (pfjetAK8NEMF[ijet]*1.)/(1.- pfjetAK8HadF[ijet]);

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
    LJet.matchAK4deepb = LJet.re_tvsb = LJet.rmu_tvsb = -100;
    
    LJets.push_back(LJet);
    
    if(LJets.size()>=njetmxAK8) break;
    
  }

  npfjetAK8 = LJets.size();
  sorted_by_pt(LJets);
  
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
    if(isMC){
      for(int ngen=0; ngen<ngenelectrons; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,genelectroneta[ngen],genelectronphi[ngen])<0.6){
	  LJets[ijet].haselectron  = true;
	  break;
	}
      }
      for(int ngen=0; ngen<ngenmuons; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,genmuoneta[ngen],genmuonphi[ngen])<0.6){
	  LJets[ijet].hasmuon = true;
	  break;
	}
      }
      
      for(int ngen=0; ngen<ngentaus; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentaueta[ngen],gentauphi[ngen])<0.6){
	  LJets[ijet].hastau = true;
	  break;
	}
      }
      
      for(int ngen=0; ngen<ngenqgs; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,genqgeta[ngen],genqgphi[ngen])<0.6){
	  LJets[ijet].hasqg = true;
	  break;
	}
      }
      
      for(int ngen=0; ngen<ngenbs; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,genbeta[ngen],genbphi[ngen])<0.6){
	  LJets[ijet].hasb = true;
	  break;
	}
      }
      for(int ngen=0; ngen<ngentops; ngen++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,gentopeta[ngen],gentopphi[ngen])<0.6){
	  LJets[ijet].hastop = true;
	  break;
	}
      }
      
      for(int ilept=0; ilept<nleptop; ilept++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,leptop4v[ilept].Rapidity(),leptop4v[ilept].Phi())<0.6){
	  LJets[ijet].hasleptop = true;
	  break;
	}
      }
      for(int ilept=0; ilept<nleptop; ilept++){
	bool match[3] = {0};
	for(int idaut=0; idaut<3; idaut++){
	  if(delta2R(LJets[ijet].y,LJets[ijet].phi,leptop4v_daught[idaut][ilept].Rapidity(),leptop4v_daught[idaut][ilept].Phi())<0.8){
	    match[idaut] = true;
	  }
	}
	if(match[0]&&match[1]&&match[2]){
	  LJets[ijet].hasleptop_alldecay = true;
	  break;
	}
      }
      
      for(int ihadt=0; ihadt<nhadtop; ihadt++){
	if(delta2R(LJets[ijet].y,LJets[ijet].phi,hadtop4v[ihadt].Rapidity(),hadtop4v[ihadt].Phi())<0.6){
	  LJets[ijet].hashadtop = true;
	  break;
	}
      }
      
      for(int ihadt=0; ihadt<nhadtop; ihadt++){
	bool match[3] = {0};
	for(int idaut=0; idaut<3; idaut++){
	  if(delta2R(LJets[ijet].y,LJets[ijet].phi,hadtop4v_daught[idaut][ihadt].Rapidity(),hadtop4v_daught[idaut][ihadt].Phi())<0.8){
	    match[idaut] = true;
	  }
	}
	if(match[0]&&match[1]&&match[2]){
	  LJets[ijet].hashadtop_alldecay = true;
	  break;
	}
      }
      
    }//if (isMC)
    
    if (LJets[ijet].elinsubpt > 0. && LJets[ijet].elinsubjpt > 0.) {
      if (delta2R(LJets[ijet].elinsubeta,LJets[ijet].elinsubphi,LJets[ijet].y,LJets[ijet].phi) < 0.8) {
	LJets[ijet].haspfelectron = true;
      }
    }
    
    int pfjetAK8mactAK4_i = -1;
    float minR = 0.4;
    
    for(unsigned kjet=0; kjet<Jets.size(); kjet++){
      if(delta2R(LJets[ijet].y,LJets[ijet].phi,Jets[kjet].eta,Jets[kjet].phi)<minR){
	minR = delta2R(LJets[ijet].y,LJets[ijet].phi,Jets[kjet].eta,Jets[kjet].phi);
	pfjetAK8mactAK4_i = kjet;
      }
    }
    
    if(pfjetAK8mactAK4_i>=0){
      LJets[ijet].matchAK4deepb = Jets[pfjetAK8mactAK4_i].btag_DeepFlav;
    }
    
    if (LJets[ijet].muinsubpt > 0. && LJets[ijet].muinsubjpt > 0.) {
      if (delta2R(LJets[ijet].muinsubeta,LJets[ijet].muinsubphi,LJets[ijet].y,LJets[ijet].phi) < 0.8) {
        LJets[ijet].haspfmuon = true;
      }
    }
        
  }//ijet
  
  int maxM(-1);
  float LMass(0);
  for(unsigned ljet=0; ljet<LJets.size(); ljet++){
    if (!(isnan(LJets[ljet].mass))) {
      if (LJets[ljet].mass > LMass) {
	LMass = LJets[ljet].mass;
	maxM = ljet;
      }
    }
  }
  
  vector <Electron> velectrons;

  for(int ie=0; ie<nelecs; ie++) {
    if (elpt[ie]<25.) continue;
    if(fabs(eleta[ie])>2.5)  continue; 
    // if(!elmvaid[ie]) continue;
    if(!elmvaid_noIso[ie]) continue;
    
    Electron velectron;

    velectron.pt = elpt[ie];
    velectron.eta = eleta[ie];
    velectron.phi = elphi[ie];
    velectron.charge = elcharge[ie];
    velectron.id = elmvaid[ie];
    velectron.Fallv2WP80 = elmvaid_Fallv2WP80[ie];
    velectron.id_noIso = elmvaid_noIso[ie];
    velectron.Fallv2WP80_noIso = elmvaid_Fallv2WP80_noIso[ie];
    velectron.p = elp[ie];
    velectron.dxy = eldxytrk[ie];
    velectron.dz = eldztrk[ie];
    velectron.ip = ((fabs(elsupcl_eta[ie])<1.4442 && fabs(eldxytrk[ie])<0.05 && fabs(eldztrk[ie])<0.1)||(fabs(elsupcl_eta[ie])>1.5660 && fabs(elsupcl_eta[ie])<2.5 && fabs(eldxytrk[ie])<0.1 && fabs(eldztrk[ie])<0.2)); 
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
     
    velectrons.push_back(velectron);
    
    if(velectrons.size() >= njetmx) break;
    
  }
  
  nelecs = velectrons.size();
  sorted_by_pt(velectrons);

  hist_npv_nopuwt->Fill(nprimi,weight);
  
  //pileup reweighting has to be updated with UL//
  if(isMC){
    
    if(npu_vert>=0 && npu_vert<100){
      puWeight = pu_rat18[npu_vert];
      puWeightUp = pu_rat18_up[npu_vert];
      puWeightDown = pu_rat18_dn[npu_vert];
    }
    if(!isnan(puWeightUp) || fabs(puWeightUp)<1.e+6){
      //weight_puwup = weight*puWeightUp;
     }
    if(!isnan(puWeightDown) || fabs(puWeightDown)<1.e+6){
      //weight_puwdown = weight*puWeightDown;
    }
    
    if(!isnan(puWeight) || fabs(puWeight)<1.e+6){
      // weight *= puWeight;  
    }
  }
  
  hist_npv->Fill(nprimi,weight);
  
  //btag weight has to be updated with UL//
  
  if((!isnan(btagwtup) || fabs(btagwtup)<1.e+6) && nbjet_cut>0){
    // weight_btagwup = weight*btagwtup;
  }
  if((!isnan(btagwtdown) || fabs(btagwtdown)<1.e+6) && nbjet_cut>0){
    //weight_btagwdown = weight*btagwtdown;
  }
  
  if((!isnan(btagwt) || fabs(btagwt)<1.e+6) && nbjet_cut>0){
    //weight *= btagwt;
    //weight_puwup *= btagwt;
    //weight_puwdown *= btagwt;
  }
  
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
     
    if(isnan(LJets[ijet].pt)) {LJets[ijet].pt = -100; }
    if(isnan(LJets[ijet].y)) { LJets[ijet].y = -100; }
    if(isnan(LJets[ijet].mass)) { LJets[ijet].mass = -100; }
    if(isnan(LJets[ijet].phi)) { LJets[ijet].phi = -100; }
    
    if(isnan(LJets[ijet].btag_DeepCSV)) { LJets[ijet].btag_DeepCSV = -100; }
    if(isnan(LJets[ijet].matchAK4deepb)) { LJets[ijet].matchAK4deepb = -100; }
    if(isnan(LJets[ijet].DeepTag_TvsQCD)) { LJets[ijet].DeepTag_TvsQCD = -100; }
    if(isnan(LJets[ijet].DeepTag_WvsQCD)) { LJets[ijet].DeepTag_WvsQCD = -100; }
    if(isnan(LJets[ijet].DeepTag_ZvsQCD)) { LJets[ijet].DeepTag_ZvsQCD = -100; }
    
    if(isnan(LJets[ijet].CHF)) { LJets[ijet].CHF = -100; }
    if(isnan(LJets[ijet].NHF)) { LJets[ijet].NHF = -100; }
    if(isnan(LJets[ijet].CEMF)) { LJets[ijet].CEMF = -100; }
    if(isnan(LJets[ijet].NEMF)) { LJets[ijet].NEMF = -100; }
    if(isnan(LJets[ijet].MUF)) { LJets[ijet].MUF = -100; }
    if(isnan(LJets[ijet].HadF)) { LJets[ijet].HadF = -100; }
    if(isnan(LJets[ijet].NHadF)) { LJets[ijet].NHadF = -100; }
    if(isnan(LJets[ijet].EMF)) { LJets[ijet].EMF = -100; }
    if(isnan(LJets[ijet].neuemfrac)) { LJets[ijet].neuemfrac = -100; }
    if(isnan(LJets[ijet].neunhadfrac)) { LJets[ijet].neunhadfrac = -100; }
    if(isnan(LJets[ijet].EEM)) { LJets[ijet].EEM = -100; }
   
    
    if(isnan(LJets[ijet].chrad)) { LJets[ijet].chrad = -100; }
    if(isnan(LJets[ijet].tau21)) { LJets[ijet].tau21 = -100; }
    if(isnan(LJets[ijet].tau32)) { LJets[ijet].tau32 = -100; }
    if(isnan(LJets[ijet].sdmass)) { LJets[ijet].sdmass = -100; }

    if(isnan(LJets[ijet].elinsubpt)) { LJets[ijet].elinsubpt = -100; }
    if(isnan(LJets[ijet].elinsubeta)) { LJets[ijet].elinsubeta = -100; }
    if(isnan(LJets[ijet].elinsubphi)) { LJets[ijet].elinsubphi = -100; }
    if(isnan(LJets[ijet].elinsubjpt)) { LJets[ijet].elinsubjpt = -100; }
    if(isnan(LJets[ijet].elinsubjeta)) { LJets[ijet].elinsubjeta = -100; }
    if(isnan(LJets[ijet].elinsubjphi)) { LJets[ijet].elinsubjphi = -100; }
    if(isnan(LJets[ijet].elinsubjmass)) { LJets[ijet].elinsubjmass = -100; }
    
    if(isnan(LJets[ijet].muinsubpt)) { LJets[ijet].muinsubpt = -100; }
    if(isnan(LJets[ijet].muinsubeta)) { LJets[ijet].muinsubeta = -100; }
    if(isnan(LJets[ijet].muinsubphi)) { LJets[ijet].muinsubphi = -100; }
    if(isnan(LJets[ijet].muinsubjpt)) { LJets[ijet].muinsubjpt = -100; }
    if(isnan(LJets[ijet].muinsubjeta)) { LJets[ijet].muinsubjeta = -100; }
    if(isnan(LJets[ijet].muinsubjphi)) { LJets[ijet].muinsubjphi = -100; }
    if(isnan(LJets[ijet].muinsubjmass)) { LJets[ijet].muinsubjmass = -100; }
    
    if(isnan(LJets[ijet].muinsubI0)) { LJets[ijet].muinsubI0 = -100; }
    if(isnan(LJets[ijet].muinsubInear)) { LJets[ijet].muinsubInear = -100; }
    if(isnan(LJets[ijet].muinsubIfar)) { LJets[ijet].muinsubIfar = -100; }
    
   
    if(isnan(LJets[ijet].sub1mass)) { LJets[ijet].sub1mass = -100; }
    if(isnan(LJets[ijet].sub1btag)) { LJets[ijet].sub1btag = -100; }
    if(isnan(LJets[ijet].sub1hadfrac)) { LJets[ijet].sub1hadfrac = -100; }
    if(isnan(LJets[ijet].sub1emfrac)) { LJets[ijet].sub1emfrac = -100; }
    if(isnan(LJets[ijet].sub2mass)) { LJets[ijet].sub2mass = -100; }
    if(isnan(LJets[ijet].sub2btag)) { LJets[ijet].sub2btag = -100; }
    if(isnan(LJets[ijet].sub2hadfrac)) { LJets[ijet].sub2hadfrac = -100; }
    if(isnan(LJets[ijet].sub2emfrac)) { LJets[ijet].sub2emfrac = -100; }

    if(isnan(LJets[ijet].subhaddiff)) { LJets[ijet].subhaddiff = -100; }
    if(isnan(LJets[ijet].subemdiff)) { LJets[ijet].subemdiff = -100; }
    if(isnan(LJets[ijet].subptdiff)) { LJets[ijet].subptdiff = -100; }
    if(isnan(LJets[ijet].subbtag)) { LJets[ijet].subbtag = -100; }

    if(isnan(LJets[ijet].haselectron)) { LJets[ijet].haselectron = -100; }
    if(isnan(LJets[ijet].haspfelectron)) { LJets[ijet].haspfelectron = -100; }
    if(isnan(LJets[ijet].hasmuon)) { LJets[ijet].hasmuon = -100; }
    if(isnan(LJets[ijet].haspfmuon)) { LJets[ijet].haspfmuon = -100; }

    if(isnan(LJets[ijet].hasleptop)) { LJets[ijet].hasleptop = -100; }
    if(isnan(LJets[ijet].hashadtop)) { LJets[ijet].hashadtop = -100; }
    if(isnan(LJets[ijet].hasqg)) { LJets[ijet].hasqg = -100; }
    if(isnan(LJets[ijet].hasb)) { LJets[ijet].hasb = -100; }
    
    
    in_pfjetAK8NHadF = -999;                                                     
    in_pfjetAK8neunhadfrac = -999;
    in_pfjetAK8subhaddiff = -999;
    in_pfjetAK8tau21 = -999;
    in_pfjetAK8chrad = -999;
    in_pfjetAK8sdmass = -999;
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
        
    if (LJets[ijet].muinsubeta != -100 && LJets[ijet].muinsubphi != -100)
      {
        if (nmuons > 0) {
          float dR_min(0.4); int nearestmu(-1);
          for(int mu=0; mu<nmuons; mu++){
            float dR = delta2R(muoneta[mu],muonphi[mu],LJets[ijet].muinsubeta,LJets[ijet].muinsubphi);
            if (dR < dR_min) {
              dR_min = dR;
              nearestmu = mu;
            }
          }
          if (nearestmu >= 0) {
	    
	    in_pfjetAK8NHadF = LJets[ijet].NHadF;
            in_pfjetAK8neunhadfrac = LJets[ijet].neunhadfrac;
            in_pfjetAK8subhaddiff = LJets[ijet].subhaddiff;
            in_pfjetAK8tau21 = LJets[ijet].tau21;
            in_pfjetAK8chrad = LJets[ijet].chrad;
            in_pfjetAK8sdmass = LJets[ijet].sdmass;
	    in_pfjetAK8matchedmuonchi = muonchi[nearestmu];
	    in_pfjetAK8matchedmuonposmatch = muonposmatch[nearestmu];
	    in_pfjetAK8matchedmuontrkink = muontrkink[nearestmu];
	    in_pfjetAK8matchedmuonsegcom = muonsegcom[nearestmu];
	    in_pfjetAK8matchedmuonhit = muonhit[nearestmu];
	    in_pfjetAK8matchedmuonmst = muonmst[nearestmu];
	    in_pfjetAK8matchedmuontrkvtx = muontrkvtx[nearestmu];
	    in_pfjetAK8matchedmuondz = muondz[nearestmu];
	    in_pfjetAK8matchedmuonpixhit = muonpixhit[nearestmu];
	    in_pfjetAK8matchedmuontrklay = muontrklay[nearestmu];
	    in_pfjetAK8matchedmuonvalfrac = muonvalfrac[nearestmu];
	    
	    TLorentzVector lep, lepbj, bj;
	    lep.SetPtEtaPhiE(LJets[ijet].muinsubpt,LJets[ijet].muinsubeta,LJets[ijet].muinsubphi,0.105658);
	    bj.SetPtEtaPhiM(LJets[ijet].muinsubjpt,LJets[ijet].muinsubjeta,LJets[ijet].muinsubjphi,LJets[ijet].muinsubjmass);
	    lepbj = lep + bj;
	    
	    in_pfjetAK8muinsubptrat = (lep.Pt())/(lepbj.Pt());
	    in_pfjetAK8muinsubmassrat = (bj.M())/(lepbj.M());
	    in_pfjetAK8muinsubinvmass = lepbj.M();
	    
	    if(isnan(in_pfjetAK8muinsubptrat)) {in_pfjetAK8muinsubptrat = -100;}
	    if(isnan(in_pfjetAK8muinsubmassrat)) {in_pfjetAK8muinsubmassrat = -100;}
	    if(isnan(in_pfjetAK8muinsubinvmass)) {in_pfjetAK8muinsubinvmass = -100;}
	    
	    in_pfjetAK8muinsubIfarbyI0 = LJets[ijet].muinsubIfar/LJets[ijet].muinsubI0;
	    in_pfjetAK8muinsubInearbyI0 = LJets[ijet].muinsubInear/LJets[ijet].muinsubI0;
	    
	    if(isnan(in_pfjetAK8muinsubIfarbyI0)) {in_pfjetAK8muinsubIfarbyI0 = -100;}
	    if(isnan(in_pfjetAK8muinsubInearbyI0)) {in_pfjetAK8muinsubInearbyI0 = -100;}
	  }
	  LJets[ijet].hasmatchmu = (nearestmu >= 0)?true:false;
	}
      }
    
    in_pfjetAK8eldxy_sv = -999;
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
    
    LJets[ijet].hasmatche = false;
    
    if (LJets[ijet].elinsubeta != -100 && LJets[ijet].elinsubphi != -100)
      {
        if (nelecs > 0) {
          float dR_min(0.4); int nearest(-1);
          for(int el=0; el<nelecs; el++){
            float dR = delta2R(eleta[el],elphi[el],LJets[ijet].elinsubeta,LJets[ijet].elinsubphi);
            if (dR < dR_min) {
              dR_min = dR;
              nearest = el;
            }
	  }
          if (nearest >= 0 && eldxy_sv[nearest] != 1000) {
	    in_pfjetAK8NHadF = LJets[ijet].NHadF;
	    in_pfjetAK8neunhadfrac = LJets[ijet].neunhadfrac;
	    in_pfjetAK8subhaddiff = LJets[ijet].subhaddiff;
	    in_pfjetAK8tau21 = LJets[ijet].tau21;
	    in_pfjetAK8chrad = LJets[ijet].chrad;
	    in_pfjetAK8sdmass = LJets[ijet].sdmass;
            in_pfjetAK8eldxy_sv = eldxy_sv[nearest];
	    in_pfjetAK8matchedelcleta = elsupcl_eta[nearest];
	    in_pfjetAK8matchedelpt = fabs(elpt[nearest]);
	    in_pfjetAK8matchedelsigmaieta = elsigmaieta[nearest];
            in_pfjetAK8matchedelsigmaiphi = elsigmaiphi[nearest];
	    in_pfjetAK8matchedelr9full = elr9full[nearest];
	    in_pfjetAK8matchedelsupcl_etaw = elsupcl_etaw[nearest];
            in_pfjetAK8matchedelsupcl_phiw = elsupcl_phiw[nearest];
	    in_pfjetAK8matchedelhcaloverecal = elhcaloverecal[nearest];
            in_pfjetAK8matchedelcloctftrkn = elcloctftrkn[nearest];
	    in_pfjetAK8matchedelcloctftrkchi2 = elcloctftrkchi2[nearest];
            in_pfjetAK8matchedele1x5bye5x5 = ele1x5bye5x5[nearest];
	    in_pfjetAK8matchedelnormchi2 = elnormchi2[nearest];
            in_pfjetAK8matchedelhitsmiss = elhitsmiss[nearest];
	    in_pfjetAK8matchedeltrkmeasure = eltrkmeasure[nearest];
	    in_pfjetAK8matchedelecloverpout = elecloverpout[nearest];
	    in_pfjetAK8matchedelecaletrkmomentum = elecaletrkmomentum[nearest];
	    in_pfjetAK8matchedeldeltaetacltrkcalo = eldeltaetacltrkcalo[nearest];
	    in_pfjetAK8matchedelsupcl_preshvsrawe = elsupcl_preshvsrawe[nearest];
	    in_pfjetAK8matchedelpfisolsumphet = elpfisolsumphet[nearest];
	    in_pfjetAK8matchedelpfisolsumchhadpt = elpfisolsumchhadpt[nearest];
	    in_pfjetAK8matchedelpfisolsumneuhadet = elpfsiolsumneuhadet[nearest];
	    in_pfjetAK8matchedeletain = eletain[nearest];
            in_pfjetAK8matchedelphiin = elphiin[nearest];
            in_pfjetAK8matchedelfbrem = elfbrem[nearest];
	    in_pfjetAK8matchedeleoverp = eleoverp[nearest];
            in_pfjetAK8matchedelhovere = elhovere[nearest];
	    
	  }
	  LJets[ijet].hasmatche = (nearest >= 0 && eldxy_sv[nearest] != 1000)?true:false;
	}
      }
    
    if (!isnan(Rho)) in_pfjetAK8matchedelRho = Rho;
    else in_pfjetAK8matchedelRho = -999;
    
    hist_2D_msd_deepak8->Fill(LJets[ijet].sdmass,LJets[ijet].DeepTag_TvsQCD);
    
    LJets[ijet].re_tvsb = reader1->EvaluateMVA("BDTG method");
    if(isnan(LJets[ijet].re_tvsb)) { LJets[ijet].re_tvsb = -100; }
    
    LJets[ijet].rmu_tvsb = reader4->EvaluateMVA("BDTG method");
    if(isnan(LJets[ijet].rmu_tvsb)) { LJets[ijet].rmu_tvsb = -100; }
    
  }//ijet
    
  vector <Electron> vfelectrons;

  for(unsigned ie=0; ie<velectrons.size(); ie++){
    if(velectrons[ie].pt<30.) continue; 
    vfelectrons.push_back(velectrons[ie]);
  }
  sorted_by_pt(vfelectrons);
  
  vector <Muon> vfmuons;
  for(unsigned imu=0; imu<vmuons.size(); imu++){
    if(vmuons[imu].pt<30.) continue;
    vfmuons.push_back(vmuons[imu]);
  }
  sorted_by_pt(vfmuons);
    
    
  int t_cand = -1;
  double remax = -200;
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
    if(LJets[ijet].re_tvsb > remax){
      remax = LJets[ijet].re_tvsb;
      t_cand = ijet;
    }
  }
  
  // top pt reweighting //
  
  if(isTT){
    
    TLorentzVector top4mom[2];
    
    int ngent = 0;
    
    for(int igen=0; igen<ngenparticles; igen++){
      
      if(abs(genpartstatus[igen])!=22) continue;
      if(!(genpartfromhard[igen])) continue;
      if(abs(genpartpdg[igen])!=6) continue;
      
      top4mom[ngent].SetPtEtaPhiM(genpartpt[igen],genparteta[igen],genpartphi[igen],genpartm[igen]);
      ngent++;
      if(ngent>=2) break;
      
    }
    
    float toppt_wt = 1;
    
    if(ngent==2){
      toppt_wt = SF_TOP(0.0615,0.0005,TMath::Min(float(500),float(top4mom[0].Pt())),TMath::Min(float(500),float(top4mom[1].Pt()))); 
    }
   
    //weight *= toppt_wt;
    // weight_puwup *= toppt_wt;
    //weight_puwdown *= toppt_wt;
    //weight_btagwup *= toppt_wt;
    //weight_btagwdown *= toppt_wt;
  }
  // top pt reweighting ends //
  
  if(isnan(weight) || weight>1.e+12) { weight = 0; }
  //if(isnan(weight_puwup) || weight_puwup>1.e+12) { weight_puwup = 0; }
  //if(isnan(weight_puwdown) || weight_puwdown>1.e+12) { weight_puwdown = 0; }
  //if(isnan(weight_btagwup) || weight_btagwup>1.e+12) { weight_btagwup = 0; }
  //if(isnan(weight_btagwdown) || weight_btagwdown>1.e+12) { weight_btagwdown = 0; }
  
  /******trigger object along with pdgid*****/
  std::vector<std::pair<int,TLorentzVector> > TrigRefObj;
  
  for (int tr=0; tr<ntrigobjs; tr++) {
    TLorentzVector trigobj;
    trigobj.SetPtEtaPhiM(trigobjpt[tr],trigobjeta[tr],trigobjphi[tr],trigobjmass[tr]);
    TrigRefObj.push_back(std::make_pair(trigobjpdgId[tr],trigobj));
  }
  
  // event selection starts
  
  if (nprimi<1) return kFALSE;  
  hist_count->Fill(1,weight);
  
  vector <Lepton> vleptons;

  for(unsigned imu=0; imu<vfmuons.size(); imu++){
    Lepton vlepton;
    vlepton.pt = vfmuons[imu].pt;
    vlepton.eta = vfmuons[imu].eta;
    vlepton.phi = vfmuons[imu].phi;
    vlepton.charge = vfmuons[imu].charge;
    vlepton.lepton_id = 1;
    vleptons.push_back(vlepton);
  }
  for(unsigned ie=0; ie<vfelectrons.size(); ie++){
    Lepton vlepton;
    vlepton.pt = vfelectrons[ie].pt;
    vlepton.eta = vfelectrons[ie].eta;
    vlepton.phi = vfelectrons[ie].phi;
    vlepton.charge = vfelectrons[ie].charge;
    vlepton.lepton_id = 2;
    vleptons.push_back(vlepton);
  }
  sorted_by_pt(vleptons);
  
  if (vleptons.size()<2) return kFALSE; //at least two leptons with pT > 30 GeV at this stage
    
  hist_count->Fill(2,weight);
  
  
  bool emu_ch = false;
  bool mumu_ch = false;
  bool ee_ch = false;
  
  if ((vleptons[0].lepton_id == 1 && vleptons[1].lepton_id == 2) || (vleptons[0].lepton_id == 2 && vleptons[1].lepton_id == 1)) emu_ch = true;
  else if (vleptons[0].lepton_id == 1 && vleptons[1].lepton_id == 1) mumu_ch = true;
  else if (vleptons[0].lepton_id == 2 && vleptons[1].lepton_id == 2) ee_ch =true;

  //Now comes trigger part//
  bool itrig_pass = false;
  bool itrigemu_pass = false;
  bool itrigmumu_pass = false;
  bool itrigee_pass = false;
  
  itrig_pass = ((hlt_Mu37Ele27==1)||(hlt_Mu27Ele37==1)||(hlt_Mu37TkMu27==1)||(hlt_DoubleEle25==1));
  if(!itrig_pass) return kFALSE; //event should at least fire a dileptonic trigger
  hist_count->Fill(3,weight);
  
#ifdef E_MU_TTBar
  if(!emu_ch) return kFALSE;
#elif defined(E_E_TTBar)
  if(!ee_ch) return kFALSE;
#elif defined(MU_MU_TTBar)
  if(!mumu_ch) return kFALSE;
#endif
  
  hist_count->Fill(4,weight);
  
#ifdef E_MU_TTBar
  itrigemu_pass = ((hlt_Mu37Ele27==1)||(hlt_Mu27Ele37==1));
  if(!itrigemu_pass) return kFALSE;
#elif defined(E_E_TTBar)
  itrigee_pass = (hlt_DoubleEle25==1);
  if(!itrigee_pass) return kFALSE;
#elif defined(MU_MU_TTBar)
  itrigmumu_pass = (hlt_Mu37TkMu27==1);
  if(!itrigmumu_pass) return kFALSE;
#endif
  
  hist_count->Fill(5,weight);
  
#ifdef E_MU_TTBar
  
  
  TLorentzVector fmucand, felcand;

  if (vleptons[0].lepton_id == 1 && vleptons[1].lepton_id == 2) {
    fmucand.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,0.105658);
    felcand.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,0.000511);
  }
  else if (vleptons[0].lepton_id == 2 && vleptons[1].lepton_id == 1) {
    fmucand.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,0.105658);
    felcand.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,0.000511);
  }
  
  /*perform trigger object matching for EMu final state*/
  bool fmuMatch = false; 
  bool felMatch = false;
  float matchN_mu(0.2), matchN_el(0.2);
  
  for (uint trv=0; trv<TrigRefObj.size(); trv++) {
    bool eltrobj(false), mutrobj(false), jettrobj(false);
    if (abs(TrigRefObj[trv].first)==13) mutrobj=true;
    else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.) jettrobj=true;
    else if (abs(TrigRefObj[trv].first)==0) eltrobj=true;
    
    if (mutrobj==true) {
      TVector3 Trvmu(TrigRefObj[trv].second.Pt(),TrigRefObj[trv].second.Eta(),TrigRefObj[trv].second.Phi());
      TVector3 fmuv(fmucand.Pt(),fmucand.Eta(),fmucand.Phi());            
      if ((fmuv-Trvmu).Mag()/fmuv.Mag() < matchN_mu) {
	matchN_mu = (fmuv-Trvmu).Mag()/fmuv.Mag();
	fmuMatch = true;
	matchN_mu = (fmuv-Trvmu).Mag()/fmuv.Mag();
      }
    }
    if (eltrobj==true) {
      TVector3 Trvel(TrigRefObj[trv].second.Pt(),TrigRefObj[trv].second.Eta(),TrigRefObj[trv].second.Phi());
      TVector3 felv(felcand.Pt(),felcand.Eta(),felcand.Phi());
      if ((felv-Trvel).Mag()/felv.Mag() < matchN_el) {
        matchN_el = (felv-Trvel).Mag()/felv.Mag();
        felMatch = true;
	matchN_el = (felv-Trvel).Mag()/felv.Mag();
      } 
    } 
  }
  if (!(fmuMatch==true && felMatch==true)) return kFALSE;
  
  //TString str;                                                                 
  //str = TString::Format("fmuMatch %u felMatch %u ievt %u",fmuMatch,felMatch,ievt);
  //if(gProofServ) gProofServ->SendAsynMessage(str);
  
  float mulpt(0.), ellpt(0.);
  bool fsttrig = false; 
  bool sndtrig = false;
  bool bothtrig = false;
  
  if (hlt_Mu37Ele27==1 && hlt_Mu27Ele37==0) {
    mulpt=40.0; ellpt=30.0;
    fsttrig=true;
  }
  else if (hlt_Mu37Ele27==0 && hlt_Mu27Ele37==1) {
    mulpt=30.0; ellpt=40.0;
    sndtrig=true;
  }
  else if (hlt_Mu37Ele27==1 && hlt_Mu27Ele37==1) {
    mulpt=30.0; ellpt=30.0;
    bothtrig=true;
  }
  if (fsttrig==1 || sndtrig==1 || bothtrig==1) {if(!(fmucand.Pt()>mulpt && felcand.Pt()>ellpt)) return kFALSE;}

  /*trigger checks for EE final state*/
#elif defined(E_E_TTBar)
  float el1pt(0.), el2pt(0.);
  if (itrigee_pass==1) {el1pt = 40.0; el2pt = 30.0;}
  
  TLorentzVector fe1cand, fe2cand;
  if (vleptons[0].lepton_id == 2 && vleptons[1].lepton_id == 2) {
    fe1cand.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,0.000511);
    fe2cand.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,0.000511);
  }
  
  /*perform trigger object matching for EE final state*/
  bool fe1Match = false;
  bool fe2Match = false;
  float matchN_e1(0.2), matchN_e2(0.2);
  
  for (uint trv=0; trv<TrigRefObj.size(); trv++) {
    bool eltrobj(false), jettrobj(false);
    if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.) jettrobj=true;
    else if (abs(TrigRefObj[trv].first)==0) eltrobj=true;

    if (eltrobj==true) {
      TVector3 Trvel(TrigRefObj[trv].second.Pt(),TrigRefObj[trv].second.Eta(),TrigRefObj[trv].second.Phi());
      TVector3 fe1v(fe1cand.Pt(),fe1cand.Eta(),fe1cand.Phi());
      TVector3 fe2v(fe2cand.Pt(),fe2cand.Eta(),fe2cand.Phi());
      if ((fe1v-Trvel).Mag()/fe1v.Mag() < matchN_e1) {
        matchN_e1 = (fe1v-Trvel).Mag()/fe1v.Mag();
        fe1Match = true;
        matchN_e1 = (fe1v-Trvel).Mag()/fe1v.Mag();
      }
      if ((fe2v-Trvel).Mag()/fe2v.Mag() < matchN_e2) {
        matchN_e2 = (fe2v-Trvel).Mag()/fe2v.Mag();
        fe2Match = true;
        matchN_e2 = (fe2v-Trvel).Mag()/fe2v.Mag();
      }
    }
  }
  
  //TString stre;
  //stre = TString::Format("fe1Match %u fe2Match %u ievt %u",fe1Match,fe2Match,ievt);
  //if(gProofServ) gProofServ->SendAsynMessage(stre);

  if (!(fe1Match==true && fe2Match==true)) return kFALSE;
  if (itrigee_pass==1) {if(!(fe1cand.Pt()>el1pt && fe2cand.Pt()>el2pt)) return kFALSE;}
  
  /*trigger checks for MuMu final state*/
#elif defined(MU_MU_TTBar)
  float mu1pt(0.), mu2pt(0.);
  if (itrigmumu_pass==1) {mu1pt = 40.0; mu2pt = 30.0;}
  
  TLorentzVector fmu1cand, fmu2cand;
  if (vleptons[0].lepton_id == 1 && vleptons[1].lepton_id == 1) {
    fmu1cand.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,0.105658);
    fmu2cand.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,0.105658);
  }
  /*perform trigger object matching for MuMu final state*/
  bool fmu1Match = false;
  bool fmu2Match = false;
  float matchN_mu1(0.2), matchN_mu2(0.2);

  for (uint trv=0; trv<TrigRefObj.size(); trv++) {
    bool mutrobj(false);
    if (abs(TrigRefObj[trv].first)==13) mutrobj=true;
    if (mutrobj==true) {
      TVector3 Trvmu(TrigRefObj[trv].second.Pt(),TrigRefObj[trv].second.Eta(),TrigRefObj[trv].second.Phi());
      TVector3 fmu1v(fmu1cand.Pt(),fmu1cand.Eta(),fmu1cand.Phi());
      TVector3 fmu2v(fmu2cand.Pt(),fmu2cand.Eta(),fmu2cand.Phi());
      if ((fmu1v-Trvmu).Mag()/fmu1v.Mag() < matchN_mu1) {
        matchN_mu1 = (fmu1v-Trvmu).Mag()/fmu1v.Mag();
        fmu1Match = true;
        matchN_mu1 = (fmu1v-Trvmu).Mag()/fmu1v.Mag();
      }
      if ((fmu2v-Trvmu).Mag()/fmu2v.Mag() < matchN_mu2) {
        matchN_mu2 = (fmu2v-Trvmu).Mag()/fmu2v.Mag();
        fmu2Match = true;
        matchN_mu2 = (fmu2v-Trvmu).Mag()/fmu2v.Mag();
      }
    }
  }
  if (!(fmu1Match==true && fmu2Match==true)) return kFALSE;

  //TString strmu;
  //strmu = TString::Format("fmu1Match %u fmu2Match %u ievt %u",fmu1Match,fmu2Match,ievt);
  //if(gProofServ) gProofServ->SendAsynMessage(strmu);
  
  if (itrigmumu_pass==1) {if(!(fmu1cand.Pt()>mu1pt && fmu2cand.Pt()>mu2pt)) return kFALSE;}
#endif
  
  hist_count->Fill(6,weight);
    
  /***** leptons should be oppositely charged****/
  if ((vleptons[0].charge*vleptons[1].charge)>0) return kFALSE;
  hist_count->Fill(7,weight);
  
  /***no other 3rd lepton other than the selected lepton set***/ 
#ifdef E_MU_TTBar  
  if (nmuons>1) return kFALSE;
  hist_count->Fill(8,weight);
  if (nelecs>1) return kFALSE;
  hist_count->Fill(9,weight);
#elif defined(E_E_TTBar)
  if (nmuons>=1) return kFALSE;
  hist_count->Fill(8,weight);
  if (nelecs>2) return kFALSE;
  hist_count->Fill(9,weight);
#elif defined(MU_MU_TTBar)
  if (nmuons>2) return kFALSE;
  hist_count->Fill(8,weight);
  if (nelecs>=1) return kFALSE;
  hist_count->Fill(9,weight);
#endif
  
  //Computation of lepton related Suman's variables//                                               
             
  TLorentzVector l1, l2;
#ifdef E_MU_TTBar
  l1 = (fmucand.Pt()>felcand.Pt()) ? fmucand : felcand;
  l2 = (fmucand.Pt()>felcand.Pt()) ? felcand : fmucand;
#elif defined(E_E_TTBar)
  l1 = fe1cand;
  l2 = fe2cand;
#elif defined(MU_MU_TTBar)
  l1 = fmu1cand;
  l2 = fmu2cand;
#endif

  if (npfjetAK4<2) return kFALSE;
  hist_count->Fill(10,weight);
  
  if (nbjetAK4<1) return kFALSE;
  hist_count->Fill(11,weight);
  
  if (npfjetAK8<1) return kFALSE;
  hist_count->Fill(12,weight);
  
  M_l1l2 = (l1+l2).M();
  //hist_new_var[0]->Fill(M_l1l2,weight);

  /***** invariant mass of lepton at least more than 20 GeV as resolved analysis cut****/
  if (M_l1l2 < 20.) return kFALSE;
  hist_count->Fill(13,weight);
  
  //Computation of selected lepton related variables (Suman's proposal)
  rat_l1pt_l2pt = l1.Pt()/l2.Pt();
  //hist_new_var[1]->Fill(rat_l1pt_l2pt,weight);
  
  deltaPhi_l1l2 = PhiInRange(l1.Phi() - l2.Phi());
  //hist_new_var[2]->Fill(deltaPhi_l1l2,weight);

  //2d iso variables for l1 and l2//                                             
  float dRl1_min(1000), dRl2_min(1000);
  int nearjet_l1(-1), nearjet_l2(-1);
  
  for(unsigned kjet=0; kjet<Jets.size(); kjet++){  
    if(delta2R(Jets[kjet].eta,Jets[kjet].phi,l1.Eta(),l1.Phi()) < dRl1_min){
      dRl1_min = delta2R(Jets[kjet].eta,Jets[kjet].phi,l1.Eta(),l1.Phi());
      nearjet_l1 = kjet;
    }
  }
  
  if(nearjet_l1>=0){
    TLorentzVector j_mom; j_mom.SetPtEtaPhiM(Jets[nearjet_l1].pt,Jets[nearjet_l1].eta,Jets[nearjet_l1].phi,Jets[nearjet_l1].mass);
    M_jl1 = (j_mom + l1).M();
    l1pt_nearjet = ((l1.Vect()).Perp(j_mom.Vect()));
  }
  //hist_new_var[3]->Fill(l1pt_nearjet,weight);

  for(unsigned kjet=0; kjet<Jets.size(); kjet++){
    if(delta2R(Jets[kjet].eta,Jets[kjet].phi,l2.Eta(),l2.Phi()) < dRl2_min){
      dRl2_min = delta2R(Jets[kjet].eta,Jets[kjet].phi,l2.Eta(),l2.Phi());
      nearjet_l2 = kjet;
    }
  }

  if(nearjet_l2>=0){
    TLorentzVector j_mom; j_mom.SetPtEtaPhiM(Jets[nearjet_l2].pt,Jets[nearjet_l2].eta,Jets[nearjet_l2].phi,Jets[nearjet_l2].mass);
    M_jl2 = (j_mom + l2).M();
    l2pt_nearjet = ((l2.Vect()).Perp(j_mom.Vect()));
  }
  //hist_new_var[4]->Fill(l2pt_nearjet,weight);

  //Computation of MET related Suman's proposed variables//
  if (PFMET != -1000 && PFMETPhi != -1000) { 
    float metx = PFMET*std::cos(PFMETPhi);
    float mety = PFMET*std::sin(PFMETPhi);
    TLorentzVector metvector;
    
    metvector.SetPxPyPzE(metx,mety,0,PFMET); //as mass and Pz components are 0, thus E = Pt  
    met_pt = PFMET;
    //hist_new_var[5]->Fill(met_pt,weight);

    met_eta = metvector.Eta(); //always be 0.
    //hist_new_var[6]->Fill(met_eta,weight);

    delta_phil1_met = PhiInRange(l1.Phi() - metvector.Phi());
    //hist_new_var[7]->Fill(delta_phil1_met,weight);

    delta_phil2_met = PhiInRange(l2.Phi() - metvector.Phi());
    //hist_new_var[8]->Fill(delta_phil2_met,weight);

    //nearest bjet of l1//
    float dRbl1_min = 1000;
    int nearbl1jet = -1;
    
    for(unsigned kjet=0; kjet<Jets.size(); kjet++){
      if(Jets[kjet].btag_DeepFlav>deep_btag_cut) {
	if(delta2R(Jets[kjet].eta,Jets[kjet].phi,l1.Eta(),l1.Phi()) < dRbl1_min){
	  dRbl1_min = delta2R(Jets[kjet].eta,Jets[kjet].phi,l1.Eta(),l1.Phi()) ;
	  nearbl1jet = kjet;
	}
      }
    }
  
    TLorentzVector nearbl1;
    TLorentzVector bl1_syst;
    if (nearbl1jet>=0) {
      nearbl1.SetPtEtaPhiM(Jets[nearbl1jet].pt,Jets[nearbl1jet].eta,Jets[nearbl1jet].phi,Jets[nearbl1jet].mass);
      bl1_syst = nearbl1 + l1;
      M_bl1 = bl1_syst.M();
      delta_phibl1_met = PhiInRange(bl1_syst.Phi() - metvector.Phi());
    }
    //Similarly do it for l2 and nearest bjet//
    float dRbl2_min = 1000;
    int nearbl2jet = -1;
    
    for(unsigned kjet=0; kjet<Jets.size(); kjet++){
      if(Jets[kjet].btag_DeepFlav>deep_btag_cut) {
        if(delta2R(Jets[kjet].eta,Jets[kjet].phi,l2.Eta(),l2.Phi()) < dRbl2_min){
          dRbl2_min = delta2R(Jets[kjet].eta,Jets[kjet].phi,l2.Eta(),l2.Phi()) ;
          nearbl2jet = kjet;
        }
      }
    }
    TLorentzVector nearbl2;
    TLorentzVector bl2_syst;
    if (nearbl2jet>=0) {
      nearbl2.SetPtEtaPhiM(Jets[nearbl2jet].pt,Jets[nearbl2jet].eta,Jets[nearbl2jet].phi,Jets[nearbl2jet].mass);
      bl2_syst = nearbl2 + l2;
      M_bl2 = bl2_syst.M();
      delta_phibl2_met = PhiInRange(bl2_syst.Phi() - metvector.Phi());
    }
    if (nearbl1jet>=0 && nearbl2jet>=0) {
      delta_phibl1bl2 = PhiInRange(bl1_syst.Phi() - bl2_syst.Phi());
    }
    rat_metpt_ak4pt = metvector.Pt()/Jets[0].pt;
    rat_metpt_ak8pt = metvector.Pt()/LJets[0].pt;
    rat_metpt_eventHT = metvector.Pt()/Event_HT;
    
    mt_of_l1met = (metvector+l1).Mt();
    mt_of_l2met = (metvector+l2).Mt();
      
    //TString str;
    //str = TString::Format("met_pt %f std::sqrt(metx*metx+ mety*mety) %f met_phi %f met_eta %f",met_pt,std::sqrt(metx*metx+ mety*mety),PFMETPhi,met_eta);            
    //if(gProofServ) gProofServ->SendAsynMessage(str);  
  }

  no_ak4jets = npfjetAK4;
  no_ak4bjets = nbjetAK4;
  no_ak8jets = npfjetAK8;

  extra_ak4j = npfjetAK4-2;

  if (extra_ak4j>=1) {

    int leadptcand(-1);
    float leadptval(-99.);
    float ptsum(0.);

    for(unsigned ijet=0; ijet<Jets.size(); ijet++){
      if (ijet == 0 || ijet == 1) continue;
      if (Jets[ijet].pt>leadptval) {
	leadptval = Jets[ijet].pt;
	leadptcand = ijet;
      }
      ptsum = ptsum + Jets[ijet].pt;
    }
    ptsum_extra_ak4 = ptsum;
  
    if (leadptcand>=0) {
      extra_ak4jqgl = Jets[leadptcand].qgl;
      extra_ak4jdeepb = Jets[leadptcand].btag_DeepFlav;
    
      if(delta2R(Jets[leadptcand].eta,Jets[leadptcand].phi,l1.Eta(),l1.Phi())<delta2R(Jets[leadptcand].eta,Jets[leadptcand].phi,l2.Eta(),l2.Phi())){
	rat_extra_ak4jpt_lpt = Jets[leadptcand].pt/l1.Pt();
      }
      else if (delta2R(Jets[leadptcand].eta,Jets[leadptcand].phi,l1.Eta(),l1.Phi())>delta2R(Jets[leadptcand].eta,Jets[leadptcand].phi,l2.Eta(),l2.Phi())){
	rat_extra_ak4jpt_lpt = Jets[leadptcand].pt/l2.Pt();
      }
    }
  }

  EventHT = Event_HT;
 
  ak81pt = LJets[0].pt;
  ak81y = LJets[0].y;
  ak81mass = LJets[0].mass;
  ak81sdmass = LJets[0].sdmass;
  ak81deep_tvsqcd = LJets[0].DeepTag_TvsQCD;
  ak81deep_wvsqcd = LJets[0].DeepTag_WvsQCD;

  if (npfjetAK8>1) {
    ak82pt = LJets[1].pt;
    ak82y = LJets[1].y;
    ak82mass = LJets[1].mass;
    ak82sdmass = LJets[1].sdmass;
    ak82deep_tvsqcd = LJets[1].DeepTag_TvsQCD;
    ak82deep_wvsqcd = LJets[1].DeepTag_WvsQCD;
  }

  deltaR_l1l2 = delta2R(l1.Eta(),l1.Phi(),l2.Eta(),l2.Phi());
  deltaR_l1b1 = delta2R(l1.Eta(),l1.Phi(),bjv[0].Eta(),bjv[0].Phi());
  if (bjv.size() >1) deltaR_l1b2 = delta2R(l1.Eta(),l1.Phi(),bjv[1].Eta(),bjv[1].Phi());
  deltaR_l2b1 = delta2R(l2.Eta(),l2.Phi(),bjv[0].Eta(),bjv[0].Phi());
  if (bjv.size() >1) deltaR_l2b2 = delta2R(l2.Eta(),l2.Phi(),bjv[1].Eta(),bjv[1].Phi());
  
  deltaR_l1j1 = delta2R(l1.Eta(),l1.Phi(),Jets[0].eta,Jets[0].phi);
  deltaR_l1j2 = delta2R(l1.Eta(),l1.Phi(),Jets[1].eta,Jets[1].phi);
  deltaR_l2j1 = delta2R(l2.Eta(),l2.Phi(),Jets[0].eta,Jets[0].phi);
  deltaR_l2j2 = delta2R(l2.Eta(),l2.Phi(),Jets[1].eta,Jets[1].phi);

  Tnewvar->Fill();
  
  if (PFMET<50.) return kFALSE; //before it was 70 GeV. 100 GeV on 17th April, 50GeV on June  
  hist_count->Fill(14,weight);
  
  hist_init[0]->Fill(nmuons,weight);
  hist_init[1]->Fill(nelecs,weight);
  hist_init[2]->Fill(PFMET,weight);
  hist_init[3]->Fill(nprimi,weight);
  hist_init[4]->Fill(npfjetAK4,weight);
  hist_init[5]->Fill(nbjetAK4,weight);
  hist_init[6]->Fill(npfjetAK8,weight);

  hist_obs[0]->Fill(LJets[0].pt,weight);
  hist_obs[1]->Fill(LJets[0].y,weight);
  hist_obs[2]->Fill(LJets[0].mass,weight);
  hist_obs[3]->Fill(LJets[0].NHadF,weight);
  hist_obs[4]->Fill(LJets[0].neunhadfrac,weight);
  hist_obs[5]->Fill(LJets[0].sdmass,weight);
  hist_obs[6]->Fill(LJets[0].chrad,weight);
  hist_obs[7]->Fill(LJets[0].subhaddiff,weight);
  
  //TString str;                                                                 
  //  str = TString::Format("NHadF %f neunhadfrac %f subhaddiff %f subhaddiff %f",LJets[0].NHadF,LJets[0].neunhadfrac,LJets[0].subhaddiff,diff_func(LJets[0].sub1hadfrac,LJets[0].sub2hadfrac));
  //if(gProofServ) gProofServ->SendAsynMessage(str);

  hist_obs[8]->Fill(LJets[0].tau21,weight);
  hist_obs[9]->Fill(LJets[0].DeepTag_TvsQCD,weight);
  hist_obs[10]->Fill(LJets[0].DeepTag_WvsQCD,weight);
  hist_obs[11]->Fill(LJets[0].DeepTag_ZvsQCD,weight);
  hist_obs[12]->Fill(LJets[0].re_tvsb,weight);
  hist_obs[13]->Fill(LJets[0].rmu_tvsb,weight);
  hist_obs[14]->Fill(LJets[0].haspfelectron,weight);
  hist_obs[15]->Fill(LJets[0].haspfmuon,weight);
  hist_obs[16]->Fill(LJets[0].hasmatche,weight);
  hist_obs[17]->Fill(LJets[0].hasmatchmu,weight);
  hist_obs[18]->Fill(delta2R(LJets[0].eta,LJets[0].phi,l1.Eta(),l1.Phi()),weight);
  hist_obs[19]->Fill(delta2R(LJets[0].eta,LJets[0].phi,l2.Eta(),l2.Phi()),weight);
  hist_obs[20]->Fill(delta2R(LJets[0].eta,LJets[0].phi,bjv[0].Eta(),bjv[0].Phi()),weight);
  if (nbjetAK4>1) hist_obs[21]->Fill(delta2R(LJets[0].eta,LJets[0].phi,bjv[1].Eta(),bjv[1].Phi()),weight); 
  if (npfjetAK8>1) { 
    hist_obs[22]->Fill(delta2R(LJets[0].eta,LJets[0].phi,LJets[1].eta,LJets[1].phi),weight);
    hist_obs[23]->Fill(LJets[1].pt,weight);
    hist_obs[24]->Fill(LJets[1].y,weight);
    hist_obs[25]->Fill(LJets[1].mass,weight);
    hist_obs[26]->Fill(LJets[1].NHadF,weight);
    hist_obs[27]->Fill(LJets[1].neunhadfrac,weight);
    hist_obs[28]->Fill(LJets[1].sdmass,weight);
    hist_obs[29]->Fill(LJets[1].chrad,weight);
    hist_obs[30]->Fill(LJets[1].subhaddiff,weight);
    hist_obs[31]->Fill(LJets[1].tau21,weight);
    hist_obs[32]->Fill(LJets[1].DeepTag_TvsQCD,weight);
    hist_obs[33]->Fill(LJets[1].DeepTag_WvsQCD,weight);
    hist_obs[34]->Fill(LJets[1].DeepTag_ZvsQCD,weight);
    hist_obs[35]->Fill(LJets[1].re_tvsb,weight);
    hist_obs[36]->Fill(LJets[1].rmu_tvsb,weight);
    hist_obs[37]->Fill(LJets[1].haspfelectron,weight);
    hist_obs[38]->Fill(LJets[1].haspfmuon,weight);
    hist_obs[39]->Fill(LJets[1].hasmatche,weight);
    hist_obs[40]->Fill(LJets[1].hasmatchmu,weight);
    hist_obs[41]->Fill(delta2R(LJets[1].eta,LJets[1].phi,l1.Eta(),l1.Phi()),weight);
    hist_obs[42]->Fill(delta2R(LJets[1].eta,LJets[1].phi,l2.Eta(),l2.Phi()),weight);
    hist_obs[43]->Fill(delta2R(LJets[1].eta,LJets[1].phi,bjv[0].Eta(),bjv[0].Phi()),weight);
    if (nbjetAK4>1) {

      TString str;                                                                                
      str = TString::Format("delta2R %f",delta2R(LJets[1].eta,LJets[1].phi,bjv[1].Eta(),bjv[1].Phi()));
      if(gProofServ) gProofServ->SendAsynMessage(str);

      hist_obs[44]->Fill(delta2R(LJets[1].eta,LJets[1].phi,bjv[1].Eta(),bjv[1].Phi()),weight);
    }
  }

#ifdef E_MU_TTBar
  hist_init[7]->Fill((fmucand + felcand).M(),weight);
  hist_init[8]->Fill(fmucand.Pt(),weight);
  hist_init[9]->Fill(fmucand.Eta(),weight);
  hist_init[10]->Fill(fmucand.Phi(),weight);
  
  hist_init[11]->Fill(felcand.Pt(),weight);
  hist_init[12]->Fill(felcand.Eta(),weight);
  hist_init[13]->Fill(felcand.Phi(),weight);
#elif defined(E_E_TTBar)
  hist_init[7]->Fill((fe1cand + fe2cand).M(),weight);
  hist_init[8]->Fill(fe1cand.Pt(),weight);
  hist_init[9]->Fill(fe1cand.Eta(),weight);
  hist_init[10]->Fill(fe1cand.Phi(),weight);

  hist_init[11]->Fill(fe2cand.Pt(),weight);
  hist_init[12]->Fill(fe2cand.Eta(),weight);
  hist_init[13]->Fill(fe2cand.Phi(),weight);
#elif defined(MU_MU_TTBar)
  hist_init[7]->Fill((fmu1cand + fmu2cand).M(),weight);
  hist_init[8]->Fill(fmu1cand.Pt(),weight);
  hist_init[9]->Fill(fmu1cand.Eta(),weight);
  hist_init[10]->Fill(fmu1cand.Phi(),weight);

  hist_init[11]->Fill(fmu2cand.Pt(),weight);
  hist_init[12]->Fill(fmu2cand.Eta(),weight);
  hist_init[13]->Fill(fmu2cand.Phi(),weight);
#endif
  hist_init[14]->Fill(bjv[0].Pt(),weight);
  hist_init[15]->Fill(bjv[0].Eta(),weight);
  hist_init[16]->Fill(bjv[0].Phi(),weight);
  /*
  hist_initpuwup[0]->Fill(nmuons,weight_puwup);
  hist_initpuwup[1]->Fill(nelecs,weight_puwup);
  hist_initpuwup[2]->Fill(PFMET,weight_puwup);
  hist_initpuwup[3]->Fill(nprimi,weight_puwup);
  hist_initpuwup[4]->Fill(npfjetAK4,weight_puwup);
  hist_initpuwup[5]->Fill(nbjetAK4,weight_puwup);
  hist_initpuwup[6]->Fill(npfjetAK8,weight_puwup);
  
#ifdef E_MU_TTBar
  hist_initpuwup[7]->Fill((fmucand + felcand).M(),weight_puwup);
  hist_initpuwup[8]->Fill(fmucand.Pt(),weight_puwup);
  hist_initpuwup[9]->Fill(fmucand.Eta(),weight_puwup);
  hist_initpuwup[10]->Fill(fmucand.Phi(),weight_puwup);

  hist_initpuwup[11]->Fill(felcand.Pt(),weight_puwup);
  hist_initpuwup[12]->Fill(felcand.Eta(),weight_puwup);
  hist_initpuwup[13]->Fill(felcand.Phi(),weight_puwup);
#elif defined(E_E_TTBar)
  hist_initpuwup[7]->Fill((fe1cand + fe2cand).M(),weight_puwup);
  hist_initpuwup[8]->Fill(fe1cand.Pt(),weight_puwup);
  hist_initpuwup[9]->Fill(fe1cand.Eta(),weight_puwup);
  hist_initpuwup[10]->Fill(fe1cand.Phi(),weight_puwup);

  hist_initpuwup[11]->Fill(fe2cand.Pt(),weight_puwup);
  hist_initpuwup[12]->Fill(fe2cand.Eta(),weight_puwup);
  hist_initpuwup[13]->Fill(fe2cand.Phi(),weight_puwup);
#elif defined(MU_MU_TTBar)
  hist_initpuwup[7]->Fill((fmu1cand + fmu2cand).M(),weight_puwup);
  hist_initpuwup[8]->Fill(fmu1cand.Pt(),weight_puwup);
  hist_initpuwup[9]->Fill(fmu1cand.Eta(),weight_puwup);
  hist_initpuwup[10]->Fill(fmu1cand.Phi(),weight_puwup);

  hist_initpuwup[11]->Fill(fmu2cand.Pt(),weight_puwup);
  hist_initpuwup[12]->Fill(fmu2cand.Eta(),weight_puwup);
  hist_initpuwup[13]->Fill(fmu2cand.Phi(),weight_puwup);
#endif
  hist_initpuwup[14]->Fill(bjv[0].Pt(),weight_puwup);
  hist_initpuwup[15]->Fill(bjv[0].Eta(),weight_puwup);
  hist_initpuwup[16]->Fill(bjv[0].Phi(),weight_puwup);

  hist_initbtagwup[0]->Fill(nmuons,weight_btagwup);
  hist_initbtagwup[1]->Fill(nelecs,weight_btagwup);
  hist_initbtagwup[2]->Fill(PFMET,weight_btagwup);
  hist_initbtagwup[3]->Fill(nprimi,weight_btagwup);
  hist_initbtagwup[4]->Fill(npfjetAK4,weight_btagwup);
  hist_initbtagwup[5]->Fill(nbjetAK4,weight_btagwup);
  hist_initbtagwup[6]->Fill(npfjetAK8,weight_btagwup);
#ifdef E_MU_TTBar
  hist_initbtagwup[7]->Fill((fmucand + felcand).M(),weight_btagwup);
  hist_initbtagwup[8]->Fill(fmucand.Pt(),weight_btagwup);
  hist_initbtagwup[9]->Fill(fmucand.Eta(),weight_btagwup);
  hist_initbtagwup[10]->Fill(fmucand.Phi(),weight_btagwup);
  hist_initbtagwup[11]->Fill(felcand.Pt(),weight_btagwup);
  hist_initbtagwup[12]->Fill(felcand.Eta(),weight_btagwup);
  hist_initbtagwup[13]->Fill(felcand.Phi(),weight_btagwup);
#elif defined(E_E_TTBar)
  hist_initbtagwup[7]->Fill((fe1cand + fe2cand).M(),weight_btagwup);
  hist_initbtagwup[8]->Fill(fe1cand.Pt(),weight_btagwup);
  hist_initbtagwup[9]->Fill(fe1cand.Eta(),weight_btagwup);
  hist_initbtagwup[10]->Fill(fe1cand.Phi(),weight_btagwup);
  hist_initbtagwup[11]->Fill(fe2cand.Pt(),weight_btagwup);
  hist_initbtagwup[12]->Fill(fe2cand.Eta(),weight_btagwup);
  hist_initbtagwup[13]->Fill(fe2cand.Phi(),weight_btagwup);
#elif defined(MU_MU_TTBar)
  hist_initbtagwup[7]->Fill((fmu1cand + fmu2cand).M(),weight_btagwup);
  hist_initbtagwup[8]->Fill(fmu1cand.Pt(),weight_btagwup);
  hist_initbtagwup[9]->Fill(fmu1cand.Eta(),weight_btagwup);
  hist_initbtagwup[10]->Fill(fmu1cand.Phi(),weight_btagwup);
  hist_initbtagwup[11]->Fill(fmu2cand.Pt(),weight_btagwup);
  hist_initbtagwup[12]->Fill(fmu2cand.Eta(),weight_btagwup);
  hist_initbtagwup[13]->Fill(fmu2cand.Phi(),weight_btagwup);
#endif
  hist_initbtagwup[14]->Fill(bjv[0].Pt(),weight_btagwup);
  hist_initbtagwup[15]->Fill(bjv[0].Eta(),weight_btagwup);
  hist_initbtagwup[16]->Fill(bjv[0].Phi(),weight_btagwup);
  
  hist_initpuwdown[0]->Fill(nmuons,weight_puwdown);
  hist_initpuwdown[1]->Fill(nelecs,weight_puwdown);
  hist_initpuwdown[2]->Fill(PFMET,weight_puwdown);
  hist_initpuwdown[3]->Fill(nprimi,weight_puwdown);
  hist_initpuwdown[4]->Fill(npfjetAK4,weight_puwdown);
  hist_initpuwdown[5]->Fill(nbjetAK4,weight_puwdown);
  hist_initpuwdown[6]->Fill(npfjetAK8,weight_puwdown);
#ifdef E_MU_TTBar
  hist_initpuwdown[7]->Fill((fmucand + felcand).M(),weight_puwdown);
  hist_initpuwdown[8]->Fill(fmucand.Pt(),weight_puwdown);
  hist_initpuwdown[9]->Fill(fmucand.Eta(),weight_puwdown);
  hist_initpuwdown[10]->Fill(fmucand.Phi(),weight_puwdown);
  hist_initpuwdown[11]->Fill(felcand.Pt(),weight_puwdown);
  hist_initpuwdown[12]->Fill(felcand.Eta(),weight_puwdown);
  hist_initpuwdown[13]->Fill(felcand.Phi(),weight_puwdown);
#elif defined(E_E_TTBar)
  hist_initpuwdown[7]->Fill((fe1cand + fe2cand).M(),weight_puwdown);
  hist_initpuwdown[8]->Fill(fe1cand.Pt(),weight_puwdown);
  hist_initpuwdown[9]->Fill(fe1cand.Eta(),weight_puwdown);
  hist_initpuwdown[10]->Fill(fe1cand.Phi(),weight_puwdown);
  hist_initpuwdown[11]->Fill(fe2cand.Pt(),weight_puwdown);
  hist_initpuwdown[12]->Fill(fe2cand.Eta(),weight_puwdown);
  hist_initpuwdown[13]->Fill(fe2cand.Phi(),weight_puwdown);
#elif defined(MU_MU_TTBar)
  hist_initpuwdown[7]->Fill((fmu1cand + fmu2cand).M(),weight_puwdown);
  hist_initpuwdown[8]->Fill(fmu1cand.Pt(),weight_puwdown);
  hist_initpuwdown[9]->Fill(fmu1cand.Eta(),weight_puwdown);
  hist_initpuwdown[10]->Fill(fmu1cand.Phi(),weight_puwdown);
  hist_initpuwdown[11]->Fill(fmu2cand.Pt(),weight_puwdown);
  hist_initpuwdown[12]->Fill(fmu2cand.Eta(),weight_puwdown);
  hist_initpuwdown[13]->Fill(fmu2cand.Phi(),weight_puwdown);
#endif
  hist_initpuwdown[14]->Fill(bjv[0].Pt(),weight_puwdown);
  hist_initpuwdown[15]->Fill(bjv[0].Eta(),weight_puwdown);
  hist_initpuwdown[16]->Fill(bjv[0].Phi(),weight_puwdown);

  hist_initbtagwdown[0]->Fill(nmuons,weight_btagwdown);
  hist_initbtagwdown[1]->Fill(nelecs,weight_btagwdown);
  hist_initbtagwdown[2]->Fill(PFMET,weight_btagwdown);
  hist_initbtagwdown[3]->Fill(nprimi,weight_btagwdown);
  hist_initbtagwdown[4]->Fill(npfjetAK4,weight_btagwdown);
  hist_initbtagwdown[5]->Fill(nbjetAK4,weight_btagwdown);
  hist_initbtagwdown[6]->Fill(npfjetAK8,weight_btagwdown);
#ifdef E_MU_TTBar  
  hist_initbtagwdown[7]->Fill((fmucand + felcand).M(),weight_btagwdown);
  hist_initbtagwdown[8]->Fill(fmucand.Pt(),weight_btagwdown);
  hist_initbtagwdown[9]->Fill(fmucand.Eta(),weight_btagwdown);
  hist_initbtagwdown[10]->Fill(fmucand.Phi(),weight_btagwdown);
  hist_initbtagwdown[11]->Fill(felcand.Pt(),weight_btagwdown);
  hist_initbtagwdown[12]->Fill(felcand.Eta(),weight_btagwdown);
  hist_initbtagwdown[13]->Fill(felcand.Phi(),weight_btagwdown);
#elif defined(E_E_TTBar)
  hist_initbtagwdown[7]->Fill((fe1cand + fe2cand).M(),weight_btagwdown);
  hist_initbtagwdown[8]->Fill(fe1cand.Pt(),weight_btagwdown);
  hist_initbtagwdown[9]->Fill(fe1cand.Eta(),weight_btagwdown);
  hist_initbtagwdown[10]->Fill(fe1cand.Phi(),weight_btagwdown);
  hist_initbtagwdown[11]->Fill(fe2cand.Pt(),weight_btagwdown);
  hist_initbtagwdown[12]->Fill(fe2cand.Eta(),weight_btagwdown);
  hist_initbtagwdown[13]->Fill(fe2cand.Phi(),weight_btagwdown);
#elif defined(MU_MU_TTBar)
  hist_initbtagwdown[7]->Fill((fmu1cand + fmu2cand).M(),weight_btagwdown);
  hist_initbtagwdown[8]->Fill(fmu1cand.Pt(),weight_btagwdown);
  hist_initbtagwdown[9]->Fill(fmu1cand.Eta(),weight_btagwdown);
  hist_initbtagwdown[10]->Fill(fmu1cand.Phi(),weight_btagwdown);
  hist_initbtagwdown[11]->Fill(fmu2cand.Pt(),weight_btagwdown);
  hist_initbtagwdown[12]->Fill(fmu2cand.Eta(),weight_btagwdown);
  hist_initbtagwdown[13]->Fill(fmu2cand.Phi(),weight_btagwdown);
#endif
  hist_initbtagwdown[14]->Fill(bjv[0].Pt(),weight_btagwdown);
  hist_initbtagwdown[15]->Fill(bjv[0].Eta(),weight_btagwdown);
  hist_initbtagwdown[16]->Fill(bjv[0].Phi(),weight_btagwdown);
  */
  if (npfjetAK8>0) {
    int npfjetAK8wmass(0);
    int ljcand(-1);
    float mass(0.);
    float maxpt(-100);
    for(unsigned ijet=0; ijet<LJets.size(); ijet++){
      if (LJets[ijet].mass > mass) {
	npfjetAK8wmass++;
	if (LJets[ijet].pt > maxpt) {
	 maxpt = LJets[ijet].pt;
	 ljcand = ijet;
	}
      }
    }
    
    if (npfjetAK8wmass>0) {
      hist_binit[0]->Fill(npfjetAK4,weight);
      hist_binit[1]->Fill(nbjetAK4,weight);
      if (npfjetAK4>0) hist_binit[2]->Fill(pfjetAK4pt[0],weight);
      if (npfjetAK4>1) hist_binit[3]->Fill(pfjetAK4pt[1],weight);
      if (nbjetAK4>0) hist_binit[4]->Fill(bjv[0].Pt(),weight);
      if (nbjetAK4>1) hist_binit[5]->Fill(bjv[1].Pt(),weight);
    }
  }


  // end //                                                                            
  
  return kTRUE;                                                                       
}
void Anal_Leptop_PROOF::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  fileOut->cd();
  fileOut->Write();
  fOutput->Add(OutFile);
  fileOut->Close();
}

void Anal_Leptop_PROOF::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
}
