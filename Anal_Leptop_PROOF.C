/*
  weight = event_weight; 
  should be in the same ntuple, event by event we are going to read it.
  
  Order lepton according to the matching of leading Ak8 jets. 
  Lepton tagging
  
 */

#define Anal_Leptop_PROOF_cxx
//#include "Anal_Leptop_PROOF.h"
#include "getobjects.h"

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
  
  OutFile = new TProofOutputFile("test_output.root");
  
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
  
  //GMA
  //Tout = new TTree("leptop","leptop");
  //Tout->Branch("weight",&weight,"weight/F");
  
  Tnewvar = new TTree("newvars","newvars");
  Tnewvar->Branch("event_pt_weight",&event_pt_weight,"event_pt_weight/F");
  Tnewvar->Branch("weight",&weight,"weight/F");
  Tnewvar->Branch("qscale",&qscale,"qscale/F");
  Tnewvar->Branch("M_l1l2",&M_l1l2,"M_l1l2/F");
  Tnewvar->Branch("rat_l2pt_l1pt",&rat_l2pt_l1pt,"rat_l2pt_l1pt/F");
  Tnewvar->Branch("deltaPhi_l1l2",&deltaPhi_l1l2,"deltaPhi_l1l2/F");
  Tnewvar->Branch("l1pt_nearjet",&l1pt_nearjet,"l1pt_nearjet/F");
  Tnewvar->Branch("l2pt_nearjet",&l2pt_nearjet,"l2pt_nearjet/F");
  Tnewvar->Branch("met_pt",&met_pt,"met_pt/F");
  Tnewvar->Branch("met_phi",&met_phi,"met_phi/F");
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
  Tnewvar->Branch("ak41pt",&ak41pt,"ak41pt/F");
  Tnewvar->Branch("ak42pt",&ak42pt,"ak42pt/F"); 
  Tnewvar->Branch("ak81pt",&ak81pt,"ak81pt/F"); 
  Tnewvar->Branch("ak81y",&ak81y,"ak81y/F");
  Tnewvar->Branch("ak81mass",&ak81mass,"ak81mass/F");
  Tnewvar->Branch("ak81sdmass",&ak81sdmass,"ak81sdmass/F");
  Tnewvar->Branch("ak81chrad",&ak81chrad,"ak81chrad/F");
  Tnewvar->Branch("ak81NHad",&ak81NHad,"ak81NHad/F");                
  Tnewvar->Branch("ak81neuhad",&ak81neuhad,"ak81neuhad/F");
  Tnewvar->Branch("ak81subhaddiff",&ak81subhaddiff,"ak81subhaddiff/F");
  Tnewvar->Branch("ak81tau21",&ak81tau21,"ak81tau21/F");
  Tnewvar->Branch("ak81tau32",&ak81tau32,"ak81tau32/F");
  Tnewvar->Branch("ak81deep_tvsqcd",&ak81deep_tvsqcd,"ak81deep_tvsqcd/F");
  Tnewvar->Branch("ak81deep_wvsqcd",&ak81deep_wvsqcd,"ak81deep_wvsqcd/F");
  Tnewvar->Branch("ak82pt",&ak82pt,"ak82pt/F");
  Tnewvar->Branch("ak82y",&ak82y,"ak82y/F");
  Tnewvar->Branch("ak82mass",&ak82mass,"ak82mass/F");
  Tnewvar->Branch("ak82sdmass",&ak82sdmass,"ak82sdmass/F");
  Tnewvar->Branch("ak82chrad",&ak82chrad,"ak82chrad/F");
  Tnewvar->Branch("ak82NHad",&ak82NHad,"ak82NHad/F");
  Tnewvar->Branch("ak82neuhad",&ak82neuhad,"ak82neuhad/F");
  Tnewvar->Branch("ak82subhaddiff",&ak82subhaddiff,"ak82subhaddiff/F");
  Tnewvar->Branch("ak82tau21",&ak82tau21,"ak82tau21/F");
  Tnewvar->Branch("ak82tau32",&ak82tau32,"ak82tau32/F");
  Tnewvar->Branch("ak82deep_tvsqcd",&ak82deep_tvsqcd,"ak82deep_tvsqcd/F");
  Tnewvar->Branch("ak82deep_wvsqcd",&ak82deep_wvsqcd,"ak82deep_wvsqcd/F");
  Tnewvar->Branch("delta_phibl1bl2",&delta_phibl1bl2,"delta_phibl1bl2/F"); 
  
  //Tnewvar->Branch("delta_phijl1jl2",&delta_phijl1jl2,"delta_phijl1jl2/F");
  Tnewvar->Branch("deltaR_l1l2",&deltaR_l1l2,"deltaR_l1l2/F");
  
  //Tnewvar->Branch("deltaR_l1b1",&deltaR_l1b1,"deltaR_l1b1/F");
  //Tnewvar->Branch("deltaR_l2b1",&deltaR_l2b1,"deltaR_l2b1/F");
  //Tnewvar->Branch("deltaR_l1b2",&deltaR_l1b2,"deltaR_l1b2/F");
  //Tnewvar->Branch("deltaR_l2b2",&deltaR_l2b2,"deltaR_l2b2/F");
  
  Tnewvar->Branch("deltaR_l1j1",&deltaR_l1j1,"deltaR_l1j1/F");
  Tnewvar->Branch("deltaR_l2j1",&deltaR_l2j1,"deltaR_l2j1/F");
  Tnewvar->Branch("deltaR_l1j2",&deltaR_l1j2,"deltaR_l1j2/F");
  Tnewvar->Branch("deltaR_l2j2",&deltaR_l2j2,"deltaR_l2j2/F");
  Tnewvar->Branch("j1_btag_sc",&j1_btag_sc,"j1_btag_sc/F"); 
  Tnewvar->Branch("j2_btag_sc",&j2_btag_sc,"j2_btag_sc/F");

  Tnewvar->Branch("j1_btag_sc_ptsort",&j1_btag_sc_ptsort,"j1_btag_sc_ptsort/F"); 
  Tnewvar->Branch("j2_btag_sc_ptsort",&j2_btag_sc_ptsort,"j2_btag_sc_ptsort/F");

  Tnewvar->Branch("genmatch_sc1",&genmatch_sc1,"genmatch_sc1/I");	
  Tnewvar->Branch("genmatch_sc2",&genmatch_sc2,"genmatch_sc2/I");
  
  Tnewvar->Branch("topmatchvar",&topmatchvar,"topmatchvar/I");	
  Tnewvar->Branch("dirgltrthr",&dirgltrthr,"dirgltrthr/D");
  Tnewvar->Branch("dirglthrmin",&dirglthrmin,"dirglthrmin/D");
    
  char name[1000];

  for (int ij=0; ij<nprvar; ij++) {
    hist_prvar[ij] = new TH1D(prvar_name[ij], prvar_name[ij], prvar_bins[ij], prvar_low[ij], prvar_high[ij]);
  }
	
  
  for (int ij=0; ij <ntypes; ij++) {
    for (int jk=0; jk<ntcount; jk++) {
      sprintf (name, "%s_cnt%i", ptprvar_name[ij], jk);
      hist_prptvar[ij][jk] = new TH2D(name, name, 60, -3.0, 3.0, nybin, ybins); 
    }
  }
  
  for (int ij=0; ij<npr_angle; ij++) {
    if (ij<2) {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 30.0, 510.0, 120, 0.0, 4.8);
    } else if (ij<8) {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 180.0, 840.0, 120, 0.0, 4.8);
    } else {
      hist_prptangle[ij] = new TH2D(pr_angle[ij], pr_angle[ij], 60, 30.0, 510.0, 120, 0.0, 2.0);
    }
  }
  
  for(int init=0; init<nhist_in; init++){
    char namein[1000]; //nameinup[1000], nameindn[1000];
    char titlein[1000];
    
    sprintf(namein,"hist_%s",initnames[init]);
    sprintf(titlein,"%s",titlenames[init]);
    hist_init[init] = new TH1D(namein,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
    hist_init[init]->Sumw2();
    
    /*
      sprintf(nameinup,"hist_%s_puup",initnames[init]);
      sprintf(nameindn,"hist_%s_pudn",initnames[init]);
      hist_init_pu_sys[init][0] = new TH1D(nameinup,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init_pu_sys[init][0]->Sumw2();
      hist_init_pu_sys[init][1] = new TH1D(nameindn,titlein,ini_nbins[init],ini_low[init],ini_up[init]);
      hist_init_pu_sys[init][1]->Sumw2();
    */
  }
	
  for(int lvar=0; lvar<nobshist; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Obs_%s",obsnames[lvar]);
    hist_obs[lvar] = new TH1D(lnamein,lnamein,obs_nbins[lvar],obs_low[lvar],obs_up[lvar]);
    hist_obs[lvar]->Sumw2();
    
  }

  for(int lvar=0; lvar<nhistbtag; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Obs_%s",obsnames_btag[lvar]);
    hist_obs_btag[lvar] = new TH1D(lnamein,lnamein,obs_nbins_btag[lvar],obs_low_btag[lvar],obs_up_btag[lvar]);
    hist_obs_btag[lvar]->Sumw2();
  }
  
  for(int lvar=0; lvar<nhistgenmatch; lvar++){
    char lnamein[1000];
    sprintf(lnamein,"Obs_%s",obsnames_genmatch[lvar]);
    hist_obs_genmatch[lvar] = new TH1D(lnamein,lnamein,obs_nbins_genmatch[lvar],obs_low_genmatch[lvar],obs_up_genmatch[lvar]);
    hist_obs_genmatch[lvar]->Sumw2();		
  }
  
  hist_2d_deltaR_vsbtagsc[0]=new TH2D("hist_2D_deltaR_vsbtagsc_1","deltaR(leading AK4 jet, leading AK8 jet) vs b tag score of leading AK4 jet",60,0,1,60,0,6);
  hist_2d_deltaR_vsbtagsc[0]->Sumw2();
  hist_2d_deltaR_vsbtagsc[1]=new TH2D("hist_2D_deltaR_vsbtagsc_2","deltaR(sub-leading AK4 jet, sub-leading AK8 jet) vs b tag score of sub-leading AK4 jet",60,0,1,60,0,6);
  hist_2d_deltaR_vsbtagsc[1]->Sumw2();
  
  hist_2d_pt_vsbtagsc[0]=new TH2D("hist_2D_ak4pt_vsbtagsc_1","pt of leading AK4 jet vs b tag score of leading AK4 jet",120,0,1000,60,0,1);
  hist_2d_pt_vsbtagsc[0]->Sumw2();
  hist_2d_pt_vsbtagsc[1]=new TH2D("hist_2D_ak4pt_vsbtagsc_2","pt of sub leading AK4 jet vs b tag score of sub-leading AK4 jet",120,0,1000,60,0,1);
  hist_2d_pt_vsbtagsc[1]->Sumw2();
  hist_2d_pt_vsbtagsc[2]=new TH2D("hist_2D_ak8pt_vsbtagsc_1","pt of leading AK8 jet vs b tag score of leading AK4 jet",120,0,2000,60,0,1);
  hist_2d_pt_vsbtagsc[2]->Sumw2();
  hist_2d_pt_vsbtagsc[3]=new TH2D("hist_2D_ak8pt_vsbtagsc_2","pt of sub leading AK8 jet vs b tag score of sub-leading AK4 jet",120,0,2000,60,0,1);
  hist_2d_pt_vsbtagsc[3]->Sumw2();
  
  hist_2d_deltaR_vspt[0]=new TH2D("hist_2d_deltaR_vspt_1","deltaR(leading AK4 jet, leading lep) vs pt of leading AK8 jet",120,200,2500,120,0,1.5);
  hist_2d_deltaR_vspt[0]->Sumw2();  
  hist_2d_deltaR_vspt[1]=new TH2D("hist_2d_deltaR_vspt_2","deltaR(sub-leading AK4 jet, sub-leading lep) vs pt of sub-leading AK8 jet",120,200,2500,120,0,1.5);
  hist_2d_deltaR_vspt[1]->Sumw2();  
  
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
  hist_npv = new TH1D(name,title,100,-0.5,99.5);//80,-0.1,79.9);
  hist_npv->Sumw2();
    
  sprintf(name,"N_PV_nopuwt");
  sprintf(title,"# of Primary Vertices");
  hist_npv_nopuwt = new TH1D(name,title,100,-0.5,99.5);
  hist_npv_nopuwt->Sumw2();
  
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
  
  /*
    reader3 = new TMVA::Reader( "BDTG_Rt" );
    reader3->AddVariable( "selpfjetAK8NHadF", &in_pfjetAK8NHadF);
    reader3->AddVariable( "selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
    reader3->AddVariable( "selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
    reader3->AddVariable( "selpfjetAK8tau21", &in_pfjetAK8tau21);
    reader3->AddVariable( "selpfjetAK8chrad", &in_pfjetAK8chrad);
    reader3->AddVariable( "selpfjetAK8sdmass", &in_pfjetAK8sdmass);
    reader3->BookMVA("BDTG method", weightfile3);
  */
  
  reader4 = new TMVA::Reader( "BDTG_Rmu" );
  reader4->AddVariable( "selpfjetAK8NHadF", &in_mupfjetAK8NHadF);
  reader4->AddVariable( "selpfjetAK8neunhadfrac", &in_mupfjetAK8neunhadfrac);
  reader4->AddVariable( "selpfjetAK8subhaddiff", &in_mupfjetAK8subhaddiff);
  reader4->AddVariable( "selpfjetAK8tau21", &in_mupfjetAK8tau21);
  reader4->AddVariable( "selpfjetAK8chrad", &in_mupfjetAK8chrad);
  reader4->AddVariable( "selpfjetAK8sdmass", &in_mupfjetAK8sdmass);

  
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
  //Tout->Fill();
  
  /*     
	 TString str;                                                            
     
	 str = TString::Format("check start evt %d ]]]",ievt);          
	 if(gProofServ) gProofServ->SendAsynMessage(str);
  */
  
  
  // event selection starts
  
  hist_prvar[1]->Fill(nprimi,weight);
  if (nprimi<1)  return kFALSE;
  hist_count->Fill(1,weight);
  hist_npv_nopuwt->Fill(nprimi,weight);
  
  //pileup weight factor needs to be updated with UL/
  
  //if(isMC){
  //if(npu_vert>=0 && npu_vert<100){
  //puWeight = pu_rat18[npu_vert];
  //puWeightUp = pu_rat18_up[npu_vert];
  //puWeightDown = pu_rat18_dn[npu_vert];
  //}
  //if(!isnan(puWeightUp) || fabs(puWeightUp)<1.e+6){
  //weight_puwup = weight*puWeightUp;                                           
               
  //}
  //if(!isnan(puWeightDown) || fabs(puWeightDown)<1.e+6){
  //weight_puwdown = weight*puWeightDown;                                        
  
  //}
  //if(!isnan(puWeight) || fabs(puWeight)<1.e+6){
  // weight *= puWeight;                                                         
             
  //}
  //}
  //hist_npv->Fill(nprimi,weight);

  // First demand that the event should at least one of the triggers considered in various topologies            
  // !CAUTION: Change/Remove this while deriving trigger efficiency 
    
  bool itrig_pass = false;
  bool itrig_onlysinglee = false;
  hist_prvar[2]->Fill(int(hlt_Mu37Ele27)*32+int(hlt_Mu27Ele37)*16+int(hlt_Mu37TkMu27)*8+ int(hlt_DoubleEle25)*4+int(hlt_Mu50)*2+int(hlt_Ele50_PFJet165), weight);

  itrig_pass = (hlt_Mu37Ele27||hlt_Mu27Ele37||hlt_Mu37TkMu27||hlt_DoubleEle25||hlt_Mu50||hlt_Ele50_PFJet165);
  itrig_onlysinglee = (!hlt_Mu37Ele27 && !hlt_Mu27Ele37 && !hlt_Mu37TkMu27 && !hlt_DoubleEle25 && !hlt_Mu50 && hlt_Ele50_PFJet165);
		
  if(!itrig_pass) return kFALSE; //event should at least fire a dileptonic/single lepton trigger            
  
  // end of basic trigger criterion //

  hist_count->Fill(2,weight);

  // Get generator-level particles //
       
  vector<GenParton> genpartons;
  vector<GenParton> LHEtops;
  vector<TopQuark> gentops;

  if(isMC){
    getPartons(genpartons);
  }
  if(isMC && isTT){
    // Get GEN-level top quarks                                                         
       
    getLHETops(LHEtops,genpartons); // before shower (get original top quarks which have decayed) --> will be usedto derive top pt reweighting                                   
              
                             
    getGENTops(gentops,genpartons); // after shower (get top quarks from its daughters) --> will tell details about the signature of ttbar events at GEN level
    
    nleptop = nhadtop = 0;
    int leptop_id_daught[2];
    
    for(auto & top: gentops){
      if(abs(top.daughter[0].pdgId)==11 || abs(top.daughter[0].pdgId)==13 || abs(top.daughter[0].pdgId)==15){
	leptop_id_daught[nleptop] = abs(top.daughter[0].pdgId);
	nleptop++;
      }
      if(abs(top.daughter[0].pdgId)==1 || abs(top.daughter[0].pdgId)==3){
	nhadtop++;
      }
    }
  
    // tagging top structure of events //
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
    
#ifdef E_MU_TTBar
    if(!(DiLeptt && EMU)) return kFALSE; //for signal EMU
    //if((DiLeptt && EMU)) return kFALSE; //for non-signal EMU TTbar
#elif defined(E_E_TTBar)
    if(!(DiLeptt && EE)) return kFALSE; //for signal EE
    //if((DiLeptt && EE)) return kFALSE; //for non-signal EE TTbar
#elif defined(MU_MU_TTBar) 
    if(!(DiLeptt && MUMU)) return kFALSE; //for signal MUMU
    //if((DiLeptt && MUMU)) return kFALSE; //for non-signal MUMU TTbar
#endif
    bool boosted = false;
    boosted = (LHEtops.size()>0 && LHEtops[0].pt>250. && LHEtops[1].pt>250);
    //if(!boosted) return kFALSE; // Use for boosted signal
  } // if(isMC && isTT)

  //Get RECO-level objects //
  //First sorted out electron and muon
  
  //Some histograms before getting electrons with analysis criteria
  for(int ie=0; ie<nelecs; ie++) {
    hist_prptvar[3][min(ntcount-1,ie)]->Fill(min(float(2.99), max(float(-2.99), eleta[ie])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(elpt[ie]))),weight);
  }
  //Here you get electrons with your criteria
  vector <Electron> velectrons;
  getelectrons(velectrons,lepton_pt_cut,absetacut);
  
  //Some histograms before getting muons with analysis criteria
  for(int mu=0; mu<nmuons; mu++){
    hist_prptvar[4][min(ntcount-1,mu)]->Fill(min(float(2.99), max(float(-2.99), muoneta[mu])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(muonpt[mu]))),weight);
  }

  //Here you get muons with your criteria
  vector <Muon> vmuons;
  getmuons(vmuons,lepton_pt_cut,absetacut);
  
  //Make lepton collection from electrons & muons (using only common variables)
  vector <Lepton> vleptons;
  getLeptons(vleptons,vmuons,velectrons,lepton_pt_cut);
    
  //int nbjet_cut = -1; Need to use when btagwt will be considered in mc event 
  //float btagwt = 1.;
  //float btagwtup = 1.;
  //float btagwtdown = 1.;
  //float btag_eff = 1;
  //int nbjetAK4 = 0;
  //int nbjetAK4_lead = 0;

  //Some histograms before getting jets with analysis criteria
  for(int ijet=0; ijet<npfjetAK4; ijet++){
    for(int ie=0; ie<(int)velectrons.size(); ie++) {                         
      double delr = delta2R(pfjetAK4eta[ijet],pfjetAK4phi[ijet], velectrons[ie].eta, velectrons[ie].phi);               
      hist_prptangle[0]->Fill(pfjetAK4pt[ijet], delr,weight);
    }
    for(int imu=0; imu<(int)vmuons.size(); imu++)                                
      {                                                                         
	double delr = delta2R(pfjetAK4eta[ijet], pfjetAK4phi[ijet], vmuons[imu].eta, vmuons[imu].phi);                  
	hist_prptangle[1]->Fill(pfjetAK4pt[ijet], delr,weight);
      }
  }
  //Here you get AK4 jets with your criteria                                     
  vector <AK4Jet> Jets;
  getAK4jets(Jets,AK4jet_pt_cut,absetacut,isMC);
  //Perform lepton-jet cleaning                                                  
  LeptonJet_cleaning(Jets,vleptons,AK4jet_pt_cut,absetacut);
  
  //Get b-tagged jets from AK4 jets                                            
  vector <AK4Jet> BJets;
  for(auto & jet: Jets){
    if(isBJet(jet,deep_btag_cut)){
      //See similar histograms as did early for other objects
      //hist_prptvar[2][min(ntcount-1,ij)]->Fill(min(2.99, max(-2.99, jet.y)), min(ybins[nybin]-1.0, max(ybins[0]+0.1, jet.pt)),weight);
      BJets.push_back(jet);
    }
  }
  
  //Fill similar histograms for largejet as did early for other objects
  for(int ijet=0; ijet<npfjetAK8; ijet++){                                       
    if(!pfjetAK8jetID[ijet]) continue;                                           
    pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet] ;                                      
    pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];
    if(isMC){                                                                    
      pfjetAK8pt[ijet] *= (1+pfjetAK8reso[ijet]) ;                               
      pfjetAK8mass[ijet] *= (1+pfjetAK8reso[ijet]) ;
    }                                                                            
    hist_prptvar[0][min(ntcount-1,ijet)]->Fill(min(float(2.99), max(float(-2.99), pfjetAK8y[ijet])), min(ybins[nybin]-1.0, max(ybins[0]+0.1, double(pfjetAK8pt[ijet]))),weight);
  }
  
  //Here you get AK8 jets with your criteria 
  vector <AK8Jet> LJets;
  getAK8jets(LJets,AK8jet_pt_cut,absetacut,isMC);
  
  // assign information from GEN-matching
  if (isMC) {
    AssignGen(LJets,genpartons);
    if(isTT){
      TopAssignment_toJet(LJets,LHEtops,gentops);
    }
  }
  
  //Get index of AK8 jet nearest to each lepton                                  
  for(auto & lep: vleptons){
    if (LJets.size()>0) lep.AK8_neighbor_index = get_nearest_AK8Jet(LJets,lep.p4);
  }
  
  
  if(LJets.size()>0) {
    for(unsigned ijet=0; ijet<Jets.size(); ijet++){
      double dr = delta2R(LJets[0].y,LJets[0].phi,Jets[ijet].y,Jets[ijet].phi);
      hist_prptangle[2]->Fill(LJets[0].pt,dr,weight);
    }
    if (LJets.size()>1 && Jets.size()>1) {
      for(unsigned ijet=0; ijet<Jets.size(); ijet++){
	double dr2 = delta2R(LJets[1].y,LJets[1].phi,Jets[ijet].y,Jets[ijet].phi);
	hist_prptangle[3]->Fill(LJets[1].pt,dr2,weight);
      }
    }
    //Get indices of nearest lepton, AK4 jet for each AK8 jet                    
    for(auto & jet: LJets){
      //Match lepton with AK8 jets
      jet.match_lepton_index = get_nearest_lepton(vleptons,jet.p4);
      // Match AK4 jets with AK8 jets //                                         
      jet.match_AK4_index  = get_nearest_AK4(Jets,jet.p4);
      if(jet.match_AK4_index>=0 && jet.match_AK4_index<int(Jets.size())){
        jet.matchAK4deepb = Jets[jet.match_AK4_index].btag_DeepFlav;
      }
    }
    //Debarati : Keeping histograms of GMA to see the distributions before selection
      for (int jk=0; jk<(int)vleptons.size(); jk++) {
      int ityp = (vleptons[jk].lepton_id==2) ? 4 : 5;
      double dr = delta2R(LJets[0].y, LJets[0].phi, vleptons[jk].eta, vleptons[jk].phi);
      hist_prptangle[ityp]->Fill(LJets[0].pt, dr,weight);
      }
      
    //Debarati : Adding GMA concept of ordering vleptons, AK4 wrt closest to leading two AK8 jets
      int ivarl(0), ivarj(0);//Keeping it in new structure to keep GMA histograms
      if (LJets[0].match_lepton_index >=0) {
	Lepton tmplep = vleptons[LJets[0].match_lepton_index];
	vleptons.erase(vleptons.begin() + LJets[0].match_lepton_index);
	vleptons.insert(vleptons.begin(), tmplep);
	ivarl = 1;
      }
      
      if (LJets[0].match_AK4_index >=0) {
	AK4Jet tmpjet = Jets[LJets[0].match_AK4_index];
	Jets.erase(Jets.begin() + LJets[0].match_AK4_index);
	Jets.insert(Jets.begin(), tmpjet);
	ivarj = 1;
	}
      if (LJets.size()>1) {
	for (int jk=ivarl; jk<(int)vleptons.size(); jk++) {
	  int ityp = (vleptons[jk].lepton_id==2) ? 6 : 7;
	  double dr = delta2R(LJets[1].y, LJets[1].phi, vleptons[jk].eta, vleptons[jk].phi);
	  hist_prptangle[ityp]->Fill(LJets[1].pt, dr,weight);
	}
	if (LJets[1].match_lepton_index >=1 && LJets[0].match_lepton_index != LJets[1].match_lepton_index) {
	  Lepton tmplep = vleptons[LJets[1].match_lepton_index];
	  vleptons.erase(vleptons.begin() + LJets[1].match_lepton_index);
	  vleptons.insert(vleptons.begin()+1, tmplep);
	  }
	
	if (LJets[1].match_AK4_index >= 1 && LJets[0].match_AK4_index != LJets[1].match_AK4_index) {
	AK4Jet tmpjet = Jets[LJets[1].match_AK4_index];
	Jets.erase(Jets.begin() + LJets[1].match_AK4_index);
	Jets.insert(Jets.begin()+1, tmpjet);
	}
      }
      // Assign electronic top tagger score //                                            
      ReadTagger(LJets,vleptons,vmuons,velectrons,reader1,reader4);
  }

  //if(isnan(weight) || weight>1.e+12) { weight = 0; }
  
  /******trigger object along with pdgid*****/
  std::vector<std::pair<int,TLorentzVector> > TrigRefObj;
  
  for (int tr=0; tr<ntrigobjs; tr++) {
    TLorentzVector trigobj;
    trigobj.SetPtEtaPhiM(trigobjpt[tr],trigobjeta[tr],trigobjphi[tr],trigobjmass[tr]);
    TrigRefObj.push_back(std::make_pair(trigobjpdgId[tr],trigobj));
  }
  
  // end of object selection //
  //// Event selection ////

  hist_prvar[3]->Fill(vleptons.size(), weight);
  //at least two leptons with pT > 30 GeV at this stage
  if (vleptons.size()<2)  return kFALSE; 
  //// Condition on number of leptons is put here only. We will need at least two leptons for trigger matching

  hist_count->Fill(3,weight);
  
  bool emu_ch = false;
  bool mumu_ch = false;
  bool ee_ch = false;
  
  if ((vleptons[0].pdgId==11 && vleptons[1].pdgId==13) || (vleptons[0].pdgId==13 && vleptons[1].pdgId==11)) emu_ch = true;
  else if (vleptons[0].pdgId == 13 && vleptons[1].pdgId == 13) mumu_ch = true;
  else if (vleptons[0].pdgId == 11 && vleptons[1].pdgId == 11) ee_ch =true;
    
  hist_prvar[4]->Fill(int(emu_ch)*4+int(mumu_ch)*2+int(ee_ch));

  // Composition of two leading leptons should be correct in respective event categories //    
  bool correct_lepton_combo = false;
#ifdef E_MU_TTBar
  if(emu_ch) {correct_lepton_combo = true;}
#elif defined(E_E_TTBar)
  if(ee_ch) {correct_lepton_combo = true;}
#elif defined(MU_MU_TTBar)
  if(mumu_ch) {correct_lepton_combo = true;}
#endif
	
  // hist_count->Fill(4,weight);

  //Now comes individual channel the triggers consideration//   
  Lepton lepcand_1, lepcand_2;
  lepcand_1 = vleptons[0];
  lepcand_2 = vleptons[1];

  vector<bool> double_hlts; vector<vector<float>> double_pt_cuts; vector<vector<int>> double_pids;
  vector<bool> single_hlts; vector<float> single_pt_cuts; vector<int> single_pids; vector<float> single_other_pt_cuts; vector<int> single_other_pids;

  // _pt_cuts <-- offline pt thresholds (should be finalized from the trigger efficiency curves                                                                                  
  // _pids <-- corresponding pdgId of the objects used in triggering
  
#ifdef E_MU_TTBar
  
  double_hlts.push_back(hlt_Mu37Ele27);
  double_pt_cuts.push_back({37+3,27+3});
  double_pids.push_back({13,11});
  double_hlts.push_back(hlt_Mu27Ele37);
  double_pt_cuts.push_back({37+3,27+3});
  double_pids.push_back({11,13});

  single_hlts.push_back(hlt_Mu50);
  single_pt_cuts.push_back(50+3);
  single_pids.push_back(13);
  single_other_pt_cuts.push_back(-99);
  single_other_pids.push_back(0);
  single_hlts.push_back(hlt_Ele50_PFJet165);
  single_pt_cuts.push_back(50+3);
  single_pids.push_back(11);
  single_other_pt_cuts.push_back(165+15);
  single_other_pids.push_back(0);

#elif defined(E_E_TTBar)
    
  double_hlts.push_back(hlt_DoubleEle25);
  double_pt_cuts.push_back({25+15,25+5});
  double_pids.push_back({11,11});

  single_hlts.push_back(hlt_Ele50_PFJet165);
  single_pt_cuts.push_back(50+3);
  single_pids.push_back(11);
  single_other_pt_cuts.push_back(165+15);
  single_other_pids.push_back(0);

#elif defined(MU_MU_TTBar)
  
  double_hlts.push_back(hlt_Mu37TkMu27);
  double_pt_cuts.push_back({37+3,27+3});
  double_pids.push_back({13,13});

  single_hlts.push_back(hlt_Mu50);
  single_pt_cuts.push_back(50+3);
  single_pids.push_back(13);
  single_other_pt_cuts.push_back(-99);
  single_other_pids.push_back(0);
  
#endif

  // remember that the other object in single lepton triggers is jet by default

  bool anytrig_pass(false);
  if(double_hlts.size()>0){
    for(unsigned ihlt=0; ihlt<double_hlts.size(); ihlt++){
      if(double_hlts[ihlt]){ anytrig_pass = true; break; }
    }
  }
  
  if(single_hlts.size()>0){
    for(unsigned ihlt=0; ihlt<single_hlts.size(); ihlt++){
      if(single_hlts[ihlt]){ anytrig_pass = true; break; }
    }
  }
  
  hist_prvar[5]->Fill(int(anytrig_pass), weight);
  
  hist_prvar[6]->Fill(min(float(359.0),vleptons[0].pt), weight);
	
  bool trig_threshold_pass(false), trig_matching_pass(false);

  vector<TH1D*> hists;
  hists.push_back(hist_init[0]);
  hists.push_back(hist_init[1]);
  //vector<TH2D*> hist2d_prptangle;
  //hist2d_prptangle.push_back(hist_prptangle[8]);
  //hist2d_prptangle.push_back(hist_prptangle[9]);

  Match_trigger(double_hlts, single_hlts,
                double_pt_cuts, single_pt_cuts,
                double_pids, single_pids,
                single_other_pt_cuts, single_other_pids,
                TrigRefObj,
                lepcand_1,lepcand_2,Jets,
                trig_threshold_pass,
                trig_matching_pass,
                hists/*,hist2d_prptangle*/
                );
  
  // end of trigger stuffs //
    
  hist_prvar[7]->Fill(lepcand_1.charge*lepcand_2.charge, weight);
  hist_prvar[8]->Fill(velectrons.size());
  hist_prvar[9]->Fill(vmuons.size());
  
  bool is_additional_muons = false;
  bool is_additional_elecs = false;

#ifdef E_MU_TTBar
  is_additional_muons = (int(vmuons.size())>1);
  is_additional_elecs = (int(velectrons.size())>1);
#elif defined(E_E_TTBar)
  is_additional_muons = (int(vmuons.size())>=1);
  is_additional_elecs = (int(velectrons.size())>2);
#elif defined(MU_MU_TTBar)
  is_additional_muons = (int(vmuons.size())>2);
  is_additional_elecs = (int(velectrons.size())>=1);
#endif

  //if (LJets.size()<2) return kFALSE;

  vector<bool> event_cuts;
  //Before this stage 3 event selection cuts applied offline => at least 1 primary vertex, at least one listed trigger fired, at least two leptons with pt > 30 GeV and they are closest (< 0.7) to leading AK8 jets.
  
  //Here follows other event cuts offline//
  
  //Event should have at least one AK8 jet \1. 
  event_cuts.push_back(LJets.size()>=2 && LJets[1].pt>300.); 
  // Composition of two leading leptons should be correct in respective event categories // \2.
  event_cuts.push_back(correct_lepton_combo);
  // Did the event fire any trigger in the topology considered? \3.                            
  event_cuts.push_back(anytrig_pass);
  // offline objects should pass kinematic thresholds in triggers \4.
  event_cuts.push_back(trig_threshold_pass);
  // offline objects should be matched to trigger objects \5.                       
  event_cuts.push_back(trig_matching_pass);
  // leptons should be oppositely charged \6.
  event_cuts.push_back(vleptons.size()>=2 && (lepcand_1.charge*lepcand_2.charge)<0);
  //Cut on the number of additional muons (pt>30 GeV) //   \7.                                
  event_cuts.push_back(!is_additional_muons);
  //Cut on the number of additional electrons (pt>30 GeV) // \8.                              
  event_cuts.push_back(!is_additional_elecs);
  //Cut on at least two AK4 jets satisfying object selection criteria// \9.
  event_cuts.push_back(Jets.size()>=2);
  //Cut on AK8 and AK4 jet matching// \10.
  event_cuts.push_back(LJets.size()>=2 && Jets.size()>=2 && delta2R(LJets[0].y,LJets[0].phi,Jets[0].y,Jets[0].phi)<0.8);
  //Cut on 2nd AK8 and AK4 jet matching// \11.
  event_cuts.push_back(LJets.size()>=2 && Jets.size()>=2 && delta2R(LJets[1].y,LJets[1].phi,Jets[1].y,Jets[1].phi)<0.8);
  //Cut on AK8 and lepton matching// \12.
  event_cuts.push_back(LJets.size()>=2 && Jets.size()>=2 && delta2R(LJets[0].y,LJets[0].phi, vleptons[0].eta, vleptons[0].phi)<0.7);
  //Cut on 2nd AK8 and lepton matching// \13.
  event_cuts.push_back(LJets.size()>=2 && Jets.size()>=2 && delta2R(LJets[1].y,LJets[1].phi, vleptons[1].eta, vleptons[1].phi)<0.7);
  //Cut on btag score of AK4 jets though (as loose cut as possible)// \14.
  event_cuts.push_back(Jets.size()>=2 && Jets[0].btag_DeepFlav>0 && Jets[1].btag_DeepFlav>0);
  // Inv mass of selected leptons should be > 20 GeV \15. 
  event_cuts.push_back(vleptons.size()>=2 && (((lepcand_1.p4+lepcand_2.p4).M())>20));
  // MET > 50 GeV \16.                                                                   
  event_cuts.push_back(PFMET>50);
  //event_cuts array should be useful to derive (n-1) cut efficiency

  //some variable distributions before using cuts on them at event selection//
  hist_prvar[10]->Fill(LJets.size(), weight);
  hist_prvar[11]->Fill(Jets.size(), weight);
  
  if(LJets.size()>0 && Jets.size()>0) {
    double lepakmatch =delta2R(LJets[0].y,LJets[0].phi,Jets[0].y,Jets[0].phi);
    hist_prvar[12]->Fill(lepakmatch, weight);
  }
  if(LJets.size()>1 && Jets.size()>1) {
    double lepakmatch2 = delta2R(LJets[1].y,LJets[1].phi,Jets[1].y,Jets[1].phi);
    hist_prvar[13]->Fill(lepakmatch2, weight);
  }
  hist_prvar[14]->Fill(BJets.size(), weight);
  hist_prvar[16]->Fill(min(float(359.0),Jets[0].pt), weight);
  
  bool event_pass = true;
  for(unsigned icut=0; icut<event_cuts.size(); icut++){
    event_pass *= event_cuts[icut];
    //  str = TString::Format("cut %u pass %o",icut+1,bool(event_cuts[icut]));                         
    //  if(gProofServ) gProofServ->SendAsynMessage(str);                                               
    if(!event_pass) break;
    if(event_pass){
    hist_count->Fill(4+icut,weight);
    }
  }

  if(!event_pass) return kFALSE;

  // end of event selection // 

  // Calculate all extra weights you need to apply (PU reweighting, b tag SF, lepton SF, top pt reweighting, ...) 
  // top pt reweighting //                                                                                                   
  
  if(isTT){
    float toppt_wt = 1;
    if(LHEtops.size()==2){
      toppt_wt = SF_TOP(0.0615,0.0005,TMath::Min(float(500),float(LHEtops[0].pt)),TMath::Min(float(500),float(LHEtops[1].pt)));
    }
  }
  
  // top pt reweighting ends //                                                                                              
  // end of extra weights                                                                                                    
  // fill basic distributions here: 
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
    hist_2D_msd_deepak8->Fill(LJets[ijet].sdmass,LJets[ijet].DeepTag_TvsQCD,weight);
  }//ijet                                                                                                                    
  for(unsigned ijet=0; ijet<Jets.size(); ijet++){
    if(isBJet(Jets[ijet],deep_btag_cut)) {
      if(abs(Jets[ijet].hadronFlavour)==5){  hist_2D_bpass_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
      if(abs(Jets[ijet].hadronFlavour)==4){  hist_2D_bpass_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
      if(abs(Jets[ijet].hadronFlavour)!=5 && abs(Jets[ijet].hadronFlavour)!=4){
        hist_2D_bpass_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight);
      }
    }

    if(abs(Jets[ijet].hadronFlavour)==5){ hist_2D_ball_flavb->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
    if(abs(Jets[ijet].hadronFlavour)==4){ hist_2D_ball_flavc->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
    if(abs(Jets[ijet].hadronFlavour)!=5 && abs(Jets[ijet].hadronFlavour)!=4){  hist_2D_ball_flavq->Fill(Jets[ijet].pt,fabs(Jets[ijet].eta),weight); }
  }
  // end of basic histograms 

  
  //Computation of lepton related Suman's variables//                                         
  TLorentzVector l1(0,0,0,0), l2(0,0,0,0);
  //l1.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,(vleptons[0].lepton_id==1) ? 0.105658 : 0.000511);
  //l2.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,(vleptons[1].lepton_id==1) ? 0.105658 : 0.000511);
  l1 = lepcand_1.p4;
  l2 = lepcand_2.p4;
  
  M_l1l2= rat_l2pt_l1pt= deltaPhi_l1l2= l1pt_nearjet= l2pt_nearjet= met_pt= met_phi= delta_phil1_met= delta_phil2_met= delta_phibl1_met= delta_phibl2_met= rat_metpt_ak4pt= rat_metpt_ak8pt= rat_metpt_eventHT= mt_of_l1met= mt_of_l2met= no_ak4jets= no_ak4bjets= no_ak8jets= EventHT= extra_ak4j= ptsum_extra_ak4= extra_ak4jqgl= extra_ak4jdeepb= rat_extra_ak4jpt_lpt= ak81pt= ak81y= ak81mass= ak81sdmass= ak81deep_tvsqcd= ak81deep_wvsqcd= ak82pt= ak82y= ak82mass= ak82sdmass= ak82deep_tvsqcd= ak82deep_wvsqcd= M_bl1= M_bl2= M_jl1= M_jl2= delta_phibl1bl2= delta_phijl1jl2= deltaR_l1l2= deltaR_l1j1= deltaR_l2j1= deltaR_l1j2= deltaR_l2j2= j1_btag_sc= j2_btag_sc = ak81NHad = ak81chrad = ak81neuhad = ak81tau21 = ak81subhaddiff = ak81tau32 = ak82NHad = ak82chrad = ak82neuhad = ak82tau21 = ak82subhaddiff = ak82tau32 = ak41pt = ak42pt = -99; 
  genmatch_sc1 = genmatch_sc2 = topmatchvar = 0; 
  
  int ixtyp=-1;
  if (vleptons.size()>=2) {
    M_l1l2 = (l1+l2).M();
    if (vleptons[0].lepton_id==2 && vleptons[1].lepton_id==2) {ixtyp=0;
    } else if (vleptons[0].lepton_id==1 && vleptons[1].lepton_id==1) {ixtyp=2;
    } else { ixtyp=1;}
    hist_prvar[17+ixtyp]->Fill(min(float(359.0), M_l1l2), weight);
  }
  //hist_new_var[0]->Fill(M_l1l2,weight);
  
  //Computation of selected lepton related variables (Suman's proposal)
  rat_l2pt_l1pt = l2.Pt()/max(1.0,l1.Pt());
  //hist_new_var[1]->Fill(rat_l1pt_l2pt,weight);
  
  deltaPhi_l1l2 = PhiInRange(l1.Phi() - l2.Phi());
  //hist_new_var[2]->Fill(deltaPhi_l1l2,weight);
	
  //2d iso variables for l1 and l2//                                             
  float dRl1_min(0.8), dRl2_min(0.8);
  int nearjet_l1(-1), nearjet_l2(-1);
	
  if(Jets.size()>0){
    TLorentzVector j_mom;
    j_mom.SetPtEtaPhiM(Jets[0].pt, Jets[0].eta, Jets[0].phi, Jets[0].mass);
    M_jl1 = (j_mom + l1).M();
    l1pt_nearjet = ((l1.Vect()).Perp(j_mom.Vect()));
    if (Jets.size()>1) { 
      //hist_new_var[3]->Fill(l1pt_nearjet,weight);
      j_mom.SetPtEtaPhiM(Jets[1].pt, Jets[1].eta, Jets[1].phi, Jets[1].mass);
      M_jl2 = (j_mom + l2).M();
      l2pt_nearjet = ((l2.Vect()).Perp(j_mom.Vect()));
    }
  }
  //hist_new_var[4]->Fill(l2pt_nearjet,weight);
  
  //Computation of MET related Suman's proposed variables//
  if (PFMET > -999) { //  && PFMETPhi != -1000) { 
    float metx = PFMET*std::cos(PFMETPhi);
    float mety = PFMET*std::sin(PFMETPhi);
    TLorentzVector metvector;
    
    metvector.SetPxPyPzE(metx,mety,0,PFMET); //as mass and Pz components are 0, thus E = Pt  
    met_pt = PFMET;
    //hist_new_var[5]->Fill(met_pt,weight);
		
    met_phi = metvector.Phi(); //always be 0.
    //hist_new_var[6]->Fill(met_eta,weight);
		
    delta_phil1_met = PhiInRange(l1.Phi() - metvector.Phi());
    //hist_new_var[7]->Fill(delta_phil1_met,weight);
		
    delta_phil2_met = PhiInRange(l2.Phi() - metvector.Phi());
    //hist_new_var[8]->Fill(delta_phil2_met,weight);
		
    TLorentzVector bl1_syst;
    TLorentzVector nearbl1;
    if (Jets.size()>0) { 
      nearbl1.SetPtEtaPhiM(Jets[0].pt, Jets[0].eta, Jets[0].phi, Jets[0].mass);
      bl1_syst = nearbl1 + l1;
      M_bl1 = bl1_syst.M();
      delta_phibl1_met = PhiInRange(bl1_syst.Phi() - metvector.Phi());
      if (Jets.size()>1) { 
	TLorentzVector nearbl2;
	TLorentzVector bl2_syst;
	nearbl2.SetPtEtaPhiM(Jets[1].pt, Jets[1].eta, Jets[1].phi, Jets[1].mass);
	bl2_syst = nearbl2 + l2;
	M_bl2 = bl2_syst.M();
	delta_phibl2_met = PhiInRange(bl2_syst.Phi() - metvector.Phi());
	
	delta_phibl1bl2 = PhiInRange(bl1_syst.Phi() - bl2_syst.Phi());
      }
    }
    if (Jets.size()>0) { 
      rat_metpt_ak4pt = metvector.Pt()/max(float(1.0),Jets[0].pt);
      rat_metpt_ak8pt = metvector.Pt()/max(float(1.0),LJets[0].pt);
    }
    //GMA defined later on    rat_metpt_eventHT = metvector.Pt()/Event_HT;
    
    mt_of_l1met = (metvector+l1).Mt();
    mt_of_l2met = (metvector+l2).Mt();
  }
  
  no_ak4jets = Jets.size();
  no_ak4bjets = BJets.size();
  no_ak8jets = LJets.size();

  int nAK4inAK8=0;
  float ptsum(0.);
  bool found=false;  /// whether we have found the leading extra ak4 jet or not
  int extra_leadak4_index=-1;
  for(int ijet=0;ijet<(int)Jets.size();ijet++)
    {
      if( delta2R(LJets[0].y,LJets[0].phi,Jets[ijet].y,Jets[ijet].phi)<0.8  || delta2R(LJets[min(1,(int)LJets.size()-1)].y,LJets[min(1,(int)LJets.size()-1)].phi,Jets[ijet].y,Jets[ijet].phi)<0.8) 
    //if( delta2R(LJets[fjet].y,LJets[fjet].phi,Jets[ijet].y,Jets[ijet].phi)<0.8 )
	{
	  nAK4inAK8++;
	  ptsum = ptsum + Jets[ijet].pt;
	}
      else if(!found)
	{
	  extra_leadak4_index=ijet;
	  found=true;
	}

    }
  extra_ak4j = Jets.size() - nAK4inAK8;
  ptsum_extra_ak4 = 0;
  for(int ijet=0;ijet<(int)Jets.size();ijet++)
    {
      ptsum_extra_ak4 = ptsum_extra_ak4 + Jets[ijet].pt;
    }
  ptsum_extra_ak4 = ptsum_extra_ak4 - ptsum;
  
    if(extra_ak4j>=1 && extra_leadak4_index >=0){
      extra_ak4jqgl = Jets[extra_leadak4_index].qgl;
      extra_ak4jdeepb = Jets[extra_leadak4_index].btag_DeepFlav;
      
      if (delta2R(Jets[extra_leadak4_index].eta, Jets[extra_leadak4_index].phi, l1.Eta(), l1.Phi()) < delta2R(Jets[extra_leadak4_index].eta, Jets[extra_leadak4_index].phi, l2.Eta(), l2.Phi())) {
	rat_extra_ak4jpt_lpt = Jets[extra_leadak4_index].pt/max(1.0,l1.Pt());
      } else {
	rat_extra_ak4jpt_lpt = Jets[extra_leadak4_index].pt/max(1.0,l2.Pt());
      }
    }
    //  EventHT = Event_HT;

    //calculate the directly global transverse thrust i.e. dirgltrthr
    dirgltrthr = 0;
    dirglthrmin = 0;
    
  if (Jets.size()>=2) {
    std::vector<TLorentzVector> allsjets_4v;
    double Pt_sum(0.);
    for(unsigned ijet=0; ijet< Jets.size(); ijet++){
      TLorentzVector sjv;
      sjv.SetPtEtaPhiM(Jets[ijet].pt,Jets[ijet].eta,Jets[ijet].phi,Jets[ijet].mass);
      allsjets_4v.push_back(sjv);
      Pt_sum += Jets[ijet].pt;
    }
		
    for(unsigned ilep=0; ilep<vleptons.size(); ilep++){
      TLorentzVector sjv;
      sjv.SetPtEtaPhiM(vleptons[ilep].pt, vleptons[ilep].eta, vleptons[ilep].phi, (vleptons[ilep].lepton_id==1) ? 0.105658 : 0.000511);
      allsjets_4v.push_back(sjv);
      Pt_sum += Jets[ilep].pt;
    }
    EventHT = Pt_sum;	
    rat_metpt_eventHT = PFMET/max(float(1.0), EventHT);
    
    std::vector<double> ThrustAxis;
    std::vector<double> Thrust;
    
    for(unsigned int j =0;j<4;j++){                                                            
      Thrust.push_back(0.);                                                                    
    }
    Thrust_calculate(allsjets_4v,Thrust);
    dirgltrthr = Thrust[3];
    
    //now comes directly global thrust minor dirglthrmin                         
    //rotate the coordinate system around the beam axis such that 
    //the thrust axis is the new y'-Axis - the projections are                   
    //simply the new y-values then                                               

    double alpha=atan2(Thrust[1],Thrust[0]);
    for(unsigned int i=0; i<allsjets_4v.size(); i++){
      dirglthrmin += fabs(-sin(alpha)*allsjets_4v[i].Px()+cos(alpha)*allsjets_4v[i].Py());
    }
    dirglthrmin = dirglthrmin/max(1.0, Pt_sum);
  }

  // ********thrust calculation ended******
  ak81pt = LJets[0].pt;
  ak81y = LJets[0].y;
  ak81mass = LJets[0].mass;
  ak81sdmass = LJets[0].sdmass;
  ak81deep_tvsqcd = LJets[0].DeepTag_TvsQCD;
  ak81deep_wvsqcd = LJets[0].DeepTag_WvsQCD;

  if (LJets.size()>1 ) {
    ak82pt = LJets[1].pt;
    ak82y = LJets[1].y;
    ak82mass = LJets[1].mass;
    ak82sdmass = LJets[1].sdmass;
    ak82deep_tvsqcd = LJets[1].DeepTag_TvsQCD;
    ak82deep_wvsqcd = LJets[1].DeepTag_WvsQCD;
  }

  deltaR_l1l2 = delta2R(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi());
  
  //deltaR_l1b1 = delta2R(l1.Eta(),l1.Phi(),bjv[0].Eta(),bjv[0].Phi());
  //if (bjv.size() >1) deltaR_l1b2 = delta2R(l1.Eta(),l1.Phi(),bjv[1].Eta(),bjv[1].Phi());
  //deltaR_l2b1 = delta2R(l2.Eta(),l2.Phi(),bjv[0].Eta(),bjv[0].Phi());
  //if (bjv.size() >1) deltaR_l2b2 = delta2R(l2.Eta(),l2.Phi(),bjv[1].Eta(),bjv[1].Phi());
  //L1 has to attached with large radius jets
	
  deltaR_l1j1 = delta2R(l1.Eta(),l1.Phi(),Jets[0].eta,Jets[0].phi);
  deltaR_l1j2 = delta2R(l1.Eta(),l1.Phi(),Jets[1].eta,Jets[1].phi);
  deltaR_l2j1 = delta2R(l2.Eta(),l2.Phi(),Jets[0].eta,Jets[0].phi);
  deltaR_l2j2 = delta2R(l2.Eta(),l2.Phi(),Jets[1].eta,Jets[1].phi);

  j1_btag_sc = Jets[0].btag_DeepFlav;
  j2_btag_sc = Jets[1].btag_DeepFlav;

  /*
    int itopf=-1;
    if(topblep[0][0] >-1 && topblep[0][0]>-1){
    if(genpartpt[topblep[0][0]] < genpartpt[topblep[1][0]])
    itopf=1;
    else itopf=0;}
    for(int it=0;it<2;it++)
    {
    if(topblep[it][0]>-1 && topblep[it][1] >-1 && topblep[it][2]>-1)
    {
    if( delta2R(LJets[0].y,LJets[0].phi,genparteta[topblep[it][2]],genpartphi[topblep[it][2]]) < 0.7)
    genmatch_sc1 += 1;

    if(delta2R(LJets[0].y,LJets[0].phi,genparteta[topblep[it][1]],genpartphi[topblep[it][1]]) < 0.8)
    genmatch_sc1 += 10;
    
    if(delta2R(LJets[0].y,LJets[0].phi,genparteta[topblep[it][0]],genpartphi[topblep[it][0]]) <  0.8)
    {
    genmatch_sc1 += 100;
    if(itopf == it) topmatchvar+=10000;
    else topmatchvar+=1000;
    }
    
    if(LJets.size()>1){
    if(delta2R(LJets[1].y,LJets[1].phi,genparteta[topblep[it][2]],genpartphi[topblep[it][2]]) < 0.7)
    genmatch_sc2 += 1;
    
    if(delta2R(LJets[1].y,LJets[1].phi,genparteta[topblep[it][1]],genpartphi[topblep[it][1]]) < 0.8)
    genmatch_sc2 += 10;
    
    if(delta2R(LJets[1].y,LJets[1].phi,genparteta[topblep[it][0]],genpartphi[topblep[it][0]]) < 0.8)
    {
    genmatch_sc2 += 100;
    if(itopf == it) topmatchvar+=100;
    else topmatchvar+=10;
    }
    }
    }
    }*/
  ak41pt=Jets[0].pt;
  ak42pt=Jets[1].pt;
  
  ak81NHad = LJets[0].NHadF;
  ak81chrad =LJets[0].chrad;
  ak81neuhad =LJets[0].neunhadfrac;
  ak81tau21 =LJets[0].tau21;
  ak81subhaddiff =LJets[0].subhaddiff;
  ak81tau32 = LJets[0].tau32;
  if(LJets.size()>1 ){
    ak82NHad =LJets[1].NHadF;
    ak82chrad =LJets[1].chrad;
    ak82neuhad =LJets[1].neunhadfrac;
    ak82tau21 =LJets[1].tau21;
    ak82subhaddiff =LJets[1].subhaddiff;
    ak82tau32 =LJets[1].tau32;
  }
  Tnewvar->Fill();
  
  if( Jets.size()>0 && LJets.size()>0)  hist_2d_pt_vsbtagsc[2]->Fill(LJets[0].pt,Jets[0].btag_DeepFlav,weight);
  if( Jets.size()>1 && LJets.size()>1)  hist_2d_pt_vsbtagsc[3]->Fill(LJets[1].pt,Jets[1].btag_DeepFlav,weight);
  
  hist_init[2]->Fill(nmuons,weight);
  hist_init[3]->Fill(nelecs,weight);
  hist_init[4]->Fill(PFMET,weight);
  hist_init[5]->Fill(nprimi,weight);
  hist_init[6]->Fill(Jets.size(),weight);
  hist_init[7]->Fill(BJets.size(),weight);
  hist_init[8]->Fill(LJets.size(),weight);

  hist_obs[0]->Fill(LJets[0].pt,weight);
  hist_obs[1]->Fill(LJets[0].y,weight);
  hist_obs[2]->Fill(LJets[0].mass,weight);
  hist_obs[3]->Fill(LJets[0].NHadF,weight);
  hist_obs[4]->Fill(LJets[0].neunhadfrac,weight);
  hist_obs[5]->Fill(LJets[0].sdmass,weight);
  hist_obs[6]->Fill(LJets[0].chrad,weight);
  hist_obs[7]->Fill(LJets[0].subhaddiff,weight);
  hist_obs[8]->Fill(LJets[0].tau21,weight);
  hist_obs[45]->Fill(LJets[0].tau32,weight);

  if(Jets[0].btag_DeepFlav > deep_btag_cut){
    hist_obs_btag[0]->Fill(LJets[0].NHadF,weight);
    hist_obs_btag[1]->Fill(LJets[0].neunhadfrac,weight);
    hist_obs_btag[2]->Fill(LJets[0].sdmass,weight);
    hist_obs_btag[3]->Fill(LJets[0].chrad,weight);
    hist_obs_btag[4]->Fill(LJets[0].subhaddiff,weight);
    hist_obs_btag[5]->Fill(LJets[0].tau21,weight);
    hist_obs_btag[6]->Fill(LJets[0].tau32,weight);
  }
  
  if(LJets[0].hasleptop_alldecay){
    hist_obs_genmatch[0]->Fill(LJets[0].NHadF,weight);
    hist_obs_genmatch[1]->Fill(LJets[0].neunhadfrac,weight);
    hist_obs_genmatch[2]->Fill(LJets[0].sdmass,weight);
    hist_obs_genmatch[3]->Fill(LJets[0].chrad,weight);
    hist_obs_genmatch[4]->Fill(LJets[0].subhaddiff,weight);
    hist_obs_genmatch[5]->Fill(LJets[0].tau21,weight);
    hist_obs_genmatch[6]->Fill(LJets[0].tau32,weight);
  }
  
 hist_2d_deltaR_vspt[0]->Fill(LJets[0].pt,delta2R(Jets[0].y,Jets[0].phi,vleptons[0].eta,vleptons[0].phi),weight);
  
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
	
  if (Jets.size()>0) hist_obs[20]->Fill(delta2R(LJets[0].eta, LJets[0].phi, Jets[0].eta, Jets[0].phi), weight);
  if (Jets.size()>1)
    {
      hist_obs[21]->Fill(delta2R(LJets[0].eta, LJets[0].phi, Jets[1].eta, Jets[1].phi), weight);
      
       if (LJets.size()>1 ) {
	 hist_2d_deltaR_vspt[1]->Fill(LJets[1].pt,delta2R(Jets[1].y,Jets[1].phi,vleptons[1].eta,vleptons[1].phi),weight);
	 hist_obs[22]->Fill(delta2R(LJets[0].y, LJets[0].phi, LJets[1].y, LJets[1].phi), weight);
	 hist_obs[23]->Fill(LJets[1].pt,weight);
	 hist_obs[24]->Fill(LJets[1].y,weight);
	 hist_obs[25]->Fill(LJets[1].mass,weight);
	 hist_obs[26]->Fill(LJets[1].NHadF,weight);
	 hist_obs[27]->Fill(LJets[1].neunhadfrac,weight);
	 hist_obs[28]->Fill(LJets[1].sdmass,weight);
	 hist_obs[29]->Fill(LJets[1].chrad,weight);
	 hist_obs[30]->Fill(LJets[1].subhaddiff,weight);
	 hist_obs[31]->Fill(LJets[1].tau21,weight);
	 hist_obs[46]->Fill(LJets[1].tau32,weight);

	 if(Jets[1].btag_DeepFlav > deep_btag_cut){
	   hist_obs_btag[7]->Fill(LJets[1].NHadF,weight);
	   hist_obs_btag[8]->Fill(LJets[1].neunhadfrac,weight);
	   hist_obs_btag[9]->Fill(LJets[1].sdmass,weight);
	   hist_obs_btag[10]->Fill(LJets[1].chrad,weight);
	   hist_obs_btag[11]->Fill(LJets[1].subhaddiff,weight);
	   hist_obs_btag[12]->Fill(LJets[1].tau21,weight);
	   hist_obs_btag[13]->Fill(LJets[1].tau32,weight);
	 }
	 
	 if(LJets[1].hasleptop_alldecay){
	   hist_obs_genmatch[7]->Fill(LJets[1].NHadF,weight);
	   hist_obs_genmatch[8]->Fill(LJets[1].neunhadfrac,weight);
	   hist_obs_genmatch[9]->Fill(LJets[1].sdmass,weight);
	   hist_obs_genmatch[10]->Fill(LJets[1].chrad,weight);
	   hist_obs_genmatch[11]->Fill(LJets[1].subhaddiff,weight);
	   hist_obs_genmatch[12]->Fill(LJets[1].tau21,weight);
	   hist_obs_genmatch[13]->Fill(LJets[1].tau32,weight);
	 }
       
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
	 hist_obs[43]->Fill(delta2R(LJets[1].y, LJets[1].phi, Jets[0].y, Jets[0].phi), weight);
	 hist_obs[44]->Fill(delta2R(LJets[1].y, LJets[1].phi, Jets[1].y, Jets[1].phi), weight);
       }
    }
#ifdef E_MU_TTBar

  TLorentzVector fmucand, felcand;
  if(lepcand_1.pdgId==13 && lepcand_2.pdgId==11) { fmucand = lepcand_1.p4; felcand = lepcand_2.p4;   }
  else if(lepcand_1.pdgId==11 && lepcand_2.pdgId==13) { fmucand = lepcand_2.p4; felcand = lepcand_1.p4;   }

  hist_init[9]->Fill((fmucand + felcand).M(),weight);
  hist_init[10]->Fill(fmucand.Pt(),weight);
  hist_init[11]->Fill(fmucand.Eta(),weight);
  hist_init[12]->Fill(fmucand.Phi(),weight);
  
  hist_init[13]->Fill(felcand.Pt(),weight);
  hist_init[14]->Fill(felcand.Eta(),weight);
  hist_init[15]->Fill(felcand.Phi(),weight);
#elif defined(E_E_TTBar)
  hist_init[9]->Fill((l1+l2).M(),weight);
  hist_init[10]->Fill(l1.Pt(),weight);
  hist_init[11]->Fill(l1.Eta(),weight);
  hist_init[12]->Fill(l1.Phi(),weight);

  hist_init[13]->Fill(l2.Pt(),weight);
  hist_init[14]->Fill(l2.Eta(),weight);
  hist_init[15]->Fill(l2.Phi(),weight);
#elif defined(MU_MU_TTBar)
  hist_init[9]->Fill((l1+l2).M(),weight);
  hist_init[10]->Fill(l1.Pt(),weight);
  hist_init[11]->Fill(l1.Eta(),weight);
  hist_init[12]->Fill(l1.Phi(),weight);

  hist_init[13]->Fill(l2.Pt(),weight);
  hist_init[14]->Fill(l2.Eta(),weight);
  hist_init[15]->Fill(l2.Phi(),weight);
#endif
  if (Jets.size()>0 && Jets[0].btag_DeepFlav>deep_btag_cut) {
    hist_init[16]->Fill(Jets[0].pt, weight);
    hist_init[17]->Fill(Jets[0].eta, weight);
    hist_init[18]->Fill(Jets[0].phi, weight);
  } else if (Jets.size()>1 && Jets[1].btag_DeepFlav>deep_btag_cut) {
    hist_init[16]->Fill(Jets[1].pt, weight);
    hist_init[17]->Fill(Jets[1].eta, weight);
    hist_init[18]->Fill(Jets[1].phi, weight);
  }

  // end //                                                                            
  //  if(gProofServ) {  str = TString::Format("check end evt %d ]]]",ievt);gProofServ->SendAsynMessage(str);	}*/
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
