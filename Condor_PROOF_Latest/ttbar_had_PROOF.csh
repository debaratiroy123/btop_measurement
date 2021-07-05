#!/bin/sh
#cd /home/chatterj/t3store/WPrime/CMSSW_8_0_29/src/Analysis_Run/Run2016G/ 
cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#eval `cmsenv`
export X509_USER_PROXY=/home/deroy/x509up_u56660
#/home/chatterj/x509up_u56530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
#slc7_amd64_gcc820
source $VO_CMS_SW_DIR/cmsset_default.sh 
eval `scramv1 runtime -sh`
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
root -l -b -q runAll.C
#./Wprime_TrigEffB.exe
