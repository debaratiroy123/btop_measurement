#!/bin/sh
cd /home/chatterj/t3store/WPrime/CMSSW_8_0_29/src/Analysis_Run/Run2016G/ 
export X509_USER_PROXY=/home/chatterj/x509up_u56530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc530
source $VO_CMS_SW_DIR/cmsset_default.sh 
eval `scramv1 runtime -sh`
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
root -l -b -q runAll.C
#./Wprime_TrigEffB.exe
