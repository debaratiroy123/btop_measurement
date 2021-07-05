#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/deroy/x509up_u56660
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 
#eval `scramv1 runtime -sh`
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
cd /home/deroy/t3store3/Muon_MuEl/SameWeight/Combination/ULCombination
root -l -b -q runAll.C

