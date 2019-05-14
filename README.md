# Commands:
//CMSSW release
cmsrel CMSSW_10_2_7
cd CMSSW_10_2_7/src
cmsenv
git clone git@github.com:sandrofonseca/BanalysisMiniAOD.git B_analysis
cd B_analysis/BAnalyzer
scram b -j8
cmsRun python/dstard0_cfg.py


Grid and Crab Commands:

//Grid Environment:
voms-proxy-init --voms cms
//To copy file to your local
xrdcp root://cmsxrootd.fnal.gov//store/path_file /local/path
//setup CRAB3
source /cvmfs/cms.cern.ch/crab3/crab.sh
//Submit and Resubmit CRAB
crab submit -c crab_file.py
//Status CRAB
crab status -d path_project/
//Resubmit failed jobs
crab resubmit -d path_project/
//Get all jobs
crab getoutput -d path_project/
//Get specific jobs
crab getoutput -d path_project/ --jobids=1,2,...
//Kill all jobs
crab kill -d path_project/
