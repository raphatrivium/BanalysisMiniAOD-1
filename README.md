# Commands:

cmsrel CMSSW_10_2_7

cd CMSSW_10_2_7/src


cmsenv


git clone git@github.com:sandrofonseca/BanalysisMiniAOD.git B_analysis


cd B_analysis/BAnalyzer


scram b -j8

cmsRun python/dstard0_cfg.py


