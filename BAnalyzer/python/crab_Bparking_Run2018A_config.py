from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = 'Bparking_Run2018A'
config.General.workArea = 'crab_Run2018A'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = '/eos/user/j/jordanm/MuonEff/CMSSW_10_2_3/src/MuonAnalysis/MuonFastFeedbackTools/python/AllJobsDATA2018_cfg.py'
config.JobType.psetName = '/afs/cern.ch/work/r/ragomesd/B_parking/CMSSW_10_2_7/src/BanalysisMiniAOD-1/BAnalyzer/python/dstard0_cfg.py'
#config.JobType.inputFiles='/eos/user/j/jordanm/MuonEff/CMSSW_10_2_3/src/MuonAnalysis'
config.JobType.outputFiles = ['D0DstarData.root']

config.Data.inputDBS = 'global'
config.Data.inputDataset = '/ParkingBPH1/Run2018A-14May2018-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2018B-PromptReco-v2/MINIAOD'

#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.runRange = '314472-325175' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/ragomesd' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag= 'SingleMu_Run2018A'
#config.Data.outputDatasetTag= 'SingleMu_Run2018B-PromptReco-v2'

#config.Site.whitelist = ['T3_US_FNALLPC']
#config.Site.storageSite     = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_US_Nebraska'
