from CRABClient.UserUtilities import config
config = config()
#config.General.requestName = 'GEN-SIM_Production_DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS_v02'
config.General.requestName = 'DR_Production_DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS_v02'
config.General.workArea = 'crab_projects'
config.JobType.pluginName = 'PrivateMC' #'Analysis'
#config.JobType.psetName = 'DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS.py'
config.JobType.psetName = 'DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_DR1.py'
#config.Data.inputDataset = 'NONE'
config.Data.inputDataset = '/CRAB_PrivateMC/ragomesd-GEN-SIM_Production_DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS_v02-ad1b94f21078c8c947f7844ad1bb0e2a/USER'
#config.Data.inputDBS = 'phys03'
config.Data.publishDBS = 'phys03'
config.Data.splitting = 'FileBased' #'EventBased'  #'FileBased'  #'EventBased'
config.Data.unitsPerJob = 2000
NJOBS = 25  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/ragomesd' # or '/store/group/<subdir>'
config.General.transferOutputs = True
config.Data.publication = True
config.Data.outputPrimaryDataset = 'CRAB_PrivateMC'
#config.Data.outputDatasetTag = 'GEN-SIM_Production_DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS_v02'
config.Data.outputDatasetTag = 'DR_Production_DStarToD0Pi_D0KPi_DStarNoPtFilter_TuneCP5_13TeV-pythia8-evtgen_GS_v02'
#config.Data.ignoreLocality = True
#config.Site.storageSite = 'T2_CH_CERNBOX'
config.Site.storageSite = 'T2_US_Nebraska'
#config.Site.whitelist = ['T2_US_Nebraska']

