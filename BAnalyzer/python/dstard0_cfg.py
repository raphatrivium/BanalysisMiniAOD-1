import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v9"

#number of events (-1 -> all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
IgnoreCompletely = cms.untracked.vstring('ProductNotFound'),
SkipEvent = cms.untracked.vstring('Error: uninitialized ProxyBase used') )

#sourcefile = 'list259431.txt'

#USE THIS BLOCK FOR CONVENTIONAL FILE OPENING!===========
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/GEN-SIM-RECO/AsymptNoPUbx25Reco_MCRUN2_74_V9-v3/00000/04B705C0-5407-E511-AB52-00A0D1EE8A20.root'	
'root://cmsxrootd.fnal.gov//store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/129EAF02-B85A-E811-AC4E-0CC47A1DF818.root' #Data mniAOD
	   )
)
#===================================================================



#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall15DR76/DStarToD0Pi_D0KPi_DStarFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/001496DA-0FD4-E511-9DA9-02163E00EB62.root'))

#)
#)


process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

process.L1T1coll=process.hltLevel1GTSeed.clone()
process.L1T1coll.L1TechTriggerSeeding = cms.bool(True)
process.L1T1coll.L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')


process.trigger = cms.EDFilter("HLTHighLevel",
TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  HLTPaths = cms.vstring('HLT_MinBias'), # provide list of HLT paths (or patterns) you want
  eventSetupPathsKey = cms.string(''),
  andOr              = cms.bool(True),
  throw              = cms.bool(False),
#  saveTags           = cms.bool(False)

)



process.analysis = cms.EDAnalyzer('DstarD0TTree',
    # Analysis
    doMC=cms.bool(True),
    doRec=cms.bool(True),
    #tracks = cms.InputTag('generalTracks'),#AOD
    tracks = cms.InputTag('packedPFCandidates'),#Minia AOD
    #recVtxs = cms.InputTag('offlinePrimaryVertices'), #AOD
    recVtxs = cms.InputTag('offlineSlimmedPrimaryVertices'), #Mini AOD
    #genParticles = cms.InputTag('genParticleCollection'),
	genParticles = cms.InputTag('genParticles'),
    ParticleFlowTag = cms.InputTag("particleFlow"),
    # Options
    comEnergy = cms.double(13000.),
#    HLTPath = cms.string("HLT_L1_BscMinBiasOR_BptxPlusORMinus"), 
    HLTPath = cms.string("HLT_MinBias"),  
    TTBIt = cms.int32(34),
    SaveROOTTree = cms.untracked.bool(True)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('D0DstarMCAnalysis.root')
)
#---------------------
#process.p = cms.Path(process.L1T1coll+process.trigger+process.primaryVertexFilter+process.noscraping+process.analysis)

process.p = cms.Path(process.analysis)
