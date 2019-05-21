import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = "101X_dataRun2_Prompt_v9"

#number of events (-1 -> all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
IgnoreCompletely = cms.untracked.vstring('ProductNotFound'),
SkipEvent = cms.untracked.vstring('Error: uninitialized ProxyBase used') )

#sourcefile = 'list259431.txt'

#USE THIS BLOCK FOR CONVENTIONAL FILE OPENING!===========
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/GEN-SIM-RECO/AsymptNoPUbx25Reco_MCRUN2_74_V9-v3/00000/04B705C0-5407-E511-AB52-00A0D1EE8A20.root'	
#'root://cmsxrootd.fnal.gov//store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/129EAF02-B85A-E811-AC4E-0CC47A1DF818.root'
#'root://cmsxrootd.fnal.gov//store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/02992D8C-205A-E811-8049-0025905C54C6.root'
#'file:129EAF02-B85A-E811-AC4E-0CC47A1DF818.root' #Data mniAOD
'file:EA334703-C759-E811-B440-008CFAE45030.root'
	   )
)
#===================================================================

process.analysis = cms.EDAnalyzer('DstarD0TTree',
    # Analysis
    doMC=cms.bool(True),
    doRec=cms.bool(True),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    #PathName = cms.untracked.string("HLT_Mu12_IP6_part0_v"), #ParkingBPH1 ######
    #PathName = cms.untracked.string("HLT_Mu7_IP4_part0_v"), #ParkingBPH1  ########
    #PathName = cms.untracked.string("HLT_Mu8_IP3_part0_v"), #ParkingBPH1  v1
    #PathName = cms.untracked.string("HLT_Mu8_IP5_part0_v"), #ParkingBPH1 ########
    #PathName = cms.untracked.string("HLT_Mu8_IP6_part0_v"), #ParkingBPH1 ########
    #PathName = cms.untracked.string("HLT_Mu9_IP0_part0_v"), #ParkingBPH1 ########
    #PathName = cms.untracked.string("HLT_Mu9_IP3_part0_v"), #ParkingBPH1 #######
    #PathName = cms.untracked.string("HLT_Mu9_IP4_part0_v"), #ParkingBPH1 #######
    #PathName = cms.untracked.string("HLT_Mu9_IP5_part0_v"), #ParkingBPH1 ########
    #PathName = cms.untracked.string("HLT_Mu9_IP5_part0_v"), #ParkingBPH1 ######## 
    PathName = cms.untracked.string("HLT_Mu9_IP6_part0_v"),  #ParkingBPH1 ----

    tracks = cms.InputTag('packedPFCandidates'),#Minia AOD
    recVtxs = cms.InputTag('offlineSlimmedPrimaryVertices'), #Mini AOD
    genParticles = cms.InputTag('genParticles'),
    ParticleFlowTag = cms.InputTag("particleFlow"),
    # Options
    comEnergy = cms.double(13000.),
    TTBIt = cms.int32(34),
    debug = cms.untracked.bool(False),
    SaveROOTTree = cms.untracked.bool(True)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('D0DstarData.root')
)

process.p = cms.Path(process.analysis)
