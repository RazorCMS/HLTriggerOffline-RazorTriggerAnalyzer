import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTING')

#load run conditions
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#get hemispheres and MET
process.load("RecoMET.METProducers.PFMET_cfi")
process.load("HLTriggerOffline.SUSYBSM.razorHemispheres_cff")

#define input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/a/anwang/work/HLT/CMSSW_7_2_1/src/HLTriggerOffline/SUSYBSM/test/samples/TTbarLepton_13_Razor7e33_1_1_dqQ.root',
#                                      'file:/afs/cern.ch/user/a/anwang/work/HLT/CMSSW_7_2_1/src/HLTriggerOffline/SUSYBSM/test/samples/TTbarLepton_13_Razor7e33_2_1_Oup.root',
#                                      'file:/afs/cern.ch/user/a/anwang/work/HLT/CMSSW_7_2_1/src/HLTriggerOffline/SUSYBSM/test/samples/TTbarLepton_13_Razor7e33_3_1_eiw.root',
#                                      'file:/afs/cern.ch/user/a/anwang/work/HLT/CMSSW_7_2_1/src/HLTriggerOffline/SUSYBSM/test/samples/TTbarLepton_13_Razor7e33_4_1_srx.root',
#                                      'file:/afs/cern.ch/user/a/anwang/work/HLT/CMSSW_7_2_1/src/HLTriggerOffline/SUSYBSM/test/samples/TTbarLepton_13_Razor7e33_5_1_IBN.root'
                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/user/amwang/11_Calo/TTbar_Lepton_New_Calo_1_1_wkN.root',
                                                              'root://xrootd.unl.edu//store/user/amwang/11_Calo/TTbar_Lepton_New_Calo_4_1_lMH.root',
                                                              'root://xrootd.unl.edu//store/user/amwang/11_Calo/TTbar_Lepton_New_Calo_2_1_KtV.root',
                                                              'root://xrootd.unl.edu//store/user/amwang/11_Calo/TTbar_Lepton_New_Calo_5_1_Gik.root',
                                                              'root://xrootd.unl.edu//store/user/amwang/11_Calo/TTbar_Lepton_New_Calo_3_1_Pmd.root'
    )
)

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razor_11_11_Calo_TTbar_Lepton.root"),
    closeFileFast = cms.untracked.bool(True)
)

#get global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

#create AK4 charged-hadron subtracted jets
process.load("CommonTools.ParticleFlow.pfNoPileUpJME_cff")
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJets = ak4PFJets.clone()
process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfNoPileUpJME', doAreaFastjet = True)

#declare analyzer module
process.razorTriggerAnalysis = cms.EDAnalyzer("RazorTriggerAnalyzer",
  trigSummary = cms.InputTag("hltTriggerSummaryAOD"),
  pfMETCollection = cms.InputTag("pfMet"),
  pfJetCollection = cms.InputTag("ak4PFJetsCHS"),
  TriggerResults = cms.InputTag('TriggerResults','','reHLT'),
  TriggerPath = cms.string('HLT_RsqMR300_Rsq0p09_MR200_v1'),
  TriggerFilter = cms.InputTag('hltRsqMR300Rsq0p09MR200', '', 'reHLT'), #the last filter in the path
  #CaloFilter = cms.InputTag('hltRsqMRNoMinRsqNoMinMRNoMinCalo', '', 'reHLT'), #filter implementing cuts on calo MR and Rsq
  CaloFilter = cms.InputTag('hltRsqMR200Rsq0p01MR100Calo', '', 'reHLT'), #filter implementing cuts on calo MR and Rsq 
  hemispheres = cms.InputTag('hemispheres')
  )

#define messagelogger (controls verbosity of the module)
process.MessageLogger = cms.Service("MessageLogger",
       destinations   = cms.untracked.vstring('detailedInfo','critical','cerr'),
       critical       = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
       detailedInfo   = cms.untracked.PSet(threshold  = cms.untracked.string('INFO') ),
       cerr           = cms.untracked.PSet(threshold  = cms.untracked.string('WARNING') )
)

process.run_module = cms.Path(process.pfNoPileUpJMESequence*process.ak4PFJets*process.ak4PFJetsCHS*cms.ignore(process.hemispheres)*process.pfMet*process.razorTriggerAnalysis)
