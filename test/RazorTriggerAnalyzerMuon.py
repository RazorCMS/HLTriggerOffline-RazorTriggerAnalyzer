import FWCore.ParameterSet.Config as cms
import os
import sys

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
process.load("HLTriggerOffline.RazorTriggerAnalyzer.razorHemispheres_cff")

#hemispheres for calo jets passing ID
process.hemispheresCaloIDPassed = process.hemispheres.clone(inputTag = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"), minJetPt = cms.double(30), maxNJ = cms.int32(20))

#define input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(30000) )


myfilelist = cms.untracked.vstring()
myfilelist.extend(['root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_102.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_103.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_104.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_105.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_106.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_109.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_111.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_112.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_113.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_114.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_115.root',
    'root://xrootd.unl.edu:///store/user/duanders/QCD_Pt-300to470_Tune4C_13TeV_pythia8/crab_triggerStudies_QCD300to470_v1/141220_091517/0000/hlt_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_HLT_117.root',  
])
process.source = cms.Source("PoolSource",
                            fileNames = myfilelist
)

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razor_nocuts_QCD300to470_test.root"),
    closeFileFast = cms.untracked.bool(True)
)

#get global tag
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')
from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V3A::All')
process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'

#create AK4 charged-hadron subtracted jets
process.load("CommonTools.ParticleFlow.pfNoPileUpJME_cff")
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJets = ak4PFJets.clone()
process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfNoPileUpJME', doAreaFastjet = True)

#declare analyzer module
process.razorTriggerAnalysis = cms.EDAnalyzer("RazorTriggerAnalyzerMuon",
  trigSummary = cms.InputTag("hltTriggerSummaryAOD"),
  pfMETCollection = cms.InputTag("pfMet"),
  caloMETCollection = cms.InputTag("caloMet"),                                              
  hltMETCollection = cms.InputTag("hltMet", '', 'reHLT'),
  hltMETJetIDCollection = cms.InputTag("hltMetCleanUsingJetID", '', 'reHLT'),
  hltPFMETCollection = cms.InputTag("hltPFMETProducer"),                                              
  muonCollection = cms.InputTag('muons'),                                              
  pfJetCollection = cms.InputTag("ak4PFJetsCHS"),
  caloJetCollection = cms.InputTag("ak4CaloJets"),                                              
  hltCaloJetCollection = cms.InputTag("hltAK4CaloJetsCorrected"),
  hltCaloJetCollectionIDPassed = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"),
  hltPFJetCollection = cms.InputTag("hltAK4PFJetsCorrected"),                                              
  TriggerResults = cms.InputTag('TriggerResults','','reHLT'),
  TriggerPath = cms.string('HLT_RsqMR300_Rsq0p09_MR200_v1'),
  TriggerFilter = cms.InputTag('hltRsqMR300Rsq0p09MR200', '', 'reHLT'), #the last filter in the path
  #CaloFilter = cms.InputTag('hltRsqMRNoMinRsqNoMinMRNoMinCalo', '', 'reHLT'), #filter implementing cuts on calo MR and Rsq
  CaloFilter = cms.InputTag('hltRsqMR200Rsq0p01MR100Calo', '', 'reHLT'), #filter implementing cuts on calo MR and Rsq 
  hemispheres = cms.InputTag('hemispheres'),
  hemispheresCaloIDPassed = cms.InputTag('hemispheresCaloIDPassed')
  )

#define messagelogger (controls verbosity of the module)
process.MessageLogger = cms.Service("MessageLogger",
       destinations   = cms.untracked.vstring('detailedInfo','critical','cerr'),
       critical       = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
       detailedInfo   = cms.untracked.PSet(threshold  = cms.untracked.string('INFO') ),
       cerr           = cms.untracked.PSet(threshold  = cms.untracked.string('WARNING') )
)

process.run_module = cms.Path(process.pfNoPileUpJMESequence*process.ak4PFJets*process.ak4PFJetsCHS*cms.ignore(process.hemispheres)*cms.ignore(process.hemispheresCaloIDPassed)*process.pfMet*process.razorTriggerAnalysis)
