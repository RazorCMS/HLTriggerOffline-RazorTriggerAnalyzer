#ifndef RazorTriggerAnalyzerMuon_H
#define RazorTriggerAnalyzerMuon_H

//ROOT
#include "TTree.h"

//event
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// MET
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

// Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

//Hemispheres
#include "HLTrigger/JetMET/interface/HLTRHemisphere.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h"

using namespace edm;

class RazorTriggerAnalyzerMuon: public EDAnalyzer{

  public:
  RazorTriggerAnalyzerMuon(const edm::ParameterSet& ps);
  static double CalcMR(TLorentzVector ja,TLorentzVector jb);
  static double CalcR(double MR, TLorentzVector ja,TLorentzVector jb, edm::Handle<edm::View<reco::MET> > met, const std::vector<math::XYZTLorentzVector>& muons);
  virtual ~RazorTriggerAnalyzerMuon();

  protected:
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup);
  void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) ;
  void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
  void endRun(edm::Run const& run, edm::EventSetup const& eSetup);

  private:

  //variables from config file
  edm::EDGetTokenT<edm::View<reco::MET> > thePfMETCollection_;
  edm::EDGetTokenT<edm::View<reco::MET> > theCaloMETCollection_;
  edm::EDGetTokenT<edm::View<reco::MET> > theHLTMETCollection_;
  edm::EDGetTokenT<edm::View<reco::MET> > theHLTMETJetIDCollection_;
  edm::EDGetTokenT<edm::View<reco::MET> > theHLTPfMETCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
  edm::EDGetTokenT<trigger::TriggerEvent> theTrigSummary_;
  edm::EDGetTokenT<reco::PFJetCollection> thePfJetCollection_;
  edm::EDGetTokenT<reco::CaloJetCollection> theCaloJetCollection_;
  edm::EDGetTokenT<reco::CaloJetCollection> theHLTCaloJetCollection_;
  edm::EDGetTokenT<reco::PFJetCollection> theHLTPFJetCollection_;
  edm::EDGetTokenT<reco::MuonCollection> theMuonCollection_;
  edm::EDGetTokenT<std::vector<math::XYZTLorentzVector> > theHemispheres_;

  std::string triggerPath_;
  edm::InputTag triggerFilter_;
  edm::InputTag caloFilter_;


  //tree for output
  TTree *outTree;

  //variables for tree
  double MR, Rsq;
  double onlineMR, onlineRsq;
  double caloMR, caloRsq;
  double pfHT, pfMET, hltPFMETProducer;
  double caloMET, hltMET, hltMETJetID;
  bool hasFired;
  bool denomFired;
  int numMuons;
  int numMuonsPassed30;
  int numCaloJetsPassed30;
  int numHLTCaloJetsPassed30;
  int numHLTPFJetsPassed30;
  double muonET;
  bool passedCaloDiJetCut;
  bool passedPFDiJetCut;
};

#endif
