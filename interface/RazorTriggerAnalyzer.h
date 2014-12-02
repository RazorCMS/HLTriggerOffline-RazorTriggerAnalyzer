#ifndef RazorTriggerAnalyzer_H
#define RazorTriggerAnalyzer_H

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

// Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

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

class RazorTriggerAnalyzer: public EDAnalyzer{

  public:
  RazorTriggerAnalyzer(const edm::ParameterSet& ps);
  static double CalcMR(TLorentzVector ja,TLorentzVector jb);
  static double CalcR(double MR, TLorentzVector ja,TLorentzVector jb, edm::Handle<edm::View<reco::MET> > met, const std::vector<math::XYZTLorentzVector>& muons);
  virtual ~RazorTriggerAnalyzer();

  protected:
  void analyze(edm::Event const& e, edm::EventSetup const& eSetup);
  void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup) ;
  void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
  void endRun(edm::Run const& run, edm::EventSetup const& eSetup);

  private:

  //variables from config file
  edm::EDGetTokenT<edm::View<reco::MET> > thePfMETCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
  edm::EDGetTokenT<trigger::TriggerEvent> theTrigSummary_;
  edm::EDGetTokenT<reco::PFJetCollection> thePfJetCollection_;
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
  double pfHT, pfMET;
  bool hasFired;
  bool denomFired;
};

#endif
