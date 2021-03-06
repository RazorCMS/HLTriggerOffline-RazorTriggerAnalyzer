#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTriggerOffline/RazorTriggerAnalyzer/interface/RazorTriggerAnalyzer.h"

RazorTriggerAnalyzer::RazorTriggerAnalyzer(const edm::ParameterSet& ps)
{
  edm::LogInfo("RazorTriggerAnalyzer") << "Constructor RazorTriggerAnalyzer::RazorTriggerAnalyzer " << std::endl;
  // Get parameters from configuration file
  theTrigSummary_ = consumes<trigger::TriggerEvent>(ps.getParameter<edm::InputTag>("trigSummary"));
  thePfMETCollection_ = consumes<edm::View<reco::MET> >(ps.getParameter<edm::InputTag>("pfMETCollection"));
  theHemispheres_ = consumes<std::vector<math::XYZTLorentzVector> >(ps.getParameter<edm::InputTag>("hemispheres"));
  triggerResults_ = consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("TriggerResults"));
  triggerPath_ = ps.getParameter<std::string>("TriggerPath");
  triggerFilter_ = ps.getParameter<edm::InputTag>("TriggerFilter");
  caloFilter_ = ps.getParameter<edm::InputTag>("CaloFilter");
  thePfJetCollection_ = consumes<reco::PFJetCollection>(ps.getParameter<edm::InputTag>("pfJetCollection"));

  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  outTree = fs->make<TTree>("TriggerInfo", "Razor variable info");

  outTree->Branch("MR", &MR, "MR/D"); //offline values
  outTree->Branch("Rsq", &Rsq, "Rsq/D");
  outTree->Branch("onlineMR", &onlineMR, "onlineMR/D");
  outTree->Branch("onlineRsq", &onlineRsq, "onlineRsq/D");
  outTree->Branch("caloMR", &caloMR, "caloMR/D");
  outTree->Branch("caloRsq", &caloRsq, "caloRsq/D");
  outTree->Branch("pfHT", &pfHT, "pfHT/D");
  outTree->Branch("pfMET", &pfMET, "pfMET/D");
  outTree->Branch("hasFired", &hasFired, "hasFired/O");
  outTree->Branch("denomFired", &denomFired, "denomFired/O");
}

RazorTriggerAnalyzer::~RazorTriggerAnalyzer()
{
   edm::LogInfo("RazorTriggerAnalyzer") << "Destructor RazorTriggerAnalyzer::~RazorTriggerAnalyzer " << std::endl;
}

void RazorTriggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg,
  edm::EventSetup const& context)
{
   edm::LogInfo("RazorTriggerAnalyzer") << "RazorTriggerAnalyzer::beginLuminosityBlock" << std::endl;
}

void RazorTriggerAnalyzer::analyze(edm::Event const& e, edm::EventSetup const& eSetup){
    
  edm::LogInfo("RazorTriggerAnalyzer") << "RazorTriggerAnalyzer::analyze" << std::endl;

  using namespace std;
  using namespace edm;
  using namespace reco;

  //reset tree variables
  MR = -999.; 
  Rsq = -999.;
  onlineMR = -999.;//online razor variables computed using PF quantities
  onlineRsq = -999.;
  caloMR = -999.;//online razor variables computed using calo quantities
  caloRsq = -999.;
  pfHT = 0;
  pfMET = -999.;
  hasFired = false;
  denomFired = false;

  // get hemispheres
  Handle< vector<math::XYZTLorentzVector> > hemispheres;
  e.getByToken (theHemispheres_,hemispheres);

  // get the MET Collection
  edm::Handle<edm::View<reco::MET> > inputMet;
  e.getByToken(thePfMETCollection_,inputMet);
  if ( !inputMet.isValid() ){
    edm::LogError ("RazorTriggerAnalyzer") << "invalid collection: PFMET" << "\n";
   return;
  }

  //get the jet collection
  edm::Handle<reco::PFJetCollection> pfJetCollection;
  e.getByToken (thePfJetCollection_,pfJetCollection);
  if ( !pfJetCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzer") << "invalid collection: PFJets" << "\n";
    return;
  }

  //check what is in the menu
  edm::Handle<edm::TriggerResults> hltresults;
  e.getByToken(triggerResults_,hltresults);
  if(!hltresults.isValid()){
    edm::LogError ("RazorTriggerAnalyzer") << "invalid collection: TriggerResults" << "\n";
    return;
  }
  
  //get the trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerSummary;
  e.getByToken(theTrigSummary_, triggerSummary);
  if(!triggerSummary.isValid()) {
    edm::LogError ("RazorTriggerAnalyzer") << "invalid collection: TriggerSummary" << "\n";
    return;
  }
  
  //get online objects
  //HLTriggerOffline/Egamma/python/TriggerTypeDefs.py contains the trigger object IDs

  //find the indices of the filters that store the online PF/Calo MR and Rsq
  size_t filterIndex = triggerSummary->filterIndex( triggerFilter_ );
  size_t caloFilterIndex = triggerSummary->filterIndex( caloFilter_ );
  
  //search for online MR and Rsq objects
  trigger::TriggerObjectCollection triggerObjects = triggerSummary->getObjects();

  if( !(filterIndex >= triggerSummary->sizeFilters()) ){
      const trigger::Keys& keys = triggerSummary->filterKeys( filterIndex );
      for( size_t j = 0; j < keys.size(); ++j ){
          trigger::TriggerObject foundObject = triggerObjects[keys[j]];
          if(foundObject.id() == 0){ //the MET object containing MR and Rsq will show up with ID = 0
              onlineMR = foundObject.px(); //razor variables stored in dummy reco::MET objects
              onlineRsq = foundObject.py();
          }
      }
  }

  //search for calo MR and Rsq objects
  if( !(caloFilterIndex >= triggerSummary->sizeFilters()) ){
      const trigger::Keys& keys = triggerSummary->filterKeys( caloFilterIndex );
      for( size_t j = 0; j < keys.size(); ++j ){
          trigger::TriggerObject foundObject = triggerObjects[keys[j]];
          if(foundObject.id() == 0){ 
	    cout << "found the object?" ; // print statement here
              caloMR = foundObject.px(); //razor variables stored in dummy reco::MET objects
              caloRsq = foundObject.py();
          }
      }
  }
  std::string denomPath = "HLT_ZeroBias_v1"; //reference trigger
  //std::string denomPath = "HLT_Ele27_eta2p1_WP85_Gsf_v1"; //trigger path used as a reference for computing turn-ons and efficiencies
  const edm::TriggerNames& trigNames = e.triggerNames(*hltresults);
  unsigned int numTriggers = trigNames.size();
  //loop over fired triggers and look for razor trigger and/or reference electron trigger
  for( unsigned int hltIndex=0; hltIndex<numTriggers; ++hltIndex ){
      if (trigNames.triggerName(hltIndex)==triggerPath_ && hltresults->wasrun(hltIndex) && hltresults->accept(hltIndex)) hasFired = true;
      if (trigNames.triggerName(hltIndex)==denomPath && hltresults->wasrun(hltIndex) && hltresults->accept(hltIndex)) denomFired = true;
  }

  //compute PF HT
  for (reco::PFJetCollection::const_iterator i_pfjet = pfJetCollection->begin(); i_pfjet != pfJetCollection->end(); ++i_pfjet){
      if (i_pfjet->pt() < 40) continue;
      if (fabs(i_pfjet->eta()) > 3.0) continue;
      pfHT += i_pfjet->pt();
  }

  pfMET = (inputMet->front()).pt();

  //this part is adapted from HLTRFilter.cc 

  // check that the input collections are available
  if (not hemispheres.isValid()){
      edm::LogError("RazorTriggerAnalyzer") << "Hemisphere object is invalid!" << "\n";
      return;
  }

  if(hasFired && denomFired){

      if(hemispheres->size() ==0){  // the Hemisphere Maker will produce an empty collection of hemispheres if the number of jets in the
          edm::LogError("RazorTriggerAnalyzer") << "Cannot calculate M_R and R^2 because there are too many jets! (trigger passed automatically without forming the hemispheres)" << endl;
          return;
      }

      //This looks unusual because of how HLTRHemisphere deals with muons -- go look there for more information.  In the current razor triggers, muons are added to the hemispheres (no special treatment), and hemispheres->size() = 2.
      if(hemispheres->size() != 0 && hemispheres->size() != 2 && hemispheres->size() != 5 && hemispheres->size() != 10){
          edm::LogError("RazorTriggerAnalyzer") << "Invalid hemisphere collection!  hemispheres->size() = " << hemispheres->size() << endl;
          return;
      }

      TLorentzVector ja(hemispheres->at(0).x(),hemispheres->at(0).y(),hemispheres->at(0).z(),hemispheres->at(0).t());
      TLorentzVector jb(hemispheres->at(1).x(),hemispheres->at(1).y(),hemispheres->at(1).z(),hemispheres->at(1).t());

      //dummy vector (this trigger does not care about muons)
      std::vector<math::XYZTLorentzVector> muonVec;

      // where we calculate M_R and R^2 

      MR = CalcMR(ja,jb);
      double R  = CalcR(MR,ja,jb,inputMet,muonVec);
      Rsq = R*R;    

      outTree->Fill();
  } 
  else if(denomFired){ //calculate M_R and R^2 for the denominator histograms

      if(hemispheres->size() ==0){  // the Hemisphere Maker will produce an empty collection of hemispheres if the number of jets in the event is larger than the threshold.  In this case we cannot compute razor variables
          return;
      }

      if(hemispheres->size() != 0 && hemispheres->size() != 2 && hemispheres->size() != 5 && hemispheres->size() != 10){
          return;
      }

      TLorentzVector ja(hemispheres->at(0).x(),hemispheres->at(0).y(),hemispheres->at(0).z(),hemispheres->at(0).t());
      TLorentzVector jb(hemispheres->at(1).x(),hemispheres->at(1).y(),hemispheres->at(1).z(),hemispheres->at(1).t());
      //dummy vector (this trigger does not care about muons)
      std::vector<math::XYZTLorentzVector> muonVec;

      MR = CalcMR(ja,jb);
      double R  = CalcR(MR,ja,jb,inputMet,muonVec);
      Rsq = R*R;    

      outTree->Fill();
  }
}

void RazorTriggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup)
{
  edm::LogInfo("RazorTriggerAnalyzer") << "RazorTriggerAnalyzer::endLuminosityBlock" << std::endl;
}


void RazorTriggerAnalyzer::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("RazorTriggerAnalyzer") << "RazorTriggerAnalyzer::endRun" << std::endl;
}

//CalcMR and CalcR borrowed from HLTRFilter.cc
double 
RazorTriggerAnalyzer::CalcMR(TLorentzVector ja, TLorentzVector jb){
  if(ja.Pt()<=0.1) return -1;

  ja.SetPtEtaPhiM(ja.Pt(),ja.Eta(),ja.Phi(),0.0);
  jb.SetPtEtaPhiM(jb.Pt(),jb.Eta(),jb.Phi(),0.0);
  
  if(ja.Pt() > jb.Pt()){
    TLorentzVector temp = ja;
    ja = jb;
    jb = temp;
  }
  
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();
  
  double MR = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
		   (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
  
  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  
  double mygamma = 1./sqrt(1.-mybeta*mybeta);
  
  //use gamma times MRstar
  return MR*mygamma;  
}

double 
  RazorTriggerAnalyzer::CalcR(double MR, TLorentzVector ja, TLorentzVector jb, edm::Handle<edm::View<reco::MET> > inputMet, const std::vector<math::XYZTLorentzVector>& muons){
  //now we can calculate MTR
  TVector3 met;
  met.SetPtEtaPhi((inputMet->front()).pt(),0.0,(inputMet->front()).phi());
  
  std::vector<math::XYZTLorentzVector>::const_iterator muonIt;
  for(muonIt = muons.begin(); muonIt!=muons.end(); muonIt++){
    TVector3 tmp;
    tmp.SetPtEtaPhi(muonIt->pt(),0,muonIt->phi());
    met-=tmp;
  }

  double MTR = sqrt(0.5*(met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect())));
  
  //filter events
  return float(MTR)/float(MR); //R
  
}

 //define this as a plug-in
DEFINE_FWK_MODULE(RazorTriggerAnalyzer);
