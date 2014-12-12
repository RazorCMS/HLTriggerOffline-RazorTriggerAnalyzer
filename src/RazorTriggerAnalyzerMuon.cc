#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTriggerOffline/RazorTriggerAnalyzer/interface/RazorTriggerAnalyzerMuon.h"

RazorTriggerAnalyzerMuon::RazorTriggerAnalyzerMuon(const edm::ParameterSet& ps)
{
  edm::LogInfo("RazorTriggerAnalyzerMuon") << "Constructor RazorTriggerAnalyzerMuon::RazorTriggerAnalyzerMuon " << std::endl;
  // Get parameters from configuration file
  theTrigSummary_ = consumes<trigger::TriggerEvent>(ps.getParameter<edm::InputTag>("trigSummary"));
  thePfMETCollection_ = consumes<edm::View<reco::MET> >(ps.getParameter<edm::InputTag>("pfMETCollection"));
  theCaloMETCollection_ = consumes<edm::View<reco::MET> >(ps.getParameter<edm::InputTag>("caloMETCollection"));
  theHemispheres_ = consumes<std::vector<math::XYZTLorentzVector> >(ps.getParameter<edm::InputTag>("hemispheres"));
  triggerResults_ = consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("TriggerResults"));
  triggerPath_ = ps.getParameter<std::string>("TriggerPath");
  triggerFilter_ = ps.getParameter<edm::InputTag>("TriggerFilter");
  caloFilter_ = ps.getParameter<edm::InputTag>("CaloFilter");
  thePfJetCollection_ = consumes<reco::PFJetCollection>(ps.getParameter<edm::InputTag>("pfJetCollection"));
  theCaloJetCollection_ = consumes<reco::CaloJetCollection>(ps.getParameter<edm::InputTag>("caloJetCollection"));
  theHLTCaloJetCollection_ = consumes<reco::CaloJetCollection>(ps.getParameter<edm::InputTag>("hltCaloJetCollection"));
  theHLTPFJetCollection_ = consumes<reco::PFJetCollection>(ps.getParameter<edm::InputTag>("hltPFJetCollection"));
  theMuonCollection_    = consumes<reco::MuonCollection>(ps.getParameter<edm::InputTag>("muonCollection"));
  theHLTMETCollection_ = consumes<edm::View<reco::MET> > (ps.getParameter<edm::InputTag>("hltMETCollection"));
  theHLTMETJetIDCollection_ = consumes<edm::View<reco::MET> > (ps.getParameter<edm::InputTag>("hltMETJetIDCollection"));
  theHLTPfMETCollection_ = consumes<edm::View<reco::MET> > (ps.getParameter<edm::InputTag>("hltPFMETCollection"));

  //theMuonCollection_ = ps.getParameter<edm::InputTag>("muonCollection");

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
  outTree->Branch("caloMET", &caloMET, "caloMET/D");
  outTree->Branch("hltMET", &hltMET, "hltMET/D");
  outTree->Branch("hltMETJetID", &hltMETJetID, "hltMETJetID/D");
  outTree->Branch("hltPFMETProducer", &hltPFMETProducer, "hltPFMETProducer/D");
  outTree->Branch("hasFired", &hasFired, "hasFired/O");
  outTree->Branch("denomFired", &denomFired, "denomFired/O");
  outTree->Branch("numMuons", &numMuons, "numMuons/I");
  outTree->Branch("numMuonsPassed30", &numMuonsPassed30, "numMuonsPassed30/I");
  outTree->Branch("numCaloJetsPassed30", &numCaloJetsPassed30, "numCaloJetsPassed30/I");
  outTree->Branch("numHLTCaloJetsPassed30", &numHLTCaloJetsPassed30, "numHLTCaloJetsPassed30/I");
  outTree->Branch("numHLTPFJetsPassed30", &numHLTPFJetsPassed30, "numHLTPFJetsPassed30/I");
  outTree->Branch("muonET", &muonET, "muonET/D");
  outTree->Branch("passedCaloDiJetCut", &passedCaloDiJetCut, "passedCaloDiJetCut/O");
  outTree->Branch("passedPFDiJetCut", &passedPFDiJetCut, "passedPFDiJetCut/O");

}

RazorTriggerAnalyzerMuon::~RazorTriggerAnalyzerMuon()
{
   edm::LogInfo("RazorTriggerAnalyzerMuon") << "Destructor RazorTriggerAnalyzerMuon::~RazorTriggerAnalyzerMuon " << std::endl;
}

void RazorTriggerAnalyzerMuon::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg,
  edm::EventSetup const& context)
{
   edm::LogInfo("RazorTriggerAnalyzerMuon") << "RazorTriggerAnalyzerMuon::beginLuminosityBlock" << std::endl;
}

void RazorTriggerAnalyzerMuon::analyze(edm::Event const& e, edm::EventSetup const& eSetup){
    
  edm::LogInfo("RazorTriggerAnalyzerMuon") << "RazorTriggerAnalyzerMuon::analyze" << std::endl;

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
  caloMET = -999.;
  hltMET = -999.;
  hltMETJetID = -999.;
  hltPFMETProducer = -999.;
  hasFired = false;
  denomFired = false;
  numMuons = 0;
  numMuonsPassed30 = 0;
  numCaloJetsPassed30 = 0;
  numHLTCaloJetsPassed30 = 0;
  numHLTPFJetsPassed30 = 0;
  muonET = 0;
  passedCaloDiJetCut = false;
  passedPFDiJetCut = false;


  // get hemispheres
  Handle< vector<math::XYZTLorentzVector> > hemispheres;
  e.getByToken (theHemispheres_,hemispheres);

  // get the MET Collection
  edm::Handle<edm::View<reco::MET> > inputMet;
  e.getByToken(thePfMETCollection_,inputMet);
  if ( !inputMet.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: PFMET" << "\n";
   return;
  }

  // get the Calo MET Collection
  edm::Handle<edm::View<reco::MET> > inputCaloMet;
  e.getByToken(theCaloMETCollection_,inputCaloMet);
  if ( !inputCaloMet.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: CaloMET" << "\n";
    return;
  }

  // get the HLT MET Collection
  edm::Handle<edm::View<reco::MET> > inputHLTMet;
  e.getByToken(theHLTMETCollection_,inputHLTMet);
  if ( !inputHLTMet.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: HLTMET" << "\n";
    return;
  }

  // get the HLT MET JetID Collection                                                                                   
  edm::Handle<edm::View<reco::MET> > inputHLTMetJetID;
  e.getByToken(theHLTMETJetIDCollection_,inputHLTMetJetID);
  if ( !inputHLTMetJetID.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: HLTMETJetID" << "\n";
    return;
  }

  // get the HLT PFMET Collection                                                                                                                            
  edm::Handle<edm::View<reco::MET> > inputHLTPFMet;
  e.getByToken(theHLTPfMETCollection_,inputHLTPFMet);
  if ( !inputHLTPFMet.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: HLTPFMET" << "\n";
    return;
  }

  //get the jet collection
  edm::Handle<reco::PFJetCollection> pfJetCollection;
  e.getByToken (thePfJetCollection_,pfJetCollection);
  if ( !pfJetCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: PFJets" << "\n";
    return;
  }

  //get calo jet collection
  edm::Handle<reco::CaloJetCollection> caloJetCollection;
  e.getByToken (theCaloJetCollection_,caloJetCollection);
  if ( !caloJetCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: CaloJets" << "\n";
    return;
  }

  //get hlt calojet collection
  edm::Handle<reco::CaloJetCollection> hltCaloJetCollection;
  e.getByToken (theHLTCaloJetCollection_,hltCaloJetCollection);
  if ( !hltCaloJetCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: HLT CaloJets" << "\n";
    return;
  }

  //get hlt pfjet collection                                                                                                                        
  edm::Handle<reco::PFJetCollection> hltPFJetCollection;
  e.getByToken (theHLTPFJetCollection_,hltPFJetCollection);
  if ( !hltPFJetCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: HLT PFJets" << "\n";
    return;
  }

  //get the muon collection
  edm::Handle<reco::MuonCollection> muonCollection;
  e.getByToken(theMuonCollection_, muonCollection);
  if ( !muonCollection.isValid() ){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: Muons" << "\n";
    return;
  }

  //check what is in the menu
  edm::Handle<edm::TriggerResults> hltresults;
  e.getByToken(triggerResults_,hltresults);
  if(!hltresults.isValid()){
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: TriggerResults" << "\n";
    return;
  }
  
  //get the trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerSummary;
  e.getByToken(theTrigSummary_, triggerSummary);
  if(!triggerSummary.isValid()) {
    edm::LogError ("RazorTriggerAnalyzerMuon") << "invalid collection: TriggerSummary" << "\n";
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
  //std::string denomPath = "HLT_ZeroBias_v1"; //reference trigger
  std::string denomPath = "HLT_Ele27_eta2p1_WP85_Gsf_v1"; //trigger path used as a reference for computing turn-ons and efficiencies
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
  caloMET = (inputCaloMet->front()).pt(); // check this
  hltMET = (inputHLTMet->front()).pt(); // check this  
  hltMETJetID = (inputHLTMetJetID->front()).pt(); // check this                                                 
  hltPFMETProducer = (inputHLTPFMet->front()).pt();

  // count muons
  for (reco::MuonCollection::const_iterator i_muon = muonCollection->begin(); i_muon != muonCollection->end(); ++i_muon){
    muonET += i_muon->pt();
    //if (fabs(i_muon->eta()) > 3.0) continue;
    numMuons += 1;
    if (i_muon->pt() < 30) continue;
    if (fabs(i_muon->eta()) > 3.0) continue;
    numMuonsPassed30 += 1;
  }

  //number of calojets that pass 30 Gev and 3.0 eta cut
  for (reco::CaloJetCollection::const_iterator i_calojet = caloJetCollection->begin(); i_calojet != caloJetCollection->end(); ++i_calojet){
    if (i_calojet->pt() < 30) continue;
    if (fabs(i_calojet->eta()) > 3.0) continue;
    numCaloJetsPassed30 += 1;
  }

  //number of hlt calojets that pass 30 Gev and 3.0 eta cut                                                                                                  
  bool lead_hltcalo = false;
  int passed60 = 0;
  bool sublead_hltcalo = false;
  for (reco::CaloJetCollection::const_iterator i_hltcalojet = hltCaloJetCollection->begin(); i_hltcalojet != hltCaloJetCollection->end(); ++i_hltcalojet){
    if (i_hltcalojet->pt() < 30) continue;
    if (fabs(i_hltcalojet->eta()) > 3.0) continue;    
    // requirements for passing dijet cut
    if (i_hltcalojet->pt() > 60) passed60++;
    if (i_hltcalojet->pt() > 70) lead_hltcalo = true;
    numHLTCaloJetsPassed30 += 1;
  }
  cout << "NumberHLTCaloJets: " << numHLTCaloJetsPassed30 << endl;
  if (passed60 > 1) sublead_hltcalo = true;
  if ((lead_hltcalo == true) && (sublead_hltcalo == true)) passedCaloDiJetCut = true;


  //number of hlt pfjets that pass 30 Gev and 3.0 eta cut                                                                                         
                                                                                                                                                    
  bool lead_hltpf = false;
  int pfpassed60 = 0;
  bool sublead_hltpf = false;
  for (reco::PFJetCollection::const_iterator i_hltpfjet = hltPFJetCollection->begin(); i_hltpfjet != hltPFJetCollection->end(); ++i_hltpfjet){
    if (i_hltpfjet->pt() < 30) continue;
    if (fabs(i_hltpfjet->eta()) > 3.0) continue;
    // requirements for passing dijet cut                                                                                                            
    if (i_hltpfjet->pt() > 60) pfpassed60++;
    if (i_hltpfjet->pt() > 70) lead_hltpf = true;
    numHLTPFJetsPassed30 += 1;
  }
  if (pfpassed60 > 1) sublead_hltpf = true;
  if ((lead_hltpf == true) && (sublead_hltpf == true)) passedPFDiJetCut = true;

  //this part is adapted from HLTRFilter.cc 

  // check that the input collections are available
  if (not hemispheres.isValid()){
      edm::LogError("RazorTriggerAnalyzerMuon") << "Hemisphere object is invalid!" << "\n";
      return;
  }

  if(hasFired && denomFired){

      if(hemispheres->size() ==0){  // the Hemisphere Maker will produce an empty collection of hemispheres if the number of jets in the
          edm::LogError("RazorTriggerAnalyzerMuon") << "Cannot calculate M_R and R^2 because there are too many jets! (trigger passed automatically without forming the hemispheres)" << endl;
          return;
      }

      //This looks unusual because of how HLTRHemisphere deals with muons -- go look there for more information.  In the current razor triggers, muons are added to the hemispheres (no special treatment), and hemispheres->size() = 2.
      if(hemispheres->size() != 0 && hemispheres->size() != 2 && hemispheres->size() != 5 && hemispheres->size() != 10){
          edm::LogError("RazorTriggerAnalyzerMuon") << "Invalid hemisphere collection!  hemispheres->size() = " << hemispheres->size() << endl;
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

void RazorTriggerAnalyzerMuon::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup)
{
  edm::LogInfo("RazorTriggerAnalyzerMuon") << "RazorTriggerAnalyzerMuon::endLuminosityBlock" << std::endl;
}


void RazorTriggerAnalyzerMuon::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("RazorTriggerAnalyzerMuon") << "RazorTriggerAnalyzerMuon::endRun" << std::endl;
}

//CalcMR and CalcR borrowed from HLTRFilter.cc
double 
RazorTriggerAnalyzerMuon::CalcMR(TLorentzVector ja, TLorentzVector jb){
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
  RazorTriggerAnalyzerMuon::CalcR(double MR, TLorentzVector ja, TLorentzVector jb, edm::Handle<edm::View<reco::MET> > inputMet, const std::vector<math::XYZTLorentzVector>& muons){
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
DEFINE_FWK_MODULE(RazorTriggerAnalyzerMuon);
