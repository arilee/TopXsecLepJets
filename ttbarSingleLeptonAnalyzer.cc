// -*- C++ -*-
//
// Package:    ttbarSingleLeptonAnalyzer
// Class:      ttbarSingleLeptonAnalyzer
// 
/**\class ttbarSingleLeptonAnalyzer ttbarSingleLeptonAnalyzer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Javier Brochero Cifuentes,512 1-001,+41227670488,
//         Created:  Tue Feb  3 09:52:55 CET 2015
// $Id$
//
//

// system include files
#include <memory>
#include <math.h> 
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm> // max

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
//
// class declaration
//

class ttbarSingleLeptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit ttbarSingleLeptonAnalyzer(const edm::ParameterSet&);
  ~ttbarSingleLeptonAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  bool IsTightMuon(cat::MuonCollection::const_iterator &i_muon_candidate);
  bool IsTightElectron(cat::ElectronCollection::const_iterator &i_electron_candidate);
  

  // ----------member data ---------------------------

  TTree *AnalysisTree = new TTree();

  unsigned int minTracks_;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int b_Event, b_Run, b_Lumi_Number;

  // PU/Vertices
  float b_PUWeight; 
  int b_nPV, b_nGoodPV;

  int b_Channel;

  // MET
  float b_MET, b_MET_phi;

  // Leptons
  float b_Lepton_px;
  float b_Lepton_py;
  float b_Lepton_pz;
  float b_Lepton_E;

  // Jets
  std::vector<float> *b_Jet_px;
  std::vector<float> *b_Jet_py;
  std::vector<float> *b_Jet_pz;
  std::vector<float> *b_Jet_E;
  // ID
  std::vector<bool> *b_Jet_LooseID;
  // Flavour
  std::vector<int> *b_Jet_partonFlavour;
  std::vector<int> *b_Jet_hadronFlavour;
  // Smearing and Shifts  
  std::vector<float> *b_Jet_smearedRes;
  std::vector<float> *b_Jet_smearedResDown;
  std::vector<float> *b_Jet_smearedResUp;
  std::vector<float> *b_Jet_shiftedEnUp;
  std::vector<float> *b_Jet_shiftedEnDown;
  // b-Jet discriminant
  std::vector<float> *b_Jet_CSV;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ttbarSingleLeptonAnalyzer::ttbarSingleLeptonAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
     

  edm::Service<TFileService> fs;
  AnalysisTree = fs->make<TTree>("AnalysisTree", "TopTree");
  
  AnalysisTree->Branch("event",      &b_Event,       "Event/I");
  AnalysisTree->Branch("run",        &b_Run,         "Run/I");
  AnalysisTree->Branch("luminumber", &b_Lumi_Number, "Lumi_Number/I");

  AnalysisTree->Branch("PUWeight", &b_PUWeight, "PUWeight/F");
  AnalysisTree->Branch("PV",       &b_nPV,      "nPV/I");
  AnalysisTree->Branch("GoodPV",   &b_nGoodPV,  "nGoodPV/I");

  AnalysisTree->Branch("channel",  &b_Channel,  "Channel/I");

  AnalysisTree->Branch("MET",     &b_MET,     "MET/F");
  AnalysisTree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  AnalysisTree->Branch("lepton_px", &b_Lepton_px, "lepton_px/F");
  AnalysisTree->Branch("lepton_py", &b_Lepton_py, "lepton_py/F");
  AnalysisTree->Branch("lepton_pz", &b_Lepton_pz, "lepton_pz/F");
  AnalysisTree->Branch("lepton_E" , &b_Lepton_E,  "lepton_E/F" );

  AnalysisTree->Branch("jet_px", "std::vector<float>", &b_Jet_px);
  AnalysisTree->Branch("jet_py", "std::vector<float>", &b_Jet_py);
  AnalysisTree->Branch("jet_pz", "std::vector<float>", &b_Jet_pz);
  AnalysisTree->Branch("jet_E" , "std::vector<float>", &b_Jet_E );

  AnalysisTree->Branch("jet_LooseID", "std::vector<bool>", &b_Jet_LooseID);
  
  AnalysisTree->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  AnalysisTree->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);
  
  AnalysisTree->Branch("jet_smearedRes",     "std::vector<float>", &b_Jet_smearedRes);
  AnalysisTree->Branch("jet_smearedResDown", "std::vector<float>", &b_Jet_smearedResDown);
  AnalysisTree->Branch("jet_smearedResUp",   "std::vector<float>", &b_Jet_smearedResUp); 
  AnalysisTree->Branch("jet_shiftedEnUp",    "std::vector<float>", &b_Jet_shiftedEnUp);  
  AnalysisTree->Branch("jet_shiftedEnDown",  "std::vector<float>", &b_Jet_shiftedEnDown);

  AnalysisTree->Branch("jet_CSV" , "std::vector<float>", &b_Jet_CSV );

}


ttbarSingleLeptonAnalyzer::~ttbarSingleLeptonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbarSingleLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  b_Jet_px = new std::vector<float>;
  b_Jet_py = new std::vector<float>;
  b_Jet_pz = new std::vector<float>;
  b_Jet_E  = new std::vector<float>;

  b_Jet_LooseID = new std::vector<bool>;
  
  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;
  
  b_Jet_smearedRes     = new std::vector<float>;
  b_Jet_smearedResDown = new std::vector<float>;
  b_Jet_smearedResUp   = new std::vector<float>;
  b_Jet_shiftedEnUp    = new std::vector<float>;
  b_Jet_shiftedEnDown  = new std::vector<float>;
  
  b_Jet_CSV  = new std::vector<float>;
  

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // PU Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  edm::Handle<double> PUWeight;
  edm::Handle<int>    nTrueInter;
  edm::Handle<int>    reco_nPV;

  iEvent.getByLabel("pileupWeight", "nTrueInteraction", nTrueInter);

  iEvent.getByLabel("pileupWeight", PUWeight);
  b_PUWeight = *PUWeight;

  iEvent.getByLabel("recoEventInfo", "pvN", reco_nPV); // Same as vtxs.size()!!

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<reco::VertexCollection> pvertex;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", pvertex);
  const reco::VertexCollection& vtxs = *(pvertex.product());
  
  int n_vtxs = 0;
  // Loop over vertices
  if (vtxs.size() != 0) {
    for (size_t i=0; i<vtxs.size(); i++) {
      
     if ( fabs(vtxs[i].z())        < 24 &&
	  vtxs[i].position().Rho() < 2  &&
	  vtxs[i].ndof()           > 4  &&
	  !(vtxs[i].isFake())              ) {
       n_vtxs ++;
     }
    }// for(vertex)
  } // if(vertex)
  
  b_nPV = vtxs.size();
  b_nGoodPV = n_vtxs;  

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Secondary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<cat::SecVertexCollection> svertex;
  iEvent.getByLabel("catSecVertex", svertex);
 
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Missing E_T
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  Handle<cat::METCollection> MET;
  iEvent.getByLabel("catMETs", MET);

  // MET-PF
  b_MET     = MET->begin()->pt();
  b_MET_phi = MET->begin()->phi();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  std::vector<cat::ElectronCollection::const_iterator> elec;

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByLabel("catElectrons", electrons); 
 
  for(cat::ElectronCollection::const_iterator i_electron = electrons->begin();
      i_electron != electrons->end();
      ++ i_electron){

    if(IsTightElectron(i_electron)) elec.push_back(i_electron);

  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  std::vector<cat::MuonCollection::const_iterator> muon;

  Handle<cat::MuonCollection> muons;
  iEvent.getByLabel("catMuons", muons); 
 
  for(cat::MuonCollection::const_iterator i_muon = muons->begin();
      i_muon != muons->end();
      ++ i_muon){

    if( IsTightMuon(i_muon) ) muon.push_back(i_muon);

  }

  //---------------------------------------------------------------------------
  //----------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;
  int ch_tag  =999;

  if(muon.size() == 1 && elec.size() == 0){
    lepton.SetPxPyPzE(muon[0]->px(), muon[0]->py(), muon[0]->pz(), muon[0]->energy());
    ch_tag = 0; //muon + jets
  }

  if(muon.size() == 0 && elec.size() == 1){
    lepton.SetPxPyPzE(elec[0]->px(), elec[0]->py(), elec[0]->pz(), elec[0]->energy());
    ch_tag = 1; //electron + jets
  }


  //---------------------------------------------------------------------------
  //----------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  if (ch_tag<2){ // Single lepton event 

    b_Channel  = ch_tag;
    
    b_Lepton_px = lepton.Px();
    b_Lepton_py = lepton.Py();
    b_Lepton_pz = lepton.Pz();
    b_Lepton_E  = lepton.E();
    
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Jets
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    Handle<cat::JetCollection> jet;
    iEvent.getByLabel("catJets", jet);  
    
    for(cat::JetCollection::const_iterator i_jet = jet->begin();
	i_jet != jet->end();
	++ i_jet){
            
      if(fabs(i_jet->eta()) < 2.4 &&
	 i_jet->pt()        > 20 ){
	
	// Basic variables
	b_Jet_px->push_back(i_jet->px());
	b_Jet_py->push_back(i_jet->py());
	b_Jet_pz->push_back(i_jet->pz());
	b_Jet_E ->push_back(i_jet->energy());

	// Jet ID (Loose)
	b_Jet_LooseID ->push_back(i_jet->LooseId());

	// Parton Flavour
	b_Jet_partonFlavour->push_back(i_jet->partonFlavour()); 
	b_Jet_hadronFlavour->push_back(i_jet->hadronFlavour());
	
	// Smeared and Shifted
	b_Jet_smearedRes     ->push_back( i_jet->smearedRes() ); 
	b_Jet_smearedResDown ->push_back(i_jet->smearedResDown());
	b_Jet_smearedResUp   ->push_back(i_jet->smearedResUp());
	b_Jet_shiftedEnUp    ->push_back(i_jet->shiftedEnUp());
	b_Jet_shiftedEnDown  ->push_back(i_jet->shiftedEnDown());

	// b-tag discriminant
	float jet_btagDis_CSV = i_jet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
	b_Jet_CSV ->push_back(jet_btagDis_CSV);
	
      }
    }
    
    
    AnalysisTree->Fill();
    
  } // if(ch_tag)

  delete b_Jet_px;
  delete b_Jet_py;
  delete b_Jet_pz;
  delete b_Jet_E;

  delete b_Jet_LooseID;
  
  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;
  
  delete b_Jet_smearedRes;
  delete b_Jet_smearedResDown;
  delete b_Jet_smearedResUp;
  delete b_Jet_shiftedEnUp;
  delete b_Jet_shiftedEnDown;

  delete b_Jet_CSV;
 

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//------------- Good Muon Selection -----------------------
bool ttbarSingleLeptonAnalyzer::IsTightMuon(cat::MuonCollection::const_iterator &i_muon_candidate)
{
  bool GoodMuon=true;
  
  // Tight cut already defined into CAT::Muon
  //GoodMuon &= (i_muon_candidate->isTightMuon());
  
  GoodMuon &= (i_muon_candidate->isPFMuon());       // PF
  GoodMuon &= (i_muon_candidate->pt()> 20);         // pT
  GoodMuon &= (fabs(i_muon_candidate->eta())< 2.4); // eta

  GoodMuon &=(i_muon_candidate->isGlobalMuon());
  GoodMuon &=(i_muon_candidate->isPFMuon());  
  GoodMuon &=(i_muon_candidate->normalizedChi2() < 10);  
  GoodMuon &=(i_muon_candidate->numberOfValidMuonHits() > 0);  
  GoodMuon &=(i_muon_candidate->numberOfMatchedStations() > 1);  
  GoodMuon &=(fabs(i_muon_candidate->dxy()) < 0.2); //mm
  GoodMuon &=(fabs(i_muon_candidate->dz()) < 0.5); //mm
  GoodMuon &=(i_muon_candidate->numberOfValidPixelHits() > 0);
  GoodMuon &=(i_muon_candidate->trackerLayersWithMeasurement() > 5);

  float PFIsoMuon=999.;
  PFIsoMuon = i_muon_candidate->chargedHadronIso(0.3) +
              std::max(0.0, i_muon_candidate->neutralHadronIso(0.3) + 
	                    i_muon_candidate->photonIso(0.3) - 
	                    0.5*i_muon_candidate->puChargedHadronIso(0.3));
    
  PFIsoMuon = PFIsoMuon/i_muon_candidate->pt();
  
  GoodMuon &=( PFIsoMuon<0.12 );

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool ttbarSingleLeptonAnalyzer::IsTightElectron(cat::ElectronCollection::const_iterator &i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate->isPF() );            // PF
  GoodElectron &= (i_electron_candidate->pt() > 20);          // pT
  GoodElectron &= (fabs(i_electron_candidate->eta()) < 2.4);  // eta
  GoodElectron &= (fabs(i_electron_candidate->eta()) < 1.4442 || 
		   fabs(i_electron_candidate->eta()) > 1.566);

  GoodElectron &= i_electron_candidate->passConversionVeto();

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate->electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium") > 0.0;

//----------------------------------------------------------------------------------------------------
//------------- The Relative Isolation is already calculated in the CAT object -----------------------
//----------------------------------------------------------------------------------------------------
  // Effective Area Parametrization
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  // Double_t AEff03 = 0.;

  // if      (fabs(i_electron_candidate->eta()) < 0.8)                                                 AEff03 = 0.1013;
  // else if (fabs(i_electron_candidate->eta()) >= 0.8   && fabs(i_electron_candidate->eta()) < 1.3)   AEff03 = 0.0988; 
  // else if (fabs(i_electron_candidate->eta()) >= 1.3   && fabs(i_electron_candidate->eta()) < 2.0)   AEff03 = 0.0572; 
  // else if (fabs(i_electron_candidate->eta()) >= 2.0   && fabs(i_electron_candidate->eta()) < 2.2)   AEff03 = 0.0842; 
  // else if (fabs(i_electron_candidate->eta()) >= 2.2)                                                AEff03 = 0.1530; 
  
  // float PFIsoElectron = ( i_electron_candidate->chargedHadronIso( 0.3 ) +
  // 			  std::max(0.0, 
  // 				   i_electron_candidate->neutralHadronIso( 0.3 ) +
  // 				   i_electron_candidate->photonIso( 0.3 ) -  
  // 				   AEff03*1
  // 				   )
  // 			  );
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------


  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;
  GoodElectron &=( i_electron_candidate->relIso( 0.3 ) < 0.12 );

  return GoodElectron;

}
// ------------ method called once each job just before starting event loop  ------------
void 
ttbarSingleLeptonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ttbarSingleLeptonAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ttbarSingleLeptonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ttbarSingleLeptonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ttbarSingleLeptonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ttbarSingleLeptonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttbarSingleLeptonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbarSingleLeptonAnalyzer);
