// -*- C++ -*-
//
// Package:    DisappearingTrackProducer
// Class:      DisappearingTrackProducer
// 
/**\class DisappearingTrackProducer DisappearingTrackProducer.cc RA2Classic/DisappearingTrackProducer/src/DisappearingTrackProducer.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Fri Apr 11 16:35:33 CEST 2014
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//
// class declaration
//

class DisappearingTrackProducer : public edm::EDProducer {
public:
	explicit DisappearingTrackProducer(const edm::ParameterSet&);
	~DisappearingTrackProducer();
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken;
  edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronsToken;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken;
  edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken;

  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEBToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEEToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsESToken;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> reducedHcalRecHitsHBToken;
  edm::EDGetTokenT<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >> reducedHcalRecHitsHFToken;
  edm::EDGetTokenT<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >> reducedHcalRecHitsHOToken;
  edm::EDGetTokenT<std::vector<reco::CaloJet>> caloJetsToken;

  double minTrackPt;
  double maxTrackEta;
  double conePtSumDR;
  double conePtSumMaxPtPercentage;
  double minTrackJetDR;
  double minTrackLeptonDR;
  double RequireNumberOfValidPixelHits;
  double RequireNumberOfValidTrackerHits;
  double track_dxy;
  double track_dz;

  double minMissingOuterHits;
  double caloEnergyDepositionMaxDR;
  double caloEnergyDepositionMaxE;
  double deadNoisyDR;     
  bool useCaloJetsInsteadOfHits;

  bool useKinematics = false;
  bool useTrackTrackerIso = false;
  bool useTrackerActivity = false;
  bool useTrackJetIso = false;
  bool useNoFakes = false;
  bool useTrackLeptonIso = false;
  bool useGapsVeto = false;
  bool useCaloEnergy = false;
  bool useMissingOuterHits = false;
  bool useNoisyCaloVeto = false;
      
  std::string doStage;
	
	// ----------member data ---------------------------
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
DisappearingTrackProducer::DisappearingTrackProducer(const edm::ParameterSet& iConfig)
{
	//register your produc	
  produces<std::vector<reco::Track> >("DisappearingTracks");
  //produces<bool>("");
	/* Examples
	 *   produces<ExampleData2>();
	 * 
	 *   //if do put with a label
	 *   produces<ExampleData2>("label");
	 * 
	 *   //if you want to put into the Run
	 *   produces<ExampleData2,InRun>();
	 */
	//now do what ever other initialization is needed
  tracksToken = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("selectedTracks"));
  electronsToken = consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("selectedElectrons"));
  muonsToken = consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("selectedMuons"));
  jetsToken = consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("selectedPFJets"));
  tracksToken = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("selectedTracks"));
  reducedEcalRecHitsEBToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsEB"));
  reducedEcalRecHitsEEToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsEE"));
  reducedEcalRecHitsESToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsES"));
  reducedHcalRecHitsHBToken = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(iConfig.getParameter<edm::InputTag>("selectedHcalRecHits"));
  reducedHcalRecHitsHFToken = consumes<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedHcalRecHits"));
  reducedHcalRecHitsHOToken = consumes<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >>(iConfig.getParameter<edm::InputTag>("selectedHcalRecHits"));
  caloJetsToken = consumes<std::vector<reco::CaloJet>>(iConfig.getParameter<edm::InputTag>("selectedCaloJets"));   

  minTrackPt = iConfig.getParameter<double>("minTrackPt");
  maxTrackEta = iConfig.getParameter<double>("maxTrackEta");
  conePtSumDR = iConfig.getParameter<double>("conePtSumDR");
  conePtSumMaxPtPercentage = iConfig.getParameter<double>("conePtSumMaxPtPercentage");
  minTrackJetDR = iConfig.getParameter<double>("minTrackJetDR");
  minTrackLeptonDR = iConfig.getParameter<double>("minTrackLeptonDR");
  RequireNumberOfValidPixelHits = iConfig.getParameter<double>("RequireNumberOfValidPixelHits");
  RequireNumberOfValidTrackerHits = iConfig.getParameter<double>("RequireNumberOfValidTrackerHits");
  track_dxy = iConfig.getParameter<double>("maxDxy");
  track_dz = iConfig.getParameter<double>("maxDz");

  minMissingOuterHits = iConfig.getParameter<double>("minMissingOuterHits");
  caloEnergyDepositionMaxDR = iConfig.getParameter<double>("caloEnergyDepositionMaxDR");
  caloEnergyDepositionMaxE = iConfig.getParameter<double>("caloEnergyDepositionMaxE");
  deadNoisyDR = iConfig.getParameter<double>("deadNoisyDR");

  useCaloJetsInsteadOfHits = iConfig.getParameter<bool>("useCaloJetsInsteadOfHits");
   
  doStage = iConfig.getParameter<std::string>("doStage");
  std::cout << "what about you guys!" << std::endl;
   
  if (doStage == "all") { useKinematics = true; useTrackTrackerIso = true; useTrackerActivity = true;
    useTrackJetIso = true; useNoFakes = true; useTrackLeptonIso = true; useGapsVeto = true;
    useCaloEnergy = true; useMissingOuterHits = true; useNoisyCaloVeto = true; }
  else if (doStage == "Kinematics") useKinematics = true;
  else if (doStage == "TrackerIso") useTrackTrackerIso = true;
  else if (doStage == "TrackerActivity") useTrackerActivity = true;
  else if (doStage == "TrackJetIso") useTrackJetIso = true;
  else if (doStage == "NoFakes") useNoFakes = true;
  else if (doStage == "TrackLeptonIso") useTrackLeptonIso = true;
  else if (doStage == "GapsVeto") useGapsVeto = true;
  else if (doStage == "CaloEnergy") useCaloEnergy = true;
  else if (doStage == "MissingOuterHits") useMissingOuterHits = true;
  else if (doStage == "NoisyCaloVeto") useNoisyCaloVeto = true;
  else std::cout << "Unknown stage label!" << std::endl;
  std::cout << "we ain't found shit!" <<  std::endl;
	
}


DisappearingTrackProducer::~DisappearingTrackProducer()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
	
}

// check if particle is isolated with respect to another particle collection:
template <typename T1, typename T2>
bool checkParticleIsolation(T1 particle, T2 pcollection, double RequiredMinDR){

  bool isolated = false;
  double minDR = 1e6;
  for( const auto& item : *pcollection){
    if (deltaR(particle, item) < minDR) minDR = deltaR(particle, item);
  }
  if (minDR > RequiredMinDR) isolated = true;
  return isolated;
}

// function for looping over one single recHit collection:
template <typename T>
double LoopOverRecHits(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double caloEnergyDepositionMaxDR){
   
  double energyDepositedPerCalorimeter = 0;      
  for( const auto& hit : *hitcollection){ 
          
    double hitEnergy = hit.energy();
    double cellEta = CaloGeomHandle.product()->getPosition(hit.detid()).eta();
    double cellPhi = CaloGeomHandle.product()->getPosition(hit.detid()).phi();
        
    if (deltaR(cellEta, cellPhi, track.eta(), track.phi()) < caloEnergyDepositionMaxDR) {
      energyDepositedPerCalorimeter += hitEnergy;
    }
  }
       
  return energyDepositedPerCalorimeter;
}

// check dead or noisy ECAL channels:
template <typename T>
bool checkNoDeadNoisyECALInTrackCone(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double DeadNoisyDR){
   
  double noDeadNoisyCellsInTrackCone = true;

  for( const auto& hit : *hitcollection){
    double cellEta = CaloGeomHandle.product()->getPosition(hit.detid()).eta();
    double cellPhi = CaloGeomHandle.product()->getPosition(hit.detid()).phi();

    std::vector<int> badFlags(2);
    badFlags[0] = EcalRecHit::kDead;
    badFlags[1] = EcalRecHit::kNoisy;

    if (hit.checkFlags(badFlags)) {
      if (deltaR(cellEta, cellPhi, track.eta(), track.phi()) < DeadNoisyDR) {
	noDeadNoisyCellsInTrackCone = false;
      }
    }
  }

  return noDeadNoisyCellsInTrackCone;
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
DisappearingTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  std::cout << "in DisappearingTrackProducer::produce" << std::endl;
  // register input collections:
  Handle<std::vector<reco::Track>> tracks; iEvent.getByToken( tracksToken, tracks);
  Handle<std::vector<reco::GsfElectron>> electrons; iEvent.getByToken( electronsToken, electrons);
  Handle<std::vector<reco::Muon>> muons; iEvent.getByToken( muonsToken, muons);
  Handle<std::vector<reco::PFJet>> jets; iEvent.getByToken( jetsToken, jets);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEB;
  iEvent.getByToken( reducedEcalRecHitsEBToken, reducedEcalRecHitsEB);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEE;
  iEvent.getByToken( reducedEcalRecHitsEEToken, reducedEcalRecHitsEE);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsES;
  iEvent.getByToken( reducedEcalRecHitsESToken, reducedEcalRecHitsES);
  Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> reducedHcalRecHitsHB;
  iEvent.getByToken( reducedHcalRecHitsHBToken, reducedHcalRecHitsHB);
  Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >> reducedHcalRecHitsHF;
  iEvent.getByToken( reducedHcalRecHitsHFToken, reducedHcalRecHitsHF);
  Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >> reducedHcalRecHitsHO;
  iEvent.getByToken( reducedHcalRecHitsHOToken, reducedHcalRecHitsHO);
  Handle<std::vector<reco::CaloJet>> caloJets;
  iEvent.getByToken( caloJetsToken, caloJets);
  edm::ESHandle<CaloGeometry> CaloGeomHandle;
  iSetup.get<CaloGeometryRecord>().get(CaloGeomHandle);


  // register output track collection:
  std::unique_ptr<std::vector<reco::Track> > candidateTracks(new std::vector<reco::Track>);

  for( const auto& track : *tracks){

    // check track kinematics:
    bool passedTrackKinematics = false;
    if (!useKinematics) passedTrackKinematics = true;
    else if ((track.pt() > minTrackPt) && (std::abs(track.eta())) < maxTrackEta) passedTrackKinematics = true;

    // check pt sum of tracker activity within a cone around current track:
    bool passedTrackTrackerIso = false;
    if (!useTrackTrackerIso) passedTrackTrackerIso = true;
    else {
      double conePtSum = 0;
      for( const auto& othertrack : *tracks){
	if (deltaR(track, othertrack) == 0) continue;
	if (deltaR(track, othertrack) < conePtSumDR) conePtSum += track.pt();
      }
      if (conePtSum / track.pt() <= conePtSumMaxPtPercentage) passedTrackTrackerIso = true;
    }

    // check track/jet isolation:
    bool passedTrackJetIso = false;
    if (!useTrackJetIso) passedTrackJetIso = true;
    else passedTrackJetIso = checkParticleIsolation(track, jets, minTrackJetDR);

    // check if fake track:
    bool passedNoFakes = false;
    if (!useNoFakes) passedNoFakes = true;
    else {
      reco::HitPattern hitpattern = track.hitPattern();
      if ( (hitpattern.numberOfValidPixelHits() >= RequireNumberOfValidPixelHits) && (hitpattern.numberOfValidTrackerHits() >= RequireNumberOfValidTrackerHits) ) {
	// no missing inner or middle hits, require consecutive pattern of hits:
	if ( hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS) == 0 ) {
	  // check track displacement:
	  if ( (std::abs(track.dxy()) < track_dxy) && (std::abs(track.dz()) < track_dz) ) {
	    passedNoFakes = true;
	  }
	}
      }
    }

    // check track / lepton isolation:
    bool passedTrackLeptonIso = false;
    if (!useTrackLeptonIso) passedTrackLeptonIso = true;
    else passedTrackLeptonIso = checkParticleIsolation(track, electrons, minTrackLeptonDR) &
	   checkParticleIsolation(track, muons, minTrackLeptonDR);

    // check gaps veto:
    bool passedGapsVeto = true;
    if ( useGapsVeto && (std::abs(track.eta()) > 0.15) && (std::abs(track.eta()) < 0.35)) passedGapsVeto = false;
    if ( useGapsVeto && (std::abs(track.eta()) > 1.55) && (std::abs(track.eta()) < 1.85)) passedGapsVeto = false;
    if ( useGapsVeto && (std::abs(track.eta()) > 1.42) && (std::abs(track.eta()) < 1.65)) passedGapsVeto = false;
      
        
    // check calo energy deposition:
    double energyDeposited = 0;
       
    if (!useCaloJetsInsteadOfHits) {        
      // loop over ecal and hcal recHit collections:
      if (reducedEcalRecHitsEB.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsEB, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
      if (reducedEcalRecHitsEE.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsEE, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
      if (reducedEcalRecHitsES.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsES, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
      if (reducedHcalRecHitsHB.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHB, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
      if (reducedHcalRecHitsHF.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHF, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
      if (reducedHcalRecHitsHO.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHO, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    } else {
      // optional cross-check, loop over calojets:
      for( const auto& jet : *caloJets){ 
	if (deltaR(jet.detectorP4().eta(), jet.detectorP4().phi(), track.eta(), track.phi()) < caloEnergyDepositionMaxDR) {
	  energyDeposited += jet.maxEInEmTowers() + jet.maxEInHadTowers();
	}
      }
    }

    bool passedCaloEnergy = false;
    if (!useCaloEnergy) passedCaloEnergy = true;
    else if (energyDeposited < caloEnergyDepositionMaxE) passedCaloEnergy = true;


    // check minimum of missing outer hits:
    bool passedMissingOuterHits = false;
    if (!useMissingOuterHits) passedMissingOuterHits = true;
    else {
      reco::HitPattern hitpattern = track.hitPattern();
      int missingOuterHits = hitpattern.stripLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS);
      if(missingOuterHits > minMissingOuterHits) passedMissingOuterHits = true;
    }


    // check dead or noisy calo cells:
    bool passedNoisyCaloVeto = false;
    if (!useNoisyCaloVeto) passedNoisyCaloVeto = true;
    else {
      passedNoisyCaloVeto = checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEB, track, CaloGeomHandle, deadNoisyDR) &
	checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEE, track, CaloGeomHandle, deadNoisyDR) &
	checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsES, track, CaloGeomHandle, deadNoisyDR);
    }

    // save new track collection:
    if (passedTrackKinematics && passedTrackTrackerIso && passedTrackJetIso &&
	passedNoFakes && passedTrackLeptonIso && passedGapsVeto &&
	passedCaloEnergy && passedMissingOuterHits && passedNoisyCaloVeto) {
      candidateTracks->push_back(track);
    }
  }
  std::cout << "Made it to here!" << std::endl;

  // save track collection in event:
  iEvent.put(std::move(candidateTracks),"DisappearingTracks");
  //std::unique_ptr<bool> newbool(new bool(true));
  //iEvent.put(std::move(newbool),"");
	
}

// ------------ method called once each job just before starting event loop  ------------
void 
DisappearingTrackProducer::beginJob()
{
  std::cout << "downright gleeful" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisappearingTrackProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DisappearingTrackProducer::beginRun(edm::Run&, edm::EventSetup const&)
{  std::cout << "running around the circles" << std::endl;
}

// ------------ method called when ending the processing of a run  ------------
void 
DisappearingTrackProducer::endRun(edm::Run&, edm::EventSetup const&)
{
  std::cout << "endrun" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DisappearingTrackProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
 std::cout << "lumi!" << std::endl;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DisappearingTrackProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisappearingTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisappearingTrackProducer);
