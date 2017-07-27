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
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
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
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


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
  edm::EDGetTokenT<reco::VertexCollection> PrimVtxTok_;

  double minTrackPt;
  double maxTrackEta;
  double coneRelIsoDR;
  double conePtSumMaxPtPercentage;
  double minTrackJetDR;
  double minTrackLeptonDR;
  double RequireNumberOfValidPixelHits;
  double RequireNumberOfValidTrackerHits;
  double maxDxy;
  double maxDz;
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
  bool doDeDx = false;
  std::string dEdxEstimator_;
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
  //register your product	
  produces<std::vector<TLorentzVector> > ("chiCands");
  produces<std::vector<double> >           ("chiCands@dxyVtx");
  produces<std::vector<double> >           ("chiCands@dzVtx");
  produces<std::vector<int> >          ("chiCands@nMissingOuterHits");
  produces<std::vector<int> >         ("chiCands@nMissingInnerHits");
  produces<std::vector<int> >        ("chiCands@nMissingMiddleHits");
  produces<std::vector<int> >         ("chiCands@nValidPixelHits");
  produces<std::vector<int> >          ("chiCands@nValidTrackerHits");
  produces<std::vector<double> >          ("chiCands@chi2perNdof");
  produces<std::vector<double> >       ("chiCands@trkRelIso");
  produces<std::vector<double> >        ("chiCands@trkMiniRelIso");
  produces<std::vector<double> >         ("chiCands@matchedCaloEnergy");
  produces<std::vector<double> >         ("chiCands@deDxHarmonic2");
  produces<std::vector<bool> >        ("chiCands@passExo16044Tag");
  produces<std::vector<bool> >        ("chiCands@passExo16044JetIso");
  produces<std::vector<bool> >        ("chiCands@passExo16044LepIso");

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
  edm::InputTag PrimVtxTag_;
  PrimVtxTag_=iConfig.getParameter<edm::InputTag>("PrimaryVertex");
  PrimVtxTok_ = consumes<reco::VertexCollection>(PrimVtxTag_);

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
  edm::InputTag estimatorTag("dedxHarmonic2");
  consumes<edm::ValueMap<reco::DeDxData>>(estimatorTag);

  minTrackPt = iConfig.getParameter<double>("minTrackPt");
  maxTrackEta = iConfig.getParameter<double>("maxTrackEta");
  coneRelIsoDR = iConfig.getParameter<double>("coneRelIsoDR");
  conePtSumMaxPtPercentage = iConfig.getParameter<double>("conePtSumMaxPtPercentage");
  minTrackJetDR = iConfig.getParameter<double>("minTrackJetDR");
  minTrackLeptonDR = iConfig.getParameter<double>("minTrackLeptonDR");
  RequireNumberOfValidPixelHits = iConfig.getParameter<double>("RequireNumberOfValidPixelHits");
  RequireNumberOfValidTrackerHits = iConfig.getParameter<double>("RequireNumberOfValidTrackerHits");
  maxDxy = iConfig.getParameter<double>("maxDxy");
  maxDz = iConfig.getParameter<double>("maxDz");

  minMissingOuterHits = iConfig.getParameter<double>("minMissingOuterHits");
  caloEnergyDepositionMaxDR = iConfig.getParameter<double>("caloEnergyDepositionMaxDR");
  caloEnergyDepositionMaxE = iConfig.getParameter<double>("caloEnergyDepositionMaxE");
  deadNoisyDR = iConfig.getParameter<double>("deadNoisyDR");

  useCaloJetsInsteadOfHits = iConfig.getParameter<bool>("useCaloJetsInsteadOfHits");
  doDeDx = iConfig.getParameter<bool>("doDeDx");
  dEdxEstimator_ = iConfig.getParameter<std::string>("dEdxEstimator");
  doStage = iConfig.getParameter<std::string>("doStage");
   
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

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxTrackHandle;
  iEvent.getByLabel(dEdxEstimator_, dEdxTrackHandle);
  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxTrackHandle.product();

  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(PrimVtxTok_, vtx_h);
  const reco::Vertex vtx = vtx_h->at(0);

  // register output track collection:
  std::unique_ptr<std::vector<TLorentzVector> > chiCands(new std::vector<TLorentzVector>);

  std::unique_ptr<std::vector<double> > chiCands_dxyVtx(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_dzVtx(new std::vector<double>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingOuterHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingInnerHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingMiddleHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nValidPixelHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nValidTrackerHits(new std::vector<int>);
  std::unique_ptr<std::vector<double> > chiCands_chi2perNdof(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trkRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trkMiniRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_matchedCaloEnergy(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_deDxHarmonic2(new std::vector<double>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044JetIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044LepIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044Tag(new std::vector<bool>);

  int itrack = -1;
  for( const auto& track : *tracks){
    itrack+=1;

    // check track kinematics:
    bool passedTrackKinematics = false;
    if (!useKinematics) passedTrackKinematics = true;
    else if ((track.pt() > minTrackPt) && (std::abs(track.eta())) < maxTrackEta) passedTrackKinematics = true;
    if (!(passedTrackKinematics)) continue;


    bool passExo16044Kinematics = track.pt() > 55 && std::abs(track.eta()) < 2.1;

    // check gaps veto:
    bool passedGapsVeto = true;
    if ( useGapsVeto && (std::abs(track.eta()) > 0.15) && (std::abs(track.eta()) < 0.35)) passedGapsVeto = false;
    if ( useGapsVeto && (std::abs(track.eta()) > 1.55) && (std::abs(track.eta()) < 1.85)) passedGapsVeto = false;
    if ( useGapsVeto && (std::abs(track.eta()) > 1.42) && (std::abs(track.eta()) < 1.65)) passedGapsVeto = false;
    if (!(passedGapsVeto)) continue;


    // check dead or noisy calo cells:
    bool passedNoisyCaloVeto = false;
    if (!useNoisyCaloVeto) passedNoisyCaloVeto = true;
    else {
      passedNoisyCaloVeto = checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEB, track, CaloGeomHandle, deadNoisyDR) &
	checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEE, track, CaloGeomHandle, deadNoisyDR) &
	checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsES, track, CaloGeomHandle, deadNoisyDR);
    }
    if (!(passedNoisyCaloVeto)) continue;


    // do vertex consistency

    bool vertex_match = std::abs(track.dxy(vtx.position())) < maxDxy && std::abs(track.dz(vtx.position())) < maxDz;
    if (!(vertex_match)) continue;


    //This completes the loose candidate selection, so now fill the products depending only on track object
    TLorentzVector tlv;
    tlv.SetPxPyPzE(track.px(),track.py(),track.pz(),track.pt());
    chiCands->push_back(tlv);
    chiCands_dxyVtx->push_back(std::abs(track.dxy(vtx.position())));
    chiCands_dzVtx->push_back(std::abs(track.dz(vtx.position())));
    reco::HitPattern hitpattern = track.hitPattern();
    chiCands_nMissingOuterHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS));
    chiCands_nMissingInnerHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS));
    chiCands_nMissingMiddleHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS));
    chiCands_nValidPixelHits->push_back(hitpattern.numberOfValidPixelHits());
    chiCands_nValidTrackerHits->push_back(hitpattern.numberOfValidTrackerHits());
    chiCands_chi2perNdof->push_back(1.0*track.chi2()/track.ndof());


    // check if fake track:
    bool passedNoFakes = false;
    if (!useNoFakes) passedNoFakes = true;
    else {
      if ( (hitpattern.numberOfValidPixelHits() >= RequireNumberOfValidPixelHits) && (hitpattern.numberOfValidTrackerHits() >= RequireNumberOfValidTrackerHits) ) {
	// no missing inner or middle hits, require consecutive pattern of hits:
	if ( hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS) == 0 &&
	     hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS) == 0) 
	  {
	    // check track displacement:
	    if ( vertex_match ) passedNoFakes = true;
	  }
      }
    }

      
        
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
	if (deltaR(jet.detectorP4().eta(), jet.detectorP4().phi(), track.eta(), track.phi()) < caloEnergyDepositionMaxDR) 
	  {
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
      int missingOuterHits = hitpattern.stripLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS);
      if(missingOuterHits > minMissingOuterHits) passedMissingOuterHits = true;
    }

    // check pt sum of tracker activity within a cone around current track:
    bool passedTrackTrackerIso = false;
    double relIso = -1;
    double miniIso = -1;
    if (!useTrackTrackerIso) passedTrackTrackerIso = true;
    else {
      double conePtSum_rel  = 0;
      double conePtSum_mini = 0;
      double miniConeDr = 0.2;
      if (track.pt()>50) miniConeDr = 10.0/track.pt();
      else if (track.pt()>200) miniConeDr = 0.05;
      for( const auto& othertrack : *tracks){
	if (deltaR(track, othertrack) == 0) continue;
	if (deltaR(track, othertrack) < coneRelIsoDR) conePtSum_rel += othertrack.pt();
	if (deltaR(track, othertrack) < miniConeDr) conePtSum_mini += othertrack.pt();
      }
      if (conePtSum_rel / track.pt() <= conePtSumMaxPtPercentage) passedTrackTrackerIso = true;
      relIso = conePtSum_rel / track.pt();
      miniIso = conePtSum_mini / track.pt();
    }
    

    // check track/jet isolation:
    bool passedTrackJetIso = false;
    if (!useTrackJetIso) passedTrackJetIso = true;
    else passedTrackJetIso = checkParticleIsolation(track, jets, minTrackJetDR);
    passedTrackJetIso = true;//hack because this is killing everything. should include below as well.

    // check track / lepton isolation:
    bool passedTrackLeptonIso = false;
    if (!useTrackLeptonIso) passedTrackLeptonIso = true;
    else passedTrackLeptonIso = checkParticleIsolation(track, electrons, minTrackLeptonDR) &
	   checkParticleIsolation(track, muons, minTrackLeptonDR);
    
    double dedx = -1;
    if (doDeDx) dedx = dEdxTrack[reco::TrackRef(tracks, itrack)].dEdx();
	
    chiCands_deDxHarmonic2->push_back(dedx);
    chiCands_trkRelIso->push_back(relIso);
    chiCands_trkMiniRelIso->push_back(miniIso);
    chiCands_matchedCaloEnergy->push_back(energyDeposited);
    chiCands_passExo16044JetIso->push_back(passedTrackJetIso);
    chiCands_passExo16044LepIso->push_back(passedTrackLeptonIso);


    // save new track collection:
    if (passExo16044Kinematics && passedTrackTrackerIso &&
	passedNoFakes && passedTrackLeptonIso && passedGapsVeto &&
	passedCaloEnergy && passedMissingOuterHits && passedNoisyCaloVeto) {//missing && passedTrackJetIso bug?
      chiCands_passExo16044Tag->push_back(true);
    }
    else chiCands_passExo16044Tag->push_back(false);
  }

  // save track collection in event:
  iEvent.put(std::move(chiCands),"chiCands");
  iEvent.put(std::move(chiCands_dxyVtx),"chiCands@dxyVtx");
  iEvent.put(std::move(chiCands_dzVtx),"chiCands@dzVtx");
  iEvent.put(std::move(chiCands_nMissingOuterHits),"chiCands@nMissingOuterHits");
  iEvent.put(std::move(chiCands_nMissingInnerHits),"chiCands@nMissingInnerHits");
  iEvent.put(std::move(chiCands_nMissingMiddleHits),"chiCands@nMissingMiddleHits");
  iEvent.put(std::move(chiCands_nValidPixelHits),"chiCands@nValidPixelHits");
  iEvent.put(std::move(chiCands_nValidTrackerHits),"chiCands@nValidTrackerHits");
  iEvent.put(std::move(chiCands_chi2perNdof),"chiCands@chi2perNdof");
  iEvent.put(std::move(chiCands_trkRelIso),"chiCands@trkRelIso");
  iEvent.put(std::move(chiCands_trkMiniRelIso),"chiCands@trkMiniRelIso");
  iEvent.put(std::move(chiCands_matchedCaloEnergy),"chiCands@matchedCaloEnergy");
  iEvent.put(std::move(chiCands_deDxHarmonic2),"chiCands@deDxHarmonic2");
  iEvent.put(std::move(chiCands_passExo16044JetIso),"chiCands@passExo16044JetIso");
  iEvent.put(std::move(chiCands_passExo16044LepIso),"chiCands@passExo16044LepIso");
  iEvent.put(std::move(chiCands_passExo16044Tag), "chiCands@passExo16044Tag");
	
} 

// ------------ method called once each job just before starting event loop  ------------
void 
DisappearingTrackProducer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisappearingTrackProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
DisappearingTrackProducer::beginRun(edm::Run&, edm::EventSetup const&)
{}

// ------------ method called when ending the processing of a run  ------------
void 
DisappearingTrackProducer::endRun(edm::Run&, edm::EventSetup const&)
{}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DisappearingTrackProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{}

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
