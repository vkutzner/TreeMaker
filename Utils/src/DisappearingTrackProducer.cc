//
// authors:  Viktor Kutzner, Sam Bein
//
// Disappearing track candidate track producer, adapted for TreeMaker
//
//         Created:  Wed, 17 May 2017 13:05:08 GMT
//


#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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
  edm::EDGetTokenT<reco::VertexCollection> PrimVtxToken;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsToken;

  double minTrackPt;
  double maxTrackEta;
  double coneRelIsoDR;
  double conePtSumMaxPtPercentage;
  double minTrackJetDR;
  double minTrackLeptonDR;
  double maxChargedPFCandSumDR;
  double maxNeutralPFCandSumDR;
  int RequireNumberOfValidPixelHits;
  int RequireNumberOfValidTrackerHits;
  double maxDxy;
  double maxDz;
  int minMissingOuterHits;
  int pixelLayersWithMeasurement;
  int trackerLayersWithMeasurement;
  bool trackQualityUndef;
  bool trackQualityLoose;
  bool trackQualityTight;
  bool trackQualityHighPurity;
  bool trackQualityConfirmed;
  bool trackQualityGoodIterative;
  bool trackQualityLooseSetWithPV;
  bool trackQualityHighPuritySetWithPV;
  bool trackQualityDiscarded;
  bool trackQualitySize;
  double caloEnergyDepositionMaxDR;
  double caloEnergyDepositionMaxE;
  double deadNoisyDR;     
  bool useCaloJetsInsteadOfHits;
  bool doDeDx = false;
  std::string dEdxEstimator;

  bool passExo16044GapsVeto;
  bool passExo16044DeadNoisyECALVeto;
  bool passPFCandVeto;
   
};


DisappearingTrackProducer::DisappearingTrackProducer(const edm::ParameterSet& iConfig)
{
  //register your product    
  produces<std::vector<TLorentzVector> > ("chiCands");
  produces<std::vector<double> >         ("chiCands@dxyVtx");
  produces<std::vector<double> >         ("chiCands@dzVtx");
  produces<std::vector<int> >            ("chiCands@nMissingOuterHits");
  produces<std::vector<int> >            ("chiCands@pixelLayersWithMeasurement");
  produces<std::vector<int> >            ("chiCands@trackerLayersWithMeasurement");
  produces<std::vector<int> >            ("chiCands@nMissingInnerHits");
  produces<std::vector<int> >            ("chiCands@nMissingMiddleHits");
  produces<std::vector<int> >            ("chiCands@nValidPixelHits");
  produces<std::vector<int> >            ("chiCands@nValidTrackerHits");
  produces<std::vector<double> >         ("chiCands@ptError");
  produces<std::vector<double> >         ("chiCands@chi2perNdof");
  produces<std::vector<double> >         ("chiCands@trkRelIso");
  produces<std::vector<double> >         ("chiCands@trkMiniRelIso");
  produces<std::vector<double> >         ("chiCands@trackLeptonIso");
  produces<std::vector<double> >         ("chiCands@trackJetIso");
  produces<std::vector<double> >         ("chiCands@matchedCaloEnergy");
  produces<std::vector<double> >         ("chiCands@matchedCaloEnergyJets");
  produces<std::vector<double> >         ("chiCands@deDxHarmonic2");
  produces<std::vector<double> >         ("chiCands@chargedPtSum");
  produces<std::vector<double> >         ("chiCands@neutralPtSum");
  produces<std::vector<double> >         ("chiCands@neutralWithoutGammaPtSum");
  produces<std::vector<bool> >           ("chiCands@passExo16044GapsVeto");
  produces<std::vector<bool> >           ("chiCands@passExo16044DeadNoisyECALVeto");
  produces<std::vector<bool> >           ("chiCands@passExo16044Tag");
  produces<std::vector<bool> >           ("chiCands@passExo16044JetIso");
  produces<std::vector<bool> >           ("chiCands@passExo16044LepIso");
  produces<std::vector<bool> >           ("chiCands@passPFCandVeto");
  produces<std::vector<bool> >           ("chiCands@trackQualityUndef");
  produces<std::vector<bool> >           ("chiCands@trackQualityLoose");
  produces<std::vector<bool> >           ("chiCands@trackQualityTight");
  produces<std::vector<bool> >           ("chiCands@trackQualityHighPurity");
  produces<std::vector<bool> >           ("chiCands@trackQualityConfirmed");
  produces<std::vector<bool> >           ("chiCands@trackQualityGoodIterative");
  produces<std::vector<bool> >           ("chiCands@trackQualityLooseSetWithPV");
  produces<std::vector<bool> >           ("chiCands@trackQualityHighPuritySetWithPV");
  produces<std::vector<bool> >           ("chiCands@trackQualityDiscarded");
  produces<std::vector<bool> >           ("chiCands@trackQualitySize");

  PrimVtxToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"));

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
  pfCandsToken = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("selectedPFCand"));

  edm::InputTag estimatorTag("dedxHarmonic2");
  consumes<edm::ValueMap<reco::DeDxData>>(estimatorTag);

  minTrackPt = iConfig.getParameter<double>("minTrackPt");
  maxTrackEta = iConfig.getParameter<double>("maxTrackEta");
  coneRelIsoDR = iConfig.getParameter<double>("coneRelIsoDR");
  conePtSumMaxPtPercentage = iConfig.getParameter<double>("conePtSumMaxPtPercentage");
  minTrackJetDR = iConfig.getParameter<double>("minTrackJetDR");
  minTrackLeptonDR = iConfig.getParameter<double>("minTrackLeptonDR");
  RequireNumberOfValidPixelHits = iConfig.getParameter<int>("RequireNumberOfValidPixelHits");
  RequireNumberOfValidTrackerHits = iConfig.getParameter<int>("RequireNumberOfValidTrackerHits");
  maxDxy = iConfig.getParameter<double>("maxDxy");
  maxDz = iConfig.getParameter<double>("maxDz");

  maxChargedPFCandSumDR = iConfig.getParameter<double>("maxChargedPFCandSumDR");
  maxNeutralPFCandSumDR = iConfig.getParameter<double>("maxNeutralPFCandSumDR");

  minMissingOuterHits = iConfig.getParameter<int>("minMissingOuterHits");
  caloEnergyDepositionMaxDR = iConfig.getParameter<double>("caloEnergyDepositionMaxDR");
  caloEnergyDepositionMaxE = iConfig.getParameter<double>("caloEnergyDepositionMaxE");
  deadNoisyDR = iConfig.getParameter<double>("deadNoisyDR");

  useCaloJetsInsteadOfHits = iConfig.getParameter<bool>("useCaloJetsInsteadOfHits");
  doDeDx = iConfig.getParameter<bool>("doDeDx");
  dEdxEstimator = iConfig.getParameter<std::string>("dEdxEstimator");
       
}


DisappearingTrackProducer::~DisappearingTrackProducer()
{}


template <typename T1, typename T2>
double getParticleIsolation(T1 particle, T2 pcollection){

  // get DR of a selected particle with respect to all other particles in a given collection

  double minDR = 1e6;
  for( const auto& item : *pcollection){
    if ((deltaR(particle, item) < minDR) && (deltaR(particle, item) != 0)) minDR = deltaR(particle, item);
  }
  return minDR;
}


template <typename T1, typename T2>
bool checkParticleIsolation(T1 particle, T2 pcollection, double RequiredMinDR){

  // checks particle isolation with respect to another particle collection,
  // returns true if two particles are not inside a cone with minimal DR.

  bool isolated = false;
  double minDR = 1e6;
  for( const auto& item : *pcollection){
    if ((deltaR(particle, item) < minDR) && (deltaR(particle, item) != 0)) minDR = deltaR(particle, item);
  }
  if (minDR > RequiredMinDR) isolated = true;

  return isolated;
}


template <typename T>
double LoopOverRecHits(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double caloEnergyDepositionMaxDR){

  // get deposited energy for a single recHit collection inside a cone of max DR around a selected track
   
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


template <typename T>
bool checkNoDeadNoisyECALInTrackCone(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double DeadNoisyDR){

  // check if selected track points to dead or noisy ECAL channels within a cone of max DR:
   
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


void DisappearingTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  // register input collections:
  Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken( tracksToken, tracks);
  Handle<std::vector<reco::GsfElectron>> electrons;
  iEvent.getByToken( electronsToken, electrons);
  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken( muonsToken, muons);
  Handle<std::vector<reco::PFJet>> jets;
  iEvent.getByToken( jetsToken, jets);
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
  Handle<std::vector<reco::PFCandidate>> pfCands;
  iEvent.getByToken( pfCandsToken, pfCands);

  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxTrackHandle;
  iEvent.getByLabel(dEdxEstimator, dEdxTrackHandle);
  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxTrackHandle.product();

  edm::ESHandle<CaloGeometry> CaloGeomHandle;
  iSetup.get<CaloGeometryRecord>().get(CaloGeomHandle);

  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(PrimVtxToken, vtx_h);
  const reco::Vertex vtx = vtx_h->at(0);

  // register output:
  std::unique_ptr<std::vector<TLorentzVector> > chiCands(new std::vector<TLorentzVector>);

  std::unique_ptr<std::vector<double> > chiCands_dxyVtx(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_dzVtx(new std::vector<double>);
  std::unique_ptr<std::vector<int> > chiCands_trackQuality(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingOuterHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_pixelLayersWithMeasurement(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_trackerLayersWithMeasurement(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingInnerHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nMissingMiddleHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nValidPixelHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > chiCands_nValidTrackerHits(new std::vector<int>);
  std::unique_ptr<std::vector<double> > chiCands_ptError(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_chi2perNdof(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trkRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trkMiniRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trackJetIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_trackLeptonIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_matchedCaloEnergy(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_matchedCaloEnergyJets(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_deDxHarmonic2(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_chargedPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_neutralPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<double> > chiCands_neutralWithoutGammaPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044GapsVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044DeadNoisyECALVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passPFCandVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044JetIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044LepIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_passExo16044Tag(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityUndef(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityLoose(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityTight(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityHighPurity(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityConfirmed(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityGoodIterative(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityLooseSetWithPV(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityHighPuritySetWithPV(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualityDiscarded(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > chiCands_trackQualitySize(new std::vector<bool>);

  int itrack = -1;
  for( const auto& track : *tracks){
    itrack+=1;

    // track kinematics:
    // bool trackKinematics = ((track.pt() > minTrackPt) && (std::abs(track.eta())) < maxTrackEta);
    bool trackKinematics = track.pt() > minTrackPt;
    if (!(trackKinematics)) continue;

    bool passExo16044Kinematics = track.pt() > 55 && std::abs(track.eta()) < 2.1;

    // vertex consistency:
    bool vertex_match = std::abs(track.dxy(vtx.position())) < maxDxy && std::abs(track.dz(vtx.position())) < maxDz;
    // if (!(vertex_match)) continue;

    //This completes the loose candidate selection.

    // check gaps veto:
    passExo16044GapsVeto = true;
    if ((std::abs(track.eta()) > 0.15) && (std::abs(track.eta()) < 0.35) &&
        (std::abs(track.eta()) > 1.55) && (std::abs(track.eta()) < 1.85) &&
        (std::abs(track.eta()) > 1.42) && (std::abs(track.eta()) < 1.65)) passExo16044GapsVeto = false;
    chiCands_passExo16044GapsVeto->push_back(passExo16044GapsVeto);

    // check dead or noisy calo cells:
    passExo16044DeadNoisyECALVeto = checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEB, track, CaloGeomHandle, deadNoisyDR) &&
                                    checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEE, track, CaloGeomHandle, deadNoisyDR) &&
                                    checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsES, track, CaloGeomHandle, deadNoisyDR);
    chiCands_passExo16044DeadNoisyECALVeto->push_back(passExo16044DeadNoisyECALVeto);

    //now fill the products depending only on track object
    TLorentzVector tlv;
    tlv.SetPxPyPzE(track.px(),track.py(),track.pz(),track.pt());
    chiCands->push_back(tlv);
    chiCands_dxyVtx->push_back(std::abs(track.dxy(vtx.position())));
    chiCands_dzVtx->push_back(std::abs(track.dz(vtx.position())));
    reco::HitPattern hitpattern = track.hitPattern();
    chiCands_pixelLayersWithMeasurement->push_back(hitpattern.pixelLayersWithMeasurement());
    chiCands_trackerLayersWithMeasurement->push_back(hitpattern.trackerLayersWithMeasurement());
    chiCands_nMissingOuterHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS));
    chiCands_nMissingInnerHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS));
    chiCands_nMissingMiddleHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS));
    chiCands_nValidPixelHits->push_back(hitpattern.numberOfValidPixelHits());
    chiCands_nValidTrackerHits->push_back(hitpattern.numberOfValidTrackerHits());
    chiCands_ptError->push_back(track.ptError());
    chiCands_chi2perNdof->push_back(1.0*track.chi2()/track.ndof());
    chiCands_trackQualityUndef->push_back(track.quality(track.undefQuality));
    chiCands_trackQualityLoose->push_back(track.quality(track.loose));
    chiCands_trackQualityTight->push_back(track.quality(track.tight));
    chiCands_trackQualityHighPurity->push_back(track.quality(track.highPurity));
    chiCands_trackQualityConfirmed->push_back(track.quality(track.confirmed));
    chiCands_trackQualityGoodIterative->push_back(track.quality(track.goodIterative));
    chiCands_trackQualityLooseSetWithPV->push_back(track.quality(track.looseSetWithPV));
    chiCands_trackQualityHighPuritySetWithPV->push_back(track.quality(track.highPuritySetWithPV));
    chiCands_trackQualityDiscarded->push_back(track.quality(track.discarded));
    chiCands_trackQualitySize->push_back(track.quality(track.qualitySize));

    // save charged/neutral pt sum within DR cone (without selected track):
    double chargedPtSum = 0;
    double neutralPtSum = 0;
    double neutralWithoutGammaPtSum = 0;
    for( const auto& pfCand : *pfCands) {
        double dR = deltaR(track, pfCand);
        if ( (dR < maxChargedPFCandSumDR) && (pfCand.charge() != 0) )
                chargedPtSum += pfCand.pt();

        if ( (dR < maxNeutralPFCandSumDR) && (pfCand.charge() == 0) ) {
                neutralPtSum += pfCand.pt();
                if (pfCand.pdgId() != 22) neutralWithoutGammaPtSum += pfCand.pt();
        }
        // PF lepton veto
        passPFCandVeto = true;
        if ( (std::abs(pfCand.pdgId()) == 11) || (std::abs(pfCand.pdgId()) == 13) || (std::abs(pfCand.pdgId()) == 15) ) {
            if (dR < 0.01) passPFCandVeto = false;
        }
    }
    chiCands_chargedPtSum->push_back(chargedPtSum);
    chiCands_neutralPtSum->push_back(neutralPtSum);
    chiCands_neutralWithoutGammaPtSum->push_back(neutralWithoutGammaPtSum);
    chiCands_passPFCandVeto->push_back(passPFCandVeto);

    // now perform the tight selection:

    // check if fake track:
    bool passExo16044NoFakes = false;
    if ( (hitpattern.numberOfValidPixelHits() >= RequireNumberOfValidPixelHits) &&
         (hitpattern.numberOfValidTrackerHits() >= RequireNumberOfValidTrackerHits) ) {
          // checked for minimal number of hits, now require no missing inner or middle hits:
          if ( hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS) == 0 &&
               hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS) == 0) {
                // check track displacement:
                if ( vertex_match ) passExo16044NoFakes = true;
          }
    }
        
    // check calo energy deposition:
    double energyDeposited = 0;
    if (reducedEcalRecHitsEB.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsEB, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    if (reducedEcalRecHitsEE.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsEE, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    if (reducedEcalRecHitsES.isValid()) energyDeposited += LoopOverRecHits(reducedEcalRecHitsES, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    if (reducedHcalRecHitsHB.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHB, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    if (reducedHcalRecHitsHF.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHF, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    if (reducedHcalRecHitsHO.isValid()) energyDeposited += LoopOverRecHits(reducedHcalRecHitsHO, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    bool passedCaloEnergy = false;
    if (energyDeposited < caloEnergyDepositionMaxE) passedCaloEnergy = true;
    chiCands_matchedCaloEnergy->push_back(energyDeposited);

    // optional cross-check, loop over calojets:
    double energyDepositedJets = 0;
    for( const auto& jet : *caloJets){ 
         if (deltaR(jet.detectorP4().eta(), jet.detectorP4().phi(), track.eta(), track.phi()) < caloEnergyDepositionMaxDR) {
             energyDepositedJets += jet.maxEInEmTowers() + jet.maxEInHadTowers();
         }
    }
    chiCands_matchedCaloEnergyJets->push_back(energyDepositedJets);

    // check minimum of missing outer hits:
    bool passedMissingOuterHits = false;
    int missingOuterHits = hitpattern.stripLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS);
    if (missingOuterHits > minMissingOuterHits) passedMissingOuterHits = true;


    // check pt sum of tracker activity within a cone around current track:
    bool passedTrackTrackerIso = false;
    double relIso = -1;
    double miniIso = -1;
    double conePtSum_rel  = 0;
    double conePtSum_mini = 0;
    double miniConeDr = 0.2;
    if (track.pt()>50) miniConeDr = 10.0/track.pt();
    else if (track.pt()>200) miniConeDr = 0.05;
    for (const auto& othertrack : *tracks){
        if (deltaR(track, othertrack) == 0) continue;
        if (deltaR(track, othertrack) < coneRelIsoDR) conePtSum_rel += othertrack.pt();
        if (deltaR(track, othertrack) < miniConeDr) conePtSum_mini += othertrack.pt();
    }
    if (conePtSum_rel / track.pt() <= conePtSumMaxPtPercentage) passedTrackTrackerIso = true;
    relIso = conePtSum_rel / track.pt();
    miniIso = conePtSum_mini / track.pt();

    chiCands_trkRelIso->push_back(relIso);
    chiCands_trkMiniRelIso->push_back(miniIso);


    // check track/jet isolation:
    double trackIsoDR = getParticleIsolation(track, jets);
    chiCands_trackJetIso->push_back(trackIsoDR);

    bool passedTrackJetIso = false;
    if (trackIsoDR > minTrackJetDR) passedTrackJetIso = true;
    chiCands_passExo16044JetIso->push_back(passedTrackJetIso);


    // check track / lepton isolation:
    double trackLeptonIso = std::min(getParticleIsolation(track, electrons), getParticleIsolation(track, muons));
    chiCands_trackLeptonIso->push_back(trackLeptonIso);

    bool passedTrackLeptonIso = false;
    if (trackLeptonIso > minTrackLeptonDR) passedTrackLeptonIso = true;
    chiCands_passExo16044LepIso->push_back(passedTrackLeptonIso);
  

    // save dE/dx:
    double dedx = -1;
    if (doDeDx) dedx = dEdxTrack[reco::TrackRef(tracks, itrack)].dEdx();
    chiCands_deDxHarmonic2->push_back(dedx);


    // check if EXO-16-044 selection criteria are satisfied:
    if (passExo16044Kinematics && passExo16044GapsVeto && passExo16044DeadNoisyECALVeto && passedTrackTrackerIso &&
        passExo16044NoFakes && passedTrackLeptonIso && passedTrackJetIso && 
        passedCaloEnergy && passedMissingOuterHits) chiCands_passExo16044Tag->push_back(true);
    else chiCands_passExo16044Tag->push_back(false);

  }

  // save track collection in event:
  iEvent.put(std::move(chiCands),                                 "chiCands");
  iEvent.put(std::move(chiCands_dxyVtx),                          "chiCands@dxyVtx");
  iEvent.put(std::move(chiCands_dzVtx),                           "chiCands@dzVtx");
  iEvent.put(std::move(chiCands_pixelLayersWithMeasurement),      "chiCands@pixelLayersWithMeasurement");
  iEvent.put(std::move(chiCands_trackerLayersWithMeasurement),    "chiCands@trackerLayersWithMeasurement");
  iEvent.put(std::move(chiCands_nMissingOuterHits),               "chiCands@nMissingOuterHits");
  iEvent.put(std::move(chiCands_nMissingInnerHits),               "chiCands@nMissingInnerHits");
  iEvent.put(std::move(chiCands_nMissingMiddleHits),              "chiCands@nMissingMiddleHits");
  iEvent.put(std::move(chiCands_nValidPixelHits),                 "chiCands@nValidPixelHits");
  iEvent.put(std::move(chiCands_nValidTrackerHits),               "chiCands@nValidTrackerHits");
  iEvent.put(std::move(chiCands_ptError),                         "chiCands@ptError");
  iEvent.put(std::move(chiCands_chi2perNdof),                     "chiCands@chi2perNdof");
  iEvent.put(std::move(chiCands_trkRelIso),                       "chiCands@trkRelIso");
  iEvent.put(std::move(chiCands_trkMiniRelIso),                   "chiCands@trkMiniRelIso");
  iEvent.put(std::move(chiCands_trackJetIso),                     "chiCands@trackJetIso");
  iEvent.put(std::move(chiCands_trackLeptonIso),                  "chiCands@trackLeptonIso");
  iEvent.put(std::move(chiCands_matchedCaloEnergy),               "chiCands@matchedCaloEnergy");
  iEvent.put(std::move(chiCands_matchedCaloEnergyJets),           "chiCands@matchedCaloEnergyJets");
  iEvent.put(std::move(chiCands_deDxHarmonic2),                   "chiCands@deDxHarmonic2");
  iEvent.put(std::move(chiCands_chargedPtSum),                    "chiCands@chargedPtSum");
  iEvent.put(std::move(chiCands_neutralPtSum),                    "chiCands@neutralPtSum");
  iEvent.put(std::move(chiCands_neutralWithoutGammaPtSum),        "chiCands@neutralWithoutGammaPtSum");
  iEvent.put(std::move(chiCands_passExo16044GapsVeto),            "chiCands@passExo16044GapsVeto");
  iEvent.put(std::move(chiCands_passExo16044DeadNoisyECALVeto),   "chiCands@passExo16044DeadNoisyECALVeto");
  iEvent.put(std::move(chiCands_passPFCandVeto),                  "chiCands@passPFCandVeto");
  iEvent.put(std::move(chiCands_passExo16044JetIso),              "chiCands@passExo16044JetIso");
  iEvent.put(std::move(chiCands_passExo16044LepIso),              "chiCands@passExo16044LepIso");
  iEvent.put(std::move(chiCands_passExo16044Tag),                 "chiCands@passExo16044Tag");
  iEvent.put(std::move(chiCands_trackQualityUndef),               "chiCands@trackQualityUndef");
  iEvent.put(std::move(chiCands_trackQualityLoose),               "chiCands@trackQualityLoose");
  iEvent.put(std::move(chiCands_trackQualityTight),               "chiCands@trackQualityTight");
  iEvent.put(std::move(chiCands_trackQualityHighPurity),          "chiCands@trackQualityHighPurity");
  iEvent.put(std::move(chiCands_trackQualityConfirmed),           "chiCands@trackQualityConfirmed");
  iEvent.put(std::move(chiCands_trackQualityGoodIterative),       "chiCands@trackQualityGoodIterative");
  iEvent.put(std::move(chiCands_trackQualityLooseSetWithPV),      "chiCands@trackQualityLooseSetWithPV");
  iEvent.put(std::move(chiCands_trackQualityHighPuritySetWithPV), "chiCands@trackQualityHighPuritySetWithPV");
  iEvent.put(std::move(chiCands_trackQualityDiscarded),           "chiCands@trackQualityDiscarded");
  iEvent.put(std::move(chiCands_trackQualitySize),                "chiCands@trackQualitySize");

}

void DisappearingTrackProducer::beginJob() {}

void DisappearingTrackProducer::endJob() {}

void DisappearingTrackProducer::beginRun(edm::Run&, edm::EventSetup const&) {}

void DisappearingTrackProducer::endRun(edm::Run&, edm::EventSetup const&) {}

void DisappearingTrackProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

void DisappearingTrackProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

void DisappearingTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisappearingTrackProducer);
