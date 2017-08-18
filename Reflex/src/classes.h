#include <vector>
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h" 

namespace {
  struct dictionary {
    std::vector<TLorentzVector> vlv;
    std::vector<pat::Jet> vpj;
    std::vector<std::vector<TLorentzVector> > vvlv;
    std::vector<std::vector<pat::Jet> > vvpj;
    std::vector<reco::Track> vtrk;
	edm::PtrVector<pat::PackedCandidate> rv2pp;
    edm::Wrapper<std::vector<TLorentzVector> > wvlv;
    edm::Wrapper<std::vector<std::vector<TLorentzVector> > > wvvlv;
	edm::Wrapper<edm::PtrVector<pat::PackedCandidate> > wrv2pp;
    edm::Wrapper<std::vector<reco::Track> > wrvtrk;
  };
}
