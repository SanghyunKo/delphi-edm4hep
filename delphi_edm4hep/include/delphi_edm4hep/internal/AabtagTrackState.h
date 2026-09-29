// AabtagTrackState.h — internal (not exported).
//
// AABTAG measures an impact parameter for the subset of tracks it can use,
// against its own primary vertex. That belongs on the track, so it is carried
// as a TrackState at AtVertex. Building it lives here because both writers
// that own a track need it: TrackingWriter in pass 1 and TrackHybridWriter in
// pass 2, where the fullDST tracks are clones and would otherwise lose it.

#pragma once

#include <edm4hep/TrackState.h>

#include <unordered_map>

namespace delphi_edm4hep::aabtag {

// PA address -> AABTAG track index (1..NTRK). Empty when AABTAG produced
// nothing usable for this event.
std::unordered_map<int, int> lpaToTrack();

// The state AABTAG measured for its track `i`, whose charge is `charge`
// (+1/-1; AABTAG's impact-parameter sign depends on it).
edm4hep::TrackState vertexState(int i, int charge);

}  // namespace delphi_edm4hep::aabtag
