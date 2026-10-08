// TrackHybridWriter — pass-2 implementation.

#include "delphi_edm4hep/Tracking/TrackHybrid.h"

#include "delphi_edm4hep/internal/AabtagTrackState.h"
#include "delphi_edm4hep/internal/PaWalk.h"

#include <edm4hep/MutableTrack.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace delphi_edm4hep::track_hybrid {

namespace {

// Clone the scalar / vector fields of a Track. Only the AtIP state is copied:
// the extrapolation states belong to the pass that decoded them, and pass 2
// appends its own from the fullDST, which carries far more of them.
// Cross-collection relations are not copied either: the fDST track gets its
// own, pointing at the fDST track elements.
void cloneTrackShallow(edm4hep::MutableTrack dst, const edm4hep::Track& src) {
  dst.setType   (src.getType());
  dst.setChi2   (src.getChi2());
  dst.setNdf    (src.getNdf());
  dst.setNholes (src.getNholes());
  for (auto n : src.getSubdetectorHitNumbers())  dst.addToSubdetectorHitNumbers(n);
  for (auto n : src.getSubdetectorHoleNumbers()) dst.addToSubdetectorHoleNumbers(n);
  for (auto ts : src.getTrackStates()) {
    if (ts.location == edm4hep::TrackState::AtIP) dst.addToTrackStates(ts);
  }
}

}  // namespace

void TrackHybridWriter::emit()
{
  edm4hep::TrackCollection trk_out;

  if (!ctx_.fdst_pa_to_sdst_track) {
    put(std::move(trk_out), "TRAC", "Tracks", Provenance::Derived);
    return;
  }
  const auto& pa_to_track = *ctx_.fdst_pa_to_sdst_track;
  const auto& sdst_tracks =
    frame_.get<edm4hep::TrackCollection>(sdstName("TRAC", "Tracks"));

  // AABTAG's impact parameters ride on the track, and these tracks are clones
  // that do not carry the pass-1 state, so they are attached again here from
  // the fullDST commons. Keyed by PA address, which the walk gives per index.
  const auto lpa_to_btag = aabtag::lpaToTrack();
  std::unordered_map<int, int> pa_to_lpa;
  pawalk::forEachPA([&](int lpa, int paIdx) { pa_to_lpa.emplace(paIdx, lpa); });

  // sDST Track index -> the fDST PAs the perigee match resolved to it.
  // Usually 1:1, but many-to-one is tolerated.
  std::vector<std::vector<int>> track_to_pas(sdst_tracks.size());
  for (int paIdx = 0; paIdx < static_cast<int>(pa_to_track.size()); ++paIdx) {
    const int t = pa_to_track[paIdx];
    if (t >= 0 && t < static_cast<int>(sdst_tracks.size())) {
      track_to_pas[t].push_back(paIdx);
    }
  }

  for (std::size_t i = 0; i < sdst_tracks.size(); ++i) {
    auto out = trk_out.create();
    cloneTrackShallow(out, sdst_tracks[i]);

    if (!ctx_.track_elements) continue;
    const auto& te = *ctx_.track_elements;
    for (int paIdx : track_to_pas[i]) {
      if (paIdx < static_cast<int>(te.pa_to_segments.size())) {
        for (const auto& seg : te.pa_to_segments[paIdx]) out.addToTracks(seg);
      }
      if (paIdx < static_cast<int>(te.pa_to_plane_hits.size())) {
        for (const auto& hit : te.pa_to_plane_hits[paIdx]) {
          out.addToTrackerHits(hit);
        }
      }
      if (ctx_.trax &&
          paIdx < static_cast<int>(ctx_.trax->pa_to_states.size())) {
        for (const auto& st : ctx_.trax->pa_to_states[paIdx]) {
          out.addToTrackStates(st);
        }
      }
      if (const auto pa = pa_to_lpa.find(paIdx); pa != pa_to_lpa.end()) {
        if (const auto b = lpa_to_btag.find(pa->second); b != lpa_to_btag.end()) {
          out.addToTrackStates(
              aabtag::vertexState(b->second, pawalk::trueCharge(pa->second)));
        }
      }
    }
  }

  put(std::move(trk_out), "TRAC", "Tracks", Provenance::Derived);
}

}  // namespace delphi_edm4hep::track_hybrid
