// MatchProvenanceWriter — pass-2 implementation.
//
// Reads sDST_TRAC_Tracks + sDST_MAIN_Particles from frame_ (loaded
// from the intermediate by the harness). Walks the fDST PA chain,
// computes a perigee match for each charged PA, builds:
//   - ctx_.fdst_pa_to_sdst_track   (paIdx -> sDST Track idx; -1 = no match)
//   - sdst_track_matched           (sDST Track idx -> bool)
// then emits the per-Particle provenance flag.

#include "delphi_edm4hep/Tracking/MatchProvenance.h"

#include "delphi_edm4hep/Helix.h"
#include "delphi_edm4hep/PerigeeMatch.h"
#include "delphi_edm4hep/internal/PaWalk.h"

#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackCollection.h>
#include <podio/UserDataCollection.h>

#include <cmath>
#include <cstdint>
#include <vector>

namespace ph = phdst;

namespace delphi_edm4hep::matchprov {

namespace {

// Build sDST Track-index -> sDST Particle-index lookup by walking the
// Particles collection's `tracks` relation once. O(N²) but N ~ 30.
// Returns a vector of size sdst_tracks.size(); -1 = no owning Particle
// (shouldn't happen in practice — every Track has one Particle).
std::vector<int> buildTrackToParticleMap(
  const edm4hep::TrackCollection& sdst_tracks,
  const edm4hep::ReconstructedParticleCollection& sdst_particles)
{
  std::vector<int> map(sdst_tracks.size(), -1);
  for (std::size_t i = 0; i < sdst_particles.size(); ++i) {
    const auto p = sdst_particles[i];
    if (p.tracks_size() == 0) continue;
    const auto p_trk = p.getTracks(0);
    for (std::size_t j = 0; j < sdst_tracks.size(); ++j) {
      if (sdst_tracks[j] == p_trk) {
        map[j] = static_cast<int>(i);
        break;
      }
    }
  }
  return map;
}

}  // namespace

void MatchProvenanceWriter::emit()
{
  // The sDST_TRAC_Tracks and sDST_MAIN_Particles collections come
  // from the intermediate frame the harness loaded. If they're
  // missing (e.g. pass-1 produced no charged tracks), emit empty
  // provenance and bail.
  const auto& sdst_tracks =
    frame_.get<edm4hep::TrackCollection>(sdstName("TRAC", "Tracks"));
  const auto& sdst_particles =
    frame_.get<edm4hep::ReconstructedParticleCollection>(sdstName("MAIN", "Particles"));

  podio::UserDataCollection<std::int32_t> prov;
  prov.resize(sdst_particles.size());
  for (std::size_t i = 0; i < sdst_particles.size(); ++i) {
    // -1 for neutral, 0 for charged-not-yet-matched. Charged that
    // get matched below are flipped to +1.
    prov[i] = (sdst_particles[i].getCharge() == 0.f) ? -1 : 0;
  }

  // Pre-compute the Track -> Particle index map so the fDST loop can
  // mark provenance in O(1) per match.
  const auto track_to_particle =
    buildTrackToParticleMap(sdst_tracks, sdst_particles);

  // Walk fDST PAs, compute perigee match for each charged PA. Build
  // both PA -> Track and PA -> Particle maps in one pass.
  std::vector<int> fdst_pa_to_track;
  std::vector<int> fdst_pa_to_particle;
  auto extend_maps_to = [&](int paIdx) {
    if (paIdx >= static_cast<int>(fdst_pa_to_track.size())) {
      fdst_pa_to_track.resize   (static_cast<std::size_t>(paIdx + 1), -1);
      fdst_pa_to_particle.resize(static_cast<std::size_t>(paIdx + 1), -1);
    }
  };
  pawalk::forEachPA([&](int lpa, int paIdx) {
    extend_maps_to(paIdx);

    const int lmain = pawalk::lphpa("MAIN", lpa);
    if (lmain <= 0) return;
    const int charge_code =
      static_cast<int>(std::lround(ph::Q(lmain + 8)));
    if (charge_code == 0) return;   // neutral — no perigee match

    const int ltrac = pawalk::lphpa("TRAC", lpa);
    if (ltrac <= 0) return;

    // DELPHI perigee (cm, rad, 1/cm), converted exactly as the sDST tracks
    // it is matched against.
    const auto h = Helix::fromPerigee(ph::Q(ltrac + 2), ph::Q(ltrac + 3),
                                      ph::Q(ltrac + 4), ph::Q(ltrac + 5),
                                      ph::Q(ltrac + 6));
    const double D0_mm    = h.params().D0;
    const double Z0_mm    = h.params().Z0;
    const double omega_pm = h.params().omega;

    // Charge sign — same encoding as Tracking.cpp (code 3 = undefined
    // -> 0, no perigee match for ambiguous charge).
    int charge_signed = 0;
    if      (charge_code == 1) charge_signed = +1;
    else if (charge_code == 2) charge_signed = -1;
    if (charge_signed == 0) return;

    const int matched = perigee_match::findMatch(
      sdst_tracks, D0_mm, Z0_mm, omega_pm, charge_signed);
    fdst_pa_to_track[paIdx] = matched;

    if (matched >= 0 && matched < static_cast<int>(track_to_particle.size())) {
      const int particle_idx = track_to_particle[matched];
      if (particle_idx >= 0 && particle_idx < static_cast<int>(prov.size())) {
        prov[particle_idx] = +1;
        fdst_pa_to_particle[paIdx] = particle_idx;
      }
    }
  });

  // Expose the per-PA -> Track / -> Particle maps for downstream
  // pass-2 writers (TOF / MTPC / TE-merge).
  ctx_.fdst_pa_to_sdst_track    = std::move(fdst_pa_to_track);
  ctx_.fdst_pa_to_sdst_particle = std::move(fdst_pa_to_particle);

  // Bank prefix is MAIN, not TRAC: this collection is parallel to
  // fDST_MAIN_Particles (charged + neutral), NOT to fDST_TRAC_Tracks
  // (charged only). The -1 = neutral value is a per-particle concept.
  put(std::move(prov), "MAIN", "Particles_MatchProvenance", Provenance::Custom);
}

}  // namespace delphi_edm4hep::matchprov
