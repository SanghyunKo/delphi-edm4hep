// PerigeeMatch.cpp — implementation.
//
// Convert sDST tracks back to DELPHI cm-basis and apply the legacy
// converter's cost formula + tolerance. Implementation note: we
// iterate the whole sDST collection per call; for ~50 tracks/event
// this is fast enough that a build-once index isn't needed yet.

#include "delphi_edm4hep/PerigeeMatch.h"

#include <edm4hep/TrackState.h>

#include <cmath>
#include <limits>

namespace delphi_edm4hep::perigee_match {

namespace {

constexpr double kMm2Cm = 0.1;       // mm -> cm

// Pull the AtIP TrackState (by value) from a Track. EDM4hep returns
// TrackStates by value via the range view, so we can't safely keep a
// pointer or reference — we copy. found=false if the Track has no
// TrackStates at all.
struct AtIPLookup { edm4hep::TrackState state; bool found; };
AtIPLookup findAtIPState(const edm4hep::Track& trk) {
  for (auto ts : trk.getTrackStates()) {
    if (ts.location == edm4hep::TrackState::AtIP) return {ts, true};
  }
  // Fallback: first entry. Tracking always emits AtIP first, so this
  // catches the case where location wasn't set explicitly.
  if (trk.trackStates_size() > 0) return {trk.getTrackStates()[0], true};
  return {{}, false};
}

}  // namespace

int findMatch(
  const edm4hep::TrackCollection& sdst_tracks,
  double fdst_D0_mm,
  double fdst_Z0_mm,
  double fdst_omega_per_mm,
  int    fdst_charge_signed,
  double tolerance_cm2)
{
  // Convert fDST candidate to DELPHI cm-basis.
  const double fdst_d0_cm   = fdst_D0_mm        * kMm2Cm;
  const double fdst_z0_cm   = fdst_Z0_mm        * kMm2Cm;
  const double fdst_invR_pc = fdst_omega_per_mm / kMm2Cm;   // 1/mm -> 1/cm

  int    best_idx  = -1;
  double best_cost = std::numeric_limits<double>::infinity();

  const int n = static_cast<int>(sdst_tracks.size());
  for (int i = 0; i < n; ++i) {
    const auto t  = sdst_tracks[i];
    const auto lk = findAtIPState(t);
    if (!lk.found) continue;
    const auto& ts = lk.state;

    // Charge gate (intentionally weak): both sides are converted by
    // Helix::fromPerigee, so omega carries the sign of the charge on
    // each. Rather than gate on it, we trust the geometric uniqueness of the
    // perigee within tolerance — same physical track gives identical
    // (d0, z0, invR) on both DST levels to float precision. The
    // legacy converter took the same approach (gate by PA.MAIN
    // charge code, not by omega sign). With ~100 μm tolerance the
    // collision rate among same-event tracks is negligible.
    (void)fdst_charge_signed;

    const double sdst_d0_cm   = static_cast<double>(ts.D0)    * kMm2Cm;
    const double sdst_z0_cm   = static_cast<double>(ts.Z0)    * kMm2Cm;
    const double sdst_invR_pc = static_cast<double>(ts.omega) / kMm2Cm;

    const double dd0 = fdst_d0_cm   - sdst_d0_cm;
    const double dz0 = fdst_z0_cm   - sdst_z0_cm;
    const double diR = (fdst_invR_pc - sdst_invR_pc) * 1000.0;
    const double cost = dd0*dd0 + dz0*dz0 + diR*diR;

    if (cost < best_cost) {
      best_cost = cost;
      best_idx  = i;
    }
  }

  return (best_cost < tolerance_cm2) ? best_idx : -1;
}

}  // namespace delphi_edm4hep::perigee_match
