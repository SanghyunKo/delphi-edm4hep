// TrackElementsWriter — implementation.

#include "delphi_edm4hep/Tracking/TrackElements.h"

#include "delphi_edm4hep/Helix.h"
#include "delphi_edm4hep/TeBank.h"
#include "delphi_edm4hep/internal/PaWalk.h"

#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/UserDataCollection.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

namespace ph = phdst;

namespace delphi_edm4hep::track_elements {

namespace {

constexpr float kCm2Mm = 10.f;

// Marks a track-state component the module did not measure. Zero cannot be
// used: it is a legal measured value.
constexpr float kNotMeasured = std::numeric_limits<float>::quiet_NaN();

// TE-basis field indices used by te_bank::Decoded::has_meas. The first two
// coordinates are (x, y) on a plane and (R*Phi, z) on a cylinder.
constexpr int kBasisX     = 0;
constexpr int kBasisY     = 1;
constexpr int kBasisRPhi  = 1;
constexpr int kBasisZ     = 2;
constexpr int kBasisTheta = 3;
constexpr int kBasisPhi   = 4;
constexpr int kBasisInvP  = 5;

// TEVF measures two coordinates and no direction, so it is a plane hit; every
// other module measures at least phi and becomes a segment.
struct TeModule { const char* name; bool plane; };
constexpr std::array<TeModule, 8> kTeModules{{
  {"TEID", false}, {"TETP", false}, {"TEOD", false}, {"TEFA", false},
  {"TEFB", false}, {"TERF", false}, {"TEST", false}, {"TEVF", true },
}};

}  // namespace

void TrackElementsWriter::emit()
{
  const std::size_t n = kTeModules.size();
  edm4hep::TrackCollection             segments;
  edm4hep::TrackerHitPlaneCollection   planes;
  podio::UserDataCollection<float>     segment_length;
  podio::UserDataCollection<float>     plane_length;

  Output out;

  // Field strength as the event recorded it, so both passes agree. Zero when
  // absent, which yields zero curvature rather than a wrong one.
  const double B = frame_.getParameter<float>(
      makeName("EVT", "BField")).value_or(0.f);

  pawalk::forEachPA([&](int lpa, int paIdx) {
    if (paIdx >= static_cast<int>(out.pa_to_segments.size())) {
      out.pa_to_segments.resize(paIdx + 1);
      out.pa_to_plane_hits.resize(paIdx + 1);
    }
    for (std::size_t m = 0; m < n; ++m) {
      const int lte = pawalk::lphpa(kTeModules[m].name, lpa);
      if (lte <= 0) continue;
      const auto mod = te_bank::decodeModule(lte, pawalk::iphreq());

      for (const auto& te : mod.elements) {
        // The label carries the detector and the reconstruction stage, so one
        // integer identifies which TE produced this object.
        const std::int32_t kind = mod.label * 10 + mod.stage;

        if (kTeModules[m].plane) {
          auto hit = planes.create();
          hit.setType(kind);
          hit.setQuality(te.descriptor);
          hit.setPosition({ te.coord[0] * kCm2Mm,
                            te.coord[1] * kCm2Mm,
                            te.coord[2] * kCm2Mm });
          out.pa_to_plane_hits[paIdx].push_back(hit);
          plane_length.push_back(te.length * kCm2Mm);
          continue;
        }

        const auto helix = Helix::fromTrackElement(
            te.coord[0], te.coord[1], te.coord[2], te.theta, te.phi,
            te.invP, te.invPt, te.is_cylindrical, te.cov, B);

        auto state = helix.toTrackState(edm4hep::TrackState::AtOther);
        // The reference point is on the fitted segment, so the helix passes
        // through it: D0 = Z0 = 0 wherever that coordinate was measured. An
        // unmeasured coordinate is a placeholder (e.g. z on an inner-detector
        // jet-chamber element), so the offset it would define is not measured
        // either. A cylinder fixes R and a plane fixes z.
        const bool transverse = te.is_cylindrical
            ? te.has_meas[kBasisRPhi]
            : te.has_meas[kBasisX] && te.has_meas[kBasisY];
        const bool longitudinal = te.is_cylindrical ? te.has_meas[kBasisZ] : true;
        state.D0 = transverse   ? 0.f : kNotMeasured;
        state.Z0 = longitudinal ? 0.f : kNotMeasured;
        if (!te.has_meas[kBasisPhi])   state.phi       = kNotMeasured;
        if (!te.has_meas[kBasisTheta]) state.tanLambda = kNotMeasured;
        if (!te.has_meas[kBasisInvP])  state.omega     = kNotMeasured;

        auto trk = segments.create();
        trk.setType(kind);
        trk.setChi2(te.chi2);
        trk.setNdf(static_cast<std::int32_t>(std::lround(te.ndf)));
        trk.addToTrackStates(state);
        out.pa_to_segments[paIdx].push_back(trk);
        segment_length.push_back(te.length * kCm2Mm);
      }
    }
  });

  // Emitted whether or not any module is on the file, so the collection set
  // does not vary between samples.
  put(std::move(segments),       "TE",   "Segments",       Provenance::Transcribed);
  put(std::move(segment_length), "TE",   "Segments_Length", Provenance::Transcribed);
  put(std::move(planes),         "TEVF", "TrackElementPlane",
      Provenance::Transcribed);
  put(std::move(plane_length),   "TEVF", "TrackElementPlane_Length",
      Provenance::Transcribed);

  ctx_.track_elements = std::move(out);
}

}  // namespace delphi_edm4hep::track_elements
