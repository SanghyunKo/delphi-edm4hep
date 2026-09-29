// TraxWriter — implementation. Runs in both passes, before the track writers.
//
// PA.TRAX (label 20) holds the track extrapolated onto a set of named
// surfaces: its own first measured point, the surface just in front of it, and
// the named detectors, ordered by increasing R in the barrel and |Z| in the
// endcap (dst_content.txt:1917-1945).
//
// Module layout, 1-based from the LPHPA address:
//   +2  number of surfaces
//   then per point, the first at +3:
//     +0  number of following words, 8 without a covariance or 23 with one
//     +1  detector id      +2  IBAR, 1 = (R, R*Phi, z) / 0 = (x, y, z)
//     +3..+5  coordinates (cm)   +6 theta   +7 phi   +8 signed 1/p
//     +9..+23 covariance
//   the next point starts one word past the last, the count word not being
//   included in its own count.
//
// Each point becomes an edm4hep::TrackState on the track built from the same
// PA, so a reader walks TRAC_Tracks -> trackStates rather than a side
// collection. The states carry no detector field, but `location` distinguishes
// the first measured point and the calorimeter crossings, and the surfaces are
// separable by radius in any case.

#include "delphi_edm4hep/Tracking/Trax.h"

#include "delphi_edm4hep/Helix.h"
#include "delphi_edm4hep/internal/PaWalk.h"

#include <edm4hep/TrackState.h>

#include <cmath>
#include <string>

namespace ph = phdst;

namespace delphi_edm4hep::trax {

namespace {

// A point record is the count word plus that many words.
constexpr int kCountWord   = 1;
constexpr int kWordsNoCov  = 8;
constexpr int kWordsWithCov = 23;
constexpr int kCovWords    = 15;

// TANAGRA detector id -> edm4hep TrackState location. The calorimeter surfaces
// are HPC(9), HAB(13), HAF(22) and EMF(26); id 0 is the track's own first
// measured point; TOF(11) and the muon chambers (14, 17, 30) have no dedicated
// location. Ids are the ones PXTRAX writes (pxdst34.car:15991-15994), named
// through TANAGRA's own table (tanagra322.car:12164-12176).
int locationForDetector(int det_id) {
  switch (det_id) {
    case 0:  return edm4hep::TrackState::AtFirstHit;
    case 9: case 13: case 22: case 26:
             return edm4hep::TrackState::AtCalorimeter;
    default: return edm4hep::TrackState::AtOther;
  }
}

// The 15 covariance words are a 5x5 lower-triangular row-major matrix over the
// free parameters: the three stored coordinates minus the one the surface
// fixes, then theta, phi and 1/p (exx.car:5035-5036).
//   cylinder: (R*Phi, z, theta, phi, 1/p)   plane: (x, y, theta, phi, 1/p)
// Scatter it into the TE basis Helix::fromTrackElement reads.
CovMatrix6 covariance(int lpt, bool cylindrical) {
  static constexpr int kCylinder[5] = {1, 2, 3, 4, 5};
  static constexpr int kPlane   [5] = {0, 1, 3, 4, 5};
  const int* idx = cylindrical ? kCylinder : kPlane;

  CovMatrix6 cov{};
  int slot = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j, ++slot) {
      cov[covOffset(idx[i], idx[j])] = ph::Q(lpt + kWordsNoCov + 1 + slot);
    }
  }
  return cov;
}

int nint(float x) { return static_cast<int>(std::lround(x)); }

}  // namespace

void TraxWriter::emit()
{
  Output out;

  // Field strength as the event recorded it, so both passes agree.
  const double B = frame_.getParameter<float>(
      makeName("EVT", "BField")).value_or(0.f);

  pawalk::forEachPA([&](int lpa, int paIdx) {
    if (paIdx >= static_cast<int>(out.pa_to_states.size())) {
      out.pa_to_states.resize(paIdx + 1);
    }
    const int ltrax = pawalk::lphpa("TRAX", lpa);
    if (ltrax <= 0) return;

    const int blen   = pawalk::iphreq();
    const int npts   = nint(ph::Q(ltrax + 2));
    const int q      = pawalk::trueCharge(lpa);

    int lpt = ltrax + 3;
    for (int ip = 0; ip < npts; ++ip) {
      const int n_words = nint(ph::Q(lpt));
      if (n_words < kWordsNoCov) break;                    // malformed
      if (lpt + n_words > ltrax + blen) break;             // past the bank

      const int  det_id      = nint(ph::Q(lpt + 1));
      const bool cylindrical = nint(ph::Q(lpt + 2)) != 0;  // IBAR
      const bool has_cov     = n_words >= kWordsWithCov;

      // The central extrapolator writes 1/p signed geometrically, but the
      // muon/calorimeter one (EXXTR: TOF, HAB, HAF, MUB, MUF, SMC) writes the
      // unsigned 1/|p| (exx.car:9858, 9902) despite the bank documentation. Its
      // covariance row still refers to the signed word, so restore that sign:
      // geometric = -charge. Harmless on the points already signed.
      const double word = ph::Q(lpt + 8);
      const double invP = q != 0 ? -q * std::fabs(word) : word;

      const auto helix = Helix::fromTrackElement(
          ph::Q(lpt + 3), ph::Q(lpt + 4), ph::Q(lpt + 5),
          ph::Q(lpt + 6), ph::Q(lpt + 7), invP,
          /*invPt=*/false,
          cylindrical,
          has_cov ? covariance(lpt, cylindrical) : CovMatrix6{},
          B);

      out.pa_to_states[paIdx].push_back(
          helix.toTrackState(locationForDetector(det_id)));

      lpt += n_words + kCountWord;
    }
  });

  ctx_.trax = std::move(out);
}

}  // namespace delphi_edm4hep::trax
