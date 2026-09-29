// Helix.h — DELPHI ↔ EDM4hep track-parameter conversion (public API).
//
// One value type that consolidates what used to be the free functions in
// BasisConversion.h (TE-basis ↔ helix, omega ↔ momentum) and
// PerigeeTrack.h (perigee weight-matrix → helix). Construct from a source
// basis via a named factory (the conversion runs inward); read out the
// EDM4hep helix params / cov / momentum, or emit an edm4hep::TrackState.
//
// The internal canonical basis is the EDM4hep helix
// (D0, phi, omega, Z0, tanLambda, time) in mm / rad / 1-over-mm, with a
// 6×6 lower-tri covariance (time row/col always 0). Analysts can use this
// header standalone to convert without running the converter.
//
//   omega [1/mm] = transverse curvature with the sign of the charge (EDM4hep).
//     DELPHI signs 1/R and 1/p geometrically: positive is counter-clockwise
//     seen from +z, opposite to the charge in its field. The perigee path sets
//     omega = -invR/10; fromTrackElement sets omega = -kOmega·B·invP, divided
//     by sin(theta) when the word is 1/p rather than 1/p_T (invPt flag).
//
// Bank *parsing* (te_bank::decode for the TE descriptor/cov, hpc::padDecode
// for PXHGET) is a separate concern and feeds the factories here — it is
// deliberately NOT part of this class.

#pragma once

#include <edm4hep/TrackState.h>
#include <edm4hep/Vector3f.h>

#include <array>

namespace delphi_edm4hep {

// 6×6 symmetric covariance, lower-triangular storage (21 floats):
//   off = i*(i+1)/2 + j   for i >= j.
// Shared by the TE-basis cov (te_bank::decode) and the helix-basis cov.
using CovMatrix6 = std::array<float, 21>;

inline constexpr int covOffset(int i, int j) {
  return (i >= j) ? (i * (i + 1) / 2 + j) : (j * (j + 1) / 2 + i);
}

// Speed-of-light curvature constant (mm·GeV/T).
inline constexpr double kOmega = 2.99792458e-4;

// EDM4hep helix parameters.
struct HelixParams {
  float D0;
  float phi;
  float omega;
  float Z0;
  float tanLambda;
  float time;
};

struct Momentum {
  double px;
  double py;
  double pz;
  double p;   // |p|
};

class Helix {
 public:
  // ---- construct from a source basis (convert inward) --------------------

  // DELPHI perigee at the IP (PA.TRAC / PA.ELTR):
  //   d0, z0 [cm], theta, phi [rad], invR [1/cm, signed].
  //   W = 15-element lower-tri weight matrix (= cov^-1); pass nullptr to
  //   skip cov. valid() is false on degenerate sin(theta) ~ 0.
  static Helix fromPerigee(float d0, float z0, float theta, float phi,
                           float invR,
                           const std::array<float, 15>* W = nullptr);

  // DELPHI TE / TRAX surface (PA.TE*, PA.TRAX): reference point
  //   (c1, c2, c3) [cm] + direction (theta, phi) [rad] + invP [1/GeV],
  //   TE-basis 6×6 cov and B [T]. invP must be signed geometrically, as
  //   DELPHI writes it (positive = counter-clockwise). D0 = Z0 = 0 at the
  //   reference point by construction.
  //   invPt: true if `invP` is 1/|p_T| (TE descriptor bit 11), false if 1/|p|;
  //   omega is formed accordingly so it equals the transverse curvature.
  //   cylindrical: true when the surface is a cylinder (TE descriptor bit 1),
  //   where the coordinates are stored as (R, R*Phi, z) rather than (x, y, z).
  //   Such an element is confined to its surface, so R is fixed and the
  //   covariance covers only (R*Phi, z): expressed in Cartesian terms its
  //   transverse block is rank one, along the azimuthal direction.
  static Helix fromTrackElement(double c1, double c2, double c3,
                                double theta, double phi, double invP,
                                bool invPt, bool cylindrical,
                                const CovMatrix6& teCov,
                                double B);

  // Already in the EDM4hep helix basis (no cov).
  static Helix fromHelix(float D0, float phi, float omega,
                         float Z0, float tanLambda);

  // ---- read out (convert outward) ----------------------------------------
  bool valid()  const { return valid_; }
  bool hasCov() const { return hasCov_; }
  const HelixParams& params() const { return p_; }
  const CovMatrix6&  cov()    const { return cov_; }

  // omega + direction + (B, q) -> momentum 4-vector.
  Momentum momentum(double B, int charge) const;

  // Emit an edm4hep::TrackState at the given location (e.g.
  // edm4hep::TrackState::AtIP / AtOther / AtFirstHit). covMatrix is filled
  // only when hasCov(); referencePoint comes from the source factory.
  edm4hep::TrackState toTrackState(int location) const;

 private:
  Helix() = default;

  HelixParams       p_{};
  CovMatrix6        cov_{};
  edm4hep::Vector3f refPoint_{0.f, 0.f, 0.f};
  bool              valid_  = true;
  bool              hasCov_ = false;
};

}  // namespace delphi_edm4hep
