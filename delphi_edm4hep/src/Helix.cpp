// Helix.cpp — implementation.
//
// The conversion kernels are transplanted verbatim from the former
// BasisConversion.cpp (TE-basis 6×6 Jacobian) and PerigeeTrack.cpp
// (perigee 5×5 Jacobian + weight-matrix inversion); the numbers are
// unchanged, only the packaging differs.

#include "delphi_edm4hep/Helix.h"

#include <edm4hep/Constants.h>   // edm4hep::TrackParams

#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include <cmath>

namespace delphi_edm4hep {

namespace {

constexpr double kCm2Mm = 10.0;

// ---- perigee (raw DELPHI) -> helix 5×5 Jacobian -------------------------
//   raw   : (d0, z0, theta, phi, invR)
//   helix : (D0, phi, omega, Z0, tanLambda)
TMatrixD perigeeJacobian(double theta_dp) {
  TMatrixD J(5, 5);
  J.Zero();
  const double s      = std::sin(theta_dp);
  const double inv_s2 = (s != 0.0) ? 1.0 / (s * s) : 0.0;
  J(0, 0) = -kCm2Mm;
  J(1, 3) = +1.0;
  J(2, 4) = -1.0 / kCm2Mm;      // omega = -invR/10, see fromPerigee
  J(3, 1) = +kCm2Mm;
  J(4, 2) = -inv_s2;
  return J;
}

// ---- TE-basis -> helix 6×6 Jacobian -------------------------------------
//   TE    : (c1, c2, c3, theta, phi, invP)  [cartesian c1=x,c2=y,c3=z]
//   helix : (D0, phi, omega, Z0, tanLambda, time)
TMatrixD teJacobian(double theta, double phi, double invP,
                    double B_tesla, bool invPt) {
  TMatrixD J(6, 6);
  J.Zero();
  const double s = std::sin(theta);
  const double c = std::cos(theta);
  if (std::fabs(s) < 1e-9) return J;
  const double inv_s2 = 1.0 / (s * s);
  const double k      = -kOmega * B_tesla;   // omega = k * invP (/ sin theta)
  J(0, 0) = -std::sin(phi) * kCm2Mm;
  J(0, 1) =  std::cos(phi) * kCm2Mm;
  // The reference point is the measured point itself. Moving it along the
  // track moves where the helix passes the reference: Z0 by -tanLambda and
  // phi by +omega per unit transverse length.
  const double omega = invPt ? k * invP : k * invP / s;
  J(1, 0) = omega * std::cos(phi) * kCm2Mm;
  J(1, 1) = omega * std::sin(phi) * kCm2Mm;
  J(3, 0) = -(c / s) * std::cos(phi) * kCm2Mm;
  J(3, 1) = -(c / s) * std::sin(phi) * kCm2Mm;
  J(1, 4) = 1.0;
  // omega = transverse curvature. The omega-row partials must match whichever
  // momentum form fromTrackElement uses: omega = k*invP (no theta dep) when the
  // bank word is 1/p_T, else omega = k*invP/sin(theta).
  if (invPt) {
    J(2, 3) = 0.0;
    J(2, 5) = k;
  } else {
    J(2, 3) = -k * c * invP * inv_s2;   // d/dtheta[k*invP/sin]
    J(2, 5) =  k / s;                    // d/d(invP)
  }
  J(3, 2) = kCm2Mm;
  J(4, 3) = -inv_s2;
  return J;
}

// Cylindrical -> Cartesian for the coordinate triple and its covariance.
// Stored as (R, R*Phi, z), so the point's azimuth is Phi = (R*Phi)/R and
//   x = R cos Phi,  y = R sin Phi,  z = z.
// R is fixed by the surface, so the free coordinates are (R*Phi, z) at TE-basis
// indices 1 and 2. Their derivatives are
//   dx/d(R*Phi) = -sin Phi,  dy/d(R*Phi) = cos Phi,  dz/dz = 1,
// the congruence applied to the covariance below; angles and 1/p pass through.
TMatrixD cylToCartJacobian(double Phi) {
  TMatrixD J(6, 6);
  J.Zero();
  J(0, 1) = -std::sin(Phi);
  J(1, 1) =  std::cos(Phi);
  J(2, 2) = 1.0;
  J(3, 3) = 1.0;   // theta
  J(4, 4) = 1.0;   // phi
  J(5, 5) = 1.0;   // 1/p
  return J;
}

TMatrixDSym covFromArray(const CovMatrix6& cov) {
  TMatrixDSym out(6);
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j <= i; ++j) {
      const double v = static_cast<double>(cov[covOffset(i, j)]);
      out(i, j) = v;
      out(j, i) = v;
    }
  return out;
}

}  // namespace

Helix Helix::fromPerigee(float d0, float z0, float theta, float phi,
                         float invR, const std::array<float, 15>* W) {
  Helix h;
  h.p_.D0    = -d0 * static_cast<float>(kCm2Mm);
  h.p_.phi   = phi;
  // DELPHI's 1/R is signed geometrically: positive is counter-clockwise seen
  // from +z, i.e. opposite to the charge in DELPHI's field (DST content, MAIN:
  // "with the sign opposite to the charge sign"). EDM4hep's omega carries the
  // sign of the charge.
  h.p_.omega = -invR / static_cast<float>(kCm2Mm);
  h.p_.Z0    = z0 * static_cast<float>(kCm2Mm);
  h.p_.time  = 0.f;

  const float s = std::sin(theta);
  if (std::fabs(s) < 1e-6f) {           // degenerate
    h.p_.tanLambda = 0.f;
    h.valid_ = false;
    return h;
  }
  h.p_.tanLambda = std::cos(theta) / s;

  if (!W) return h;                      // params only, no cov

  // Build the 5×5 weight matrix (lower-tri, mirror), invert -> raw cov.
  TMatrixDSym Wm(5);
  {
    int k = 0;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j <= i; ++j) { Wm(i, j) = (*W)[k]; Wm(j, i) = (*W)[k]; ++k; }
  }
  double det = 0.0;
  TMatrixDSym C_dp(Wm);
  C_dp.Invert(&det);
  if (det == 0.0 || !std::isfinite(det)) {
    h.hasCov_ = true;                    // zero cov, but a valid track
    return h;
  }
  const TMatrixD J  = perigeeJacobian(theta);
  const TMatrixD JC (J, TMatrixD::kMult, C_dp);
  const TMatrixD Jt (TMatrixD::kTransposed, J);
  const TMatrixD Cf (JC, TMatrixD::kMult, Jt);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j <= i; ++j)
      h.cov_[covOffset(i, j)] =
        static_cast<float>(0.5 * (Cf(i, j) + Cf(j, i)));
  h.hasCov_ = true;
  return h;
}

Helix Helix::fromTrackElement(double c1, double c2, double c3,
                              double theta, double phi, double invP,
                              bool invPt, bool cylindrical,
                              const CovMatrix6& teCov,
                              double B) {
  Helix h;

  // On a cylinder the stored triple is (R, R*Phi, z); convert it, and its
  // covariance, into the Cartesian basis the helix Jacobian expects. Phi is
  // undefined on the axis, so guard R = 0.
  double x = c1, y = c2;
  CovMatrix6 cov = teCov;
  if (cylindrical) {
    // The stored second coordinate is the arc length R*Phi, so dividing by
    // the radius recovers the azimuth of the point itself. This is the
    // position angle, unrelated to `phi`, which is the track direction.
    const double Phi = (c1 != 0.0) ? c2 / c1 : 0.0;
    x = c1 * std::cos(Phi);
    y = c1 * std::sin(Phi);

    // Congruence C' = Jc * C * Jc^T, same construction as the helix step
    // below: multiply, transpose, multiply, then symmetrise to kill the
    // rounding asymmetry the two products leave behind.
    const TMatrixD    Jc = cylToCartJacobian(Phi);
    const TMatrixDSym Cc = covFromArray(teCov);
    const TMatrixD    JC (Jc, TMatrixD::kMult, Cc);
    const TMatrixD    Jt (TMatrixD::kTransposed, Jc);
    const TMatrixD    Cx (JC, TMatrixD::kMult, Jt);
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j <= i; ++j)
        cov[covOffset(i, j)] = static_cast<float>(0.5 * (Cx(i, j) + Cx(j, i)));
  }

  h.refPoint_ = { static_cast<float>(x  * kCm2Mm),
                  static_cast<float>(y  * kCm2Mm),
                  static_cast<float>(c3 * kCm2Mm) };
  h.p_.D0   = 0.f;
  h.p_.Z0   = 0.f;
  h.p_.time = 0.f;
  h.p_.phi  = static_cast<float>(phi);
  h.hasCov_ = true;

  const double s = std::sin(theta);
  if (std::fabs(s) < 1e-9) {
    h.p_.tanLambda = 0.f;
    h.p_.omega     = 0.f;
    h.valid_ = false;
    return h;                            // cov_ stays zero
  }
  h.p_.tanLambda = static_cast<float>(std::cos(theta) / s);
  // invP is signed geometrically, as DELPHI writes it (positive =
  // counter-clockwise), so omega = -kOmega*B*invP carries the sign of the
  // charge; divided by sin(theta) when the word is 1/p rather than 1/p_T.
  h.p_.omega     = static_cast<float>(
      invPt ? -kOmega * B * invP
            : -kOmega * B * invP / s);

  // Cov push-forward C_helix = J · C_te · J^T.
  const TMatrixD    J  = teJacobian(theta, phi, invP, B, invPt);
  const TMatrixDSym Ct = covFromArray(cov);
  const TMatrixD    JC (J, TMatrixD::kMult, Ct);
  const TMatrixD    Jt (TMatrixD::kTransposed, J);
  const TMatrixD    Cf (JC, TMatrixD::kMult, Jt);
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j <= i; ++j)
      h.cov_[covOffset(i, j)] =
        static_cast<float>(0.5 * (Cf(i, j) + Cf(j, i)));
  return h;
}

Helix Helix::fromHelix(float D0, float phi, float omega,
                       float Z0, float tanLambda) {
  Helix h;
  h.p_ = { D0, phi, omega, Z0, tanLambda, 0.f };
  return h;
}

Momentum Helix::momentum(double B, int charge) const {
  if (charge == 0 || std::fabs(B) < 1e-9 ||
      std::fabs(p_.omega) < 1e-15) {
    return Momentum{0.0, 0.0, 0.0, 0.0};
  }
  const double pT = kOmega * std::fabs(static_cast<double>(charge) * B)
                  / std::fabs(static_cast<double>(p_.omega));
  const double px = pT * std::cos(p_.phi);
  const double py = pT * std::sin(p_.phi);
  const double pz = pT * static_cast<double>(p_.tanLambda);
  return Momentum{px, py, pz, std::sqrt(px*px + py*py + pz*pz)};
}

edm4hep::TrackState Helix::toTrackState(int location) const {
  edm4hep::TrackState ts;
  ts.location       = location;
  ts.D0             = p_.D0;
  ts.phi            = p_.phi;
  ts.omega          = p_.omega;
  ts.Z0             = p_.Z0;
  ts.tanLambda      = p_.tanLambda;
  ts.time           = p_.time;
  ts.referencePoint = refPoint_;
  if (hasCov_) {
    using P = edm4hep::TrackParams;
    static constexpr P kP[6] = {
      P::d0, P::phi, P::omega, P::z0, P::tanLambda, P::time
    };
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b <= a; ++b) {
        const float v = cov_[covOffset(a, b)];
        ts.covMatrix.setValue(v, kP[a], kP[b]);
        if (a != b) ts.covMatrix.setValue(v, kP[b], kP[a]);
      }
  }
  return ts;
}

}  // namespace delphi_edm4hep
