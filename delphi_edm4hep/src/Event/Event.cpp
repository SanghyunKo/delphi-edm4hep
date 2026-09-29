// EventWriter — implementation.
//
// Reads per-event scalars from the SKELANA / PHDST commons populated by PSBEG
// and writes them into the Frame as named parameters under "<source>_EVT_*".
// No bank walking — pure common-block reads.
//
// The parameters fall into three groups by origin, declared separately below:
// values copied from the DST pilot record, values SKELANA counts or sums for
// itself, and the beam spot, which comes from an external database rather
// than from the DST at all.

#include "delphi_edm4hep/Event/Event.h"
#include "delphi_edm4hep/internal/BeamSpotStatus.h"

#include "phdst/functions.hpp"  // IPHPIC
#include "phdst/phciii.hpp"
#include "phdst/uxcom.hpp"      // IQ
#include "phdst/uxlink.hpp"     // LDTOP
#include "skelana/pscbsp.hpp"
#include "skelana/pscevt.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

// PHGEN BPILOT subroutine (current event's solenoid B in Tesla,
// curvature-to-momentum factor in GeV/cm). No header wrapper.
extern "C" void bpilot_(float* btesla, float* bgevcm);

// DSTQID (dstana) returns the DST processing identifier as "YYLN": two-digit
// year, DELANA processing letter, short/mini DST number. gfortran passes the
// character length as a hidden trailing argument.
extern "C" void dstqid_(char* tag, std::size_t len);

namespace ph = phdst;
namespace sk = skelana;

namespace {
std::string processingTag() {
  char buf[4] = {};
  dstqid_(buf, sizeof(buf));
  std::string tag(buf, sizeof(buf));
  tag.erase(tag.find_last_not_of(' ') + 1);
  return tag;
}

// Refuse the event, naming it. A writer that throws aborts the job through the
// harness's callback guard, which marks the partial output unpublishable.
[[noreturn]] void refuse(const std::string& what) {
  throw std::runtime_error(what + " (run " + std::to_string(ph::IIIRUN)
                           + ", event " + std::to_string(ph::IIIEVT) + ')');
}

// The beam energy has a failure mode SKELANA does not report: it substitutes a
// value and carries on, and the substitute reaches the output looking like a
// measurement. (The beam spot's is handled in internal/BeamSpotStatus.h.)
void checkEventIsUsable() {
  // Without the DANA pilot blocklet the centre-of-mass energy is not read but
  // set to 91.250 GeV (skelana.car:2743-2748). That substitution cannot be
  // recognised from the energy itself — 91.250 is also the true stored value
  // for LEP1 simulation — so test for the blocklet.
  if (ph::IPHPIC("DANA", 0) <= 0) {
    refuse("DANA pilot blocklet absent, so the centre-of-mass energy would be "
           "SKELANA's 91.250 GeV fallback rather than a stored value");
  }

  // EBEAM = ECMAS/2 is both a track cut and a denominator in SKELANA; at zero
  // it rejects every charged track.
  if (sk::ECMAS <= 0.f) {
    refuse("centre-of-mass energy is " + std::to_string(sk::ECMAS) + " GeV");
  }
}
}  // namespace

namespace delphi_edm4hep::event {

void EventWriter::emit() {
  checkEventIsUsable();

  constexpr float kCm2Mm = 10.0f;

  // Pilot-record words, written when the DST was produced.
  auto stored = parameters("EVT", Provenance::Transcribed);
  stored("runNumber",   ph::IIIRUN);
  stored("eventNumber", ph::IIIEVT);
  stored("fileSeq",     ph::IIFILE);
  stored("date",        ph::IIIDAT);   // yymmdd
  stored("time",        ph::IIITIM);   // hhmmss
  stored("fillNumber",  ph::IIFILL);
  stored("experiment",  ph::IIIEXP);
  stored("dstVersion",  sk::ISVER);

  // Era identifiers. The processing tag selects the calibration SKELANA uses —
  // its first two characters pick the RICH refractive index — and the PXDST
  // version selects between algorithm generations.
  stored("dstProcessingTag", processingTag());
  stored("pxdstVersion",     ph::LDTOP > 0 ? ph::IQ(ph::LDTOP + 3) : 0);

  // Centre-of-mass energy from the DANA pilot blocklet. SKELANA substitutes
  // 91.250 GeV when that blocklet is absent, which at LEP2 would be wrong by
  // more than a factor of two; check BeamSpotErrorCode and the energy itself
  // before relying on it.
  stored("ECMS", sk::ECMAS);

  // Per-event solenoid field and the curvature-to-momentum factor, both from
  // the pilot record. pT [GeV/c] = BFieldGevPerCm / |omega_DELPHI [1/cm]|.
  float btesla = 0.f, bgevcm = 0.f;
  bpilot_(&btesla, &bgevcm);
  stored("BField",         btesla);
  stored("BFieldGevPerCm", bgevcm);

  // Multiplicities and summed energies counted by SKELANA over the PA chain
  // under cuts hard-coded in PSHEVT — a momentum, track-length, impact-
  // parameter and polar-angle selection independent of IFLCUT and of the
  // LVLOCK track selection. Two events with the same nCharged can therefore
  // disagree with a count made from MAIN_Particles and LVLOCK.
  auto counted = parameters("EVT", Provenance::Derived);
  counted("hadronicTagTeam4", sk::IHAD4);
  counted("nChargedTeam4",    sk::NCTR4);
  counted("nCharged",         sk::NCTRK);
  counted("nNeutral",         sk::NNTRK);
  counted("EChargedTotal",    sk::ECHAR);
  counted("ENeutralEM",       sk::EMNEU);
  counted("ENeutralHad",      sk::EHNEU);

  // Beam spot from the per-processing database under $DELPHI_DAT, selected by
  // the DSTQID tag and looked up per run and file — not from the DST. For
  // simulation it is instead the generated interaction point plus a Gaussian
  // smear, so it follows the event rather than describing a beam. Without a
  // usable entry (BeamSpotErrorCode 2) position and widths are NaN.
  const bool usable = beamspot::positionUsable(sk::IERRBS);
  auto bs = [&](float cm) {
    return usable ? cm * kCm2Mm : std::numeric_limits<float>::quiet_NaN();
  };
  auto beamspot = parameters("EVT", Provenance::Derived);
  beamspot("BeamSpotX",      bs(sk::XYZBS(1)));
  beamspot("BeamSpotY",      bs(sk::XYZBS(2)));
  beamspot("BeamSpotZ",      bs(sk::XYZBS(3)));
  beamspot("BeamSpotSigmaX", bs(sk::DXYZBS(1)));
  beamspot("BeamSpotSigmaY", bs(sk::DXYZBS(2)));
  beamspot("BeamSpotSigmaZ", bs(sk::DXYZBS(3)));
  beamspot("BeamSpotErrorCode", sk::IERRBS);
}

}  // namespace delphi_edm4hep::event
