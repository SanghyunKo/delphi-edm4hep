// Truth domain — TruthGenWriter + TruthRecoLinkWriter.
//
// PSCLUJ (filled by PSHLUJ on sDST or PSFLUJ on fDST per bank-presence
// gating) feeds TruthGenWriter. PSCTBL exact tables (IPAST then ISTLU)
// feed TruthRecoLinkWriter, replacing the legacy helix-NN match.

#include "delphi_edm4hep/Truth/Truth.h"

#include "phdst/uxcom.hpp"
#include "phdst/uxlink.hpp"
#include "skelana/pscluj.hpp"
#include "skelana/psctbl.hpp"

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/RecoMCParticleLinkCollection.h>

#include <array>
#include <cmath>
#include <optional>
#include <vector>

// JETSET helpers. Charges via LUCHGE in units of e/3 (exhaustive PDG
// coverage; replaces the 16-entry whitelist of the legacy code).
extern "C" {
  void pshluj_();
  void psfluj_();
  int  luchge_(int* kf);
}

namespace ph = phdst;
namespace sk = skelana;

namespace delphi_edm4hep::truth {

namespace {

inline float charge_from_pdg(int pdg) {
  return static_cast<float>(luchge_(&pdg)) / 3.0f;
}

// Gate the PSCLUJ unpacker on bank presence in the LDTOP chain, invoke
// the sDST or fDST routine, return sk::NP. (full-DST SH chain at
// LDTOP-3, sDST condensed LU at LDTOP-28/-29.)
int unpack_lujets() {
  sk::NP = 0;
  if (ph::LDTOP <= 0) return 0;
  const bool short_sim = (ph::IQ(ph::LDTOP - 2) > 28)
                      && (ph::LQ(ph::LDTOP - 28) != 0)
                      && (ph::LQ(ph::LDTOP - 29) != 0);
  const bool full_sim  = (ph::LQ(ph::LDTOP - 3) != 0);
  if      (short_sim) pshluj_();
  else if (full_sim)  psfluj_();
  return sk::NP;
}

constexpr double kCm2Mm = 10.0;

// Per-event MC truth interaction point (sim-PV), in mm. Handles BOTH sim
// layouts: shortDST (LPVS at LDTOP-28, first sim-PV x/y/z at Q(ip+4..6) cm) and
// fullDST (LDTOP-3 -> LSH -> LST -> LPV, x/y/z at Q(lpv+5..7) cm); same decode
// as delphi-raw-nanoaod fillSimPV. Returns nullopt when no sim banks are present
// (real data / missing bank). Used as the production vertex of generator-only
// entries, which DELSIM gives no vertex of their own.
std::optional<std::array<double, 3>> read_sim_pv_mm() {
  if (ph::LDTOP <= 0) return std::nullopt;
  // shortDST sim layout: LPVS at LDTOP-28, x/y/z at +4/+5/+6 (cm). Guard the
  // down-link index exactly as unpack_lujets() does -- LDTOP-28 is only a valid
  // link when the LDTOP bank has >28 down-links -- so this is safe on real data
  // (and on fullDST, where it falls through to the LSH-chain walk below).
  if (ph::IQ(ph::LDTOP - 2) > 28) {
    const int lpvs = ph::LQ(ph::LDTOP - 28);
    if (lpvs > 0) {
      const int npvs = ph::IQ(lpvs + 1);
      if (npvs >= 1) {
        const int ip = lpvs + 1 + npvs;
        return std::array<double, 3>{ph::Q(ip + 4) * kCm2Mm,
                                     ph::Q(ip + 5) * kCm2Mm,
                                     ph::Q(ip + 6) * kCm2Mm};
      }
    }
  }
  // fullDST sim layout (our DELSIM SDST takes this path): walk the LSH chain at
  // LDTOP-3, first SH with Q(lsh+1)!=0 -> LST at LSH-4 -> LPV at LST+1; x/y/z at
  // LPV+5/+6/+7 (cm). Same decode as raw-nanoaod fillSimPV (PSFLUJ, skelana.car ~L7088).
  for (int lsh = ph::LQ(ph::LDTOP - 3); lsh > 0; lsh = ph::LQ(lsh)) {
    if (std::lround(ph::Q(lsh + 1)) == 0) continue;
    const int lst = ph::LQ(lsh - 4);
    if (lst <= 0) continue;
    const int lpv = ph::LQ(lst + 1);
    if (lpv <= 0) continue;
    return std::array<double, 3>{ph::Q(lpv + 5) * kCm2Mm,
                                 ph::Q(lpv + 6) * kCm2Mm,
                                 ph::Q(lpv + 7) * kCm2Mm};
  }
  return std::nullopt;
}

}  // namespace

// ---------------------------------------------------------------------------
void TruthGenWriter::emit() {
  GenParticleResult result;

  const int nGen = unpack_lujets();
  edm4hep::MCParticleCollection mc;
  result.handles.reserve(static_cast<std::size_t>(nGen));

  // Production vertices come from DELSIM's vertex banks, in cm and already in
  // the frame of the reconstruction (PSFLUJ, skelana.car:7125-7127; PSHLUJ,
  // skelana.car:5390-5392). Entries DELSIM did not track (generator-only) have
  // VP = 0; they are placed at the simulated primary vertex. Gated on nGen>=1
  // so real data (no LUJETS) never reads sim banks.
  const auto spv = nGen >= 1 ? read_sim_pv_mm() : std::nullopt;

  // First pass: per-LU-index handle creation. KP(i,1)=status, KP(i,2)=PDG,
  // PP(i,1..3)=p, PP(i,5)=mass, VP(i,1..3)=production vertex (cm).
  for (int i = 1; i <= nGen; ++i) {
    auto mp = mc.create();
    mp.setPDG(sk::KP(i, 2));
    mp.setGeneratorStatus(static_cast<std::int16_t>(sk::KP(i, 1)));
    mp.setMomentum({sk::PP(i, 1), sk::PP(i, 2), sk::PP(i, 3)});
    mp.setMass(sk::PP(i, 5));
    const bool tracked = sk::VP(i, 1) != 0.f || sk::VP(i, 2) != 0.f ||
                         sk::VP(i, 3) != 0.f;
    if (tracked) {
      mp.setVertex({static_cast<double>(sk::VP(i, 1)) * kCm2Mm,
                    static_cast<double>(sk::VP(i, 2)) * kCm2Mm,
                    static_cast<double>(sk::VP(i, 3)) * kCm2Mm});
    } else if (spv) {
      mp.setVertex({(*spv)[0], (*spv)[1], (*spv)[2]});
    }
    mp.setCharge(charge_from_pdg(sk::KP(i, 2)));
    result.handles.push_back(mp);
  }

  // Second pass: parent/daughter graph (parent edges only; daughters
  // are recoverable by reverse traversal).
  // Endpoint (decay vertex): a decay's real daughters sit at the (possibly
  // displaced) decay vertex, but LUJETS documentation-copy daughters (status 21)
  // sit at the PARENT's production vertex. Pick the MOST-displaced daughter so a
  // coincident doc-copy can't clobber a real displaced vertex -- the previous
  // "last-daughter-wins" silently dropped ~40% of displaced B/D endpoints. A
  // particle whose only daughters are coincident gets a zero-flight endpoint at
  // its own vertex (correct); stable particles (no daughters) keep the default.
  std::vector<double> best_disp2(static_cast<std::size_t>(nGen), -1.0);
  for (int i = 1; i <= nGen; ++i) {
    const int parent_lu = sk::KP(i, 3);
    if (parent_lu >= 1 && parent_lu <= nGen) {
      const std::size_t child = static_cast<std::size_t>(i - 1);
      const std::size_t pidx  = static_cast<std::size_t>(parent_lu - 1);
      result.handles[child].addToParents(result.handles[pidx]);
      result.handles[pidx].addToDaughters(result.handles[child]);
      const auto pv = result.handles[pidx].getVertex();
      const auto cv = result.handles[child].getVertex();
      const double dx = cv.x - pv.x, dy = cv.y - pv.y, dz = cv.z - pv.z;
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 > best_disp2[pidx]) {
        result.handles[pidx].setEndpoint(cv);
        best_disp2[pidx] = d2;
      }
    }
  }

  put(std::move(mc), "LUJ", "GenParticles", Provenance::Derived);
  ctx_.gen_truth = std::move(result);
}

// ---------------------------------------------------------------------------
void TruthRecoLinkWriter::emit() {
  edm4hep::RecoMCParticleLinkCollection links;

  // Need both upstream writers' outputs.
  if (!ctx_.gen_truth || !ctx_.tracking) {
    put(std::move(links), "TBL", "RecoToGen", Provenance::Transcribed);   // emit empty + return
    return;
  }
  const auto& gen      = *ctx_.gen_truth;
  const auto& tracking = *ctx_.tracking;

  // PSCTBL.NPA = # PA particles = VECP entries. For each VECP index j:
  //   ist = IPAST(j)                 (ST index; 0 if no MC ancestor)
  //   ilu = ISTLU(ist)               (LU index)
  //   particle_idx = vecp_to_particle[j]   (-1 if Tracking dropped it)
  const int nPA   = sk::NPA();
  const int nGen  = static_cast<int>(gen.handles.size());
  const auto& v2p = tracking.vecp_to_particle;

  for (int j = 1; j <= nPA; ++j) {
    if (j >= static_cast<int>(v2p.size())) break;
    const int particle_idx = v2p[j];
    if (particle_idx < 0) continue;
    const int ist = sk::IPAST(j);
    if (ist <= 0) continue;
    const int ilu = sk::ISTLU(ist);
    if (ilu <= 0 || ilu > nGen) continue;

    auto link = links.create();
    link.setFrom(tracking.particle_handles[particle_idx]);
    link.setTo  (gen.handles[ilu - 1]);
    link.setWeight(1.0f);
  }

  put(std::move(links), "TBL", "RecoToGen", Provenance::Transcribed);
}

}  // namespace delphi_edm4hep::truth
