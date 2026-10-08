// EmcaWriter — implementation. Runs in both passes.
//
// PA.EMCA bank layout (per shower, header at show_ofs):
//   show_ofs+1   len  (length in words of this shower record)
//   show_ofs+2   idet (9 = HPC barrel, 26 = FEMC endcap)
//   show_ofs+3   E    (GeV)
//   show_ofs+4..+6  X, Y, Z (cm) — shower-level centroid
//   show_ofs+10  nclu (per-pad or per-layer count)
//   show_ofs+11  nwdcl (words per cluster: 3 = HPC, 2 = FEMC)
//   show_ofs+11+k*nwdcl..   k-th per-cluster record
//
// HPC clusters use the 3-word PXHGET packing decoded by
// internal::hpc::padDecode → (E, layer, X, Y, Z, sigma_z).
// FEMC layers are (energy, packed=layer*1000+nhits).

#include "delphi_edm4hep/Calorimeter/Emca.h"

#include "delphi_edm4hep/internal/HpcPadDecoder.h"
#include "delphi_edm4hep/internal/PaWalk.h"

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>

#include <cmath>
#include <cstdint>
#include <limits>

namespace ph = phdst;

namespace delphi_edm4hep::emca {

namespace {

constexpr float kNotMeasured = std::numeric_limits<float>::quiet_NaN();

constexpr float kCm2Mm = 10.f;
constexpr int   kIdetHPC  = 9;
constexpr int   kIdetFEMC = 26;
constexpr int   kNwdHPC   = 3;
constexpr int   kNwdFEMC  = 2;

}  // namespace

void EmcaWriter::emit()
{
  edm4hep::CalorimeterHitCollection hpc_col;
  edm4hep::CalorimeterHitCollection femc_col;

  Output out;

  auto record = [&](int paIdx, int ns, const edm4hep::CalorimeterHit& hit) {
    if (paIdx < 0) return;
    if (paIdx >= static_cast<int>(out.pa_shower_hits.size())) {
      out.pa_shower_hits.resize(static_cast<std::size_t>(paIdx) + 1);
    }
    auto& showers = out.pa_shower_hits[static_cast<std::size_t>(paIdx)];
    if (ns >= static_cast<int>(showers.size())) {
      showers.resize(static_cast<std::size_t>(ns) + 1);
    }
    showers[static_cast<std::size_t>(ns)].push_back(hit);
  };

  pawalk::forEachPA([&](int lpa, int paIdx) {
    const int lemca = pawalk::lphpa("EMCA", lpa);
    if (lemca <= 0) return;
    const int nsh = static_cast<int>(std::lround(ph::Q(lemca + 2)));
    if (nsh <= 0) return;

    int show_ofs = lemca + 2;   // shower 1 starts here
    for (int ns = 0; ns < nsh; ++ns) {
      const int len   = static_cast<int>(std::lround(ph::Q(show_ofs + 1)));
      const int idet  = static_cast<int>(std::lround(ph::Q(show_ofs + 2)));
      const int nclu  = static_cast<int>(std::lround(ph::Q(show_ofs + 10)));
      const int nwdcl = static_cast<int>(std::lround(ph::Q(show_ofs + 11)));

      if (idet == kIdetHPC && nclu > 0 && nclu < 200 && nwdcl == kNwdHPC) {
        for (int nc = 0; nc < nclu; ++nc) {
          const int base = show_ofs + 11 + kNwdHPC * nc;
          const std::int32_t i1 = ph::IQ(base + 1);
          const std::int32_t i2 = ph::IQ(base + 2);
          const std::int32_t i3 = ph::IQ(base + 3);
          const float w1 = ph::Q(base + 1);
          const float w2 = ph::Q(base + 2);
          const float w3 = ph::Q(base + 3);
          float pe, pxc, pyc, pzc, psz; int player;
          hpc::padDecode(i1, i2, i3, w1, w2, w3,
                         pe, player, pxc, pyc, pzc, psz);
          if (pe <= 0.f) continue;
          auto hit = hpc_col.create();
          hit.setEnergy(pe);
          hit.setEnergyError(kNotMeasured);   // PXHGET gives no energy error
          hit.setTime(0.f);          // not available in PXHGET
          hit.setType(player);       // layer 1..10
          hit.setPosition({
            pxc * kCm2Mm,
            pyc * kCm2Mm,
            pzc * kCm2Mm,
          });
          record(paIdx, ns, hit);
        }
      } else if (idet == kIdetFEMC && nclu > 0 && nclu <= 10 && nwdcl == kNwdFEMC) {
        // FEMC: per-layer (energy, packed = layer*1000 + nhits).
        // Position is the shower-level centroid (FEMC bank doesn't
        // store per-layer X/Y); type carries the layer index, cellID
        // carries nhits.
        const float showX = ph::Q(show_ofs + 4);
        const float showY = ph::Q(show_ofs + 5);
        const float showZ = ph::Q(show_ofs + 6);
        for (int nc = 0; nc < nclu; ++nc) {
          const int base = show_ofs + 11 + kNwdFEMC * nc;
          const float layer_e = ph::Q(base + 1);
          const int   packed  = static_cast<int>(std::lround(ph::Q(base + 2)));
          const int   layer   = packed / 1000;
          const int   nhits   = packed % 1000;
          auto hit = femc_col.create();
          hit.setEnergy(layer_e);
          hit.setEnergyError(0.f);
          hit.setTime(0.f);
          hit.setType(layer);
          hit.setCellID(static_cast<std::uint64_t>(nhits));
          hit.setPosition({
            showX * kCm2Mm,
            showY * kCm2Mm,
            showZ * kCm2Mm,
          });
          record(paIdx, ns, hit);
        }
      }

      if (len <= 0) break;
      show_ofs += len;
    }
  });

  ctx_.emca = std::move(out);
  put(std::move(hpc_col),  "EMCA", "HPCClusters", Provenance::Transcribed);
  put(std::move(femc_col), "EMCA", "FEMCLayers", Provenance::Transcribed);
}

}  // namespace delphi_edm4hep::emca
