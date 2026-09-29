// BeamSpotStatus — implementation of the database check.

#include "delphi_edm4hep/internal/BeamSpotStatus.h"

#include "phdst/phciii.hpp"
#include "skelana/pscbsp.hpp"

#include <stdexcept>
#include <string>

namespace delphi_edm4hep::beamspot {

void requireDatabase() {
  namespace sk = skelana;
  if (sk::IERRBS == kOwnEntry && sk::DXYZBS(1) == 0.f &&
      sk::DXYZBS(2) == 0.f && sk::DXYZBS(3) == 0.f) {
    throw std::runtime_error(
        "beam-spot database was not read: error code 0 with zero widths (run " +
        std::to_string(phdst::IIIRUN) + ", event " +
        std::to_string(phdst::IIIEVT) + ')');
  }
}

}  // namespace delphi_edm4hep::beamspot
