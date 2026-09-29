// BeamSpotStatus.h -- whether SKELANA's beam spot describes the current event.
//
// PSBEAM fills PSCBSP through VDBSPT from the per-processing database under
// $DELPHI_DAT (e.g. 95C2.DB), looked up by run and file. IERRBS reports the
// outcome (vdbeam.car:582-651):
//   0  the entry for this run and file;
//   1  no entry for this file; the preceding entry is used instead;
//   2  no usable entry. The position is then whatever row the neighbour search
//      stopped on -- another run's beam spot, or an out-of-range read that can
//      give (0,0,0) -- so it describes nothing.
// If the database file cannot be opened, VDBSPT returns on the first event
// without setting IERRBS (vdbeam.car:1167-1171) and reports 2 afterwards, so a
// missing file would otherwise pass as a run of code-2 events.

#pragma once

namespace delphi_edm4hep::beamspot {

enum ErrorCode : int { kOwnEntry = 0, kPrecedingEntry = 1, kNoEntry = 2 };

// True when the position and widths describe a measured beam spot and may be
// written; otherwise they are written as NaN.
constexpr bool positionUsable(int errorCode) {
  return errorCode == kOwnEntry || errorCode == kPrecedingEntry;
}

// Refuse the event when the database was never read. A successful lookup
// always has non-zero widths, because PSBEAM folds the beam size into them
// (SETBS, vdbeam.car:656-672), so code 0 with all widths zero can only be the
// unopened-file path. Both passes call this first thing for every event.
void requireDatabase();

}  // namespace delphi_edm4hep::beamspot
