// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DataFormatsZDC/BCData.h"
#include "DataFormatsZDC/ChannelData.h"
#include <bitset>

using namespace o2::zdc;

void BCData::print(uint32_t triggerMask) const
{
  printf("Orbit %9u bc %4u nch %2d pos %d\n", ir.orbit, ir.bc, ref.getEntries(), ref.getFirstEntry());
  printf("Read:");
  for (int ic = 0; ic < NDigiChannels; ic++) {
    if (ic % NChPerModule == 0) {
      if (ic == 0) {
        printf(" %d[", ic / NChPerModule);
      } else {
        printf("] %d[", ic / NChPerModule);
      }
    }
    if (channels & (0x1 << ic)) {
      printf("R");
    } else {
      printf(" ");
    }
  }
  printf("]\nHits:");
  for (int ic = 0; ic < NDigiChannels; ic++) {
    if (ic % NChPerModule == 0) {
      if (ic == 0) {
        printf(" %d[", ic / NChPerModule);
      } else {
        printf("] %d[", ic / NChPerModule);
      }
    }
    bool is_hit = triggers & (0x1 << ic);
    bool is_trig = triggerMask & (0x1 << ic);
    if (is_trig) {
      if (is_hit) {
        printf("T");
      } else {
        printf(".");
      }
    } else {
      if (is_hit) {
        printf("H");
      } else {
        printf(" ");
      }
    }
  }
  printf("]\nAUTO:");
  for (int i = 0; i < NModules; i++) {
    std::bitset<10> bb(moduleTriggers[i]);
    printf(" %d %s%s%s%s%s", i, bb[8] ? "3" : "-", bb[7] ? "2" : "-", bb[6] ? "1" : "-", bb[5] ? "0" : "-", bb[4] ? "M" : "-");
  }
  printf("\nALIT:");
  for (int i = 0; i < NModules; i++) {
    std::bitset<10> bb(moduleTriggers[i]);
    printf(" %d %s%s%s%s ", i, bb[3] ? "3" : "-", bb[2] ? "2" : "-", bb[1] ? "1" : "-", bb[0] ? "0" : "-");
  }
  printf("\n");
}

gsl::span<const ChannelData> BCData::getBunchChannelData(const gsl::span<const ChannelData> tfdata) const
{
  // extract the span of channel data for this bunch from the whole TF data
  return ref.getEntries() ? gsl::span<const ChannelData>(&tfdata[ref.getFirstEntry()], ref.getEntries()) : gsl::span<const ChannelData>();
}
