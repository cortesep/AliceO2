// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// DataFormats/Detectors/ZDC/src/RawEventData.cxx
#include "DataFormatsZDC/RawEventData.h"

using namespace o2::zdc;

//ClassImp(EventData);

//______________________________________________________________________________
void EventChData::reset()
{
  static constexpr int payloadSize = NWPerGBTW * sizeof(UInt_t);
  memset((void*)&w[0][0], 0, payloadSize);
}

//______________________________________________________________________________
void EventData::print() const
{
  for (Int_t im = 0; im < o2::zdc::NModules; im++) {
    for (Int_t ic = 0; ic < o2::zdc::NChPerModule; ic++) {
      if (data[im][ic].f.fixed_0 == Id_w0 && data[im][ic].f.fixed_1 == Id_w1 && data[im][ic].f.fixed_2 == Id_w2) {
        // Not empty event
        auto f = data[im][ic].f;
        // Word 0
        printf("%04x %08x %08x ", data[im][ic].w[0][2], data[im][ic].w[0][1], data[im][ic].w[0][0]);
        printf("orbit %-9u bc %-4u hits %-4u offset %+6i Board %2u Ch %1u\n", f.orbit, f.bc, f.hits, f.offset, f.ch, f.board);
        // Word 1
        printf("%04x %08x %08x ", data[im][ic].w[1][2], data[im][ic].w[1][1], data[im][ic].w[1][0]);
        printf("     %s %s %s %s 0-5 ", f.Alice_0 ? "A0" : "  ", f.Alice_1 ? "A1" : "  ", f.Alice_2 ? "A2" : "  ", f.Alice_3 ? "A3" : "  ");
        printf(" %5d %5d %5d %5d %5d %5d EC=%u\n", f.s00, f.s01, f.s02, f.s03, f.s04, f.s05, f.error);
        // Word 2
        printf("%04x %08x %08x ", data[im][ic].w[2][2], data[im][ic].w[2][1], data[im][ic].w[2][0]);
        printf("%s %s %s %s %s %s 6-b ", f.Hit ? "H" : " ", f.Auto_m ? "TM" : "  ", f.Auto_0 ? "T0" : "  ", f.Auto_1 ? "T1" : "  ", f.Auto_2 ? "T2" : "  ", f.Auto_3 ? "T3" : "  ");
        printf(" %5d %5d %5d %5d %5d %5d\n", f.s06, f.s07, f.s08, f.s09, f.s10, f.s11);
      } else if (data[im][ic].f.fixed_0 == 0 && data[im][ic].f.fixed_1 == 0 && data[im][ic].f.fixed_2 == 0) {
        // Empty channel
      } else {
        // Wrong data format. Insert warning.
      }
    }
  }
}

//______________________________________________________________________________
void EventData::reset()
{
  // Clear GBT words
  static constexpr int payloadSize = NModules * NChPerModule * NWPerGBTW * sizeof(UInt_t);
  memset((void*)&data[0][0], 0, payloadSize);
  // Clear decoded samples
  static constexpr int sampleSize = NModules * NChPerModule * NTimeBinsPerBC * sizeof(Short_t);
  memset((void*)&s[0][0][0], 0, sampleSize);
}

//______________________________________________________________________________
void EventData::decode()
{
  // Convert raw samples to signed integers
  UShort_t us[NTimeBinsPerBC];
  for (Int_t im = 0; im < o2::zdc::NModules; im++) {
    for (Int_t ic = 0; ic < o2::zdc::NChPerModule; ic++) {
      if (data[im][ic].f.fixed_0 == Id_w0 && data[im][ic].f.fixed_1 == Id_w1 && data[im][ic].f.fixed_2 == Id_w2) {
        auto& chd = data[im][ic];
        us[0] = chd.f.s00;
        us[1] = chd.f.s01;
        us[2] = chd.f.s02;
        us[3] = chd.f.s03;
        us[4] = chd.f.s04;
        us[5] = chd.f.s05;
        us[6] = chd.f.s06;
        us[7] = chd.f.s07;
        us[8] = chd.f.s08;
        us[9] = chd.f.s09;
        us[10] = chd.f.s10;
        us[11] = chd.f.s11;
        for (Int_t i = 0; i < NTimeBinsPerBC; i++) {
          if (us[i] > ADCMax) {
            s[im][ic][i] = us[i] - ADCRange;
          } else {
            s[im][ic][i] = us[i];
          }
        }
      }
    }
  }
}

//______________________________________________________________________________
void EventData::decodeCh(Int_t im, Int_t ic)
{
  // Convert raw samples to signed integers
  UShort_t us[NTimeBinsPerBC];
  if (im >= 0 && im < o2::zdc::NModules && ic >= 0 && ic < o2::zdc::NChPerModule) {
    if (data[im][ic].f.fixed_0 == Id_w0 && data[im][ic].f.fixed_1 == Id_w1 && data[im][ic].f.fixed_2 == Id_w2) {
      auto& chd = data[im][ic];
      us[0] = chd.f.s00;
      us[1] = chd.f.s01;
      us[2] = chd.f.s02;
      us[3] = chd.f.s03;
      us[4] = chd.f.s04;
      us[5] = chd.f.s05;
      us[6] = chd.f.s06;
      us[7] = chd.f.s07;
      us[8] = chd.f.s08;
      us[9] = chd.f.s09;
      us[10] = chd.f.s10;
      us[11] = chd.f.s11;
      for (Int_t i = 0; i < NTimeBinsPerBC; i++) {
        if (us[i] > ADCMax) {
          s[im][ic][i] = us[i] - ADCRange;
        } else {
          s[im][ic][i] = us[i];
        }
      }
    }
  }
}
