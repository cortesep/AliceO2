// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
//file RawReaderZDC.h class  for RAW data reading

#ifndef ALICEO2_FIT_RAWREADERZDC_H_
#define ALICEO2_FIT_RAWREADERZDC_H_
#include <iostream>
#include <vector>
#include <Rtypes.h>
#include "ZDCRaw/RawReaderZDCBase.h"
#include "/DataFormatsZDC/RawEventData.h"
#include "DataFormatsZDC/ChannelData.h"
#include "DataFormatsZDC/BCData.h"
#include "DataFormatsZDC/PedestalData.h"

#include "Framework/ProcessingContext.h"
#include "Framework/DataAllocator.h"
#include "Framework/OutputSpec.h"
#include <gsl/span>

namespace o2
{
namespace zdc
{
class RawReaderZDC : public RawReaderZDCBaseNorm
{
 public:
  RawReaderZDC(bool dumpData) : mDumpData(dumpData) {}
  RawReaderZDC(const RawReaderZDC&) = default;

  RawReaderZDC() = default;
  ~RawReaderZDC() = default;
  void clear()
  {
    mVecDigits.clear();
    mVecChannelData.clear();
  }
  void accumulateDigits()
  {
    getDigits(mVecDigits, mVecChannelData);
    LOG(INFO) << "Number of Digits: " << mVecDigits.size();
    LOG(INFO) << "Number of ChannelData: " << mVecChannelData.size();
    if (mDumpData) {
      DigitBlockZDC::print(mVecDigits, mVecChannelData);
    }
  }
  static void prepareOutputSpec(std::vector<o2::framework::OutputSpec>& outputSpec)
  {
    outputSpec.emplace_back("ZDC", "DIGITSBC", 0, Lifetime::Timeframe);
    outputSpec.emplace_back("ZDC", "DIGITSCH", 0, Lifetime::Timeframe);
    outputSpec.emplace_back("ZDC", "DIGITSPD", 0, Lifetime::Timeframe);
  }
  void makeSnapshot(o2::framework::ProcessingContext& pc)
  {
    pc.outputs().snapshot(o2::framework::Output{o2::header::gDataOriginZDC, "DIGITSBC", 0, o2::framework::Lifetime::Timeframe}, mDigitsBC);
    pc.outputs().snapshot(o2::framework::Output{o2::header::gDataOriginZDC, "DIGITSCH", 0, o2::framework::Lifetime::Timeframe}, mDigitsCh);
    pc.outputs().snapshot(o2::framework::Output{o2::header::gDataOriginZDC, "DIGITSCH", 0, o2::framework::Lifetime::Timeframe}, mPedestalData);
  }
  bool mDumpData;
  std::vector<o2::zdc::BCData> mDigitsBC;
  std::vector<o2::zdc::ChannelData> mDigitsCh;
  std::vector<o2::zdc::PedestalData> mPedestalData;
};
} // namespace zdc
} // namespace o2

#endif
