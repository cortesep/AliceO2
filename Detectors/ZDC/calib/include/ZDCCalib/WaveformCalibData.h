// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ZDC_WAVEFORMCALIB_DATA_H
#define ZDC_WAVEFORMCALIB_DATA_H

#include "ZDCBase/Constants.h"
#include "ZDCCalib/WaveformCalibConfig.h"
#include <array>

/// \file WaveformCalibData.h
/// \brief Waveform calibration intermediate data
/// \author pietro.cortese@cern.ch

namespace o2
{
namespace zdc
{

struct WaveformCalibChData {
  static constexpr int NBT = WaveformCalibConfig::NBT;
  static constexpr int NW = WaveformCalibConfig::NW;

  int mFirstValid = 0;
  int mLastValid = 0;
  uint32_t mEntries = 0;
  std::array<float, NW> mData = {0};

  WaveformCalibChData& operator+=(const WaveformCalibChData& other);
  int getEntries() const;
  int getFirstValid() const;
  int getLastValid() const;
  void clear();
  void setN(int n);
  ClassDefNV(WaveformCalibChData, 1);
};

struct WaveformCalibData {
  static constexpr int NBB = WaveformCalibConfig::NBB;
  static constexpr int NBA = WaveformCalibConfig::NBA;
  static constexpr int NBT = WaveformCalibConfig::NBT;
  static constexpr int NW = WaveformCalibConfig::NW;

  uint64_t mCTimeBeg = 0; /// Time of processed time frame
  uint64_t mCTimeEnd = 0; /// Time of processed time frame
  int mN = 0;             /// Number of bunches in waveform
  int mPeak = 0;          /// Peak position

  std::array<WaveformCalibChData, NChannels> mWave;
  WaveformCalibData& operator+=(const WaveformCalibData& other);
  inline void setFirstValid(int isig, int ipos)
  {
    if (ipos > mWave[isig].mFirstValid) {
      mWave[isig].mFirstValid = ipos;
    }
  }
  inline void setLastValid(int isig, int ipos)
  {
    if (ipos < mWave[isig].mLastValid) {
      mWave[isig].mLastValid = ipos;
    }
  }
  inline void addEntry(int isig)
  {
    mWave[isig].mEntries++;
  }
  int getEntries(int is) const;
  int getFirstValid(int is) const;
  int getLastValid(int is) const;
  void print() const;
  void clear();
  void setCreationTime(uint64_t ctime);
  void setN(int n);
  int write(const std::string fn);
  ClassDefNV(WaveformCalibData, 1);
};

} // namespace zdc
} // namespace o2

#endif
