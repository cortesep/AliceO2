// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_ZDC_RECOPARAM_H_
#define O2_ZDC_RECOPARAM_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "ZDCBase/Constants.h"

namespace o2
{
namespace zdc
{
// parameters of ZDC reconstruction

struct ZDCRecoParam : public o2::conf::ConfigurableParamHelper<ZDCRecoParam> {
  Int_t tsh[NTDCChannels]={4,4,4,4,4,4,4,4,4,4}; // Trigger shift
  Int_t tth[NTDCChannels]={8,8,8,8,8,8,8,8,8,8}; // Trigger threshold
  Int_t tmod[NTDCChannels]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; // Position of TDC channel in raw data
  Int_t tch[NTDCChannels]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; // Position of TDC channel in raw data
  O2ParamDef(ZDCRecoParam, "ZDCRecoParam");
};
} // namespace zdc
} // namespace o2

#endif
