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

/// @file   EnergyCalibSpec.cxx
/// @brief  ZDC reconstruction
/// @author pietro.cortese@cern.ch

#include <iostream>
#include <vector>
#include <string>
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "Framework/Logger.h"
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsZDC/BCData.h"
#include "DataFormatsZDC/ChannelData.h"
#include "DataFormatsZDC/OrbitData.h"
#include "DataFormatsZDC/RecEvent.h"
#include "ZDCBase/ModuleConfig.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "ZDCReconstruction/RecoConfigZDC.h"
#include "ZDCReconstruction/ZDCEnergyParam.h"
#include "ZDCCalib/InterCalib.h"
#include "ZDCWorkflow/EnergyCalibSpec.h"

using namespace o2::framework;

namespace o2
{
namespace zdc
{

EnergyCalibSpec::EnergyCalibSpec()
{
  mTimer.Stop();
  mTimer.Reset();
}

EnergyCalibSpec::EnergyCalibSpec(const int verbosity)
  : mVerbosity(verbosity)
{
  mTimer.Stop();
  mTimer.Reset();
}

void EnergyCalibSpec::init(o2::framework::InitContext& ic)
{
  mVerbosity = ic.options().get<int>("verbosity-level");
}

void EnergyCalibSpec::updateTimeDependentParams(ProcessingContext& pc)
{
  // we call these methods just to trigger finaliseCCDB callback
  pc.inputs().get<o2::zdc::ZDCEnergyParam*>("energycalib");
}

void EnergyCalibSpec::run(ProcessingContext& pc)
{
  updateTimeDependentParams(pc);
  if (!mInitialized) {
    mInitialized = true;
    std::string loadedConfFiles = "Loaded ZDC configuration files:";

    // Energy calibration
    auto energyParam = pc.inputs().get<o2::zdc::ZDCEnergyParam*>("energycalib");
    if (!energyParam) {
      LOG(fatal) << "Missing ZDCEnergyParam calibration object";
      return;
    } else {
      loadedConfFiles += " ZDCEnergyParam";
      if (mVerbosity > DbgZero) {
        LOG(info) << "Loaded Energy calibration ZDCEnergyParam";
        energyParam->print();
      }
    }

    LOG(info) << loadedConfFiles;
    mTimer.CpuTime();
    mTimer.Start(false);
  }

  auto bcrec = pc.inputs().get<gsl::span<o2::zdc::BCRecData>>("bcrec");
  auto energy = pc.inputs().get<gsl::span<o2::zdc::ZDCEnergy>>("energy");
  auto tdc = pc.inputs().get<gsl::span<o2::zdc::ZDCTDCData>>("tdc");
  auto info = pc.inputs().get<gsl::span<uint16_t>>("info");
  InterCalib work;
  work.init();
  work.process(bcrec, energy, tdc, info);
}

void EnergyCalibSpec::endOfStream(EndOfStreamContext& ec)
{
  mTimer.Stop();
  LOGF(info, "ZDC TDC calibration total timing: Cpu: %.3e Real: %.3e s in %d slots",
       mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

framework::DataProcessorSpec getEnergyCalibSpec(const int verbosity = 0)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("bcrec", "ZDC", "BCREC", 0, Lifetime::Timeframe);
  inputs.emplace_back("energy", "ZDC", "ENERGY", 0, Lifetime::Timeframe);
  inputs.emplace_back("tdc", "ZDC", "TDCDATA", 0, Lifetime::Timeframe);
  inputs.emplace_back("info", "ZDC", "INFO", 0, Lifetime::Timeframe);
  inputs.emplace_back("energycalib", "ZDC", "ENERGYCALIB", 0, Lifetime::Condition, o2::framework::ccdbParamSpec(fmt::format("{}", o2::zdc::CCDBPathEnergyCalib.data())));

  std::vector<OutputSpec> outputs;

  return DataProcessorSpec{
    "zdc-energy-calib",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<EnergyCalibSpec>(verbosity)},
    o2::framework::Options{{"verbosity-level", o2::framework::VariantType::Int, 0, {"Verbosity level"}}}};
}

} // namespace zdc
} // namespace o2
