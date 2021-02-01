// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "Framework/WorkflowSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataSpecUtils.h"
#include "Framework/ControlService.h"
#include "Framework/Logger.h"
#include "Framework/ConfigParamSpec.h"
#include "DPLUtils/DPLRawParser.h"
#include "Headers/DataHeader.h"
#include "DataFormatsZDC/RawEventData.h"
#include "ZDCSimulation/Digits2Raw.h"
#include "ZDCRaw/RawReconstruction.h"
#include <vector>
#include <sstream>

using namespace o2::framework;

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(
    ConfigParamSpec{
      "input-spec", VariantType::String, "A:ZDC/RAWDATA", {"selection string input specs"}});
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& config)
{
  WorkflowSpec workflow;
  workflow.emplace_back(DataProcessorSpec{
    "zdc-raw-reco",
    select(config.options().get<std::string>("input-spec").c_str()),
    Outputs{},
    AlgorithmSpec{[](InitContext& setup) {
        auto loglevel = setup.options().get<int>("log-level");
	std::string ccdb_url = setup.options().get<std::string>("ccdb-url");
        return adaptStateless([loglevel,ccdb_url](InputRecord& inputs, DataAllocator& outputs) {
          long timeStamp = 0;
	  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
	  mgr.setURL(ccdb_url);
	  if (timeStamp == mgr.getTimestamp()) {
	    LOG(FATAL) << "Cannot get imestamp";
	  }
	  mgr.setTimestamp(timeStamp);
	  auto moduleConfig = mgr.get<o2::zdc::ModuleConfig>(o2::zdc::CCDBPathConfigModule);
	  if (!moduleConfig) {
	    LOG(FATAL) << "Cannot module configuratio for timestamp " << timeStamp;
	    return;
	  }
	  LOG(INFO) << "Loaded module configuration for timestamp " << timeStamp;
	  o2::zdc::RawReconstruction zdc_rr;
	  zdc_rr.setModuleConfig(moduleConfig);
          zdc_rr.init();
          zdc_rr.setVerbosity(loglevel);
          DPLRawParser parser(inputs);
          o2::header::DataHeader const* lastDataHeader = nullptr;
          std::stringstream rdhprintout;
          for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
            // retrieving RDH v4
            auto const* rdh = it.get_if<o2::header::RAWDataHeaderV4>();
            // retrieving the raw pointer of the page
            auto const* raw = it.raw();
            // retrieving payload pointer of the page
            auto const* payload = it.data();
            // size of payload
            size_t payloadSize = it.size();
            // offset of payload in the raw page
            size_t offset = it.offset();
            // Note: the following code is only for printing out raw page information
            const auto* dh = it.o2DataHeader();
            if (loglevel > 0) {
              if (dh != lastDataHeader) {
                // print the DataHeader information only for the first part or if we have high verbosity
                if (loglevel > 1 || dh->splitPayloadIndex == 0) {
                  rdhprintout << dh->dataOrigin.as<std::string>() << "/"
                              << dh->dataDescription.as<std::string>() << "/"
                              << dh->subSpecification << "  ";
                  // at high verbosity print part number, otherwise only the total number of parts
                  if (loglevel > 1) {
                    rdhprintout << "part " + std::to_string(dh->splitPayloadIndex) + " of " + std::to_string(dh->splitPayloadParts);
                  } else {
                    rdhprintout << " " + std::to_string(dh->splitPayloadParts) + " part(s)";
                  }
                  rdhprintout << " payload size " << dh->payloadSize << std::endl;
                }
                if (!rdhprintout.str().empty()) {
                  LOG(INFO) << rdhprintout.str();
                  rdhprintout.str(std::string());
                }
              }
              if (payload != nullptr) {
                for (Int_t ip = 0; ip < payloadSize; ip += 16) {
                  zdc_rr.processWord((const UInt_t*)&payload[ip]);
                }
              }
            }
            lastDataHeader = dh;
          }
          if (loglevel > 0) {
            LOG(INFO) << rdhprintout.str();
          }
          zdc_rr.write();
        }); }},
    Options{
      {"log-level", VariantType::Int, 1, {"Logging level [0-2]"}},
      {"ccdb-url", VariantType::String, "http://ccdb-test.cern.ch:8080", {"Path to CCDB"}}}});
  return workflow;
}
