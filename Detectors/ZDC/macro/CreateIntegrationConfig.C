// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "FairLogger.h"
#include "CCDB/CcdbApi.h"
#include "ZDCRaw/ZDCIntegrationParam.h"
#include "ZDCBase/Constants.h"
#include <string>
#include <TFile.h>
#include <map>

#endif

using namespace o2::zdc;
using namespace std;

void CreateIntegrationConfig(long tmin = 0, long tmax = -1,
                             std::string ccdbHost = "http://ccdb-test.cern.ch:8080")
{

  ZDCIntegrationParam conf;

  // Integration limits for all signals
  // Values should be in range 0-11
  conf.setIntegration(IdZNAC, 5, 9, 5, 9);
  conf.setIntegration(IdZNA1, 5, 9, 5, 9);
  conf.setIntegration(IdZNA2, 5, 9, 5, 9);
  conf.setIntegration(IdZNA3, 5, 9, 5, 9);
  conf.setIntegration(IdZNA4, 5, 9, 5, 9);
  conf.setIntegration(IdZNASum, 5, 9, 5, 9);
  //
  conf.setIntegration(IdZPAC, 5, 9, 5, 9);
  conf.setIntegration(IdZPA1, 5, 9, 5, 9);
  conf.setIntegration(IdZPA2, 5, 9, 5, 9);
  conf.setIntegration(IdZPA3, 5, 9, 5, 9);
  conf.setIntegration(IdZPA4, 5, 9, 5, 9);
  conf.setIntegration(IdZPASum, 5, 9, 5, 9);
  //
  conf.setIntegration(IdZEM1, 5, 9, 5, 9);
  conf.setIntegration(IdZEM2, 5, 9, 5, 9);
  //
  conf.setIntegration(IdZNCC, 5, 9, 5, 9);
  conf.setIntegration(IdZNC1, 5, 9, 5, 9);
  conf.setIntegration(IdZNC2, 5, 9, 5, 9);
  conf.setIntegration(IdZNC3, 5, 9, 5, 9);
  conf.setIntegration(IdZNC4, 5, 9, 5, 9);
  conf.setIntegration(IdZNCSum, 5, 9, 5, 9);
  //
  conf.setIntegration(IdZPCC, 5, 9, 5, 9);
  conf.setIntegration(IdZPC1, 5, 9, 5, 9);
  conf.setIntegration(IdZPC2, 5, 9, 5, 9);
  conf.setIntegration(IdZPC3, 5, 9, 5, 9);
  conf.setIntegration(IdZPC4, 5, 9, 5, 9);
  conf.setIntegration(IdZPCSum, 5, 9, 5, 9);

  conf.print();

  o2::ccdb::CcdbApi api;
  map<string, string> metadata; // can be empty
  api.init(ccdbHost.c_str());   // or http://localhost:8080 for a local installation
  // store abitrary user object in strongly typed manner
  api.storeAsTFileAny(&conf, CCDBPathTDCCalib, metadata, tmin, tmax);

  // return conf;
}
