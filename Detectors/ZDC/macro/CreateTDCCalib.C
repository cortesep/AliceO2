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
#include "ZDCRaw/ZDCTDCParam.h"
#include "ZDCBase/Constants.h"
#include <string>
#include <TFile.h>
#include <map>

#endif

using namespace o2::zdc;
using namespace std;

void CreateTDCCalib(long tmin = 0, long tmax = -1,
                    std::string ccdbHost = "http://ccdb-test.cern.ch:8080")
{

  ZDCTDCParam conf;

  int modID;

  // Recentering of TDCs
  conf.setShift(TDCZNAC, 14);
  conf.setShift(TDCZNAS, 14);
  conf.setShift(TDCZPAC, 14);
  conf.setShift(TDCZPAS, 14);
  conf.setShift(TDCZEM1, 14);
  conf.setShift(TDCZEM2, 14);
  conf.setShift(TDCZNCC, 14);
  conf.setShift(TDCZNCS, 14);
  conf.setShift(TDCZPCC, 14);
  conf.setShift(TDCZPCS, 14);

  conf.print();

  o2::ccdb::CcdbApi api;
  map<string, string> metadata; // can be empty
  api.init(ccdbHost.c_str());   // or http://localhost:8080 for a local installation
  // store abitrary user object in strongly typed manner
  api.storeAsTFileAny(&conf, CCDBPathTDCCalib, metadata, tmin, tmax);

  // return conf;
}
