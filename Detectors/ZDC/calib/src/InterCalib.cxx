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

#include <TROOT.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TMinuit.h>
#include <THnBase.h>
#include <THnSparse.h>
#include <TAxis.h>
#include "ZDCCalib/InterCalib.h"
#include "ZDCReconstruction/ZDCEnergyParam.h"
#include "ZDCReconstruction/ZDCTowerParam.h"
#include "Framework/Logger.h"

using namespace o2::zdc;

double InterCalib::add[InterCalib::NPAR][InterCalib::NPAR] = {0};
std::mutex InterCalib::mtx;

int InterCalib::init()
{
  std::string hn[NH] = {"ZNA", "ZPA", "ZNC", "ZPC", "ZEM"};
  if (mInterCalibConfig == nullptr) {
    LOG(fatal) << "o2::zdc::InterCalib: missing configuration object";
  }
  clear();
  for (int32_t ih = 0; ih < (2 * NH); ih++) {
    if (h[ih]) {
      delete h[ih];
    }
  }
  for (int32_t ih = 0; ih < NH; ih++) {
    if (hc[ih]) {
      delete hc[ih];
    }
  }
  int ih;
  // clang-format off
  ih = 0; h[ih] = new TH1F("hZNAS","ZNA sum",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 1; h[ih] = new TH1F("hZPAS","ZPA sum",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 2; h[ih] = new TH1F("hZNCS","ZNC sum",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 3; h[ih] = new TH1F("hZPCS","ZPC sum",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 4; h[ih] = new TH1F("hZEM2","ZEM2"   ,mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 0; h[NH+ih] = new TH1F("hZNAC","ZNA TC",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 1; h[NH+ih] = new TH1F("hZPAC","ZPA TC",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 2; h[NH+ih] = new TH1F("hZNCC","ZNC TC",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 3; h[NH+ih] = new TH1F("hZPCC","ZPC TC",mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 4; h[NH+ih] = new TH1F("hZEM1","ZEM1",  mInterCalibConfig->nb1[ih],mInterCalibConfig->amin1[ih],mInterCalibConfig->amax1[ih]);
  ih = 0; hc[ih] = new TH2F("cZNA","ZNA;TC;SUM",mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih],mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih]);
  ih = 1; hc[ih] = new TH2F("cZPA","ZPA;TC;SUM",mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih],mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih]);
  ih = 2; hc[ih] = new TH2F("cZNC","ZNC;TC;SUM",mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih],mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih]);
  ih = 3; hc[ih] = new TH2F("cZPC","ZPC;TC;SUM",mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih],mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih]);
  ih = 4; hc[ih] = new TH2F("cZEM","ZEM;ZEM1;ZEM2",mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih],mInterCalibConfig->nb2[ih],mInterCalibConfig->amin2[ih],mInterCalibConfig->amax2[ih]);
  // clang-format on
  mInitDone = true;
  return 0;
}

int InterCalib::process(const gsl::span<const o2::zdc::BCRecData>& RecBC,
                        const gsl::span<const o2::zdc::ZDCEnergy>& Energy,
                        const gsl::span<const o2::zdc::ZDCTDCData>& TDCData,
                        const gsl::span<const uint16_t>& Info)
{
  if (!mInitDone) {
    init();
  }
  LOG(info) << "o2::zdc::InterCalib processing " << RecBC.size() << " b.c.";
  o2::zdc::RecEventFlat ev;
  ev.init(RecBC, Energy, TDCData, Info);
  while (ev.next()) {
    if (ev.getNInfo() > 0) {
      auto& decodedInfo = ev.getDecodedInfo();
      for (uint16_t info : decodedInfo) {
        uint8_t ch = (info >> 10) & 0x1f;
        uint16_t code = info & 0x03ff;
        // hmsg->Fill(ch, code);
      }
      ev.print();
      // Need clean data (no messages)
      // We are sure there is no pile-up in any channel (too restrictive?)
      continue;
    }
    if (ev.getNEnergy() > 0 && ev.mCurB.triggers == 0) {
      printf("%9u.%04u Untriggered bunch\n", ev.mCurB.ir.orbit, ev.mCurB.ir.bc);
      // Skip!
      continue;
    }
    //     // Trigger bits are not propagated!!!
    //     heznac->Fill(ev.EZNAC());
    //     auto tdcid = o2::zdc::TDCZNAC;
    //     auto nhit = ev.NtdcV(tdcid);
    //     if (ev.NtdcA(tdcid) != nhit) {
    //       fprintf(stderr, "Mismatch in TDC data\n");
    //       continue;
    //     }
    //     if (nhit > 0) {
    //       double bc_d = uint32_t(ev.ir.bc / 100);
    //       double bc_m = uint32_t(ev.ir.bc % 100);
    //       hbznac->Fill(bc_m, -bc_d);
    //       for (int ihit = 0; ihit < nhit; ihit++) {
    //         htznac->Fill(ev.tdcV(tdcid, ihit), ev.tdcA(tdcid, ihit));
    //       }
    //     }
  }
  return 0;
}

int InterCalib::process(const char* hname, int ic)
{
  // Run 2 ZDC calibration is based on multi-dimensional histograms
  // with dimensions:
  // TC, T1, T2, T3, T4, Trigger
  // ic is the number of the selected trigger class
  THnSparse* hs = (THnSparse*)gROOT->FindObject(hname);
  if (hs == 0) {
    printf("Not found: %s\n", hname);
    return -1;
  }
  if (!hs->InheritsFrom(THnSparse::Class())) {
    printf("Not a THnSparse: %s\n", hname);
    hs->IsA()->Print();
    return -1;
  }
  TString hn = hname;
  int ih = -1;
  if (hn.EqualTo("hZNA")) {
    ih = 0;
  } else if (hn.EqualTo("hZPA")) {
    ih = 1;
  } else if (hn.EqualTo("hZNC")) {
    ih = 2;
  } else if (hn.EqualTo("hZPC")) {
    ih = 3;
  } else if (hn.EqualTo("hZEM")) {
    ih = 4;
  } else {
    printf("Not recognized histogram name: %s\n", hname);
    return -1;
  }
  clear(ih);
  const int32_t dim = 6;
  double x[dim];
  int32_t bins[dim];
  int64_t nb = hs->GetNbins();
  int64_t nn = 0;
  printf("Histogram %s has %ld bins\n", hname, nb);
  for (int64_t i = 0; i < nb; i++) {
    double cont = hs->GetBinContent(i, bins);
    if (cont <= 0) {
      continue;
    }
    for (int32_t d = 0; d < dim; ++d) {
      x[d] = hs->GetAxis(d)->GetBinCenter(bins[d]);
    }
    if (TMath::Nint(x[5] - ic) == 0 && x[0] > mInterCalibConfig->cutLow[ih] && x[0] < mInterCalibConfig->cutHigh[ih]) {
      nn++;
      cumulate(ih, x[0], x[1], x[2], x[3], x[4], cont);
    }
  }
  mini(ih);
  printf("Processed RUN2 data for %s ih = %d\n", hname, ih);
  printf("Trigger class selection %d and %d cuts (%g:%g): %ld\n", ic, nn);
  return 0;
}

void InterCalib::clear(int ih)
{
  int ihstart = 0;
  int ihstop = NH;
  if (ih >= 0 && ih < NH) {
    ihstart = ih;
    ihstop = ih + 1;
  }
  for (int32_t ii = ihstart; ii < ihstop; ii++) {
    for (int32_t i = 0; i < NPAR; i++) {
      for (int32_t j = 0; j < NPAR; j++) {
        sum[ii][i][j] = 0;
      }
    }
  }
}

void InterCalib::cumulate(int ih, double tc, double t1, double t2, double t3, double t4, double w = 1)
{
  // TODO: add histogram
  // TODO: store data to redraw histograms
  if (tc < mInterCalibConfig->cutLow[ih] || tc > mInterCalibConfig->cutHigh[ih]) {
    return;
  }
  double val[NPAR] = {0, 0, 0, 0, 0, 1};
  val[0] = tc;
  val[1] = t1;
  val[2] = t2;
  val[3] = t3;
  val[4] = t4;
  for (int32_t i = 0; i < 6; i++) {
    for (int32_t j = i; j < 6; j++) {
      sum[ih][i][j] += val[i] * val[j] * w;
    }
  }
  h[ih]->Fill(val[1] + val[2] + val[3] + val[4]);
  h[ih + NH]->Fill(val[0]);
  hc[ih]->Fill(val[0], val[1] + val[2] + val[3] + val[4]);
}

void InterCalib::fcn(int& npar, double* gin, double& chi, double* par, int iflag)
{
  chi = 0;
  for (int32_t i = 0; i < NPAR; i++) {
    for (int32_t j = 0; j < NPAR; j++) {
      chi += (i == 0 ? par[i] : -par[i]) * (j == 0 ? par[j] : -par[j]) * add[i][j];
    }
  }
  // printf("%g\n", chi);
}

int InterCalib::mini(int ih)
{
  // Copy to static object and symmetrize matrix
  // We use a static function and therefore we can do only one minimization at a time
  mtx.lock();
  for (int32_t i = 0; i < NPAR; i++) {
    for (int32_t j = 0; j < NPAR; j++) {
      if (j < i) {
        add[i][j] = sum[ih][j][i];
      } else {
        add[i][j] = sum[ih][i][j];
      }
    }
  }
  double arglist[10];
  int ierflg = 0;
  double l_bnd = 0.2;
  double u_bnd = 5.;
  double start = 1.0;
  double step = 0.1;
  TMinuit minuit(NPAR);
  minuit.SetFCN(fcn);
  minuit.mnparm(0, "c0", 1., 0., 1., 1., ierflg);
  minuit.mnparm(1, "c1", start, step, l_bnd, u_bnd, ierflg);
  if (ih == 4) {
    // Only two ZEM calorimeters: equalize response
    l_bnd = 0;
    u_bnd = 0;
    start = 0;
    step = 0;
  }
  minuit.mnparm(2, "c2", start, step, l_bnd, u_bnd, ierflg);
  minuit.mnparm(3, "c3", start, step, l_bnd, u_bnd, ierflg);
  minuit.mnparm(4, "c4", start, step, l_bnd, u_bnd, ierflg);
  l_bnd = -20;
  u_bnd = 20;
  start = 0;
  step = 1;
  step = 0;
  minuit.mnparm(5, "offset", start, step, l_bnd, u_bnd, ierflg);
  minuit.mnexcm("MIGRAD", arglist, 0, ierflg);
  for (Int_t i = 0; i < NPAR; i++) {
    minuit.GetParameter(i, par[ih][i], err[ih][i]);
  }
  mtx.unlock();
  return 0;
}
