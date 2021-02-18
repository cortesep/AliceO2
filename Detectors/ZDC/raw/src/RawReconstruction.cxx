#include <TROOT.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <limits>
#include "ZDCRaw/RawReconstruction.h"
#include "ZDCSimulation/Digits2Raw.h"
#include "FairLogger.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

using namespace o2::zdc;

void RawReconstruction::setStat(TH1* h)
{
  TString hn = h->GetName();
  h->Draw();
  gPad->Update();
  TPaveStats* st = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  st->SetFillStyle(0);
  st->SetBorderSize(1);
  if (hn.BeginsWith("hp")) {
    st->SetOptStat(111111);
    st->SetX1NDC(0.1);
    st->SetX2NDC(0.3);
    st->SetY1NDC(0.640);
    st->SetY2NDC(0.9);
  } else if (hn.BeginsWith("hc")) {
    st->SetOptStat(1111);
    st->SetX1NDC(0.799);
    st->SetX2NDC(0.999);
    st->SetY1NDC(0.829);
    st->SetY2NDC(0.999);
  } else if (hn.BeginsWith("hs") || hn.BeginsWith("hb")) {
    st->SetOptStat(11);
    st->SetX1NDC(0.799);
    st->SetX2NDC(0.9995);
    st->SetY1NDC(0.904);
    st->SetY2NDC(0.999);
  }
}

void RawReconstruction::init()
{
  auto& sopt = ZDCSimParam::Instance();
  mIsContinuous = sopt.continuous;
  mNBCAHead = mIsContinuous ? sopt.nBCAheadCont : sopt.nBCAheadTrig;

  if (!mModuleConfig) {
    LOG(FATAL) << "Missing ModuleConfig configuration object";
    return;
  }

  // Prepare tapered sinc function
  // tsc/TSN =3.75 (~ 4) and TSL*TSN*sqrt(2)/tsc >> 1 (n. of sigma)
  const Double_t tsc = 750;
  Int_t n = TSL * TSN;
  for (Int_t tsi = 0; tsi <= n; tsi++) {
    Double_t arg1 = TMath::Pi() * Double_t(tsi) / Double_t(TSN);
    Double_t fs = 1;
    if (arg1 != 0)
      fs = TMath::Sin(arg1) / arg1;
    Double_t arg2 = Double_t(tsi) / tsc;
    Double_t fg = TMath::Exp(-arg2 * arg2);
    mTS[n + tsi] = fs * fg;
    mTS[n - tsi] = mTS[n + tsi]; // Function is even
  }

  // Open debug file
  mDbg = new TFile("ZDCReco.root", "recreate");
  mTDbg = new TTree("zdcr", "ZDCReco");
  mTDbg->Branch("zdcr", "RecEvent", &mRec);
  // Update reconstruction parameters
  //auto& ropt=ZDCRecoParam::Instance();
  o2::zdc::ZDCRecoParam& ropt = const_cast<o2::zdc::ZDCRecoParam&>(ZDCRecoParam::Instance());

  // Fill maps
  for (Int_t itdc = 0; itdc < o2::zdc::NTDCChannels; itdc++) {
    // If the reconstruction parameters were not manually set
    if (ropt.tmod[itdc] < 0 || ropt.tch[itdc] < 0) {
      Int_t isig = TDCSignal[itdc];
      for (Int_t im = 0; im < NModules; im++) {
        for (UInt_t ic = 0; ic < NChPerModule; ic++) {
          if (mModuleConfig->modules[im].channelID[ic] == isig && mModuleConfig->modules[im].readChannel[ic]) {
            //ropt.updateFromString(TString::Format("ZDCRecoParam.tmod[%d]=%d;",itdc,im));
            //ropt.updateFromString(TString::Format("ZDCRecoParam.tch[%d]=%d;",itdc,ic));
            ropt.tmod[itdc] = im;
            ropt.tch[itdc] = ic;
            goto next_itdc;
          }
        }
      }
    }
  next_itdc:;
    LOG(INFO) << "TDC " << itdc << " mod " << ropt.tmod[itdc] << " ch " << ropt.tch[itdc];
  }

  // TDC calibration
  for (Int_t itdc = 0; itdc < o2::zdc::NTDCChannels; itdc++) {
    float fval = ropt.tdc_shift[itdc];
    // If the reconstruction parameters were not manually set
    if (fval < 0) {
      // Check if calibration object is present
      if (!mTDCParam) {
        LOG(FATAL) << "TDC " << itdc << " missing configuration object and no manual override";
      } else {
        fval = mTDCParam->getShift(itdc) * 96.;
      }
    }
    auto val = std::nearbyint(fval);
    if (val < kMinShort) {
      LOG(FATAL) << "Shift for TDC " << itdc << " " << val << " is out of range";
    }
    if (val > kMaxShort) {
      LOG(FATAL) << "Shift for TDC " << itdc << " " << val << " is out of range";
    }
    tdc_shift[itdc] = val;
    LOG(INFO) << itdc << " " << ChannelNames[TDCSignal[itdc]] << " shift= " << tdc_shift[itdc] << " i.s. = " << val / 96. << " ns";
  }

  // Fill maps channel maps for integration
  for (Int_t ich = 0; ich < NChannels; ich++) {
    // If the reconstruction parameters were not manually set
    if (ropt.amod[ich] < 0 || ropt.ach[ich] < 0) {
      for (Int_t im = 0; im < NModules; im++) {
        for (UInt_t ic = 0; ic < NChPerModule; ic++) {
          if (mModuleConfig->modules[im].channelID[ic] == ich && mModuleConfig->modules[im].readChannel[ic]) {
            ropt.amod[ich] = im;
            ropt.ach[ich] = ic;
            goto next_ich;
          }
        }
      }
    }
  next_ich:;
    LOG(INFO) << "ADC " << ich << "(" << ChannelNames[ich] << ") mod " << ropt.amod[ich] << " ch " << ropt.ach[ich];
  }

  // Integration ranges
  for (Int_t ich = 0; ich < NChannels; ich++) {
    // If the reconstruction parameters were not manually set
    if (ropt.beg_int[ich] < 0 || ropt.end_int[ich] < 0) {
      if (!mIntParam) {
        LOG(FATAL) << "Integration for signal " << ich << " missing configuration object and no manual override";
      } else {
        ropt.beg_int[ich] = mIntParam->beg_int[ich];
        ropt.end_int[ich] = mIntParam->end_int[ich];
      }
    }
    if (ropt.beg_ped_int[ich] < 0 || ropt.end_ped_int[ich] < 0) {
      if (!mIntParam) {
        LOG(FATAL) << "Integration for pedestal " << ich << " missing configuration object and no manual override";
      } else {
        ropt.beg_ped_int[ich] = mIntParam->beg_ped_int[ich];
        ropt.end_ped_int[ich] = mIntParam->end_ped_int[ich];
      }
    }
    LOG(INFO) << ChannelNames[ich] << " integration: signal=[" << ropt.beg_int[ich] << "-" << ropt.end_int[ich] << "] pedestal=[" << ropt.beg_ped_int[ich] << "-" << ropt.end_ped_int[ich] <<"]";
  }
  // Filling debug histograms
  gROOT->SetBatch();
  int nbx = (sopt.nBCAheadTrig + 1) * NTimeBinsPerBC;
  Double_t xmin = -sopt.nBCAheadTrig * NTimeBinsPerBC - 0.5;
  Double_t xmax = NTimeBinsPerBC - 0.5;
  for (uint32_t i = 0; i < NDigiChannels; i++) {
    uint32_t imod = i / NChPerModule;
    uint32_t ich = i % NChPerModule;
    if (mBaseline[i]) {
      mBaseline[i]->Reset();
    } else {
      TString hname = TString::Format("hp%d%d", imod, ich);
      TString htit = TString::Format("Baseline mod. %d ch. %d;Average orbit baseline", imod, ich);
      //mBaseline[i]=new TH1F(hname,htit,ADCRange,ADCMin-0.5,ADCMax+0.5);
      mBaseline[i] = new TH1F(hname, htit, 16378, -0.125, ADCMax + 0.125);
    }
    if (mCounts[i]) {
      mCounts[i]->Reset();
    } else {
      TString hname = TString::Format("hc%d%d", imod, ich);
      TString htit = TString::Format("Counts mod. %d ch. %d; Orbit hits", imod, ich);
      mCounts[i] = new TH1F(hname, htit, o2::constants::lhc::LHCMaxBunches + 1, -0.5, o2::constants::lhc::LHCMaxBunches + 0.5);
    }
    if (mSignal[i]) {
      mSignal[i]->Reset();
    } else {
      TString hname = TString::Format("hs%d%d", imod, ich);
      TString htit = TString::Format("Signal mod. %d ch. %d; Sample; ADC", imod, ich);
      mSignal[i] = new TH2F(hname, htit, nbx, xmin, xmax, ADCRange, ADCMin - 0.5, ADCMax + 0.5);
    }
    if (mBunch[i]) {
      mBunch[i]->Reset();
    } else {
      TString hname = TString::Format("hb%d%d", imod, ich);
      TString htit = TString::Format("Bunch mod. %d ch. %d; Sample; ADC", imod, ich);
      mBunch[i] = new TH2F(hname, htit, 100, -0.5, 99.5, 36, -35.5, 0.5);
    }
  }
  // Word id not present in payload
  mCh.f.fixed_0 = Id_wn;
  mCh.f.fixed_1 = Id_wn;
  mCh.f.fixed_2 = Id_wn;
}

int RawReconstruction::write()
{
  // This function should call processOrbit to make sure
  // all available data has been processed before closing
  // the output file
  Int_t risp = processOrbit(true);
  mDbg->cd();
  for (UInt_t i = 0; i < NDigiChannels; i++) {
    if (mBunch[i] && mBunch[i]->GetEntries() > 0) {
      setStat(mBunch[i]);
      mBunch[i]->Write();
    }
  }
  for (UInt_t i = 0; i < NDigiChannels; i++) {
    if (mBaseline[i] && mBaseline[i]->GetEntries() > 0) {
      setStat(mBaseline[i]);
      mBaseline[i]->Write();
    }
  }
  for (UInt_t i = 0; i < NDigiChannels; i++) {
    if (mCounts[i] && mCounts[i]->GetEntries() > 0) {
      setStat(mCounts[i]);
      mCounts[i]->Write();
    }
  }
  for (UInt_t i = 0; i < NDigiChannels; i++) {
    if (mSignal[i] && mSignal[i]->GetEntries() > 0) {
      setStat(mSignal[i]);
      mSignal[i]->Write();
    }
  }
  mDbg->Write();
  mDbg->Close();
  return risp;
}

inline int RawReconstruction::getHPos(uint32_t board, uint32_t ch)
{
  int ih = board * 4 + ch;
  if (ih < NDigiChannels) {
    return ih;
  } else {
    LOG(ERROR) << "Wrong ih " << ih << " board " << board << " ch " << ch;
    return -1;
  }
}

int RawReconstruction::processWord(const UInt_t* word)
{
  if (word == nullptr) {
    LOG(ERROR) << "NULL pointer";
    return 1;
  }
  if ((word[0] & 0x3) == Id_w0) {
    mCh.reset();
    for (Int_t iw = 0; iw < NWPerGBTW; iw++) {
      mCh.w[0][iw] = word[iw];
    }
  } else if ((word[0] & 0x3) == Id_w1) {
    if (mCh.f.fixed_0 == Id_w0) {
      for (Int_t iw = 0; iw < NWPerGBTW; iw++) {
        mCh.w[1][iw] = word[iw];
      }
    } else {
      LOG(ERROR) << "Wrong word sequence";
      mCh.f.fixed_0 = Id_wn;
      mCh.f.fixed_1 = Id_wn;
      mCh.f.fixed_2 = Id_wn;
    }
  } else if ((word[0] & 0x3) == Id_w2) {
    if (mCh.f.fixed_0 == Id_w0 && mCh.f.fixed_1 == Id_w1) {
      for (Int_t iw = 0; iw < NWPerGBTW; iw++) {
        mCh.w[2][iw] = word[iw];
      }
      // At this stage all the payload for the current channel is available
      processCh();
    } else {
      LOG(ERROR) << "Wrong word sequence";
    }
    // Invalidate processed information
    mCh.f.fixed_0 = Id_wn;
    mCh.f.fixed_1 = Id_wn;
    mCh.f.fixed_2 = Id_wn;
  } else {
    // Word not present in payload
    LOG(FATAL) << "Event format error";
    return 1;
  }
  return 0;
}

int RawReconstruction::processCh()
{
  if (mData.empty()) {
    // It should be empty only at first use
    mData.emplace_back();
    auto& data = mData.back();
    data.reset();
    mReco.emplace_back();
    auto& reco = mReco.back();
    reco.ir.orbit = mCh.f.orbit;
    reco.ir.bc = mCh.f.bc;
    LOG(INFO) << "New entry: orbit " << reco.ir.orbit << " bc " << reco.ir.bc;
  } else {
    auto& last_reco = mReco.back();
    // Detect event change and orbit change
    if (last_reco.ir.orbit != mCh.f.orbit || last_reco.ir.bc != mCh.f.bc) {
      // If orbit changed process last orbit
      if (last_reco.ir.orbit != mCh.f.orbit) {
        LOG(INFO) << "Orbit change " << last_reco.ir.orbit << " -> " << mCh.f.orbit;
        processOrbit();
      }
      // Then add new event
      mData.emplace_back();
      auto& data = mData.back();
      data.reset();
      mReco.emplace_back();
      auto& reco = mReco.back();
      reco.ir.orbit = mCh.f.orbit;
      reco.ir.bc = mCh.f.bc;
      LOG(INFO) << "Add entry: orbit " << reco.ir.orbit << " bc " << reco.ir.bc;
    }
  }
  // Assign channel data
  auto& data = mData.back();
  uint32_t board = mCh.f.board;
  uint32_t ch = mCh.f.ch;
  // Check if it is in allowed range
  if (board >= NModules || ch >= NChPerModule) {
    LOG(ERROR) << "board " << board << " ch " << ch << " not in allowed range";
    return 1;
  } else {
    auto& chd = data.data[board][ch];
    for (Int_t i = 0; i < NWPerBc; i++) {
      for (Int_t j = 0; j < NWPerGBTW; j++) {
        chd.w[i][j] = mCh.w[i][j];
      }
    }
    // Convert raw samples to signed integers
    UShort_t us[NTimeBinsPerBC];
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
        data.s[board][ch][i] = us[i] - ADCRange;
      } else {
        data.s[board][ch][i] = us[i];
      }
    }
  }
  //auto& reco = mReco.back();
  //LOG(INFO) << "Processed: orbit " << reco.ir.orbit << " bc " << reco.ir.bc << " board " << board << " ch " << ch;
  return 0;
}

void RawReconstruction::printIR()
{
  std::deque<RecEvent>::iterator it = mReco.begin();
  while (it != mReco.end()) {
    LOG(INFO) << (*it).ir.orbit << "." << (*it).ir.bc;
    it++;
  }
}

//int16_t
int RawReconstruction::processOrbit(bool islast)
{
  // This function starts reconstruction of a full orbit
  // Additional data from previous orbit may be present at the beginning of payload
  auto& last = mReco.back();
  auto last_orbit = last.ir.orbit;
  auto last_bc = last.ir.bc;
  auto& first = mReco.front();
  auto first_orbit = first.ir.orbit;
  auto first_bc = first.ir.bc;
  LOG(INFO) << "processOrbit(" << islast << ") " << first_orbit << "." << first_bc << " to " << last_orbit << "." << last_bc;
  // Check if orbit number is not increasing
  if (last_orbit < first_orbit) {
    LOG(FATAL) << "Orbit number is not increasing " << first_orbit << " -> " << last_orbit;
    mReco.clear();
    mData.clear();
    return 1;
  }

  // Check if there is a gap in the orbit
  Int_t nb = 0;
  if ((last_orbit - first_orbit) > 1) {
    std::deque<RecEvent>::iterator it = mReco.begin();
    while (it != mReco.end()) {
      if ((*it).ir.orbit == first_orbit) {
        nb++;
      } else {
        break;
      }
    }
    LOG(WARNING) << "Orbit gap " << first_orbit << " -> " << last_orbit << " erasing " << nb << " bunches";
    mReco.erase(mReco.begin(), mReco.begin() + nb);
    mData.erase(mData.begin(), mData.begin() + nb);
    return processOrbit();
  }

  // Identify last bunch in orbit to get pedestal information
  mLast = -1;
  for (std::deque<RecEvent>::reverse_iterator rit = mReco.rbegin(); rit != mReco.rend(); ++rit) {
    auto bc = (*rit).ir.bc;
    auto orbit = (*rit).ir.orbit;
    if (bc == 3563 && orbit == last_orbit) {
      int index = rit - mReco.rbegin();
      mLast = mReco.size() - index - 1;
      // Should always break at first iteration
      break;
    }
  }
  // Prepare pedestal for current orbit
  for (Int_t im = 0; im < NModules; im++) {
    for (Int_t ic = 0; ic < NChPerModule; ic++) {
      if (mLast >= 0) {
        float pedestal = mData[mLast].data[im][ic].f.offset;
        pedestal = (mData[mLast].data[im][ic].f.offset - 32768.) / 8.;
        mOffset[im][ic] = pedestal;
      } else {
        // TODO: find more meaningful missing value
        mOffset[im][ic] = 0;
      }
    }
  }

  // Find consecutive bc data belonging to current orbit
  std::deque<RecEvent>::iterator it = mReco.begin();
  auto* bcd_last = &((*it).ir);
  Int_t seq_beg = 0;
  Int_t seq_end = 0;
  Int_t iti = 0;
  for (; it != mReco.end(); it++, iti++) {
    auto* bcd = &((*it).ir);
    int64_t bcdiff = bcd->differenceInBC(*bcd_last);
    // Find a gap
    if (bcdiff > 1) {
      // Process consecutive BCs
      if (seq_beg == seq_end) {
        LOG(INFO) << "Lonely bunch " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc;
      } else {
        LOG(INFO) << "Processing " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc << " - " << mReco[seq_end].ir.orbit << "." << mReco[seq_end].ir.bc;
        reconstruct(seq_beg, seq_end);
      }
      seq_beg = iti;
      seq_end = iti;
    } else {
      seq_end = iti;
    }
    bcd_last = bcd;
  }
  if (seq_beg != seq_end) {
    // No recipe to process a lonely bunch
    LOG(INFO) << "Processing " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc << " - " << mReco[seq_end].ir.orbit << "." << mReco[seq_end].ir.bc;
    reconstruct(seq_beg, seq_end);
  }
  // Cleanup
  if (islast) {
    LOG(INFO) << "End of reconstruction and cleanup";
    mReco.clear();
    mData.clear();
  } else {
    if (last_bc != 3563) {
      // There is a gap at the end of the processed orbit (bunch 3563 not acquired)
      mReco.clear();
      mData.clear();
      LOG(INFO) << "Erasing orbit " << last_orbit;
    } else {
      // Cleanup current orbit but do not delete last bunches
      int nb = mReco.size();
      int nkeep = 0;
      int last_bc = 3563;
      for (std::deque<RecEvent>::reverse_iterator rit = mReco.rbegin(); rit != mReco.rend() || nkeep == 3; ++rit) {
        int bc = (*rit).ir.bc;
        if (bc == 3563) {
          nkeep++;
        } else {
          if ((last_bc - bc) > 1) {
            break;
          } else {
            nkeep++;
          }
        }
        LOG(INFO) << "Erasing orbit " << last_orbit << " keeping bc from " << last_bc;
        mReco.erase(mReco.begin(), mReco.begin() + nb - nkeep);
        mData.erase(mData.begin(), mData.begin() + nb - nkeep);
      }
    }
  }
  printIR();
  return 0;
}

int RawReconstruction::process(const EventData ev)
{
  for (Int_t im = 0; im < NModules; im++) {
    for (Int_t ic = 0; ic < NChPerModule; ic++) {
      if (ev.data[im][ic].f.fixed_0 == Id_w0 && ev.data[im][ic].f.fixed_1 == Id_w1 && ev.data[im][ic].f.fixed_2 == Id_w2) {
        //        process(ev.data[im][ic]);
      } else if (ev.data[im][ic].f.fixed_0 == 0 && ev.data[im][ic].f.fixed_1 == 0 && ev.data[im][ic].f.fixed_2 == 0) {
        // Empty channel
      } else {
        LOG(ERROR) << "Data format error";
      }
    }
  }
}

void RawReconstruction::processTrigger(int itdc, int ibeg, int iend)
{
  LOG(INFO) << __func__ << " itdc=" << itdc << " " << ibeg << "-" << iend;
  auto& ropt = ZDCRecoParam::Instance();
  Int_t nbun = iend - ibeg + 1;
  Int_t maxs2 = NTimeBinsPerBC * nbun - 1;
  Int_t im = ropt.tmod[itdc];
  Int_t ic = ropt.tch[itdc];
  Int_t shift = ropt.tsh[itdc];
  Int_t thr = ropt.tth[itdc];
  // Redundant check if the TDC channel is connected
  if (im >= 0 && ic >= 0) {
    Int_t is1 = 0, is2 = 1;
    Int_t isfired[3] = {0};
    for (;;) {
      // Shift data
      for (Int_t i = 1; i < 3; i++) {
        isfired[i] = isfired[i - 1];
      }
      Int_t b1 = ibeg + is1 / NTimeBinsPerBC;
      Int_t b2 = ibeg + is2 / NTimeBinsPerBC;
      Int_t s1 = is1 % NTimeBinsPerBC;
      Int_t s2 = is2 % NTimeBinsPerBC;
      int diff = mData[b1].s[im][ic][s1] - mData[b2].s[im][ic][s2];
      // Triple trigger condition
      if (diff > thr) {
        isfired[0] = 1;
        if (isfired[1] == 1 && isfired[2] == 1) {
          // Fired bit is assigned to the second sample, i.e. to the one that can identify the
          // signal peak position
          mReco[b2].fired[itdc][s2] = 1;
          LOG(INFO) << itdc << " Fired @ " << mReco[b2].ir.orbit << "." << mReco[b2].ir.bc << " " << b2 - ibeg << "." << s2 << "=" << NTimeBinsPerBC * (b2 - ibeg) + s2;
        }
      }
      if (is2 >= shift) {
        is1++;
      }
      if (is2 < maxs2) {
        is2++;
      }
      if (is1 == maxs2) {
        break;
      }
    }
    interpolate(itdc, ibeg, iend);
  }
}

void RawReconstruction::reconstruct(int ibeg, int iend)
{
  auto& ropt = ZDCRecoParam::Instance();
  // Apply differential discrimination with triple condition
  for (Int_t itdc = 0; itdc < NTDCChannels; itdc++) {
    Int_t im = ropt.tmod[itdc];
    Int_t ic = ropt.tch[itdc];
    // Check if the TDC channel is connected
    if (im >= 0 && ic >= 0) {
      // Check if channel has valid data for consecutive bunches in current bunch range
      // N.B. there are events recorded from ibeg-iend but we are not sure if it is the
      // case for every TDC channel
      int istart = -1, istop = -1;
      // Loop allows for gaps in the data sequence for each TDC channel
      for (int ibun = ibeg; ibun <= iend; ibun++) {
        auto& ch = mData[ibun].data[im][ic];
        if (ch.f.fixed_0 == Id_w0 && ch.f.fixed_1 == Id_w1 && ch.f.fixed_2 == Id_w2) {
          if (istart < 0) {
            istart = ibun;
          }
          istop = ibun;
        } else {
          // A gap is detected gap
          if (istart >= 0 && (istop - istart) > 0) {
            // Need data for at least two consecutive bunch crossings
            processTrigger(itdc, istart, istop);
          }
          istart = -1;
          istop = -1;
        }
      }
      // Check if there are consecutive bunch crossings at the end of group
      if (istart >= 0 && (istop - istart) > 0) {
        processTrigger(itdc, istart, istop);
      }
    }
  }
  // Reconstruct integrated charges and fill output tree
  // TODO: compare average pedestal with estimation from current event
  // TODO: failover in case of discrepancy
  for (Int_t ibun = ibeg; ibun <= iend; ibun++) {
    mRec = mReco[ibun];
    for (Int_t itdc = 0; itdc < NTDCChannels; itdc++) {
      printf("%d %u.%u %d ", ibun, mReco[ibun].ir.orbit, mReco[ibun].ir.bc, itdc);
      for (Int_t isam = 0; isam < NTimeBinsPerBC; isam++) {
        printf("%d", mRec.fired[itdc][isam]);
      }
      printf("\n");
      for (Int_t i = 0; i < MaxTDCValues; i++) {
        mRec.TdcChannels[itdc][i] = kMinShort;
        mRec.TdcAmplitudes[itdc][i] = -999;
      }
      Int_t i = 0;
      mRec.pattern[itdc] = 0;
      for (int16_t val : mReco[ibun].tdcChannels[itdc]) {
        LOG(INFO) << "TdcChannels[" << itdc << "][" << i << "]=" << val;
        mRec.TdcChannels[itdc][i] = val;
        // There is a TDC value in the search zone around main-main position
        if (std::abs(mRec.TdcChannels[itdc][i]) < ropt.tdc_search[itdc]) {
          mRec.pattern[itdc] = 1;
        }
        i++;
      }
      i = 0;
      for (float val : mReco[ibun].tdcAmplitudes[itdc]) {
        //mRec.tdcAmplitudes[itdc].push_back(val);
        //LOG(INFO) << itdc << " valt=" << val;
        LOG(INFO) << "TdcAmplitudes[" << itdc << "][" << i << "]=" << val;
        mRec.TdcAmplitudes[itdc][i] = val;
        i++;
      }
    }
    printf("%d PATTERN: ", ibun);
    for (Int_t itdc = 0; itdc < NTDCChannels; itdc++) {
      printf("%d", mRec.pattern[itdc]);
    }
    printf("\n");

    // Check if coincidence of common PM and sum of towers is satisfied
    bool fired[NChannels] = {0};
    // Side A
    if ((mRec.pattern[TDCZNAC] || ropt.bitset[TDCZNAC]) && (mRec.pattern[TDCZNAS] || ropt.bitset[TDCZNAS])) {
      for (Int_t ich = IdZNAC; ich <= IdZNASum; ich++) {
        fired[ich] = true;
      }
    }
    if ((mRec.pattern[TDCZPAC] || ropt.bitset[TDCZPAC]) && (mRec.pattern[TDCZPAS] || ropt.bitset[TDCZPAS])) {
      for (Int_t ich = IdZPAC; ich <= IdZPASum; ich++) {
        fired[ich] = true;
      }
    }
    // ZEM1 and ZEM2 are not in coincidence
    fired[IdZEM1] = mRec.pattern[TDCZEM1];
    fired[IdZEM2] = mRec.pattern[TDCZEM2];
    // Side C
    if ((mRec.pattern[TDCZNCC] || ropt.bitset[TDCZNCC]) && (mRec.pattern[TDCZNCS] || ropt.bitset[TDCZNCS])) {
      for (Int_t ich = IdZNCC; ich <= IdZNCSum; ich++) {
        fired[ich] = true;
      }
    }
    if ((mRec.pattern[TDCZPCC] || ropt.bitset[TDCZPCC]) && (mRec.pattern[TDCZPCS] || ropt.bitset[TDCZPCS])) {
      for (Int_t ich = IdZPCC; ich <= IdZPCSum; ich++) {
        fired[ich] = true;
      }
    }

    // Access samples from raw data
    for (Int_t i = 0; i < NChannelsZEM; i++) {
      mRec.energyZEM[i] = -999;
    }
    for (Int_t i = 0; i < NChannelsZN; i++) {
      mRec.energyZNA[i] = -999;
    }
    for (Int_t i = 0; i < NChannelsZN; i++) {
      mRec.energyZNC[i] = -999;
    }
    for (Int_t i = 0; i < NChannelsZP; i++) {
      mRec.energyZPA[i] = -999;
    }
    for (Int_t i = 0; i < NChannelsZP; i++) {
      mRec.energyZPC[i] = -999;
    }
    auto& data = mData[ibun];
    printf("%d FIRED: ", ibun);
    for (Int_t ich = 0; ich < NChannels; ich++) {
      printf("%d", fired[ich]);
    }
    printf("\n");
    for (Int_t ich = 0; ich < NChannels; ich++) {
      // Check if the corresponding TDC is fired
      if (fired[ich]) {
        // Check if channels are present in payload
        Int_t im = ropt.amod[ich];
        Int_t ic = ropt.ach[ich];
        // Check if the ADC channel is connected
        if (im >= 0 && ic >= 0) {
          // Check if the ADC has payload
          auto& ch = data.data[im][ic];
          if (ch.f.fixed_0 == Id_w0 && ch.f.fixed_1 == Id_w1 && ch.f.fixed_2 == Id_w2) {
            float sum = 0;
            for (Int_t is = ropt.beg_int[ich]; is <= ropt.end_int[ich]; is++) {
              // TODO: fallback if offset is missing
              // TODO: fallback if channel has pile-up
              sum += (mOffset[im][ic] - float(data.s[im][ic][is]));
            }
            printf("CH %d %s: %f\n", ich, ChannelNames[ich].data(), sum);
            if (ich == IdZNAC) {
              mRec.energyZNA[0] = sum;
            }
            if (ich == IdZNA1) {
              mRec.energyZNA[1] = sum;
            }
            if (ich == IdZNA2) {
              mRec.energyZNA[2] = sum;
            }
            if (ich == IdZNA3) {
              mRec.energyZNA[3] = sum;
            }
            if (ich == IdZNA4) {
              mRec.energyZNA[4] = sum;
            }
            if (ich == IdZNASum) {
              mRec.energyZNA[5] = sum;
            }
            if (ich == IdZPAC) {
              mRec.energyZPA[0] = sum;
            }
            if (ich == IdZPA1) {
              mRec.energyZPA[1] = sum;
            }
            if (ich == IdZPA2) {
              mRec.energyZPA[2] = sum;
            }
            if (ich == IdZPA3) {
              mRec.energyZPA[3] = sum;
            }
            if (ich == IdZPA4) {
              mRec.energyZPA[4] = sum;
            }
            if (ich == IdZPASum) {
              mRec.energyZPA[5] = sum;
            }
            if (ich == IdZEM1) {
              mRec.energyZEM[0] = sum;
            }
            if (ich == IdZEM2) {
              mRec.energyZEM[1] = sum;
            }
            if (ich == IdZNCC) {
              mRec.energyZNC[0] = sum;
            }
            if (ich == IdZNC1) {
              mRec.energyZNC[1] = sum;
            }
            if (ich == IdZNC2) {
              mRec.energyZNC[2] = sum;
            }
            if (ich == IdZNC3) {
              mRec.energyZNC[3] = sum;
            }
            if (ich == IdZNC4) {
              mRec.energyZNC[4] = sum;
            }
            if (ich == IdZNCSum) {
              mRec.energyZNC[5] = sum;
            }
            if (ich == IdZPCC) {
              mRec.energyZPC[0] = sum;
            }
            if (ich == IdZPC1) {
              mRec.energyZPC[1] = sum;
            }
            if (ich == IdZPC2) {
              mRec.energyZPC[2] = sum;
            }
            if (ich == IdZPC3) {
              mRec.energyZPC[3] = sum;
            }
            if (ich == IdZPC4) {
              mRec.energyZPC[4] = sum;
            }
            if (ich == IdZPCSum) {
              mRec.energyZPC[5] = sum;
            }
          }
        }
      }
    }
    // TODO: energy calibration
    mTDbg->Fill();
  }
}

void RawReconstruction::interpolate(int itdc, int ibeg, int iend)
{
  // TODO: the bunch 3563 (last in orbit) should be re-reconstructed if it is adjacent to b.c. 0
  constexpr int MaxTimeBin = NTimeBinsPerBC - 1;  //< number of samples per BC
  constexpr int tsnh = TSN / 2;                   // Half number of points in interpolation
  constexpr int nsbun = TSN * NTimeBinsPerBC;     // Total number of interpolated points per bunch crossing
  Int_t nbun = iend - ibeg + 1;                   // Number of adjacent bunches
  Int_t nsam = nbun * NTimeBinsPerBC;             // Number of acquired samples
  Int_t ntot = nsam * TSN;                        // Total number of points in the interpolated arrays
  Int_t nint = (nbun * NTimeBinsPerBC - 1) * TSN; // Total points in the interpolation region (-1)
  constexpr int nsp = 5;                          // Number of points to be searched

  auto& ropt = ZDCRecoParam::Instance();
  Int_t imod = ropt.tmod[itdc]; // Module corresponding to TDC channel
  Int_t ich = ropt.tch[itdc];   // Hardware channel corresponding to TDC channel

  // Check if the TDC channel is connected
  if (imod >= 0 && ich >= 0) {
    Double_t first_sample = mData[ibeg].s[imod][ich][0];
    Double_t last_sample = mData[iend].s[imod][ich][NTimeBinsPerBC - 1];
    // Constant extrapolation at the beginning and at the end of the array
    // Assign value of first sample
    for (Int_t i = 0; i < tsnh; i++) {
      mReco[ibeg].inter[itdc][i] = first_sample;
    }
    // Assign value of last sample
    for (Int_t i = ntot - tsnh; i < ntot; i++) {
      Int_t isam = i % nsbun;
      mReco[iend].inter[itdc][isam] = last_sample;
    }
    // Interpolation between acquired points (n.b. loop from 0 to nint)
    for (Int_t i = 0; i < nint; i++) {
      // Identification of the point to be assigned (need to add tsnh to identify the point)
      Int_t ibun = ibeg + (i + tsnh) / nsbun;
      Int_t isam = (i + tsnh) % nsbun;
      Int_t im = i % TSN;
      if (im == 0) {
        // This is an acquired point
        Int_t ip = (i / TSN) % NTimeBinsPerBC;
        Int_t ib = ibeg + (i / TSN) / NTimeBinsPerBC;
        if (ib != ibun) {
          LOG(FATAL) << "ib=" << ib << " ibun=" << ibun;
          return;
        }
        mReco[ibun].inter[itdc][isam] = mData[ib].s[imod][ich][ip];
      } else {
        // Do the actual interpolation
        Double_t y = 0;
        Int_t ip = i / TSN;
        Double_t sum = 0;
        for (Int_t is = TSN - im, ii = ip - TSL + 1; is < NTS; is += TSN, ii++) {
          // Default is first point in the array
          Double_t yy = first_sample;
          if (ii > 0) {
            if (ii < nsam) {
              Int_t ip = ii % NTimeBinsPerBC;
              Int_t ib = ibeg + ii / NTimeBinsPerBC;
              yy = mData[ib].s[imod][ich][ip];
            } else {
              // Last acquired point
              yy = last_sample;
            }
          }
          sum += mTS[is];
          y += yy * mTS[is];
        }
        y = y / sum;
        mReco[ibun].inter[itdc][isam] = y;
      }
    }
    // Looking for a local maximum in a searching zone
    float amp = std::numeric_limits<float>::infinity(); // Amplitude to be stored
    Int_t isam_amp = 0;                                 // Sample at maximum amplitude (relative to beginning of group)
    Int_t ip_old = -1, ip_cur = -1, ib_cur = -1;        // Current and old points
    bool is_searchable = false;                         // Flag for point in the search zone for maximum amplitude
    bool was_searchable = false;                        // Flag for point in the search zone for maximum amplitude
    Int_t ib[nsp] = {-1, -1, -1, -1, -1};
    Int_t ip[nsp] = {-1, -1, -1, -1, -1};
    // Points at the extremes are constant therefore no local maximum
    // can occur in these two regions
    for (Int_t ibun = ibeg; ibun <= iend; ibun++) {
      printf(" ");
      for (Int_t isam = 0; isam < NTimeBinsPerBC; isam++) {
        printf("%d", mReco[ibun].fired[itdc][isam]);
      }
    }
    printf("\n");
    for (Int_t i = 0; i < nint; i++) {
      Int_t isam = i + tsnh;
      // Check if trigger is fired for this point
      // For the moment we don't take into account possible extensions of the search zone
      // ip_cur can span several bunches and is used just to identify transitions
      ip_cur = isam / TSN;
      // Speed up computation
      if (ip_cur != ip_old) {
        ip_old = ip_cur;
        for (Int_t j = 0; j < 5; j++) {
          ib[j] = -1;
          ip[j] = -1;
        }
        // There are three possible triple conditions that involve current point (middle is current point)
        ip[2] = ip_cur % NTimeBinsPerBC;
        ib[2] = ibeg + ip_cur / NTimeBinsPerBC;
        ib_cur = ib[2];
        if (ip[2] > 0) {
          ip[1] = ip[2] - 1;
          ib[1] = ib[2];
        } else if (ip[2] == 0) {
          if (ib[2] > ibeg) {
            ib[1] = ib[2] - 1;
            ip[1] = MaxTimeBin;
          }
        }
        if (ip[1] > 0) {
          ip[0] = ip[1] - 1;
          ib[0] = ib[1];
        } else if (ip[1] == 0) {
          if (ib[1] > ibeg) {
            ib[0] = ib[1] - 1;
            ip[0] = MaxTimeBin;
          }
        }
        if (ip[2] < MaxTimeBin) {
          ip[3] = ip[2] + 1;
          ib[3] = ib[2];
        } else if (ip[2] == MaxTimeBin) {
          if (ib[2] < iend) {
            ib[3] = ib[2] + 1;
            ip[3] = 0;
          }
        }
        if (ip[3] < MaxTimeBin) {
          ip[4] = ip[3] + 1;
          ib[4] = ib[3];
        } else if (ip[3] == MaxTimeBin) {
          if (ib[3] < iend) {
            ib[4] = ib[3] + 1;
            ip[4] = 0;
          }
        }
        // meet the threshold condition
        was_searchable = is_searchable;
        // Search conditions with list of allowed patterns
	// No need to double check ib[?] and ip[?] because either we assign both or none
	uint16_t triggered=0x0000;
        for (Int_t j = 0; j < 5; j++) {
          if (ib[j] >= 0 && mReco[ib[j]].fired[itdc][ip[j]] > 0) {
            triggered |= (0x1<<j);
          }
        }
        // Reject conditions:
        // 00000
        // 10001
	// One among 10000 and 00001
        // Accept conditions:
        constexpr uint16_t accept[14] = {
//          0x01, // 00001 extend search zone before maximum
          0x02, // 00010
          0x04, // 00100
          0x08, // 01000
          0x10, // 10000 extend after
          0x03, // 00011
          0x06, // 00110
          0x0c, // 01100
          0x18, // 11000
          0x07, // 00111
          0x0e, // 01110
          0x1c, // 11100
          0x0f, // 01111
          0x1e, // 11110
          0x1f  // 11111
	  };
        // All other are not correct (->reject)
        is_searchable = 0;
	if(triggered!=0){
	  for(Int_t j=0; j<14; j++){
	    if(triggered==accept[j]){
              is_searchable = 1;
	      break;
	    }
	  }
	}

//         // Search condition with logical conditions
// 	// No need to double check ib[?] and ip[?] because either we assign both or none
//         if (ib[2] >= 0 && mReco[ib[2]].fired[itdc][ip[2]] > 0) {
//           is_searchable = 1;
//         } else {
//           is_searchable = 0;
//           bool triggered[5] = {false, false, false, false, false};
//           for (Int_t j = 0; j < 5; j++) {
//             // No need to double check ib[?] and ip[?] because either we assign both or none
//             if (ib[j] >= 0 && mReco[ib[j]].fired[itdc][ip[j]] > 0) {
//               triggered[j] = 1;
//             }
//           }
//           // Removed (triggered[0]&&(!triggered[1])&&(!triggered[3])&&triggered[4])
//           // This accepts several bit sequences that are not possible
//           if ((triggered[1] && (!triggered[3])) || ((!triggered[0]) && (!triggered[1]) && (triggered[3] || triggered[4]))) {
//             is_searchable = 1;
//           } else {
//             is_searchable = 0;
//           }
//         }
      }
      // We do not restrict search zone around expected main-main collision
      // because we would like to be able to identify pile-up from collisions
      // with satellites (buggy)

      // If we exit from searching zone
      if (was_searchable && !is_searchable) {
        if (amp <= ADCMax) {
          // Store identified peak
          Int_t ibun = ibeg + isam_amp / nsbun;
          Int_t tdc = isam_amp % nsbun;
          amp = mOffset[imod][ich] - amp;
          assignTDC(ibun, ibeg, iend, itdc, tdc, amp);
        }
        amp = std::numeric_limits<float>::infinity();
        isam_amp = 0;
        was_searchable = 0;
      }
      if (is_searchable) {
        Int_t mysam = isam % nsbun;
        if (mReco[ib_cur].inter[itdc][mysam] < amp) {
          amp = mReco[ib_cur].inter[itdc][mysam];
          isam_amp = isam;
        }
      }
    }
    // Trigger flag still present at the of the scan
    if (is_searchable) {
      // Add last identified peak
      if (amp <= ADCMax) {
        // Store identified peak
        Int_t ibun = ibeg + isam_amp / nsbun;
        Int_t tdc = isam_amp % nsbun;
        float famp = mOffset[imod][ich] - Float_t(amp);
        assignTDC(ibun, ibeg, iend, itdc, tdc, famp);
      }
    }
  } // Connected TDC
}

inline void RawReconstruction::assignTDC(int ibun, int ibeg, int iend, int itdc, int tdc, float amp)
{
  constexpr int nsbun = TSN * NTimeBinsPerBC; // Total number of interpolated points per bunch crossing
  constexpr int tdc_max = nsbun / 2;
  constexpr int tdc_min = -tdc_max;

  // Apply tdc shift correction
  Int_t tdc_cor = tdc - tdc_shift[itdc];
  // Correct bunch assignment
  if (tdc_cor < tdc_min && ibun >= ibeg) {
    // Assign to preceding bunch
    ibun = ibun - 1;
    tdc_cor = tdc_cor + nsbun;
  } else if (tdc_cor >= tdc_max && ibun < iend) {
    // Assign to following bunch
    ibun = ibun + 1;
    tdc_cor = tdc_cor - nsbun;
  }
  if (tdc_cor < kMinShort) {
    LOG(ERROR) << "TDC " << itdc << " " << tdc_cor << " is out of range";
    tdc_cor = kMinShort;
  }
  if (tdc_cor > kMaxShort) {
    LOG(ERROR) << "TDC " << itdc << " " << tdc_cor << " is out of range";
    tdc_cor = kMaxShort;
  }
  // Assign to correct bunch
  mReco[ibun].tdcChannels[itdc].push_back(tdc_cor);
  mReco[ibun].tdcAmplitudes[itdc].push_back(amp);
  LOG(INFO) << mReco[ibun].ir.orbit << "." << mReco[ibun].ir.bc << " "
            << "ibun=" << ibun << " itdc=" << itdc << " tdc=" << tdc << " tdc_cor=" << tdc_cor << " amp=" << amp;
}
