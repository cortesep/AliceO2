#include <TROOT.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveStats.h>
#include "ZDCRaw/RawReconstruction.h"
#include "ZDCSimulation/Digits2Raw.h"
#include "FairLogger.h"
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
  const Int_t tsn=200; // tsn-1 = number of interpolated point between each pair
  const Int_t tsl=6; // number of zeros on the right (and on the left) of central peak
  const Int_t tst=2*tsl*tsn+1;
  if(tst != mNTS){
    LOG(FATAL) << "Error in the definition of interpolating kernel";
    return;
  }
  const Double_t tsc=750;
  // tsc/tsn ~ 4 (3.75) and tsl*tsn*sqrt(2)/tsc >> 1 (n di sigma)
  Int_t n=tsl*tsn;
  for(Int_t tsi=0; tsi<=n; tsi++){
    Double_t arg1=TMath::Pi()*Double_t(tsi)/Double_t(tsn);
    Double_t fs=1;
    if(arg1!=0)fs=TMath::Sin(arg1)/arg1;
    Double_t arg2=Double_t(tsi)/tsc;
    Double_t fg=TMath::Exp(-arg2*arg2);
    mTS[n+tsi]=fs*fg;
    mTS[n-tsi]=mTS[n+tsi]; // Function is even
  }

  // Update reconstruction parameters
  //auto& ropt=ZDCRecoParam::Instance();
  o2::zdc::ZDCRecoParam& ropt = const_cast<o2::zdc::ZDCRecoParam&>(ZDCRecoParam::Instance());
  
  // Fill maps
  for(Int_t itdc=0; itdc<o2::zdc::NTDCChannels; itdc++){
    // If the reconstruction parameters were not manually set
    if(ropt.tmod[itdc]<0 || ropt.tch[itdc]<0){
      Int_t isig=TDCSignal[itdc];
      for (Int_t im = 0; im < NModules; im++) {
	for (UInt_t ic = 0; ic < NChPerModule; ic++) {
	  if (mModuleConfig->modules[im].channelID[ic]==isig && mModuleConfig->modules[im].readChannel[ic]) {
	    //ropt.updateFromString(TString::Format("ZDCRecoParam.tmod[%d]=%d;",itdc,im));
	    //ropt.updateFromString(TString::Format("ZDCRecoParam.tch[%d]=%d;",itdc,ic));
	    ropt.tmod[itdc]=im;
	    ropt.tch[itdc]=ic;
	    goto next_itdc;
	  }
	}
      }
    }
    next_itdc:;
    LOG(INFO) << "TDC " << itdc << " mod " << ropt.tmod[itdc] << " ch " << ropt.tch[itdc];
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
  TFile* f = new TFile("ZDCReco.root", "recreate");
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
  f->Close();
  return processOrbit(true);
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
    auto chd = data.data[board][ch];
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
  auto& reco = mReco.back();
  //LOG(INFO) << "Processed: orbit " << reco.ir.orbit << " bc " << reco.ir.bc << " board " << board << " ch " << ch;
  return 0;
}

void RawReconstruction::printIR(){
  std::deque<RecEvent>::iterator it = mReco.begin();
  while (it != mReco.end()) {
    LOG(INFO) << (*it).ir.orbit << "." << (*it).ir.bc;
    it++;
  }
}

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

  // Find consecutive bc data belonging to current orbit
  std::deque<RecEvent>::iterator it = mReco.begin();
  auto *bcd_last = &((*it).ir);
  Int_t seq_beg=0;
  Int_t seq_end=0;
  Int_t iti=0;
  for (; it != mReco.end(); it++, iti++){
    auto *bcd = &((*it).ir);
    int64_t bcdiff=bcd->differenceInBC(*bcd_last);
    // Find a gap
    if(bcdiff>1){
      // Process consecutive BCs
      if(seq_beg==seq_end){
	LOG(INFO) << "Lonely bunch " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc;
      }else{
	LOG(INFO) << "Processing " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc << " - " << mReco[seq_end].ir.orbit << "." << mReco[seq_end].ir.bc;
	reconstruct(seq_beg,seq_end);
      }
      seq_beg=iti;
      seq_end=iti;
    }else{
      seq_end=iti;
    }
    bcd_last=bcd;
  }
  if(seq_beg != seq_end){
    // No recipe to process a lonely bunch
    LOG(INFO) << "Processing " << mReco[seq_beg].ir.orbit << "." << mReco[seq_beg].ir.bc << " - " << mReco[seq_end].ir.orbit << "." << mReco[seq_end].ir.bc;
    reconstruct(seq_beg,seq_end);
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

void RawReconstruction::reconstruct(int ibeg, int iend)
{
  auto& ropt=ZDCRecoParam::Instance();
  // Apply differential discrimination with triple condition
  Int_t nbun=iend-ibeg+1;
  Int_t maxs2=NTimeBinsPerBC*nbun-1;
  for(Int_t itdc=0; itdc<NTDCChannels; itdc++){
    Int_t im=ropt.tmod[itdc];
    Int_t ic=ropt.tch[itdc];
    Int_t shift = ropt.tsh[itdc];
    Int_t thr=ropt.tth[itdc];
    // Check if the TDC channel is connected
    if(im>=0 && ic>=0){
      Int_t is1=0, is2=1;
      Int_t isfired[3]={0};
      for(;;){
	// Shift data
	for(Int_t i=1; i<3; i++)isfired[i]=isfired[i-1];
	Int_t b1=ibeg+is1/NTimeBinsPerBC;
	Int_t b2=ibeg+is2/NTimeBinsPerBC;
	Int_t s1=is1%NTimeBinsPerBC;
	Int_t s2=is2%NTimeBinsPerBC;
	int diff=mData[b1].s[im][ic][s1]-mData[b2].s[im][ic][s2];
	if(diff>thr){
	  isfired[0]=1;
	  if(isfired[1]==1&&isfired[2]==1){
	    mReco[b2].fired[itdc][s2]=1;
	    LOG(INFO) << "Fired @ " << b2-ibeg << "." << s2;
	  }
	}
	if(is2>=shift)is1++;
	if(is2<maxs2)is2++;
	if(is1==maxs2)break;
      }
    }
    break;
  }
}

void RawReconstruction::replayTrigger(uint32_t orbit)
{
  auto& ropt=ZDCRecoParam::Instance();
  auto& last = mReco.back();
  auto last_orbit = last.ir.orbit;
  auto last_bc = last.ir.bc;
  auto& first = mReco.front();
  auto first_orbit = first.ir.orbit;
  auto first_bc = first.ir.bc;
  LOG(INFO) << "Processing from " << first_orbit << "." << first_bc << " to " << last_orbit << "." << last_bc;

  
  // Apply differential discrimination with triple condition
//   for(Int_t itdc=0; itdc<NTDCChannels; itdc++){
//     Int_t im=ropt.tmod[itdc];
//     Int_t ic=ropt.tch[itdc];
//     // Check if the TDC channel is connected
//     if(im>=0 && ic>=0){
//     }
//   }
  /*
  Int_t maxa = ns - 2 - shift;
  Int_t maxf = ns - shift;
  Int_t isarmed = 0;
  Double_t* y = gr[tch[itch]]->GetY();
  for (Int_t is = 0; is < maxf; is++) {
    if (isarmed == 0 && is < maxa) {
      if ((y[0 + is] - y[shift + is]) > thr && (y[1 + is] - y[1 + shift + is]) > thr && (y[2 + is] - y[2 + shift + is]) > thr) {
        printf("%s armed @ is=%d [%g%+g=%g][%g%+g=%g][%g%+g=%g]\n", o2::zdc::channelName(tch[itch]), is, y[0 + is], -y[shift + 0 + is], y[0 + is] - y[shift + 0 + is],
               y[1 + is], -y[shift + 1 + is], y[1 + is] - y[shift + 1 + is], y[2 + is], -y[shift + 2 + is], y[2 + is] - y[shift + 2 + is]);
        isarmed = 1;
        narmed[itch]++;
        trig[itch][is] = 1;
      }
    } else if (isarmed == 1) {
      if ((y[0 + is] - y[shift + is]) <= 0) {
        isarmed = 0;
        nfired[itch]++;
        Double_t sfired = ((shift + is) + (0 + is)) / 2.;
        Int_t isfired = TMath::Nint(sfired);
        printf("%s fired @ is=%d ymax @ %g=%g [%g%+g=%g]\n", o2::zdc::channelName(tch[itch]), is, sfired, y[isfired], y[0 + is], -y[shift + is], y[0 + is] - y[shift + is]);
        trig[itch][isfired] = 2;
      }
    }
  }
  */
}
