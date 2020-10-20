#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMath.h>
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "ZDCBase/Constants.h"
#include "ZDCBase/ModuleConfig.h"
#include "ZDCSimulation/Digitizer.h"
#include "DataFormatsZDC/BCData.h"
#include "DataFormatsZDC/ChannelData.h"
#include "ZDCSimulation/Digits2Raw.h"
#include "ZDCSimulation/MCLabel.h"
#include "ZDCSimulation/ZDCSimParam.h"
#include "FairLogger.h"

using namespace o2::zdc;

//ClassImp(Digits2Raw);
void Digits2Raw::processDigits(const std::string& outDir, const std::string& fileDigitsName)
{
  auto& sopt = ZDCSimParam::Instance();
  mIsContinuous = sopt.continuous;

  if (!mModuleConfig) {
    LOG(ERROR) << "Missing ModuleConfig configuration object";
    return;
  }

  if (!mSimCondition) {
    LOG(ERROR) << "Missing SimCondition configuration object";
    return;
  }

  if(mNEmpty<0){
    LOG(ERROR) << "Bunch crossing map is not initialized";
    return;
  }

  if(mNEmpty==0){
    LOG(WARNING) << "Bunch crossing map has zero clean empty bunches";
  }

  setTriggerMask();

  std::string outd = outDir;
  if (outd.back() != '/') {
    outd += '/';
  }

  std::unique_ptr<TFile> digiFile(TFile::Open(fileDigitsName.c_str()));
  if (!digiFile || digiFile->IsZombie()) {
    LOG(ERROR) << "Failed to open input digits file " << fileDigitsName;
    return;
  }

  TTree* digiTree = (TTree*)digiFile->Get("o2sim");
  if (!digiTree) {
    LOG(ERROR) << "Failed to get digits tree";
    return;
  }

  o2::dataformats::MCTruthContainer<o2::zdc::MCLabel>* labelsPtr = nullptr;

  if (digiTree->GetBranch("ZDCDigitBC")) {
    digiTree->SetBranchAddress("ZDCDigitBC", &mzdcBCDataPtr);
  } else {
    LOG(ERROR) << "Branch ZDCDigitBC is missing";
    return;
  }

  if (digiTree->GetBranch("ZDCDigitCh")) {
    digiTree->SetBranchAddress("ZDCDigitCh", &mzdcChDataPtr);
  } else {
    LOG(ERROR) << "Branch ZDCDigitCh is missing";
    return;
  }

  if (digiTree->GetBranch("ZDCDigitLabels")) {
    digiTree->SetBranchAddress("ZDCDigitLabels", &labelsPtr);
    LOG(INFO) << "Branch ZDCDigitLabels is connected";
  } else {
    LOG(INFO) << "Branch ZDCDigitLabels is missing";
  }

  for (int ient = 0; ient < digiTree->GetEntries(); ient++) {
    digiTree->GetEntry(ient);
    int nbc = mzdcBCData.size();
    LOG(INFO) << "Entry " << ient << " : " << nbc << " BCs stored";
    for (int ibc = 0; ibc < nbc; ibc++) {
      // Detect orbit change and insert last bunch if needed
      convertDigits(ibc);
      writeDigits();
    }
    // If last event is not last bunch insert it
  }
  digiFile->Close();
}

void Digits2Raw::setTriggerMask()
{
  mTriggerMask = 0;
  mPrintTriggerMask = "";
  for (Int_t im = 0; im < NModules; im++) {
    if (im > 0)
      mPrintTriggerMask += " ";
    mPrintTriggerMask += std::to_string(im);
    mPrintTriggerMask += "[";
    for (UInt_t ic = 0; ic < NChPerModule; ic++) {
      if (mModuleConfig->modules[im].trigChannel[ic]) {
        UInt_t tmask = 0x1 << (im * NChPerModule + ic);
        mTriggerMask = mTriggerMask | tmask;
        mPrintTriggerMask += "T";
      } else {
        mPrintTriggerMask += " ";
      }
    }
    mPrintTriggerMask += "]";
    UInt_t mytmask = mTriggerMask >> (im * NChPerModule);
    printf("Trigger mask for module %d 0123 %s%s%s%s\n", im,
           mytmask & 0x1 ? "T" : "N",
           mytmask & 0x2 ? "T" : "N",
           mytmask & 0x4 ? "T" : "N",
           mytmask & 0x8 ? "T" : "N");
  }
  printf("trigger_mask=0x%08x %s\n", mTriggerMask, mPrintTriggerMask.data());
}

// void Digits2Raw::insertLastBunch(int ibc){
//   // Number of bunch crossings stored
//   int nbc = mzdcBCData.size();
//   
// }

void Digits2Raw::convertDigits(int ibc)
{
  // Number of bunch crossings stored
  int nbc = mzdcBCData.size();

  // Orbit and bunch crossing identifiers
  const auto& bcd = mzdcBCData[ibc];
  UShort_t bc = bcd.ir.bc;
  UInt_t orbit = bcd.ir.orbit;

  // Reset scalers at orbit change
  if (orbit != mLastOrbit) {
    for (Int_t im = 0; im < NModules; im++) {
      for (Int_t ic = 0; ic < NChPerModule; ic++){
        mScalers[im][ic] = 0;
	mSumPed[im][ic] = 0;
	mPed[im][ic] = 0;
      }
    }
    mLastOrbit = orbit;
    mLastNEmpty = 0;
  }

  // Compute or update baseline reference
  if(mEmpty[bc]>0 && mEmpty[bc]!=mLastNEmpty){
    for (Int_t im = 0; im < NModules; im++) {
      for (Int_t ic = 0; ic < NChPerModule; ic++){
        // Identify connected channel
	auto id=mModuleConfig->modules[im].channelID[ic];
	auto base_m=mSimCondition->channels[id].pedestal;      // Average pedestal
	auto base_s=mSimCondition->channels[id].pedestalFluct; // Baseline oscillations
	auto base_n=mSimCondition->channels[id].pedestalNoise; // Electronic noise
	Double_t deltan=mEmpty[bc]-mLastNEmpty;
	// We assume to have a fluctuation every two bunch crossings
	// Will need to tune this parameter
	mSumPed[im][ic]+=12.*deltan*gRandom->Gaus(base_m,base_s*TMath::Sqrt(deltan/2.)); 
	mSumPed[im][ic]+=12.*deltan*gRandom->Gaus(0,base_n*TMath::Sqrt(12.*deltan));
	Double_t myped=TMath::Nint(8.*mSumPed[im][ic]/Double_t(mEmpty[bc])/12.+32768);
	if(myped<0)myped=0;
	if(myped>65535)myped=65535;
	mPed[im][ic]=myped;
      }
    }
    mLastNEmpty=mEmpty[bc];
  }

  // Increment scalers and reset output structure
  for (UInt_t im = 0; im < NModules; im++) {
    for (UInt_t ic = 0; ic < NChPerModule; ic++) {
      // Fixed words
      mZDC.data[im][ic].w[0][0] = id_w0;
      mZDC.data[im][ic].w[0][1] = 0;
      mZDC.data[im][ic].w[0][2] = 0;
      mZDC.data[im][ic].w[1][0] = id_w1;
      mZDC.data[im][ic].w[1][1] = 0;
      mZDC.data[im][ic].w[1][2] = 0;
      mZDC.data[im][ic].w[2][0] = id_w2;
      mZDC.data[im][ic].w[2][1] = 0;
      mZDC.data[im][ic].w[2][2] = 0;
      // Module and channel numbers
      mZDC.data[im][ic].f.board = im;
      mZDC.data[im][ic].f.ch = ic;
      // Orbit and bunch crossing
      mZDC.data[im][ic].f.orbit = orbit;
      mZDC.data[im][ic].f.bc = bc;
      // If channel is hit in current bunch crossing
      if (bcd.triggers & (0x1 << (im * NChPerModule + ic))) {
        mScalers[im][ic]++;          // increment scalers
        mZDC.data[im][ic].f.Hit = 1; // flag bunch crossing
      }
      mZDC.data[im][ic].f.hits = mScalers[im][ic];
      mZDC.data[im][ic].f.offset = mPed[im][ic];
    }
  }

  // Compute autotrigger bits and assigna ALICE trigger bits
  // triggers refer to the HW trigger conditions (32 possible channels)

  // Autotrigger, current bunch crossing
  UInt_t triggers_0 = bcd.triggers;

  // ALICE current bunch crossing
  if (bcd.ext_triggers)
    for (UInt_t im = 0; im < NModules; im++)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Alice_0 = 1;

  // Next bunch crossings (ALICE and autotrigger)
  UInt_t triggers_1 = 0, triggers_2 = 0, triggers_3 = 0, triggers_m = 0;
  for (Int_t is = 1; is < 4; is++) {
    Int_t ibc_peek = ibc + is;
    if (ibc_peek >= nbc)
      break;
    const auto& bcd_peek = mzdcBCData[ibc_peek];
    UShort_t bc_peek = bcd_peek.ir.bc;
    UInt_t orbit_peek = bcd_peek.ir.orbit;
    if (bcd_peek.triggers) {
      if (orbit_peek == orbit) {
        if ((bc_peek - bc) == 1) {
          triggers_1 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_1 = 1;
        } else if ((bc_peek - bc) == 2) {
          triggers_2 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_2 = 1;
        } else if ((bc_peek - bc) == 3) {
          triggers_3 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_3 = 1;
          break;
        }
      } else if (orbit_peek == (orbit + 1)) {
        if ((bc_peek + 3564 - bc) == 1) {
          triggers_1 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_1 = 1;
        } else if ((bc_peek + 3564 - bc) == 2) {
          triggers_2 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_2 = 1;
        } else if ((bc_peek + 3564 - bc) == 3) {
          triggers_3 = bcd_peek.triggers;
          if (bcd_peek.ext_triggers)
            for (UInt_t im = 0; im < NModules; im++)
              for (UInt_t ic = 0; ic < NChPerModule; ic++)
                mZDC.data[im][ic].f.Alice_3 = 1;
          break;
        }
      } else {
        break;
      }
    }
  }

  // Previous bunch crossing just for autotrigger
  {
    Int_t ibc_peek = ibc - 1;
    if (ibc_peek >= 0) {
      const auto& bcd_peek = mzdcBCData[ibc - 1];
      UShort_t bc_peek = bcd_peek.ir.bc;
      UInt_t orbit_peek = bcd_peek.ir.orbit;
      if (bcd_peek.triggers) {
        if (orbit_peek == orbit) {
          if ((bc - bc_peek) == 1)
            triggers_m = bcd_peek.triggers;
        } else if (orbit_peek == (orbit - 1)) {
          if (bc == 0 && bc_peek == 3563)
            triggers_m = bcd_peek.triggers;
        }
      }
    }
  }

  // Assign trigger bits in payload
  for (Int_t im = 0; im < NModules; im++) {
    UInt_t tmask = (0xf << (im * NChPerModule)) & mTriggerMask;
    if (triggers_m & tmask)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Auto_m = 1;
    if (triggers_0 & tmask)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Auto_0 = 1;
    if (triggers_1 & tmask)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Auto_1 = 1;
    if (triggers_2 & tmask)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Auto_2 = 1;
    if (triggers_3 & tmask)
      for (UInt_t ic = 0; ic < NChPerModule; ic++)
        mZDC.data[im][ic].f.Auto_3 = 1;
  }

  bcd.print();
  printf("Mask: %s\n", mPrintTriggerMask.data());
  int chEnt = bcd.ref.getFirstEntry();
  for (int ic = 0; ic < bcd.ref.getEntries(); ic++) {
    const auto& chd = mzdcChData[chEnt++];
    chd.print();
    UShort_t bc = bcd.ir.bc;
    UInt_t orbit = bcd.ir.orbit;
    for (Int_t is = 0; is < o2::zdc::NTimeBinsPerBC; is++) {
      int16_t myval = chd.data[is];
    }
    // Look for channel ID in
    for (Int_t im = 0; im < NModules; im++) {
      for (UInt_t ic = 0; ic < NChPerModule; ic++) {
        if (mModuleConfig->modules[im].channelID[ic] == chd.id &&
            mModuleConfig->modules[im].readChannel[ic]) {
          Int_t is = 0;
          mZDC.data[im][ic].f.s00 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s01 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s02 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s03 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s04 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s05 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s06 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s07 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s08 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s09 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s10 = chd.data[is];
          is++;
          mZDC.data[im][ic].f.s11 = chd.data[is];
          break;
        }
      }
    }
  }
}

void Digits2Raw::writeDigits()
{
  // TODO baseline estimation
  // TODO treatment of last bunch crossing of each orbit

  for (UInt_t im = 0; im < o2::zdc::NModules; im++) {
    // Check if module has been filled with data
    // N.B. All channels are initialized if module is supposed to be readout
    // Trigger bits are the same for all the channels connected to a module
    bool TM = mZDC.data[im][0].f.Auto_m;
    bool T0 = mZDC.data[im][0].f.Auto_0;
    bool T1 = mZDC.data[im][0].f.Auto_1;
    bool T2 = mZDC.data[im][0].f.Auto_2;
    bool T3 = mZDC.data[im][0].f.Auto_3;
    bool A0 = mZDC.data[im][0].f.Alice_0;
    bool A1 = mZDC.data[im][0].f.Alice_1;
    bool A2 = mZDC.data[im][0].f.Alice_2;
    bool A3 = mZDC.data[im][0].f.Alice_3;
    bool tcond_continuous = T0 || T1;
    bool tcond_triggered = A0 || A1 || (A2 && (T0 || TM)) || (A3 && T0);
    bool tcond_last = mZDC.data[im][0].f.bc == 3563;
    if (mVerbosity > 1) {
      if (tcond_continuous)
        printf("M%d %s T0=%d || T1=%d -> %d\n", im, mIsContinuous ? "Continuous mode " : "Triggered mode  ", T0, T1, tcond_continuous);
      if (tcond_triggered)
        printf("M%d A0=%d || A1=%d || (A2=%d && (T0=%d || TM=%d))=%d || (A3=%d && T0=%d )=%d -> %d\n", im, A0, A1, A2, T0, TM, A2 && (T0 || TM), A3, T0, A3 && T0, tcond_triggered);
      if (mZDC.data[im][0].f.bc == 3563)
        printf("M%d Is last BC\n", im);
    }
    // Condition to write GBT data
    if (tcond_triggered || (mIsContinuous && tcond_continuous) || (mZDC.data[im][0].f.bc == 3563)) {
      for (UInt_t ic = 0; ic < o2::zdc::NChPerModule; ic++) {
        if (mModuleConfig->modules[im].readChannel[ic]) {
          for (Int_t iw = 0; iw < o2::zdc::NWPerBc; iw++) {
            if (mVerbosity > 1)
              print_gbt_word(&mZDC.data[im][ic].w[iw][0]);
	    // Check link for channel and insert data
          }
        }
      }
    }
  }
}

void Digits2Raw::print_gbt_word(UInt_t* word)
{
  unsigned __int128 val = word[2];
  val = val << 32;
  val = val | word[1];
  val = val << 32;
  val = val | word[0];
  static UInt_t last_orbit = 0, last_bc = 0;
  // First word
  struct zdc_payload_v0_w0 {
    unsigned fixed_0 : 2;
    unsigned board : 4;
    unsigned ch : 2;
    unsigned offset : 16;
    unsigned hits : 12;
    unsigned bc : 12;
    UInt_t orbit;
    unsigned empty_0 : 16;
    UInt_t empty_1;
  };
  // Second word
  struct zdc_payload_v0_w1 {
    unsigned fixed_1 : 2;
    unsigned error : 2;
    unsigned Alice_0 : 1;
    unsigned Alice_1 : 1;
    unsigned Alice_2 : 1;
    unsigned Alice_3 : 1;
    unsigned s00 : 12;
    unsigned s01 : 12;
    unsigned s02 : 12;
    unsigned s03 : 12;
    unsigned s04 : 12;
    unsigned s05 : 12;
    unsigned empty_2 : 16;
    UInt_t empty_3;
  };
  // Third word
  struct zdc_payload_v0_w2 {
    unsigned fixed_2 : 2;
    unsigned Hit : 1;
    unsigned Auto_m : 1;
    unsigned Auto_0 : 1;
    unsigned Auto_1 : 1;
    unsigned Auto_2 : 1;
    unsigned Auto_3 : 1;
    unsigned s06 : 12;
    unsigned s07 : 12;
    unsigned s08 : 12;
    unsigned s09 : 12;
    unsigned s10 : 12;
    unsigned s11 : 12;
    unsigned empty_4 : 16;
    UInt_t empty_5;
  };

  ULong64_t lsb = val;
  ULong64_t msb = val >> 64;
  UInt_t a = word[0];
  UInt_t b = word[1];
  UInt_t c = word[2];
  //UInt_t d=(msb>>32)&0xffffffff;
  //printf("\n%llx %llx ",lsb,msb);
  //printf("\n%8x %8x %8x %8x ",d,c,b,a);
  if ((a & 0x3) == 0) {
    union {
      unsigned __int128 w;
      struct zdc_payload_v0_w0 f;
    } w0;
    UInt_t myorbit = (val >> 48) & 0xffffffff;
    UInt_t mybc = (val >> 36) & 0xfff;
    if (myorbit != last_orbit || mybc != last_bc) {
      printf("Orbit %9u BC %4u\n", myorbit, mybc);
      last_orbit = myorbit;
      last_bc = mybc;
    }
    printf("%04x %08x %08x ", c, b, a);
    UInt_t hits = (val >> 24) & 0xfff;
    Int_t offset = (lsb >> 8) & 0xffff-32768;
    Float_t foffset=offset/8.;
    UInt_t board = (lsb >> 2) & 0xf;
    UInt_t ch = (lsb >> 6) & 0x3;
    //printf("orbit %9u bc %4u hits %4u offset %+6i Board %2u Ch %1u", myorbit, mybc, hits, offset, board, ch);
    printf("orbit %9u bc %4u hits %4u offset %+8.3f Board %2u Ch %1u", myorbit, mybc, hits, foffset, board, ch);
  } else if ((a & 0x3) == 1) {
    printf("%04x %08x %08x ", c, b, a);
    union {
      unsigned __int128 w;
      struct zdc_payload_v0_w1 f;
    } w1;
    printf("     %s %s %s %s ", a & 0x10 ? "A0" : "  ", a & 0x20 ? "A1" : "  ", a & 0x40 ? "A2" : "  ", a & 0x80 ? "A3" : "  ");
    printf("0-5 ");
    Short_t s[6];
    val = val >> 8;
    for (Int_t i = 0; i < 6; i++) {
      s[i] = val & 0xfff;
      if (s[i] > 2047)
        s[i] = s[i] - 4096;
      val = val >> 12;
    }
    printf(" %5d %5d %5d %5d %5d %5d", s[0], s[1], s[2], s[3], s[4], s[5]);
  } else if ((a & 0x3) == 2) {
    printf("%04x %08x %08x ", c, b, a);
    union {
      unsigned __int128 w;
      struct zdc_payload_v0_w2 f;
    } w2;
    printf("%s %s %s %s %s %s ", a & 0x4 ? "H" : " ", a & 0x8 ? "TM" : "  ", a & 0x10 ? "T0" : "  ", a & 0x20 ? "T1" : "  ", a & 0x40 ? "T2" : "  ", a & 0x80 ? "T3" : "  ");
    printf("6-b ");
    Short_t s[6];
    val = val >> 8;
    for (Int_t i = 0; i < 6; i++) {
      s[i] = val & 0xfff;
      if (s[i] > 2047)
        s[i] = s[i] - 4096;
      val = val >> 12;
    }
    printf(" %5d %5d %5d %5d %5d %5d", s[0], s[1], s[2], s[3], s[4], s[5]);
  } else if ((a & 0x3) == 3) {
    printf("%04x %08x %08x ", c, b, a);
    printf("HB ");
  }
  printf("\n");
}

void Digits2Raw::emptyBunches(std::bitset<3564>& bunchPattern){
  const int LHCMaxBunches=o2::constants::lhc::LHCMaxBunches;
  mNEmpty=0;
  for(Int_t ib=0; ib<LHCMaxBunches; ib++){
    Int_t mb=(ib+31)%LHCMaxBunches; // beam gas from back of calorimeter
    Int_t m1=(ib+1)%LHCMaxBunches; // previous bunch
    Int_t cb=ib; // current bunch crossing
    Int_t p1=(ib-1)%LHCMaxBunches; // colliding + 1
    Int_t p2=(ib+1)%LHCMaxBunches; // colliding + 2
    Int_t p3=(ib+1)%LHCMaxBunches; // colliding + 3
    if(bunchPattern[mb] || bunchPattern[m1] || bunchPattern[cb] || bunchPattern[p1] || bunchPattern[p2] || bunchPattern[p3]){
      mEmpty[ib]=mNEmpty;
    }else{
      mNEmpty++;
      mEmpty[ib]=mNEmpty;
    }
  }
  LOG(INFO) << "There are " << mNEmpty << " clean empty bunches";
}
