#include <deque>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include "ZDCBase/Constants.h"
#include "ZDCSimulation/ZDCSimParam.h"
#include "ZDCRaw/ZDCRecoParam.h"
#include "ZDCRaw/ZDCTDCParam.h"
#include "ZDCRaw/ZDCIntegrationParam.h"
#include "ZDCBase/ModuleConfig.h"
#include "DataFormatsZDC/RawEventData.h"
#include "DataFormatsZDC/RecEvent.h"
#ifndef ALICEO2_ZDC_RAWRECONSTRUCTION_H_
#define ALICEO2_ZDC_RAWRECONSTRUCTION_H_
namespace o2
{
namespace zdc
{
class RawReconstruction
{
 public:
  RawReconstruction() = default;
  void init();
  int process(const EventData ev);
  int processCh();
  int processOrbit(bool islast = false);
  int processWord(const UInt_t* word);
  int reconstructOrbit(UInt_t orbit);
  int getHPos(uint32_t board, uint32_t ch);
  int write();
  void setVerbosity(int v)
  {
    mVerbosity = v;
  }
  int getVerbosity() const { return mVerbosity; }

  void setModuleConfig(const ModuleConfig* moduleConfig) { mModuleConfig = moduleConfig; };
  const ModuleConfig* getModuleConfig() { return mModuleConfig; };
  void setTDCParam(const ZDCTDCParam* param) { mTDCParam = param; };
  const ZDCTDCParam* getTDCParam() { return mTDCParam; };
  void setIntegrationParam(const ZDCIntegrationParam* param) { mIntParam = param; };
  const ZDCIntegrationParam* getIntegrationParam() { return mIntParam; };

 private:
  void setStat(TH1* h);
  int mVerbosity = 1;
  // Debug histograms
  TH1* mBaseline[NDigiChannels] = {0}; /// Average baseline in orbit
  TH1* mCounts[NDigiChannels] = {0};   /// Number of bunches with hit in orbit
  TH2* mSignal[NDigiChannels] = {0};   /// Signal line shape
  TH2* mBunch[NDigiChannels] = {0};    /// Bunches with hit
  EventChData mCh;                     /// Temporary store of channel data
  bool mIsContinuous = true;           /// continuous (self-triggered) or externally-triggered readout
  int mNBCAHead = 0;                   /// when storing triggered BC, store also mNBCAHead BCs
  bool mFirst = 1;
  uint32_t mPrevOrbit;
  uint32_t mCurOrbit;
  uint32_t mEvOrbit;
  std::deque<EventData> mData;
  std::deque<RecEvent> mReco;
  void printIR(); // Dump interaction record

  const ModuleConfig* mModuleConfig = nullptr;    /// Trigger/readout configuration object
  const ZDCTDCParam* mTDCParam = nullptr;         /// TDC calibration object
  const ZDCIntegrationParam* mIntParam = nullptr; /// Configuration of integration

  void reconstruct(int ibeg, int iend);
  void processTrigger(int itdc, int ibeg, int iend);
  void interpolate(int itdc, int ibeg, int iend);
  inline void assignTDC(int ibun, int ibeg, int iend, int itdc, int tdc, float amp);
  Double_t mTS[NTS];                           /// Tapered sinc function
  TFile* mDbg = nullptr;                       /// Debug output
  TTree* mTDbg = nullptr;                      /// Debug tree
  RecEvent mRec;                               /// Debug reconstruction event
  int mLast = -1;                              /// Index of last bunch crossing in orbit
  float mOffset[NModules][NChPerModule] = {0}; /// Offset estimated on FEE for current orbit
  int16_t tdc_shift[NTDCChannels] = {0};       /// TDC correction (units of 1/96 ns)
  ClassDefNV(RawReconstruction, 1);
};
} // namespace zdc
} // namespace o2

#endif
