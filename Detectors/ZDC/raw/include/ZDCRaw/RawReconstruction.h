#include <deque>
#include <TH1.h>
#include <TH2.h>
#include "ZDCBase/Constants.h"
#include "ZDCSimulation/ZDCSimParam.h"
#include "ZDCRaw/ZDCRecoParam.h"
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

 private:
  void setStat(TH1* h);
  int mVerbosity = 1;
  TH1* mBaseline[NDigiChannels] = {0};
  TH1* mCounts[NDigiChannels] = {0};
  TH2* mSignal[NDigiChannels] = {0};
  TH2* mBunch[NDigiChannels] = {0};
  EventChData mCh;           // Temporary store of channel data
  bool mIsContinuous = true; // continuous (self-triggered) or externally-triggered readout
  int mNBCAHead = 0;         // when storing triggered BC, store also mNBCAHead BCs
  bool mFirst = 1;
  uint32_t mPrevOrbit;
  uint32_t mCurOrbit;
  uint32_t mEvOrbit;
  std::deque<EventData> mData;
  std::deque<RecEvent> mReco;
  void printIR(); // Dump interaction record
  const ModuleConfig* mModuleConfig = nullptr;                          /// Trigger/readout configuration object
  void reconstruct(int ibeg, int iend);
  void replayTrigger(uint32_t orbit);
  const static int mNTS=2401; /// Tapered sinc function array size
  Double_t mTS[mNTS]; /// Tapered sinc function

  ClassDefNV(RawReconstruction, 1);
};
} // namespace zdc
} // namespace o2

#endif
