// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "AnalysisDataModel/PID/PIDResponse.h"
#include "AnalysisDataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TOFSpectraTask {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  static constexpr std::string_view hp[Np] = {"p/El", "p/Mu", "p/Pi", "p/Ka", "p/Pr", "p/De", "p/Tr", "p/He", "p/Al"};
  static constexpr std::string_view hpt[Np] = {"pt/El", "pt/Mu", "pt/Pi", "pt/Ka", "pt/Pr", "pt/De", "pt/Tr", "pt/He", "pt/Al"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    histos.add("p/Unselected", "Unselected;#it{p} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    histos.add("pt/Unselected", "Unselected;#it{p}_{T} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    for (int i = 0; i < Np; i++) {
      histos.add(hp[i].data(), Form("%s;#it{p} (GeV/#it{c})", pT[i]), kTH1F, {{100, 0, 20}});
      histos.add(hpt[i].data(), Form("%s;#it{p}_{T} (GeV/#it{c})", pT[i]), kTH1F, {{100, 0, 20}});
    }
    histos.add("electronbeta/hp_El", ";#it{p} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    histos.add("electronbeta/hpt_El", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {{100, 0, 20}});
    histos.add("electronbeta/hlength_El", ";Track Length (cm);Tracks", kTH1D, {{100, 0, 1000}});
    histos.add("electronbeta/htime_El", ";TOF Time (ns);Tracks", kTH1D, {{1000, 0, 600}});
    histos.add("electronbeta/hp_beta_El", ";#it{p} (GeV/#it{c});#beta - #beta_{e};Tracks", kTH2D, {{100, 0, 20}, {100, -0.01, 0.01}});
    histos.add("electronbeta/hp_betasigma_El", ";#it{p} (GeV/#it{c});(#beta - #beta_{e})/#sigma;Tracks", kTH2D, {{100, 0, 20}, {100, -5, 5}});
  }

  template <std::size_t i, typename T>
  void fillParticleHistos(const T& track, const float nsigma[])
  {
    if (abs(nsigma[i]) > nsigmacut.value) {
      return;
    }
    histos.fill(HIST(hp[i]), track.p());
    histos.fill(HIST(hpt[i]), track.pt());
  }

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacut{"nsigmacut", 3, "Value of the Nsigma cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::isGlobalTrack == (uint8_t) true) && (aod::track::tofSignal > 0.f);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::pidRespTOF, aod::pidRespTOFbeta, aod::TrackSelection>>;
  void process(TrackCandidates::iterator const& track)
  {
    const float nsigma[Np] = {track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(),
                              track.tofNSigmaKa(), track.tofNSigmaPr(), track.tofNSigmaDe(),
                              track.tofNSigmaTr(), track.tofNSigmaHe(), track.tofNSigmaAl()};
    histos.fill(HIST("p/Unselected"), track.p());
    histos.fill(HIST("pt/Unselected"), track.pt());

    fillParticleHistos<0>(track, nsigma);
    fillParticleHistos<1>(track, nsigma);
    fillParticleHistos<2>(track, nsigma);
    fillParticleHistos<3>(track, nsigma);
    fillParticleHistos<4>(track, nsigma);
    fillParticleHistos<5>(track, nsigma);
    fillParticleHistos<6>(track, nsigma);
    fillParticleHistos<7>(track, nsigma);
    fillParticleHistos<8>(track, nsigma);

    //
    if (TMath::Abs(track.separationbetael() < 1.f)) {
      histos.fill(HIST("electronbeta/hp_El"), track.p());
      histos.fill(HIST("electronbeta/hpt_El"), track.pt());
      histos.fill(HIST("electronbeta/hlength_El"), track.length());
      histos.fill(HIST("electronbeta/htime_El"), track.tofSignal() / 1000);
      histos.fill(HIST("electronbeta/hp_beta_El"), track.p(), track.diffbetael());
      histos.fill(HIST("electronbeta/hp_betasigma_El"), track.p(), track.separationbetael());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  WorkflowSpec workflow{adaptAnalysisTask<TOFSpectraTask>("tofspectra-task")};
  return workflow;
}
