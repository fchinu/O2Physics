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

/// \file multiplicityCombination.cxx
/// \brief Task for combining multiplicity information from different detectors
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Universita and INFN Torino

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

using namespace o2;
using namespace o2::framework;

enum CentralityDetectors {
  FV0A,
  FT0A,
  FT0C,
  FDDA,
  FDDC,
  NTracksPV,
  NTracksGlobal,
  kNDetectors
};

struct multiplicityCombination {

  Configurable<bool> zEqualize{"zEqualize", false, "Do z equalization"};
  Configurable<std::string> pathZEq{"pathZEq", "zEq.root", "Path to .root file with weights for Z equalization"};
  Configurable<bool> useWeights{"useWeights", false, "Reweights the multiplicities"};
  Configurable<std::string> pathWeights{"pathWeights", "weights.root", "Path to .root file with weights to apply"};

  ConfigurableAxis binsMultFV0A{"binsMultFV0A", {500, 0., 20000}, "bins for FV0A multiplicity"};
  ConfigurableAxis binsMultFT0A{"binsMultFT0A", {500, 0., 7000}, "bins for FT0A multiplicity"};
  ConfigurableAxis binsMultFT0C{"binsMultFT0C", {500, 0., 2000}, "bins for FT0C multiplicity"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {500, 0., 9000}, "bins for FT0M multiplicity"};
  ConfigurableAxis binsMultFDDA{"binsMultFDDA", {500, 0., 30000}, "bins for FDDA multiplicity"};
  ConfigurableAxis binsMultFDDC{"binsMultFDDC", {500, 0., 20000}, "bins for FDDC multiplicity"};
  ConfigurableAxis binsMultNTPV{"binsMultNTPV", {100, 0., 100}, "bins for NTracksPV multiplicity"};
  ConfigurableAxis binsMultNTracksGlobal{"binsMultNTracksGlobal", {60, 0., 60}, "bins for NTracksGlobal multiplicity"};
  ConfigurableAxis binsMultCombined{"binsMultCombined", {500, 0., 40}, "bins for Combined multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {200, -25, 25}, "bins for z vertex"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using JoinedMults = soa::Join<aod::FV0Mults, aod::FT0Mults, aod::FDDMults, aod::MultsGlobal, aod::PVMults, aod::MultsExtra>;

  std::array<TH1F*, CentralityDetectors::kNDetectors> zEq;
  TH1F* weights;
  static constexpr std::array<std::string, CentralityDetectors::kNDetectors> detectorNames = {"FV0A", "FT0A", "FT0C", "FDDA", "FDDC", "NTPV", "NTracksGlobal"};

  void init(InitContext&)
  {

    AxisSpec axisMultFV0A = {binsMultFV0A, "Multiplicity FV0A"};
    AxisSpec axisMultFT0A = {binsMultFT0A, "Multiplicity FT0A"};
    AxisSpec axisMultFT0C = {binsMultFT0C, "Multiplicity FT0C"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "Multiplicity FT0M"};
    AxisSpec axisMultFDDA = {binsMultFDDA, "Multiplicity FDDA"};
    AxisSpec axisMultFDDC = {binsMultFDDC, "Multiplicity FDDC"};
    AxisSpec axisMultNTPV = {binsMultNTPV, "Multiplicity NTPV"};
    AxisSpec axisMultNTracksGlobal = {binsMultNTracksGlobal, "Multiplicity NTracksGlobal"};
    AxisSpec axisMultCombined = {binsMultCombined, "Combined multiplicity"};
    AxisSpec axisZVtx = {binsZVtx, "z_{vtx} [cm]"};

    std::array<AxisSpec, CentralityDetectors::kNDetectors> axesEstimators = {axisMultFV0A, axisMultFT0A, axisMultFT0C, axisMultFDDA, axisMultFDDC, axisMultNTPV, axisMultNTracksGlobal};

    for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
      histos.add<TH2>(("MultVsZVtx/hMult" + detectorNames[i] + "VsZVtx").c_str(), ("hMult" + detectorNames[i] + "VsZVtx").c_str(), kTH2F, {axisZVtx, axesEstimators[i]});
      histos.add<TH2>(("MultVsNTracksPV/hMult" + detectorNames[i] + "VsNTPV").c_str(), ("hMult" + detectorNames[i] + "VsNTPV").c_str(), kTH2F, {axisMultNTPV, axesEstimators[i]});
      histos.add<TH2>(("MultsVsGlobalTracks/hMult" + detectorNames[i] + "VsGlobal").c_str(), ("hMult" + detectorNames[i] + "VsGlobal").c_str(), kTH2F, {axisMultNTracksGlobal, axesEstimators[i]});
    }

    histos.add<TH2>("MultVsNTracksPV/hMultFT0MVsNTPV", "hMultFT0MVsNTPV", kTH2D, {axisMultNTPV, axisMultFT0M});
    histos.add<TH2>("MultVsNTracksPV/hMultCombinedVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
    histos.add<TH2>("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
    histos.add<TH2>("MultsVsGlobalTracks/hMultFT0MVsGlobal", "hMultFT0MVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultFT0M});
    histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedVsGlobal", "hMultCombinedVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultCombined});
    histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedNoFDDCVsGlobal", "hMultCombinedVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultCombined});

    if (zEqualize) {
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        auto zEqFile = TFile::Open(pathZEq.value.c_str(), "READ");
        if (!zEqFile) {
          LOGP(fatal, "Could not open file {}", pathZEq.value.c_str());
        }
        zEq[i] = (TH1F*)(zEqFile->Get(("hWeights" + detectorNames[i]).c_str()));
        if (!zEq[i]) {
          LOGP(fatal, "Could not find histogram {} in file {}", ("hWeights" + detectorNames[i]).c_str(), pathZEq.value.c_str());
        }
        zEq[i]->SetDirectory(0);
        zEqFile->Close();


        histos.add<TH1>(("ZEqWeights/ZEqWeights" + detectorNames[i]).c_str(), ("ZEqWeights" + detectorNames[i]).c_str(), kTH1F, {axisZVtx});
        histos.add<TH2>(("MultVsZVtxZEq/hMult" + detectorNames[i] + "zEqualised").c_str(), ("hMult" + detectorNames[i] + "zEqualised").c_str(), kTH2F, {axisZVtx, {500, 0, 5}});

      }

      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        LOG(info) << "Filling zEq histograms for detector " << detectorNames[index] << zEq[index]->GetNbinsX();
        for (int iZvtx = 1; iZvtx <= zEq[index]->GetNbinsX(); iZvtx++) {
          float zPV = zEq[index]->GetBinCenter(iZvtx);

          histos.fill(HIST("ZEqWeights/ZEqWeights") + HIST(detectorNames[index]), zPV, zEq[index]->GetBinContent(iZvtx));
        }
      });
    }


    if (useWeights) {
      auto weightFile = TFile::Open(pathWeights.value.c_str(), "READ");
      if (!weightFile) {
        LOGP(fatal, "Could not open file {}", pathWeights.value.c_str());
      }
      weights = (TH1F*)(weightFile->Get("hWeights"));
      if (!weights) {
        LOGP(fatal, "Could not find histogram {} in file {}", "hWeights", pathWeights.value.c_str());
      }
      weights->SetDirectory(0);
      weightFile->Close();
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        histos.add<TH2>(("MultVsZVtxReweighted/hMult" + detectorNames[i] + "VsZVtxReweighted").c_str(), ("hMult" + detectorNames[i] + "VsZVtxReweighted").c_str(), kTH2F, {axisZVtx, axesEstimators[i]});
      }
    }
  }

  void doZequalisation(JoinedMults const& collisions)
  {

    for (auto coll : collisions) {
      float zPV = coll.multPVz();
      histos.fill(HIST("MultVsZVtx/hMultFV0AVsZVtx"), zPV, coll.multFV0A());
      histos.fill(HIST("MultVsZVtx/hMultFT0AVsZVtx"), zPV, coll.multFT0A());
      histos.fill(HIST("MultVsZVtx/hMultFT0CVsZVtx"), zPV, coll.multFT0C());
      histos.fill(HIST("MultVsZVtx/hMultFDDAVsZVtx"), zPV, coll.multFDDA());
      histos.fill(HIST("MultVsZVtx/hMultFDDCVsZVtx"), zPV, coll.multFDDC());
      histos.fill(HIST("MultVsZVtx/hMultNTPVVsZVtx"), zPV, coll.multNTracksPVeta1());
      histos.fill(HIST("MultVsZVtx/hMultNTracksGlobalVsZVtx"), zPV, coll.multNTracksGlobal());
    }

    static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
      constexpr int index = i.value;
      std::shared_ptr<TH2> hMult = histos.get<TH2>(HIST("MultVsZVtx/hMult") + HIST(detectorNames[index]) + HIST("VsZVtx"));
      for (int iZvtx = 1; iZvtx <= hMult->GetNbinsX(); iZvtx++) {
        std::shared_ptr<TH1> hWeight = histos.get<TH1>(HIST("ZEqWeights/ZEqWeights") + HIST(detectorNames[index]));
        std::shared_ptr<TH1> hProjection(hMult->ProjectionY(("hWeightsZEq" + detectorNames[index]).c_str(), iZvtx, iZvtx));
        hProjection->GetXaxis()->SetRange(2, hProjection->GetNbinsX()); // exclude the bin with 0 multiplicity
        if (hProjection->GetEntries() == 0) {
          hWeight->SetBinContent(iZvtx, 0);
        } else {
          hProjection->ComputeIntegral();
          double nq = 1;
          double quantiles[1]{0.99};
          double mult99perc[1]{0};
          int nQuantiles = hProjection->GetQuantiles(nq, mult99perc, quantiles);
          if (mult99perc[0] == 0) {
            LOGP(warning, "Could not find 99th percentile for detector {}", detectorNames[index]);
          }
          if (nQuantiles != nq) {
            hWeight->SetBinContent(iZvtx, 0);
          } else {
            hWeight->SetBinContent(iZvtx, 1 / mult99perc[0]);
          }
        }
      }
    });
  }

  void process(JoinedMults::iterator const& coll)
  {
    float zPV = coll.multPVz(); 

    histos.fill(HIST("MultVsZVtx/hMultFV0AVsZVtx"), zPV, coll.multFV0A());
    histos.fill(HIST("MultVsZVtx/hMultFT0AVsZVtx"), zPV, coll.multFT0A());
    histos.fill(HIST("MultVsZVtx/hMultFT0CVsZVtx"), zPV, coll.multFT0C());
    histos.fill(HIST("MultVsZVtx/hMultFDDAVsZVtx"), zPV, coll.multFDDA());
    histos.fill(HIST("MultVsZVtx/hMultFDDCVsZVtx"), zPV, coll.multFDDC());
    histos.fill(HIST("MultVsZVtx/hMultNTPVVsZVtx"), zPV, coll.multNTracksPVeta1());
    histos.fill(HIST("MultVsZVtx/hMultNTracksGlobalVsZVtx"), zPV, coll.multNTracksGlobal());

    std::vector<double> weigthsOfCurrentColl(CentralityDetectors::kNDetectors, 1),
      zEqWeights(CentralityDetectors::kNDetectors, 1),
      zEqCentrality(CentralityDetectors::kNDetectors, 1),
      reweightedCentr(CentralityDetectors::kNDetectors, 1);

    if (zEqualize) {
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        zEqWeights[i] *= zEq[i]->GetBinContent(zEq[i]->FindBin(zPV));
      }

      zEqCentrality = {
        coll.multFV0A() * zEqWeights[CentralityDetectors::FV0A],
        coll.multFT0A() * zEqWeights[CentralityDetectors::FT0A],
        coll.multFT0C() * zEqWeights[CentralityDetectors::FT0C],
        coll.multFDDA() * zEqWeights[CentralityDetectors::FDDA],
        coll.multFDDC() * zEqWeights[CentralityDetectors::FDDC],
        coll.multNTracksPVeta1() * zEqWeights[CentralityDetectors::NTracksPV],
        coll.multNTracksGlobal() * zEqWeights[CentralityDetectors::NTracksGlobal]};

      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        histos.fill(HIST("MultVsZVtxZEq/hMult") + HIST(detectorNames[index]) + HIST("zEqualised"), zPV, zEqCentrality[index]);
      });
    }
    
    if (useWeights) {

      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        weigthsOfCurrentColl[i] *= weights->GetBinContent(i+1);
      }
      
      weigthsOfCurrentColl = {
        coll.multFV0A() > 0 ? weigthsOfCurrentColl[CentralityDetectors::FV0A] : 0,
        coll.multFT0A() > 0 ? weigthsOfCurrentColl[CentralityDetectors::FT0A] : 0,
        coll.multFT0C() > 0 ? weigthsOfCurrentColl[CentralityDetectors::FT0C] : 0,
        coll.multFDDA() > 0 ? weigthsOfCurrentColl[CentralityDetectors::FDDA] : 0,
        coll.multFDDC() > 0 ? weigthsOfCurrentColl[CentralityDetectors::FDDC] : 0,
        coll.multNTracksPVeta1() > 0 ? weigthsOfCurrentColl[CentralityDetectors::NTracksPV] : 0,
        coll.multNTracksGlobal() > 0 ? weigthsOfCurrentColl[CentralityDetectors::NTracksGlobal] : 0};

      reweightedCentr = {
        coll.multFV0A() * weigthsOfCurrentColl[CentralityDetectors::FV0A],
        coll.multFT0A() * weigthsOfCurrentColl[CentralityDetectors::FT0A],
        coll.multFT0C() * weigthsOfCurrentColl[CentralityDetectors::FT0C],
        coll.multFDDA() * weigthsOfCurrentColl[CentralityDetectors::FDDA],
        coll.multFDDC() * weigthsOfCurrentColl[CentralityDetectors::FDDC],
        coll.multNTracksPVeta1() * weigthsOfCurrentColl[CentralityDetectors::NTracksPV],
        coll.multNTracksGlobal() * weigthsOfCurrentColl[CentralityDetectors::NTracksGlobal]};

      //TODO: add histograms for weights without zeq

      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        histos.fill(HIST("MultVsZVtxReweighted/hMult") + HIST(detectorNames[index]) + HIST("VsZVtxReweighted"), zPV, reweightedCentr[index]);
      });

      if (zEqualize) {
        for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
          reweightedCentr[i] = zEqCentrality[i] * weigthsOfCurrentColl[i];
        }
      }
    }

    double combinedMult = zEqCentrality[CentralityDetectors::FV0A] * weigthsOfCurrentColl[CentralityDetectors::FV0A] + zEqCentrality[CentralityDetectors::FT0A] * weigthsOfCurrentColl[CentralityDetectors::FT0A] + zEqCentrality[CentralityDetectors::FT0C] * weigthsOfCurrentColl[CentralityDetectors::FT0C] + zEqCentrality[CentralityDetectors::FDDA] * weigthsOfCurrentColl[CentralityDetectors::FDDA];

    double totalWeights = weigthsOfCurrentColl[CentralityDetectors::FV0A] + weigthsOfCurrentColl[CentralityDetectors::FT0A] + weigthsOfCurrentColl[CentralityDetectors::FT0C] + weigthsOfCurrentColl[CentralityDetectors::FDDA];

    histos.fill(HIST("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV"), coll.multNTracksPVeta1(), combinedMult / totalWeights);
    histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedNoFDDCVsGlobal"), coll.multNTracksGlobal(), combinedMult / totalWeights);

    combinedMult += zEqCentrality[CentralityDetectors::FDDC] * weigthsOfCurrentColl[CentralityDetectors::FDDC];
    totalWeights += weigthsOfCurrentColl[CentralityDetectors::FDDC];

    histos.fill(HIST("MultVsNTracksPV/hMultCombinedVsNTPV"), coll.multNTracksPVeta1(), combinedMult / totalWeights);
    histos.fill(HIST("MultVsNTracksPV/hMultFV0AVsNTPV"), coll.multNTracksPVeta1(), coll.multFV0A());
    histos.fill(HIST("MultVsNTracksPV/hMultFT0AVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0A());
    histos.fill(HIST("MultVsNTracksPV/hMultFT0CVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0C());
    histos.fill(HIST("MultVsNTracksPV/hMultFT0MVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0M());
    histos.fill(HIST("MultVsNTracksPV/hMultFDDAVsNTPV"), coll.multNTracksPVeta1(), coll.multFDDA());
    histos.fill(HIST("MultVsNTracksPV/hMultFDDCVsNTPV"), coll.multNTracksPVeta1(), coll.multFDDC());
    histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedVsGlobal"), coll.multNTracksGlobal(), combinedMult / totalWeights);
    histos.fill(HIST("MultsVsGlobalTracks/hMultFV0AVsGlobal"), coll.multNTracksGlobal(), coll.multFV0A());
    histos.fill(HIST("MultsVsGlobalTracks/hMultFT0AVsGlobal"), coll.multNTracksGlobal(), coll.multFT0A());
    histos.fill(HIST("MultsVsGlobalTracks/hMultFT0CVsGlobal"), coll.multNTracksGlobal(), coll.multFT0C());
    histos.fill(HIST("MultsVsGlobalTracks/hMultFT0MVsGlobal"), coll.multNTracksGlobal(), coll.multFT0M());
    histos.fill(HIST("MultsVsGlobalTracks/hMultFDDAVsGlobal"), coll.multNTracksGlobal(), coll.multFDDA());
    histos.fill(HIST("MultsVsGlobalTracks/hMultFDDCVsGlobal"), coll.multNTracksGlobal(), coll.multFDDC());
  }
 PROCESS_SWITCH(multiplicityCombination, process, "main process function", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityCombination>(cfgc)};
}
