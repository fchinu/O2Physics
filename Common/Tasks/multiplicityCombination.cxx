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
#include "Common/Tools/Multiplicity/multCalibrator.h"
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
  FT0M,
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
  Configurable<std::string> weightsHistName{"weightsHistName", "hWeights", "Name of the weights histogram"};
  Configurable<bool> produceDistributionsWithPercentiles{"produceDistributionsWithPercentiles", false, "Produce distributions with percentiles"};
  Configurable<std::string> pathCalibration{"pathCalibration", "calibration.root", "Path to .root file with multiplicity calibration"};

  ConfigurableAxis binsMultFV0A{"binsMultFV0A", {500, 0., 20000}, "bins for FV0A multiplicity"};
  ConfigurableAxis binsMultFT0A{"binsMultFT0A", {500, 0., 7000}, "bins for FT0A multiplicity"};
  ConfigurableAxis binsMultFT0C{"binsMultFT0C", {500, 0., 2000}, "bins for FT0C multiplicity"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {500, 0., 9000}, "bins for FT0M multiplicity"};
  ConfigurableAxis binsMultFDDA{"binsMultFDDA", {500, 0., 30000}, "bins for FDDA multiplicity"};
  ConfigurableAxis binsMultFDDC{"binsMultFDDC", {500, 0., 20000}, "bins for FDDC multiplicity"};
  ConfigurableAxis binsMultNTPV{"binsMultNTPV", {100, 0., 100}, "bins for NTracksPV multiplicity"};
  ConfigurableAxis binsMultNTracksGlobal{"binsMultNTracksGlobal", {60, 0., 60}, "bins for NTracksGlobal multiplicity"};
  ConfigurableAxis binsMultCombined{"binsMultCombined", {500, 0., 40}, "bins for Combined multiplicity"};
  ConfigurableAxis binsMultZEq{"binsMultZEq", {500, 0., 1}, "bins for Z equalized multiplicity"};
  ConfigurableAxis binsZVtx{"binsZVtx", {200, -25, 25}, "bins for z vertex"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using JoinedMults = soa::Join<aod::FV0Mults, aod::FT0Mults, aod::FDDMults, aod::MultsGlobal, aod::PVMults, aod::MultsExtra>;

  std::array<TH1F*, CentralityDetectors::kNDetectors> zEq;
  std::array<TH1F*, CentralityDetectors::kNDetectors + 1> hCalibrations; // +1 for combined
  std::array<std::array<std::shared_ptr<TH2>, CentralityDetectors::kNDetectors>, CentralityDetectors::kNDetectors> hCorrelations;
  TH1F* weights;

  static constexpr std::array<std::string, CentralityDetectors::kNDetectors> detectorNames = {"FV0A", "FT0A", "FT0C", "FT0M", "FDDA", "FDDC", "NTPV", "NTracksGlobal"};

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

    std::array<AxisSpec, CentralityDetectors::kNDetectors> axesEstimators = {axisMultFV0A, axisMultFT0A, axisMultFT0C, axisMultFT0M, axisMultFDDA, axisMultFDDC, axisMultNTPV, axisMultNTracksGlobal};

    for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
      histos.add<TH2>(("MultVsZVtx/hMult" + detectorNames[i] + "VsZVtx").c_str(), ("hMult" + detectorNames[i] + "VsZVtx").c_str(), kTH2F, {axisZVtx, axesEstimators[i]});
      histos.add<TH1>(("MultDistributions/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH1F, {axesEstimators[i]});
      histos.add<TH2>(("MultVsNTracksPV/hMult" + detectorNames[i] + "VsNTPV").c_str(), ("hMult" + detectorNames[i] + "VsNTPV").c_str(), kTH2F, {axisMultNTPV, binsMultZEq});
      histos.add<TH2>(("MultsVsGlobalTracks/hMult" + detectorNames[i] + "VsGlobal").c_str(), ("hMult" + detectorNames[i] + "VsGlobal").c_str(), kTH2F, {axisMultNTracksGlobal, binsMultZEq});
    }

    histos.add<TH2>("MultVsNTracksPV/hMultCombinedVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
    histos.add<TH2>("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
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

        histos.add<TH1>(("MultDistributions/zEq/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH1F, {binsMultZEq});
        histos.add<TH1>(("ZEqWeights/ZEqWeights" + detectorNames[i]).c_str(), ("ZEqWeights" + detectorNames[i]).c_str(), kTH1F, {axisZVtx});
        histos.add<TH2>(("MultVsZVtxZEq/hMult" + detectorNames[i] + "zEqualised").c_str(), ("hMult" + detectorNames[i] + "zEqualised").c_str(), kTH2F, {axisZVtx, binsMultZEq});

        for (int j = 0; j < CentralityDetectors::kNDetectors; j++) {
          hCorrelations[i][j] = histos.add<TH2>(("MultDistributions/zEq/hCorr" + detectorNames[i] + "_" + detectorNames[j]).c_str(), ("hCorr" + detectorNames[i] + "_" + detectorNames[j]).c_str(), kTH2F, {{1000, 0, 2}, {1000, 0, 2}});
        }
      }

      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
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
      weights = (TH1F*)(weightFile->Get(weightsHistName.value.c_str()));
      if (!weights) {
        LOGP(fatal, "Could not find histogram {} in file {}", weightsHistName.value.c_str(), pathWeights.value.c_str());
      }
      weights->SetDirectory(0);
      weightFile->Close();
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        histos.add<TH2>(("MultVsZVtxReweighted/hMult" + detectorNames[i] + "VsZVtxReweighted").c_str(), ("hMult" + detectorNames[i] + "VsZVtxReweighted").c_str(), kTH2F, {axisZVtx, binsMultZEq});
        histos.add<TH1>(("MultDistributions/Reweighted/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH1F, {binsMultZEq});
      }
      histos.add<TH1>("MultDistributions/Reweighted/hMultCombinedNoFDDC", "hMultCombinedNoFDDC", kTH1F, {binsMultZEq});
      histos.add<TH1>("MultDistributions/Reweighted/hMultCombined", "hMultCombined", kTH1F, {binsMultZEq});
    }

    if (produceDistributionsWithPercentiles) {
      auto calibFile = TFile::Open(pathCalibration.value.c_str(), "READ");
      if (!calibFile) {
        LOGP(fatal, "Could not open file {}", pathCalibration.value.c_str());
      }
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        hCalibrations[i] = (TH1F*)(calibFile->Get(("hCalib" + detectorNames[i]).c_str()));
        if (!hCalibrations[i]) {
          LOGP(fatal, "Could not find histogram {} in file {}", ("hCalib" + detectorNames[i]).c_str(), pathCalibration.value.c_str());
        }
        histos.add<TH2>(("MultDistributionsCalibratedVsNTracksPV/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH2F, {axisMultNTPV, {100, 0, 100}});
        histos.add<TH2>(("MultDistributionsCalibratedVsGlobal/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH2F, {axisMultNTracksGlobal, {100, 0, 100}});
      }
      hCalibrations[CentralityDetectors::kNDetectors] = (TH1F*)(calibFile->Get("hCalibCombined"));
      calibFile->Close();
      histos.add<TH2>("MultDistributionsCalibratedVsNTracksPV/hMultCombined", "hMultCombined", kTH2F, {axisMultNTPV, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsNTracksPV/hMultCombinedNoFDDC", "hMultCombinedNoFDDC", kTH2F, {axisMultNTPV, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsGlobal/hMultCombined", "hMultCombined", kTH2F, {axisMultNTracksGlobal, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsGlobal/hMultCombinedNoFDDC", "hMultCombinedNoFDDC", kTH2F, {axisMultNTracksGlobal, {100, 0, 100}});
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

    std::vector<double> multiplicity{coll.multFV0A(), coll.multFT0A(), coll.multFT0C(), coll.multFT0M(), coll.multFDDA(), coll.multFDDC(), static_cast<double>(coll.multNTracksPVeta1()), static_cast<double>(coll.multNTracksGlobal())};

    static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
      constexpr int index = i.value;
      histos.fill(HIST("MultVsZVtx/hMult") + HIST(detectorNames[index]) + HIST("VsZVtx"), zPV, multiplicity[index]);
      histos.fill(HIST("MultDistributions/hMult") + HIST(detectorNames[index]), multiplicity[index]);
    });

    std::vector<double> weigthsOfCurrentColl(CentralityDetectors::kNDetectors, 1),
      zEqCentrality(multiplicity),
      reweightedCentr(multiplicity);

    if (zEqualize) {
      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        zEqCentrality[index] *= zEq[index]->GetBinContent(zEq[index]->FindBin(zPV));
        histos.fill(HIST("MultVsZVtxZEq/hMult") + HIST(detectorNames[index]) + HIST("zEqualised"), zPV, zEqCentrality[index]);
        histos.fill(HIST("MultDistributions/zEq/hMult") + HIST(detectorNames[index]), zEqCentrality[index]);
      });
    }

    if (useWeights) {

      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        weigthsOfCurrentColl[index] = multiplicity[index] > 0 ? weights->GetBinContent(index + 1) : 0;
        reweightedCentr[index] *= weigthsOfCurrentColl[index];
        if (zEqualize) {
          reweightedCentr[index] *= zEq[index]->GetBinContent(zEq[index]->FindBin(zPV));
        }
        histos.fill(HIST("MultVsZVtxReweighted/hMult") + HIST(detectorNames[index]) + HIST("VsZVtxReweighted"), zPV, reweightedCentr[index]);
        histos.fill(HIST("MultDistributions/Reweighted/hMult") + HIST(detectorNames[index]), reweightedCentr[index]);
        for (int j = 0; j < CentralityDetectors::kNDetectors; j++) {
          hCorrelations[i][j]->Fill(zEqCentrality[index], zEqCentrality[j]);
        }
        //static_for<i + 1, CentralityDetectors::kNDetectors - 1>([&, index](auto j) {
        //  constexpr int jndex = j.value;
        //  histos.fill(HIST("MultDistributions/Reweighted/hCorr") + HIST(detectorNames[index]) + HIST("_") + HIST(detectorNames[jndex]), multiplicity[index], multiplicity[jndex]);
        //});
      });
      //TODO: add histograms for weights without zeq
    }

    if (produceDistributionsWithPercentiles) {
      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMult") + HIST(detectorNames[index]), coll.multNTracksPVeta1(), hCalibrations[index]->GetBinContent(hCalibrations[index]->FindBin(zEqCentrality[index])));
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMult") + HIST(detectorNames[index]), coll.multNTracksGlobal(), hCalibrations[index]->GetBinContent(hCalibrations[index]->FindBin(zEqCentrality[index])));
      });
    }

    double summedMult = reweightedCentr[CentralityDetectors::FV0A] + reweightedCentr[CentralityDetectors::FT0A] + reweightedCentr[CentralityDetectors::FT0C] + reweightedCentr[CentralityDetectors::FDDA];
    double totalWeights = weigthsOfCurrentColl[CentralityDetectors::FV0A] + weigthsOfCurrentColl[CentralityDetectors::FT0A] + weigthsOfCurrentColl[CentralityDetectors::FT0C] + weigthsOfCurrentColl[CentralityDetectors::FDDA];
    double combinedMult = summedMult / totalWeights;
    double combinedMultPercentile = hCalibrations[CentralityDetectors::kNDetectors]->GetBinContent(hCalibrations[CentralityDetectors::kNDetectors]->FindBin(combinedMult));

    if (summedMult > 0) {

      if (produceDistributionsWithPercentiles) {
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMultCombinedNoFDDC"), coll.multNTracksPVeta1(), combinedMultPercentile);
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMultCombinedNoFDDC"), coll.multNTracksGlobal(), combinedMultPercentile);
      }

      histos.fill(HIST("MultDistributions/Reweighted/hMultCombinedNoFDDC"), combinedMult);
      histos.fill(HIST("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV"), coll.multNTracksPVeta1(), combinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedNoFDDCVsGlobal"), coll.multNTracksGlobal(), combinedMult);
    }

    summedMult += zEqCentrality[CentralityDetectors::FDDC] * weigthsOfCurrentColl[CentralityDetectors::FDDC];
    totalWeights += weigthsOfCurrentColl[CentralityDetectors::FDDC];
    combinedMult = summedMult / totalWeights;

    if (summedMult > 0) {

      if (produceDistributionsWithPercentiles) {
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMultCombined"), coll.multNTracksPVeta1(), combinedMultPercentile);
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMultCombined"), coll.multNTracksGlobal(), combinedMultPercentile);
      }

      histos.fill(HIST("MultDistributions/Reweighted/hMultCombined"), combinedMult);
      histos.fill(HIST("MultVsNTracksPV/hMultCombinedVsNTPV"), coll.multNTracksPVeta1(), combinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedVsGlobal"), coll.multNTracksGlobal(), combinedMult);
    }

    static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
      constexpr int index = i.value;
      histos.fill(HIST("MultVsNTracksPV/hMult") + HIST(detectorNames[index]) + HIST("VsNTPV"), coll.multNTracksPVeta1(), zEqCentrality[index]);
      histos.fill(HIST("MultsVsGlobalTracks/hMult") + HIST(detectorNames[index]) + HIST("VsGlobal"), coll.multNTracksGlobal(), zEqCentrality[index]);
    });

    if (combinedMult > 4) {
      LOGP(warning, "Mult: {}, zPv: {}, FT0A: {}, FT0C: {}, FDDA: {}, FDDC: {}, NTPV: {}, NTracksGlobal: {}",
           combinedMult, zPV, coll.multFT0A(), coll.multFT0C(), coll.multFDDA(), coll.multFDDC(), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
    }
  }
  PROCESS_SWITCH(multiplicityCombination, process, "main process function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityCombination>(cfgc)};
}
