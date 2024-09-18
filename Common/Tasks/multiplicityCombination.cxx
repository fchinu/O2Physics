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

  std::array<TH1F*, CentralityDetectors::kNDetectors> zEqWeightsVsZvtx;
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
      histos.add<TH2>(("MultVsNTracksPV/hMult" + detectorNames[i] + "VsNTPV").c_str(), ("hMult" + detectorNames[i] + "VsNTPV").c_str(), kTH2F, {{binsMultZEq, "N_{tracks}^{PV}"}, {binsMultZEq, (detectorNames[i]+" multiplicity").c_str()}});
      histos.add<TH2>(("MultsVsGlobalTracks/hMult" + detectorNames[i] + "VsGlobal").c_str(), ("hMult" + detectorNames[i] + "VsGlobal").c_str(), kTH2F, {{binsMultZEq, "N_{tracks}^{global}"}, {binsMultZEq, (detectorNames[i]+" multiplicity").c_str()}});
    }

    histos.add<TH2>("MultVsNTracksPV/hMultCombinedVsNTPV", "hMultCombinedVsNTPV", kTH2D, {{binsMultZEq, "N_{tracks}^{PV}"}, axisMultCombined});
    histos.add<TH2>("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV", "hMultCombinedVsNTPV", kTH2D, {{binsMultZEq, "N_{tracks}^{PV}"}, axisMultCombined});
    histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedVsGlobal", "hMultCombinedVsGlobal", kTH2D, {{binsMultZEq, "N_{tracks}^{global}"}, axisMultCombined});
    histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedNoFDDCVsGlobal", "hMultCombinedVsGlobal", kTH2D, {{binsMultZEq, "N_{tracks}^{global}"}, axisMultCombined});

    if (zEqualize) {
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        auto zEqFile = TFile::Open(pathZEq.value.c_str(), "READ");
        if (!zEqFile) {
          LOGP(fatal, "Could not open file {}", pathZEq.value.c_str());
        }
        zEqWeightsVsZvtx[i] = (TH1F*)(zEqFile->Get(("hWeights" + detectorNames[i]).c_str()));
        if (!zEqWeightsVsZvtx[i]) {
          LOGP(fatal, "Could not find histogram {} in file {}", ("hWeights" + detectorNames[i]).c_str(), pathZEq.value.c_str());
        }
        zEqWeightsVsZvtx[i]->SetDirectory(0);
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
        for (int iZvtx = 1; iZvtx <= zEqWeightsVsZvtx[index]->GetNbinsX(); iZvtx++) {
          float zPV = zEqWeightsVsZvtx[index]->GetBinCenter(iZvtx);

          histos.fill(HIST("ZEqWeights/ZEqWeights") + HIST(detectorNames[index]), zPV, zEqWeightsVsZvtx[index]->GetBinContent(iZvtx));
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
        histos.add<TH2>(("MultDistributionsCalibratedVsNTracksPV/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH2F, {{binsMultZEq, "N_{tracks}^{PV}"}, {100, 0, 100}});
        histos.add<TH2>(("MultDistributionsCalibratedVsGlobal/hMult" + detectorNames[i]).c_str(), ("hMult" + detectorNames[i]).c_str(), kTH2F, {{binsMultZEq, "N_{tracks}^{global}"}, {100, 0, 100}});
      }
      hCalibrations[CentralityDetectors::kNDetectors] = (TH1F*)(calibFile->Get("hCalibCombined"));
      calibFile->Close();
      histos.add<TH2>("MultDistributionsCalibratedVsNTracksPV/hMultCombined", "hMultCombined", kTH2F, {{binsMultZEq, "N_{tracks}^{PV}"}, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsNTracksPV/hMultCombinedNoFDDC", "hMultCombinedNoFDDC", kTH2F, {{binsMultZEq, "N_{tracks}^{PV}"}, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsGlobal/hMultCombined", "hMultCombined", kTH2F, {{binsMultZEq, "N_{tracks}^{global}"}, {100, 0, 100}});
      histos.add<TH2>("MultDistributionsCalibratedVsGlobal/hMultCombinedNoFDDC", "hMultCombinedNoFDDC", kTH2F, {{binsMultZEq, "N_{tracks}^{global}"}, {100, 0, 100}});
    }
  }


  /**
   * @brief Calculates the equalized multiplicities for each detector in the collision.
   * 
   * This function takes a collision of type T and calculates the equalized multiplicities
   * for each detector in the collision. The equalized multiplicities are then used to fill
   * various histograms for analysis purposes.
   * 
   * @tparam T The type of the collision.
   * @param coll The collision used to retrieve the midrapidity multiplicity.
   * @param zEqMultiplicity The vector to store the equalized multiplicities.
   */
  template <typename T>
  void getZEqMultiplicity(const T& coll, std::vector<double>& zEqMultiplicity)
  {
    float zPV = coll.multPVz();
    static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
      constexpr int index = i.value;
      zEqMultiplicity[index] *= zEqWeightsVsZvtx[index]->GetBinContent(zEqWeightsVsZvtx[index]->FindBin(zPV));
      // Firstly fill the histogram with the equalized multiplicities
      histos.fill(HIST("MultVsZVtxZEq/hMult") + HIST(detectorNames[index]) + HIST("zEqualised"), zPV, zEqMultiplicity[index]);
      histos.fill(HIST("MultDistributions/zEq/hMult") + HIST(detectorNames[index]), zEqMultiplicity[index]);

      // Fill the histograms for the equalized multiplicities vs the number of tracks
      histos.fill(HIST("MultVsNTracksPV/hMult") + HIST(detectorNames[index]) + HIST("VsNTPV"), zEqMultiplicity[CentralityDetectors::NTracksPV], zEqMultiplicity[index]);
      histos.fill(HIST("MultsVsGlobalTracks/hMult") + HIST(detectorNames[index]) + HIST("VsGlobal"), zEqMultiplicity[CentralityDetectors::NTracksGlobal], zEqMultiplicity[index]);
    });

    // Fill the correlation histograms for the equalized multiplicities once all the detectors have been processed
    for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
      for (int j = 0; j < CentralityDetectors::kNDetectors; j++) {
        hCorrelations[i][j]->Fill(zEqMultiplicity[i], zEqMultiplicity[j]);
      }
    }
  }

  /**
   * @brief Calculates the reweighted multiplicity for a given collision.
   * 
   * This function takes a collision `coll` and calculates the reweighted multiplicity
   * based on the values in the `multiplicity` vector and the weights stored in the `weightsOfCurrentColl` vector.
   * The reweighted multiplicity is stored in the `multiplicity` vector.
   * 
   * @tparam T The type of the collision.
   * @param coll The collision used to retrieve the midrapidity multiplicity.
   * @param multiplicity The vector containing the multiplicity values.
   * @param weigthsOfCurrentColl The vector containing the weights of the current collision.
   */
  template <typename T>
  void getReweightedMulitplicity(const T& coll, std::vector<double>& multiplicity, std::vector<double>& weigthsOfCurrentColl)
  {
    float zPV = coll.multPVz();
    static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
      constexpr int index = i.value;
      weigthsOfCurrentColl[index] = multiplicity[index] > 0 ? weights->GetBinContent(index + 1) : 0;
      multiplicity[index] *= weigthsOfCurrentColl[index];
      histos.fill(HIST("MultVsZVtxReweighted/hMult") + HIST(detectorNames[index]) + HIST("VsZVtxReweighted"), zPV, multiplicity[index]);
      histos.fill(HIST("MultDistributions/Reweighted/hMult") + HIST(detectorNames[index]), multiplicity[index]);
    });
  }

  /**
   * @brief Computes the combined multiplicity based on the given collision, multiplicities, and weights.
   * 
   * This function calculates the combined multiplicity by summing the multiplicities from different centrality detectors
   * and dividing it by the total weights. It also supports producing distributions with percentiles if specified.
   * 
   * @tparam T The type of the collision.
   * @param coll The collision used to retrieve the midrapidity multiplicity.
   * @param multiplicity The vector of multiplicities from different centrality detectors.
   * @param weights The vector of weights corresponding to the multiplicities.
   * @param produceDistributionsWithPercentiles Flag indicating whether to produce distributions with percentiles.
   */
  template <typename T>
  void computeCombinedMultiplicity(const T& coll, const std::vector<double>& multiplicity, const std::vector<double>& zEqMultiplicity, const std::vector<double>& weights, bool produceDistributionsWithPercentiles)
  {
    // Firstly calculate the combined multiplicity without FDDC
    double summedMult = multiplicity[CentralityDetectors::FV0A] + multiplicity[CentralityDetectors::FT0A] + multiplicity[CentralityDetectors::FT0C] + multiplicity[CentralityDetectors::FDDA];
    double totalWeights = weights[CentralityDetectors::FV0A] + weights[CentralityDetectors::FT0A] + weights[CentralityDetectors::FT0C] + weights[CentralityDetectors::FDDA];
    double combinedMult = summedMult / totalWeights;

    if (summedMult > 0) {

      if (produceDistributionsWithPercentiles) {
        double combinedMultPercentile = hCalibrations[CentralityDetectors::kNDetectors]->GetBinContent(hCalibrations[CentralityDetectors::kNDetectors]->FindBin(combinedMult));
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMultCombinedNoFDDC"), zEqMultiplicity[CentralityDetectors::NTracksPV], combinedMultPercentile);
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMultCombinedNoFDDC"), zEqMultiplicity[CentralityDetectors::NTracksGlobal], combinedMultPercentile);
      }

      histos.fill(HIST("MultDistributions/Reweighted/hMultCombinedNoFDDC"), combinedMult);
      histos.fill(HIST("MultVsNTracksPV/hMultCombinedNoFDDCVsNTPV"), zEqMultiplicity[CentralityDetectors::NTracksPV], combinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedNoFDDCVsGlobal"), zEqMultiplicity[CentralityDetectors::NTracksGlobal], combinedMult);
    }

    // Now calculate the combined multiplicity with FDDC
    summedMult += multiplicity[CentralityDetectors::FDDC] * weights[CentralityDetectors::FDDC];
    totalWeights += weights[CentralityDetectors::FDDC];
    combinedMult = summedMult / totalWeights;

    if (summedMult > 0) {
      if (produceDistributionsWithPercentiles) {
      double combinedMultPercentile = hCalibrations[CentralityDetectors::kNDetectors]->GetBinContent(hCalibrations[CentralityDetectors::kNDetectors]->FindBin(combinedMult));
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMultCombined"), zEqMultiplicity[CentralityDetectors::NTracksPV], combinedMultPercentile);
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMultCombined"), zEqMultiplicity[CentralityDetectors::NTracksGlobal], combinedMultPercentile);
      }
      histos.fill(HIST("MultDistributions/Reweighted/hMultCombined"), combinedMult);
      histos.fill(HIST("MultVsNTracksPV/hMultCombinedVsNTPV"), zEqMultiplicity[CentralityDetectors::NTracksPV], combinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedVsGlobal"), zEqMultiplicity[CentralityDetectors::NTracksGlobal], combinedMult);
    }
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

    std::vector<double> weigthsOfCurrentColl(CentralityDetectors::kNDetectors, 1);
    std::vector<double> zEqMultiplicity(multiplicity);

    if (zEqualize) {
      getZEqMultiplicity(coll, zEqMultiplicity);
    }

    std::vector<double> reweightedMultiplicity(zEqMultiplicity);

    if (useWeights) {
      getReweightedMulitplicity(coll, reweightedMultiplicity, weigthsOfCurrentColl);
    }

    if (produceDistributionsWithPercentiles) {
      static_for<0, CentralityDetectors::kNDetectors - 1>([&](auto i) {
        constexpr int index = i.value;
        histos.fill(HIST("MultDistributionsCalibratedVsNTracksPV/hMult") + HIST(detectorNames[index]), zEqMultiplicity[CentralityDetectors::NTracksPV], hCalibrations[index]->GetBinContent(hCalibrations[index]->FindBin(zEqMultiplicity[index])));
        histos.fill(HIST("MultDistributionsCalibratedVsGlobal/hMult") + HIST(detectorNames[index]), zEqMultiplicity[CentralityDetectors::NTracksGlobal], hCalibrations[index]->GetBinContent(hCalibrations[index]->FindBin(zEqMultiplicity[index])));
      });
    }

    computeCombinedMultiplicity(coll, reweightedMultiplicity, zEqMultiplicity, weigthsOfCurrentColl, produceDistributionsWithPercentiles);

  }
  PROCESS_SWITCH(multiplicityCombination, process, "main process function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityCombination>(cfgc)};
}
