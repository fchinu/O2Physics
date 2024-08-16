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

  Configurable<bool> useWeights{"useWeights", false, "Reweights the multiplicities"};
  Configurable<std::vector<std::string>> pathsWeights{"pathsWeights", {"weights.root","weights.root","weights.root","weights.root","weights.root","weights.root","weights.root"}, "Path to the weights file"};
  Configurable<std::vector<std::string>> histoNames{"histoNames", {"weigths", "weigths", "weigths", "weigths", "weigths", "weigths", "weigths"}, "Name of the histogram in the weights file"};

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
  
  std::array<TH1F*, CentralityDetectors::kNDetectors> weights;
  static constexpr std::array<std::string, CentralityDetectors::kNDetectors> detectorNames = {"FV0A", "FT0A", "FT0C", "FDDA", "FDDC", "NTPV", "NTracksGlobal"};  
  
  void init(InitContext&) {


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
      histos.add<TH2>(("MultVsZVtx/hMult"+detectorNames[i]+"VsZVtx").c_str(), ("hMult"+detectorNames[i]+"VsZVtx").c_str(), kTH2F, {axisZVtx, axesEstimators[i]});
    }

    if (useWeights) {
      if (pathsWeights->size() != CentralityDetectors::kNDetectors || histoNames->size() != CentralityDetectors::kNDetectors) {
        LOGP(fatal, "Number of paths and histogram names must be equal to the number of detectors");
      }
      
      for (int i = 0; i < CentralityDetectors::kNDetectors; i++) {
        auto weightFile = TFile::Open((pathsWeights->at(i)).c_str(), "READ");
        if (!weightFile) {
          LOGP(fatal, "Could not open file {}", (pathsWeights->at(i)).c_str());
        }
        weights[i] = (TH1F*)(weightFile->Get((histoNames->at(i)).c_str()));
        weights[i]->SetDirectory(0);
        if (!weights[i]) {
          LOGP(fatal, "Could not find histogram {} in file {}", (histoNames->at(i)).c_str(), (pathsWeights->at(i)).c_str());
        }
        weightFile->Close();
        histos.add<TH2>(("MultVsZVtxReweighted/hMult"+detectorNames[i]+"VsZVtxReweighted").c_str(), ("hMult"+detectorNames[i]+"VsZVtxReweighted").c_str(), kTH2F, {axisZVtx, {500,0,40}});
        histos.add<TH2>(("MultVsNTracksPV/hMult"+detectorNames[i]+"VsNTPV").c_str(), ("hMult"+detectorNames[i]+"VsNTPV").c_str(), kTH2F, {axisMultNTPV, axesEstimators[i]});
        histos.add<TH2>(("MultsVsGlobalTracks/hMult"+detectorNames[i]+"VsGlobal").c_str(), ("hMult"+detectorNames[i]+"VsGlobal").c_str(), kTH2F, {axisMultNTracksGlobal, axesEstimators[i]});
      }

      histos.add<TH2>("MultVsNTracksPV/hMultFT0MVsNTPV", "hMultFT0MVsNTPV", kTH2D, {axisMultNTPV, axisMultFT0M});
      histos.add<TH2>("MultVsNTracksPV/hMultCombinedVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
      histos.add<TH2>("MultVsNTracksPV/hMultCombinedNoFDDVsNTPV", "hMultCombinedVsNTPV", kTH2D, {axisMultNTPV, axisMultCombined});
      histos.add<TH2>("MultsVsGlobalTracks/hMultFT0MVsGlobal", "hMultFT0MVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultFT0M});
      histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedVsGlobal", "hMultCombinedVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultCombined});
      histos.add<TH2>("MultsVsGlobalTracks/hMultCombinedNoFDDVsGlobal", "hMultCombinedVsGlobal", kTH2D, {axisMultNTracksGlobal, axisMultCombined});
      histos.add<TH2>("hMultFT0AVsFV0A", "hMultFT0CVsFV0A", kTH2D, {axisMultFV0A, axisMultFT0A});
    }

  }

  void process(JoinedMults::iterator const& coll)
  {
    float zPV = coll.multPVz();
 
    if (useWeights) {
      std::vector<double> weigthsOfCurrentColl{
                                                weights[CentralityDetectors::FV0A]->GetBinContent(weights[CentralityDetectors::FV0A]->FindBin(zPV)),
                                                weights[CentralityDetectors::FT0A]->GetBinContent(weights[CentralityDetectors::FT0A]->FindBin(zPV)),
                                                weights[CentralityDetectors::FT0C]->GetBinContent(weights[CentralityDetectors::FT0C]->FindBin(zPV)),
                                                weights[CentralityDetectors::FDDA]->GetBinContent(weights[CentralityDetectors::FDDA]->FindBin(zPV)),
                                                weights[CentralityDetectors::FDDC]->GetBinContent(weights[CentralityDetectors::FDDC]->FindBin(zPV)),
                                                weights[CentralityDetectors::NTracksPV]->GetBinContent(weights[CentralityDetectors::NTracksPV]->FindBin(zPV)),
                                                weights[CentralityDetectors::NTracksGlobal]->GetBinContent(weights[CentralityDetectors::NTracksGlobal]->FindBin(zPV))
      };

      std::vector<double> reweightedCentr = {
                                              coll.multFV0A()*weigthsOfCurrentColl[CentralityDetectors::FV0A],
                                              coll.multFT0A()*weigthsOfCurrentColl[CentralityDetectors::FT0A],
                                              coll.multFT0C()*weigthsOfCurrentColl[CentralityDetectors::FT0C],
                                              coll.multFDDA()*weigthsOfCurrentColl[CentralityDetectors::FDDA],
                                              coll.multFDDC()*weigthsOfCurrentColl[CentralityDetectors::FDDC],
                                              coll.multNTracksPVeta1()*weigthsOfCurrentColl[CentralityDetectors::NTracksPV],
                                              coll.multNTracksGlobal()*weigthsOfCurrentColl[CentralityDetectors::NTracksGlobal]
      };

      static_for<0, CentralityDetectors::kNDetectors-1>([&](auto i) {
        constexpr int index = i.value;
        histos.fill(HIST("MultVsZVtxReweighted/hMult")+HIST(detectorNames[index])+HIST("VsZVtxReweighted"), zPV, reweightedCentr[index]);
      });


      double CombinedMult =   coll.multFV0A()*weigthsOfCurrentColl[CentralityDetectors::FV0A] + 
                              coll.multFT0A()*weigthsOfCurrentColl[CentralityDetectors::FT0A] + 
                              coll.multFT0C()*weigthsOfCurrentColl[CentralityDetectors::FT0C];

      histos.fill(HIST("MultVsNTracksPV/hMultCombinedNoFDDVsNTPV"), coll.multNTracksPVeta1(), CombinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedNoFDDVsGlobal"), coll.multNTracksGlobal(), CombinedMult);

      CombinedMult += coll.multFDDA()*weigthsOfCurrentColl[CentralityDetectors::FDDA] + coll.multFDDC()*weigthsOfCurrentColl[CentralityDetectors::FDDC];

      histos.fill(HIST("MultVsNTracksPV/hMultCombinedVsNTPV"), coll.multNTracksPVeta1(), CombinedMult);
      histos.fill(HIST("MultVsNTracksPV/hMultFV0AVsNTPV"), coll.multNTracksPVeta1(), coll.multFV0A());
      histos.fill(HIST("MultVsNTracksPV/hMultFT0AVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0A());
      histos.fill(HIST("MultVsNTracksPV/hMultFT0CVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0C());
      histos.fill(HIST("MultVsNTracksPV/hMultFT0MVsNTPV"), coll.multNTracksPVeta1(), coll.multFT0M());
      histos.fill(HIST("MultVsNTracksPV/hMultFDDAVsNTPV"), coll.multNTracksPVeta1(), coll.multFDDA());
      histos.fill(HIST("MultVsNTracksPV/hMultFDDCVsNTPV"), coll.multNTracksPVeta1(), coll.multFDDC());
      histos.fill(HIST("MultsVsGlobalTracks/hMultCombinedVsGlobal"), coll.multNTracksGlobal(), CombinedMult);
      histos.fill(HIST("MultsVsGlobalTracks/hMultFV0AVsGlobal"), coll.multNTracksGlobal(), coll.multFV0A());
      histos.fill(HIST("MultsVsGlobalTracks/hMultFT0AVsGlobal"), coll.multNTracksGlobal(), coll.multFT0A());
      histos.fill(HIST("MultsVsGlobalTracks/hMultFT0CVsGlobal"), coll.multNTracksGlobal(), coll.multFT0C());
      histos.fill(HIST("MultsVsGlobalTracks/hMultFT0MVsGlobal"), coll.multNTracksGlobal(), coll.multFT0M());
      histos.fill(HIST("MultsVsGlobalTracks/hMultFDDAVsGlobal"), coll.multNTracksGlobal(), coll.multFDDA());
      histos.fill(HIST("MultsVsGlobalTracks/hMultFDDCVsGlobal"), coll.multNTracksGlobal(), coll.multFDDC());


    } 

    histos.fill(HIST("MultVsZVtx/hMultFV0AVsZVtx"), zPV, coll.multFV0A());
    histos.fill(HIST("MultVsZVtx/hMultFT0AVsZVtx"), zPV, coll.multFT0A());
    histos.fill(HIST("MultVsZVtx/hMultFT0CVsZVtx"), zPV, coll.multFT0C());
    histos.fill(HIST("MultVsZVtx/hMultFDDAVsZVtx"), zPV, coll.multFDDA());
    histos.fill(HIST("MultVsZVtx/hMultFDDCVsZVtx"), zPV, coll.multFDDC());
    histos.fill(HIST("MultVsZVtx/hMultNTPVVsZVtx"), zPV, coll.multNTracksPVeta1());
    histos.fill(HIST("MultVsZVtx/hMultNTracksGlobalVsZVtx"), zPV, coll.multNTracksGlobal());
    histos.fill(HIST("hMultFT0AVsFV0A"), coll.multFV0A(), coll.multFT0A());

    

  }
  PROCESS_SWITCH(multiplicityCombination, process, "main process function", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<multiplicityCombination>(cfgc)};
}
