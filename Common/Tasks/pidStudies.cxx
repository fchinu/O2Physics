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

/// \file pidStudies.cxx
/// \brief task for studies of PID performance
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università and INFN Torino
/// \author Stefano Politanò <stefano.politano@cern.ch>, INFN Torino
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct pidStudies {
  HistogramRegistry registry{"registry", {}};
  static constexpr int8_t nEstimators = 8;
  static constexpr std::string estimatorsNames[nEstimators] = {"FV0A", "FT0A", "FT0C", "FT0M", "FDDA", "FDDC", "FDDM", "NTPV"};

  std::vector<unsigned> consideredParticles = {
    11,  // e
    13,  // mu
    211, // pi
    321, // K
    2122 // p
  };

  ConfigurableAxis axisFV0A = {"axisFV0A", {100, 0., 20000.}, "axis for FV0A estimator"};
  ConfigurableAxis axisFT0A = {"axisFT0A", {100, 0., 10000.}, "axis for FT0A estimator"};
  ConfigurableAxis axisFT0C = {"axisFT0C", {100, 0., 5000.}, "axis for FT0C estimator"};
  ConfigurableAxis axisFT0M = {"axisFT0M", {100, 0., 10000.}, "axis for FT0M estimator"};
  ConfigurableAxis axisFDDA = {"axisFDDA", {100, 0., 20000.}, "axis for FDDA estimator"};
  ConfigurableAxis axisFDDC = {"axisFDDC", {100, 0., 5000.}, "axis for FDDC estimator"};
  ConfigurableAxis axisFDDM = {"axisFDDM", {100, 0., 20000.}, "axis for FDDM estimator"};
  ConfigurableAxis axisNTPV = {"axisNTPV", {100, 0., 100.}, "axis for NTPV estimator"};
  ConfigurableAxis axisdNdEta = {"axisdNdEta", {100, 0., 100.}, "axis for dN/deta"};

  std::vector<ConfigurableAxis*> estimatorsAxes = {&axisFV0A, &axisFT0A, &axisFT0C, &axisFT0M, &axisFDDA, &axisFDDC, &axisFDDM, &axisNTPV};

  Preslice<aod::McParticles> particlesPerCollision = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  using CollisionsWithMult = soa::Join<aod::Collisions, aod::PVMults, aod::MultZeqs, aod::EvSels, aod::McCollisionLabels>;

  void init(InitContext&)
  {
    for (int8_t i = 0; i < nEstimators; i++) {
      registry.add<TH2>(("etaPFive/" + estimatorsNames[i] + "VsdNdeta").c_str(), (estimatorsNames[i] + "VsdNdeta;" + estimatorsNames[i] + ";<dN_{ch}/d#eta>").c_str(), HistType::kTH2F, {*(estimatorsAxes[i]), axisdNdEta});
      registry.add<TH2>(("etaOne/" + estimatorsNames[i] + "VsdNdeta").c_str(), (estimatorsNames[i] + "VsdNdeta;" + estimatorsNames[i] + ";<dN_{ch}/d#eta>").c_str(), HistType::kTH2F, {*(estimatorsAxes[i]), axisdNdEta});
    }
  }

  void process(CollisionsWithMult const& collisions,
               aod::McCollisions const& mcCollisions,
               aod::McParticles const& particles,
               soa::Join<aod::BCs, aod::Timestamps> const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidStudies>(cfgc)};
}
