#ifndef TRUTHSMEAR_H
#define TRUTHSMEAR_H

#include <vector>
#include <string>

#include "TRandom3.h"

class TruthEvent;
class AnalysisEvent;
class UpgradePerformanceFunctions;

class TruthSmear
{
 public:
  TruthSmear(std::vector<std::string>& options);
  TruthEvent *smearEvent(AnalysisEvent *event);

 private:
  bool smearElectrons;
  bool smearMuons;
  bool smearTaus;
  bool smearPhotons;
  bool smearJets;
  bool smearMET;
  bool addPileupJets;
  UpgradePerformanceFunctions *m_upgrade;
  TRandom3 m_random;
};


#endif
