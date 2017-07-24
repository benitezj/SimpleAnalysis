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
  TruthSmear(std::vector<std::string>&);
  TruthEvent *smearEvent(AnalysisEvent*);

 private:
  bool        smearElectrons;
  bool        smearMuons;
  bool        smearTaus;
  bool        smearPhotons;
  bool        smearJets;
  bool        smearMET;
  bool        addPileupJets;
  bool        useHGTD0;
  bool	      useHGTDbtag;
  bool        useTrackConfirm;
  std::string puEffScheme;
  float        puEff;
  UpgradePerformanceFunctions *m_upgrade;
  TRandom3 m_random;
};


#endif
