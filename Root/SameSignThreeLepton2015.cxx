#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(SameSignThreeLepton2015)

void SameSignThreeLepton2015::Init()
{
  addRegions({"SR0b3j","SR0b5j","SR1b","SR3b"});
  addRegions({"SR1b3L","SR1b3LSS"});
}

void SameSignThreeLepton2015::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons = filterCrack(event->getElectrons(10, 2.47, ELooseLH|ED0Sigma5));
  auto baselineMuons = event->getMuons(10, 2.5, MuMedium);
  auto jets          = event->getJets(20., 2.8, JVT50Jet);
  auto metVec        = event->getMET();

  // Standard SUSY overlap removal      
  jets               = overlapRemoval(jets, baselineElectrons, 0.2, NOT(BTag80MV2c20));
  jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks);
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, jets, 0.4);

  auto signalElectrons = filterObjects(baselineElectrons, 10, 2.0, ETightLH|EZ05mm|EIsoFixedCutTight);
  auto signalMuons     = filterObjects(baselineMuons, 10, 2.5, MuD0Sigma3|MuZ05mm|MuIsoFixedCutTightTrackOnly);
  
  auto signalLeptons = signalElectrons+signalMuons;

  int nLeptons = signalLeptons.size();
  if (nLeptons<2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];
  if (lep1.Pt()<20) return;
  if (lep0.charge()!=lep1.charge() && nLeptons<3) return;

  // Variables we'll be cutting on
  double met   = metVec.Et();
  int nJets50  = countObjects(jets, 50.);
  int nBJets20 = countObjects(jets, 20., 2.5, BTag70MV2c20);
  double meff  = sumObjectsPt(jets) + sumObjectsPt(signalLeptons) + met; // CHECK: meff uses 20 GeV jets?

  if (nLeptons>=3 && nBJets20==0 && nJets50>=3 && met>200 && meff>550) accept("SR0b3j");
  if (nLeptons>=2 && nBJets20==0 && nJets50>=5 && met>125 && meff>650) accept("SR0b5j");
  if (nLeptons>=2 && nBJets20>=1 && nJets50>=4 && met>150 && meff>550) accept("SR1b");
  if (nLeptons>=2 && nBJets20>=3 &&               met>125 && meff>650) accept("SR3b");
  if (nLeptons>=3 && nBJets20>=1 && nJets50>=4 && met>150 && meff>550) accept("SR1b3L");
  if (nLeptons>=3 && nBJets20>=1 && nJets50>=4 && met>150 && meff>550) {
    if ((lep0.charge()==lep1.charge() && lep0.charge()==signalLeptons[2].charge()) || nLeptons>=4) accept("SR1b3LSS");
  }
  return;
}
