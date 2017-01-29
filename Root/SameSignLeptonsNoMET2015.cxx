#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(SameSignLeptonsNoMET2015)

void SameSignLeptonsNoMET2015::Init()
{
  addRegions({"SR4jSS","SR4j0bSS","SR4j0b3L","SR6jSS","SR6j0bSS","SR8j0bSS","SR8jSS","SR8j1bSS","SR8j3bSS"});
}

void SameSignLeptonsNoMET2015::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons = filterCrack(event->getElectrons(10, 2.47, ELooseLH|ED0Sigma5));
  auto baselineMuons = event->getMuons(10, 2.5, MuMedium);
  auto jets          = event->getJets(20., 2.8, JVT50Jet);

  // CHECK: No veto of events with bad jets?

  // Standard SUSY overlap removal      
  jets               = overlapRemoval(jets, baselineElectrons, 0.2, NOT(BTag80MV2c20));
  jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks);
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, jets, 0.4);

  // CHECK: Not sure isolation is actually the one listed
  auto signalElectrons = filterObjects(baselineElectrons, 20, 2.0, ETightLH|EZ05mm|EIsoGradientLoose);
  auto signalMuons     = filterObjects(baselineMuons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  
  AnalysisObjects signalLeptons = signalElectrons+signalMuons;

  int nLeptons = signalLeptons.size();
  if (nLeptons<2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];
  
  if (lep0.charge()!=lep1.charge() && nLeptons<3) return;
 
  // Variables we'll be cutting on
  int nJets20  = countObjects(jets, 20.);
  int nBJets20 = countObjects(jets, 20., 2.5, BTag70MV2c20);

  if (nLeptons>=2 && nJets20>=4 )                    accept("SR4jSS");
  if (nLeptons>=2 && nJets20>=4 && nBJets20==0 )     accept("SR4j0bSS");
  if (nLeptons>=3 && nJets20>=4 && nBJets20==0 )     accept("SR4j0b3L");

  if (nLeptons>=2 && nJets20>=6 )                    accept("SR6jSS");
  if (nLeptons>=2 && nJets20>=6 && nBJets20==0 )     accept("SR6j0bSS");
  if (nLeptons>=2 && nJets20>=8 && nBJets20==0 )     accept("SR8j0bSS");

  if (nLeptons>=2 && nJets20>=8 )                    accept("SR8jSS");
  if (nLeptons>=2 && nJets20>=8 && nBJets20>=1 )     accept("SR8j1bSS");
  if (nLeptons>=2 && nJets20>=8 && nBJets20>=3 )     accept("SR8j3bSS");
  
  return;
}

