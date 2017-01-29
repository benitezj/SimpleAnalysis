#include "SimpleAnalysis/NtupleMaker.h"

// This is a special instance of the analysis class for writing out slimmed ntuples

void NtupleMaker::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(_minElecPt,4.2,0); // Filter on pT, eta and "ID"
  auto muons      = event->getMuons(_minMuonPt,4.2,0);
  auto taus       = event->getTaus(_minTauPt,4.2,0);
  auto photons    = event->getPhotons(_minPhotonPt,4.2,0);
  auto jets       = event->getJets(_minJetPt,5.2,0);
  auto fatjets    = event->getFatJets(_minFatJetPt,5.2,0);
  auto met        = event->getMET();

  ntupVar("el",electrons);
  ntupVar("mu",muons);
  ntupVar("tau",taus);
  ntupVar("ph",photons);
  ntupVar("jet",jets,true);
  ntupVar("fatjet",fatjets,true);
  ntupVar("met",met);
  return;
}
