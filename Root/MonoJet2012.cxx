#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(MonoJet2012)

void MonoJet2012::Init()
{
  for(int ii=1;ii<=9;ii++) {
    addRegion(("SR"+std::to_string(ii)).c_str());
  }
}

void MonoJet2012::ProcessEvent(AnalysisEvent *event)
{

  auto electrons  = event->getElectrons(7,2.47);
  auto muons      = event->getMuons(7,2.5);
  auto jets       = event->getJets(30.,4.5);
  auto metVec     = event->getMET();
  float met = metVec.Pt();

  // Count electrons and muons
  int baselineMuons = muons.size();
  int baselineElectrons = electrons.size();

  // Apply lepton veto
  bool passLeptonVeto = (baselineMuons+baselineElectrons)==0;
  if (!passLeptonVeto) return;
  if (jets.size()==0) return;
  if ((jets[0].Pt()<120.)||abs(jets[0].Eta())>2.0) return;

  float dphimin=minDphi(metVec,jets);
  if (dphimin<0.5) return;
  for(int ii=1;ii<=9;ii++) {
    float metCut=100.+50.*ii;
    if (ii>6) metCut+=50*(ii-6);
    if (met<metCut) break;
    if ( (jets[0].Pt()/met>0.5) && (dphimin>1.0) ) {
      accept(("SR"+std::to_string(ii)).c_str());
    }
  }
  return;
}
