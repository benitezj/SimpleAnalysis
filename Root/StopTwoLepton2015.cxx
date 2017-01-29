#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StopTwoLepton2015)

void StopTwoLepton2015::Init()
{
  addRegions({"SRDF","SRSF","CRT","CRSF"});
}

void StopTwoLepton2015::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons = filterCrack(event->getElectrons(10, 2.47, ELooseLH|ED0Sigma5));
  auto baselineMuons = event->getMuons(10, 2.5, MuMedium|MuQoPSignificance);
  auto jets          = event->getJets(20., 2.8, JVT50Jet);
  auto metVec        = event->getMET();

  // SUSY overlap removal      
  jets               = overlapRemoval(jets, baselineElectrons, 0.2);
  jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks);
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, jets, 0.4);

  //FIXME: check isolation
  auto signalElectrons = filterObjects(baselineElectrons, 10, 2.0, EMediumLH|EZ05mm|EIsoGradientLoose);
  auto signalMuons     = filterObjects(baselineMuons, 10, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  
  AnalysisObjects signalLeptons = signalElectrons+signalMuons;

  int nLeptons = signalLeptons.size();
  if (nLeptons!=2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];
  if (lep0.Pt()<25) return;
  if (lep1.Pt()<15) return;
  if (lep0.charge()!=lep1.charge()) return;

  // Variables we'll be cutting on
  double met   = metVec.Et();
  double meff  = sumObjectsPt(jets, 2, 50) + sumObjectsPt(signalLeptons) + met; 
  double mll   = (lep0 + lep1).M();
  double mt2   = calcMT2(lep0,lep1,metVec);
  double R1    = met/meff;
  double pTbll = (metVec + lep0 + lep1).Pt();
  
  if (mll<20) return;
  if (lep0.type()!=lep1.type() && mt2>145 && R1>0.3) accept("SRDF");
  if (lep0.type()==lep1.type() && mt2>145 && R1>0.3 && (mll<71 || mll>111)) accept("SRSF");
  if (lep0.type()!=lep1.type() && mt2>60 && mt2<110 && pTbll>30 && R1>0.4) accept("CRT");
  if (lep0.type()==lep1.type() && mll>71 && mll<111 && mt2>110 && R1>0.3) accept("CRSF");

  return;
}
