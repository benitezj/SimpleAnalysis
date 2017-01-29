#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ZJets2015)

void ZJets2015::Init()
{
  addRegions({"SRZ"});
  addRegions({"CRZ","CRFS","CRT"});
}

void ZJets2015::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons = event->getElectrons(10, 2.47, ELooseLH|ED0Sigma5);
  auto baselineMuons = event->getMuons(10, 2.4, MuMedium);
  auto jets          = event->getJets(20., 2.8, JVT50Jet);
  auto metVec        = event->getMET();

  // CHECK: Veto of events with bad jets before or after JVT cut?
  if (countObjects(jets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // overlap removal      
  auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return 0.04 + 10/muon.Pt(); };
  jets               = overlapRemoval(jets, baselineElectrons, 0.2, NOT(BTag77MV2c20));
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, 0.4);
  auto candBjets     = filterObjects(jets, 20, 2.5, BTag77MV2c20);
  baselineMuons      = overlapRemoval(baselineMuons, candBjets, 0.2);
  jets               = overlapRemoval(jets, baselineMuons, 0.2);
  baselineMuons      = overlapRemoval(baselineMuons, jets, radiusCalcMuon);

  baselineElectrons  = overlapRemoval(baselineElectrons, baselineMuons, 0.01);

  auto signalElectrons = filterObjects(baselineElectrons, 25, 2.47, EMediumLH|EZ05mm|EIsoGradientLoose);
  auto signalMuons     = filterObjects(baselineMuons, 25, 2.4, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  auto signalJets      = filterObjects(jets, 30, 2.5);

  AnalysisObjects signalLeptons = signalElectrons+signalMuons;

  if (signalLeptons.size()<2 || signalJets.size()<2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];
  if (lep0.charge()==lep1.charge()) return;
  if (lep0.Pt()<50) return;
  float mll=(lep0 + lep1).M();
  bool SF = lep0.type()==lep1.type();

  float dphiMin = minDphi(metVec, signalJets, 2);
  if (dphiMin<0.4) return;

  double HT = sumObjectsPt(signalJets) + lep0.Pt() + lep1.Pt();
  if (HT<600) return;

  bool atZ     = mll>81 && mll<101;
  bool atWideZ = mll>61 && mll<121;
  float met = metVec.Et();

  if (met>225 && atZ && SF) accept("SRZ");
  if (met<60  && atZ && SF) accept("CRZ");
  if (met>225 && atWideZ && !SF) accept("CRFS");
  if (met>225 && !atZ && SF) accept("CRT");
  
  return;
}
