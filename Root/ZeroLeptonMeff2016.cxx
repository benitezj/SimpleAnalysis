#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ZeroLeptonMeff2016)

void ZeroLeptonMeff2016::Init()
{
  // Meff based SR's
  addRegions({"SR_Meff_2j_1200","SR_Meff_2j_1600","SR_Meff_2j_2000","SR_Meff_2j_2400","SR_Meff_2j_2800","SR_Meff_2j_3600","SR_Meff_2j_2100","SR_Meff_3j_1300"});
  addRegions({"SR_Meff_4j_1000","SR_Meff_4j_1400","SR_Meff_4j_1800","SR_Meff_4j_2200","SR_Meff_4j_2600","SR_Meff_4j_3000","SR_Meff_5j_1700"});
  addRegions({"SR_Meff_5j_1600","SR_Meff_5j_2000","SR_Meff_5j_2600","SR_Meff_6j_1200","SR_Meff_6j_1800","SR_Meff_6j_2200","SR_Meff_6j_2600"});
  addRegions({"SR_Meff_2jB_1600","SR_Meff_2jB_2400"});
}

void ZeroLeptonMeff2016::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto muons      = event->getMuons(7, 2.7, MuMedium);
  auto jets       = event->getJets(50., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  // Reject events with bad jets
  if (countObjects(jets, 50, 2.8, NOT(LooseBadJet))!=0) return;
  //if (jets.size()>0 && jets[0].Pt()>100 && jets[0].pass(NOT(TightBadJet))) return;
  //if (jets.size()>1 && jets[1].Pt()>100 && jets[1].pass(NOT(TightBadJet))) return;

  // Standard SUSY overlap removal
  jets       = overlapRemoval(jets, electrons, 0.2);
  electrons  = overlapRemoval(electrons, jets, 0.4);
  muons      = overlapRemoval(muons, jets, 0.4);

  // not doing overlap removal between electrons and muons with identical track
  // FIXME: not doing overlap removal between electrons within DeltaR=0.5

  // auto signalMuons     = filterObjects(muons, 25, 2.7, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  // auto signalElectrons = filterObjects(electrons, 25, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  // auto goodJets        = filterObjects(jets, 50);
  // auto bjets           = filterObjects(goodJets, 50, 2.5, BTag77MV2c20);

  // baseline lepton veto
  auto leptons=electrons + muons;
  if (leptons.size() > 0) {
    return;
  }

  // preselection
  float meffIncl    = sumObjectsPt(jets) + met;
  if(met < 250) return;
  if(jets.size() < 2) return;
  if(jets[0].Pt() < 200.) return;
  if(jets[1].Pt() < 50.) return;
  if(meffIncl < 800) return;

  int Njets = jets.size();

  // Meff based analysis regions and selections
  float meff[7];
  for(int nJet=2; nJet<=6; nJet++)
    meff[nJet] = sumObjectsPt(jets, nJet) + met;
  float metSig      = met/sqrt(meffIncl - met);
  float dphiMin3    = minDphi(metVec, jets, 3);
  float dphiMinRest = minDphi(metVec, jets);
  float Ap          = aplanarity(jets);

  // Meff-2j SR's
  if(Njets >= 2){
    if(met > 250. && jets[0].Pt() > 250. && jets[1].Pt() > 250. && fabs(jets[0].Eta()) < 0.8 && fabs(jets[1].Eta()) < 0.8 &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 14 && meffIncl > 1200)
      accept("SR_Meff_2j_1200");
    if(met > 250. && jets[0].Pt() > 300. && jets[1].Pt() > 300. && fabs(jets[0].Eta()) < 1.2 && fabs(jets[1].Eta()) < 1.2 &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 && meffIncl > 1600)
      accept("SR_Meff_2j_1600");
    if(met > 250. && jets[0].Pt() > 350. && jets[1].Pt() > 350. && fabs(jets[0].Eta()) < 1.2 && fabs(jets[1].Eta()) < 1.2 &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 && meffIncl > 2000)
      accept("SR_Meff_2j_2000");
    if(met > 250. && jets[0].Pt() > 350. && jets[1].Pt() > 350. && fabs(jets[0].Eta()) < 1.2 && fabs(jets[1].Eta()) < 1.2 &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 && meffIncl > 2400)
      accept("SR_Meff_2j_2400");
    if(met > 250. && jets[0].Pt() > 350. && jets[1].Pt() > 350. && fabs(jets[0].Eta()) < 1.2 && fabs(jets[1].Eta()) < 1.2 &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 && meffIncl > 2800)
      accept("SR_Meff_2j_2800");
    if(met > 250. && jets[0].Pt() > 350. && jets[1].Pt() > 350. &&
       dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 && meffIncl > 3600)
      accept("SR_Meff_2j_3600");
    if(met > 250. && jets[0].Pt() > 600. && jets[1].Pt() > 50. &&
       dphiMin3 > 0.4 && dphiMinRest > 0.2 && metSig > 26 && meffIncl > 2100)
      accept("SR_Meff_2j_2100");
  }
  // Meff-3j SR's
  if(Njets >= 3){
    if(met > 250. && jets[0].Pt() > 700. && jets[1].Pt() > 50. && jets[2].Pt() > 50. &&
       dphiMin3 > 0.4 && dphiMinRest > 0.2 && metSig > 18 && meffIncl > 1300)
      accept("SR_Meff_3j_1300");
  }
  // Meff-4j SR's
  if(Njets >= 4){
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 100. && fabs(jets[0].Eta()) < 1.2 && fabs(jets[1].Eta()) < 1.2 &&
       fabs(jets[2].Eta()) < 1.2 && fabs(jets[3].Eta()) < 1.2 && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.3 &&
       Ap > 0.04 && meffIncl > 1000)
      accept("SR_Meff_4j_1000");
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 100. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. &&
       fabs(jets[2].Eta()) < 2. && fabs(jets[3].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.25 &&
       Ap > 0.04 && meffIncl > 1400)
      accept("SR_Meff_4j_1400");
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 100. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. &&
       fabs(jets[2].Eta()) < 2. && fabs(jets[3].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.25 &&
       Ap > 0.04 && meffIncl > 1800)
      accept("SR_Meff_4j_1800");
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 100. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. &&
       fabs(jets[2].Eta()) < 2. && fabs(jets[3].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.25 &&
       Ap > 0.04 && meffIncl > 2200)
      accept("SR_Meff_4j_2200");
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 150. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. &&
       fabs(jets[2].Eta()) < 2. && fabs(jets[3].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.2 &&
       Ap > 0.04 && meffIncl > 2600)
      accept("SR_Meff_4j_2600");
    if(met > 250. && jets[0].Pt() > 200. && jets[3].Pt() > 150. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. &&
       fabs(jets[2].Eta()) < 2. && fabs(jets[3].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && met/meff[4] > 0.2 &&
       Ap > 0.04 && meffIncl > 3000)
      accept("SR_Meff_4j_3000");
  }
  // Meff-5j SR's
  if(Njets >= 5){
    if(met > 250. && jets[0].Pt() > 700. && jets[4].Pt() > 50. && dphiMin3 > 0.4 && dphiMinRest > 0.2 && met/meff[5] > 0.3 &&
       meffIncl > 1700)
      accept("SR_Meff_5j_1700");
    if(met > 250. && jets[0].Pt() > 200. && jets[4].Pt() > 50. && dphiMin3 > 0.4 && dphiMinRest > 0.2 && met/meff[5] > 0.15 &&
       Ap > 0.08 && meffIncl > 1600)
      accept("SR_Meff_5j_1600");
    if(met > 250. && jets[0].Pt() > 200. && jets[4].Pt() > 50. && dphiMin3 > 0.4 && dphiMinRest > 0.4 && metSig > 15 &&
       meffIncl > 2000)
      accept("SR_Meff_5j_2000");
    if(met > 250. && jets[0].Pt() > 200. && jets[4].Pt() > 50. && dphiMin3 > 0.8 && dphiMinRest > 0.4 && metSig > 18 &&
       meffIncl > 2600)
      accept("SR_Meff_5j_2600");
  }
  // Meff-6j SR's
  if(Njets >= 6){
    if(met > 250. && jets[0].Pt() > 200. && jets[5].Pt() > 50. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. && fabs(jets[2].Eta()) < 2. &&
       fabs(jets[3].Eta()) < 2. && fabs(jets[4].Eta()) < 2. && fabs(jets[5].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.2 &&
       met/meff[6] > 0.25 && meffIncl > 1200)
      accept("SR_Meff_6j_1200");
    if(met > 250. && jets[0].Pt() > 200. && jets[5].Pt() > 100. && fabs(jets[0].Eta()) < 2. && fabs(jets[1].Eta()) < 2. && fabs(jets[2].Eta()) < 2. &&
       fabs(jets[3].Eta()) < 2. && fabs(jets[4].Eta()) < 2. && fabs(jets[5].Eta()) < 2. && dphiMin3 > 0.4 && dphiMinRest > 0.2 &&
       met/meff[6] > 0.2 && Ap > 0.04 && meffIncl > 1800)
      accept("SR_Meff_6j_1800");
    if(met > 250. && jets[0].Pt() > 200. && jets[5].Pt() > 100. && dphiMin3 > 0.4 && dphiMinRest > 0.2 &&
       met/meff[6] > 0.2 && Ap > 0.08 && meffIncl > 2200)
      accept("SR_Meff_6j_2200");
    if(met > 250. && jets[0].Pt() > 200. && jets[5].Pt() > 100. && dphiMin3 > 0.4 && dphiMinRest > 0.2 &&
       met/meff[6] > 0.15 && Ap > 0.08 && meffIncl > 2600)
      accept("SR_Meff_6j_2600");
  }
  // Meff boosted boson regions
  auto fatjets = reclusterJets(jets, 1., 50.);
  if(fatjets.size() >= 2){
    if(met < 250 && fatjets[0].Pt() > 200. && fatjets[1].Pt() > 200. && fatjets[0].M() > 60. && fatjets[0].M() < 110. &&
       fatjets[1].M() > 60. && fatjets[1].M() < 110. && dphiMin3 > 0.6 && dphiMinRest > 0.4 && metSig > 20.){
      if(meffIncl > 1600.)
	accept("SR_Meff_2jB_1600");
      if(meffIncl > 1400.)
	accept("SR_Meff_2jB_2400");
    }
  }

  return;
}
