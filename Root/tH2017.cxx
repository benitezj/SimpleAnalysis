#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(tH2017)

void tH2017::Init()
{
  addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{
  auto hardElectrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto softElectrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto hardMuons      = event->getMuons(6,2.5, MuMedium);
  auto softMuons      = event->getMuons(6,2.5, MuMedium);
  auto preJets       = event->getJets(20., 2.8); 
  auto metVec     = event->getMET();
  float met = metVec.Pt();

  // Reject events with bad jets
  if (countObjects(preJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // hard lepton overlap removal and signal lepton/jet definitions
  auto hardJets            = overlapRemoval(preJets, hardElectrons, 0.2, NOT(BTag77MV2c10));
       hardElectrons       = overlapRemoval(hardElectrons, hardJets, 0.4);
       hardElectrons       = overlapRemoval(hardElectrons, hardMuons, 0.01);
       hardJets            = overlapRemoval(hardJets, hardMuons, 0.4, LessThan3Tracks);
       hardMuons           = overlapRemoval(hardMuons, hardJets, 0.4);
       hardJets            = filterObjects(hardJets, 30, 2.8, JVT50Jet);
  auto hardBJets           = filterObjects(hardJets, 30, 2.5, BTag77MV2c10);
  auto hardSignalElectrons = filterObjects(hardElectrons, 35, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto hardSignalMuons     = filterObjects(hardMuons, 35, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);

  // Hard Lepton signal and control regions
  auto hardLeptons=hardElectrons+hardMuons;
  auto hardSignalLeptons=hardSignalElectrons+hardSignalMuons;
  if (hardLeptons.size()==1 && hardSignalLeptons.size()==1) {
    //float mt   = calcMT(hardSignalLeptons[0], metVec);
 
    // b-Tag signal regions
    if (hardBJets.size()==0) accept("SR1lep0b");
    if (hardBJets.size()==1) accept("SR1lep1b");
    if (hardBJets.size()==2) accept("SR1lep2b");
    if (hardBJets.size()==3) accept("SR1lep3b");
    if (hardBJets.size()==4) accept("SR1lep4b");

  }


  return;
}
