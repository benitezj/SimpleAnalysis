#include "SimpleAnalysis/AnalysisClass.h"
#include <unordered_set>

DefineAnalysis(EwkTwoLepton2012)

void EwkTwoLepton2012::Init()
{
  addRegions({"Zjets","mT2a","mT2b","mT2c","WWa","WWb","WWc"});
}

static AnalysisObjects RemoveLowMassPairs (const AnalysisObjects &cands,float massCut) {
  std::unordered_set<int> toRemove;
  for(unsigned int ii=0;ii<cands.size();ii++) {
    for(unsigned int jj=ii+1;jj<cands.size();jj++) {
      if ( (cands[ii]+cands[jj]).M()<massCut) {
	toRemove.insert(ii);
	toRemove.insert(jj);
      }
    }
  }
  AnalysisObjects goodCands;
  for(unsigned int ii=0;ii<cands.size();ii++) {
    if (toRemove.count(ii)==0) goodCands.push_back(cands[ii]);
  }
  return goodCands;
}

void EwkTwoLepton2012::ProcessEvent(AnalysisEvent *event)
{
  static int entry=-1;
  entry++;

  auto baselineElectrons  = event->getElectrons(10, 2.47);
  auto baselineMuons      = event->getMuons(10, 2.5);
  auto baselineTaus       = event->getTaus(20, 2.47); 
  auto baselineJets       = event->getJets(20., 4.5);
  auto metVec     = event->getMET();

  // Extended SUSY overlap removal      
  baselineElectrons  = overlapRemoval(baselineElectrons, baselineElectrons, 0.05);
  baselineJets       = overlapRemoval(baselineJets, baselineElectrons, 0.2);
  baselineTaus       = overlapRemoval(baselineTaus, baselineElectrons, 0.2);
  baselineTaus       = overlapRemoval(baselineTaus, baselineMuons, 0.2);
  baselineElectrons  = overlapRemoval(baselineElectrons, baselineJets, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, baselineJets, 0.4);
  
  AnalysisObjects copyElectrons = baselineElectrons;
  baselineElectrons  = overlapRemoval(baselineElectrons, baselineMuons, 0.01);
  baselineMuons      = overlapRemoval(baselineMuons, copyElectrons, 0.01);
  baselineMuons      = overlapRemoval(baselineMuons, baselineMuons, 0.01);

  baselineMuons      = RemoveLowMassPairs(baselineMuons, 12);
  baselineElectrons  = RemoveLowMassPairs(baselineElectrons, 12);

  AnalysisObjects signalLeptons = baselineElectrons+baselineMuons;
  if (signalLeptons.size() != 2) return;
  //if (baselineTaus.size() != 0) return; //FIXME: should be enabled, but was disabled in orig.
  if (signalLeptons[0].Pt() < 35 || signalLeptons[1].Pt() < 20) return;

  bool isEM = (baselineElectrons.size() == 1);

  float mll = (signalLeptons[0]+signalLeptons[1]).M();
  if (mll < 20) return;
  AnalysisObjects centralJets=filterObjects(baselineJets,20,2.4);

  int njets_B20 = countObjects(centralJets, 20., 2.4,true);
  int njets_L20 = countObjects(centralJets, 20., 2.4)-njets_B20;
  int njets_C30 = countObjects(centralJets, 30., 2.4);
  int njets_F30 = countObjects(baselineJets, 30., 4.5)-njets_C30;
  if (njets_B20+njets_F30 != 0) return;

  //calculate metRel
  AnalysisObjects allCentral=centralJets+signalLeptons;
  float dphimin=minDphi(metVec,allCentral);
  float metrel=metVec.Pt();
  if (dphimin < 3.14159265358979323846/2.) metrel *= sin(dphimin);

  float deltaMZ = fabs(mll-91.2);  
  float mt2     = calcMT2(signalLeptons[0],signalLeptons[1],metVec);
  float drll    = signalLeptons[0].DeltaR(signalLeptons[1]);
  float ptll    = (signalLeptons[0]+signalLeptons[1]).Pt();;
  float mjj     = 0;
  if (baselineJets.size() >= 2) mjj = (baselineJets[0]+baselineJets[1]).M();

  if (!isEM && njets_L20>=2 && baselineJets[1].Pt()>45 && deltaMZ<10 && mjj>50 && mjj<100 &&  
      metrel>80 && ptll>80 && drll>0.3 && drll<1.5) accept("Zjets");

  if (njets_L20 !=0) return;
  if (!isEM && deltaMZ<10) return;

  if (mt2>90)  accept("mT2a");
  if (mt2>120) accept("mT2b");
  if (mt2>150) accept("mT2c");
  
  if (metrel>80 && ptll>80 && mll<120) accept("WWa");
  if (mt2>90 && mll<170)               accept("WWb");
  if (mt2>100)                         accept("WWc");
}
