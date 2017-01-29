#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(MultiJets2015)

void MultiJets2015::Init()
{
  addRegions({"SR8j50","SR8j50_1b","SR8j50_2b","SR9j50","SR9j50_1b","SR9j50_2b",
                       "SR10j50","SR10j50_1b","SR10j50_2b"});  
  addRegions({"SR7j80","SR7j80_1b","SR7j80_2b","SR8j80","SR8j80_1b","SR8j80_2b"});

  addRegions({"CR6j50","CR6j50_1b","CR6j50_2b"});
  addRegions({"CR5j80","CR5j80_1b","CR5j80_2b"});

  addRegions({"CR7j50_0b","CR7j50_1b","CR8j50_0b","CR8j50_1b","CR9j50_0b","CR9j50_1b"});
  addRegions({"CR6j80_0b","CR6j80_1b","CR7j80_0b","CR7j80_1b"});
}

void MultiJets2015::ProcessEvent(AnalysisEvent *event)
{

  auto electrons  = event->getElectrons(10, 2.47, ELooseLH);
  auto muons      = event->getMuons(10, 2.5, MuMedium);
  auto jets       = event->getJets(20., 2.8, JVT50Jet);
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  // Reject events with bad jets
  if (countObjects(jets, 20, 2.8, NOT(LooseBadJet))!=0) return;
  // FIXME: not doing event cleaning for JVT-miss jets near MET or inefficient hadronic calo

  // Standard SUSY overlap removal      
  jets       = overlapRemoval(jets, electrons, 0.2);
  jets       = overlapRemoval(jets, muons, 0.4, LessThan3Tracks);
  electrons  = overlapRemoval(electrons, jets, 0.4);
  muons      = overlapRemoval(muons, jets, 0.4);
  //not doing overlap removal between electrons and muons with identical track

  // Leptons for control regions
  auto signalMuons     = filterObjects(muons, 10, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  auto signalElectrons = filterObjects(electrons, 10, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto signalLeptons   = signalMuons+signalElectrons;

  // Count objects
  int nBaselineLeptons = muons.size() + electrons.size();

  // Count jets
  int n50   = countObjects(jets, 50., 2.0);
  int n80   = countObjects(jets, 80., 2.0);
  int nbtag = countObjects(jets, 40., 2.5, BTag70MV2c20);

  // Define HT and MetSig
  Double_t SUSY_MetSig=0.;
  Double_t SUSY_HT=sumObjectsPt(jets, 999, 40);
  if (SUSY_HT>0.) SUSY_MetSig=met/sqrt(SUSY_HT);

  if(SUSY_MetSig>4 && nBaselineLeptons==0) { // apply this cut to all 
    if (n50>=8)              accept("SR8j50");
    if (n50>=8 && nbtag>=1)  accept("SR8j50_1b");
    if (n50>=8 && nbtag>=2)  accept("SR8j50_2b");

    if (n50>=9)              accept("SR9j50");
    if (n50>=9 && nbtag>=1)  accept("SR9j50_1b");
    if (n50>=9 && nbtag>=2)  accept("SR9j50_2b");

    if (n50>=10)             accept("SR10j50");
    if (n50>=10 && nbtag>=1) accept("SR10j50_1b");
    if (n50>=10 && nbtag>=2) accept("SR10j50_2b");

    if (n80>=7)              accept("SR7j80");
    if (n80>=7 && nbtag>=1)  accept("SR7j80_1b");
    if (n80>=7 && nbtag>=2)  accept("SR7j80_2b");
    
    if (n80>= 8)             accept("SR8j80");
    if (n80>= 8 && nbtag>=1) accept("SR8j80_1b");
    if (n80>= 8 && nbtag>=2) accept("SR8j80_2b");

    // multi-jet control regions
    if (n50==6)              accept("CR6j50");
    if (n50==6 && nbtag>=1)  accept("CR6j50_1b");
    if (n50==8 && nbtag>=2)  accept("CR6j50_2b");

    if (n80==5)              accept("CR5j80");
    if (n80==5 && nbtag>=1)  accept("CR5j80_1b");
    if (n80==5 && nbtag>=2)  accept("CR5j80_2b");

  }
  if (SUSY_MetSig>3 && signalLeptons.size()==1 && signalLeptons[0].Pt()>20) {
    int n50CR=n50 + (signalLeptons[0].Pt()>50);
    int n80CR=n80 + (signalLeptons[0].Pt()>80);
    float mt=calcMT(signalLeptons[0], metVec);
    std::string type="_0b";
    if (nbtag>0) type="_1b";
    if (mt<120) {
      if (n50CR>=7) accept("CR7j50"+type);
      if (n50CR>=8) accept("CR8j50"+type);
      if (n50CR>=9) accept("CR9j50"+type);
      if (n80CR>=6) accept("CR6j80"+type);
      if (n80CR>=7) accept("CR7j80"+type);
    }
  }

  return;
}
