#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(OneLeptonMultiJets2015)

void OneLeptonMultiJets2015::Init()
{
  addRegions({"preSelection","8j0b","8j4b","9j0b","9j4b","10j0b","10j4b"});
}

void OneLeptonMultiJets2015::ProcessEvent(AnalysisEvent *event)
{
  auto hardElectrons = event->getElectrons(35,2.47);
  auto hardMuons     = event->getMuons(35,2.4);
  auto preJets       = event->getJets(40.,2.5);

  // overlap removal and signal lepton/jet definitions
  auto goodJets       = overlapRemoval(preJets,hardElectrons,0.2);
  auto bjets          = filterObjects(goodJets,40.,2.5,true);

  hardMuons = overlapRemoval(hardMuons,goodJets,0.4);
  hardElectrons = overlapRemoval(hardElectrons,goodJets,0.4);

  // Signal regions
  auto hardLeptons=hardElectrons+hardMuons;
  if ( hardLeptons.size()>=1 && goodJets.size()>5 && goodJets[0].Pt()>100 && goodJets[1].Pt()>75 ) {
    accept("preSelection");
    if ( bjets.size()>=4 ) {
      if ( goodJets.size()>=8 )  accept("8j4b");
      if ( goodJets.size()>=9 )  accept("9j4b");
      if ( goodJets.size()>=10 ) accept("10j4b");
    } else if ( bjets.size()==0 ) {
      if ( goodJets.size()>=8 )  accept("8j0b");
      if ( goodJets.size()>=9 )  accept("9j0b");
      if ( goodJets.size()>=10 ) accept("10j0b");
    }
  }
  return;
}
