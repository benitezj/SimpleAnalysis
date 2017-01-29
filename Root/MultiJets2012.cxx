#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(MultiJets2012)

void MultiJets2012::Init()
{
  addRegions({"SR_8j50_0b","SR_8j50_1b","SR_8j50_2b","SR_9j50_0b","SR_9j50_1b","SR_9j50_2b","SR_10j50"});  
  addRegions({"SR_7j80_0b","SR_7j80_1b","SR_7j80_2b","SR_8j80_0b","SR_8j80_1b","SR_8j80_2b"});
  addRegions({"SR_8j50_MJ340","SR_8j50_MJ420","SR_9j50_MJ340","SR_9j50_MJ420","SR_10j50_MJ340","SR_10j50_MJ420"});
}

void MultiJets2012::ProcessEvent(AnalysisEvent *event)
{

  auto electrons  = event->getElectrons(10,2.47);
  auto muons      = event->getMuons(10,2.4);
  auto jets_eta28 = event->getJets(20.,2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();
  //  double metPhi   = metVec.Phi();

  // Standard SUSY overlap removal      
  jets_eta28 = overlapRemoval(jets_eta28,electrons,0.2);
  electrons  = overlapRemoval(electrons,jets_eta28,0.4);
  muons      = overlapRemoval(muons,jets_eta28,0.4);
   
  // Count electrons and muons
  int baselineMuons = muons.size();
  int baselineElectrons = electrons.size();

  // Define  central jets
  auto jets=filterObjects(jets_eta28,40.,2.0);

  // Count jets
  int Count50Jet_eta28 = countObjects(jets_eta28,50.);
  int Count50Jet = countObjects(jets,50.);
  int Count80Jet = countObjects(jets,80.);

  // New b-tagging count
  int btag_eta25 = countObjects(jets_eta28,40.,2.5,true);

  // Apply lepton veto
  bool passLeptonVeto = (baselineMuons+baselineElectrons)==0;

  // Define HT and MetSig
  Double_t SUSY_MetSig=0.;
  Double_t SUSY_HT=sumObjectsPt(jets_eta28,999,40);
  if (SUSY_HT>0.) SUSY_MetSig=met/sqrt(SUSY_HT);

  auto fatJets = reclusterJets(jets_eta28,1.0,100);
  fatJets = filterObjects(fatJets,100,1.5);
  double m_sigma_fat=0.;
  for (const auto& jet : fatJets) m_sigma_fat += jet.M();

  // Count50Jet is number of jets with pt > 50, |eta| < 2.0 and no overlap
  // btag_eta25 is number of b-tagged jets with pt > 40 and |eta|<2.5

  if(SUSY_MetSig > 4 && passLeptonVeto) { // apply this cut to all 
    if (Count50Jet == 8 && btag_eta25 == 0) accept("SR_8j50_0b");
    if (Count50Jet == 8 && btag_eta25 == 1) accept("SR_8j50_1b");
    if (Count50Jet == 8 && btag_eta25 >= 2) accept("SR_8j50_2b");

    if (Count50Jet == 9 && btag_eta25 == 0) accept("SR_9j50_0b");
    if (Count50Jet == 9 && btag_eta25 == 1) accept("SR_9j50_1b");
    if (Count50Jet == 9 && btag_eta25 >= 2) accept("SR_9j50_2b");

    if (Count50Jet >=10) accept("SR_10j50");

    // Now for jets with pT > 80 GeV
    // N.B Count80Jet us number of jets in event with pt > 80 and |eta|<2.0
    if (Count80Jet == 7 && btag_eta25 == 0) accept("SR_7j80_0b");
    if (Count80Jet == 7 && btag_eta25 == 1) accept("SR_7j80_1b");
    if (Count80Jet == 7 && btag_eta25 >= 2) accept("SR_7j80_2b");
    
    if (Count80Jet >= 8 && btag_eta25 == 0) accept("SR_8j80_0b");
    if (Count80Jet >= 8 && btag_eta25 == 1) accept("SR_8j80_1b");
    if (Count80Jet >= 8 && btag_eta25 >= 2) accept("SR_8j80_2b");
      
    // "fat" Jet collection
    if (Count50Jet_eta28 >= 8) {
      if(m_sigma_fat > 340) accept("SR_8j50_MJ340");
      if(m_sigma_fat > 420) accept("SR_8j50_MJ420");
    }
    if (Count50Jet_eta28 >= 9) {
      if(m_sigma_fat > 340) accept("SR_9j50_MJ340");
      if(m_sigma_fat > 420) accept("SR_9j50_MJ420");
    }
    if (Count50Jet_eta28 >= 10) {
      if(m_sigma_fat > 340) accept("SR_10j50_MJ340");
      if(m_sigma_fat > 420) accept("SR_10j50_MJ420");
    }
  }

  return;
}
