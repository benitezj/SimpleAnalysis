#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(SameSignThreeLepton2016)

void SameSignThreeLepton2016::Init()
{
  addRegions({"Rpc2L2bS","Rpc2L2bH"});
  addRegions({"Rpc2Lsoft1b","Rpc2Lsoft2b"});
  addRegions({"Rpc2L0bS","Rpc2L0bH"});
  addRegions({"Rpc3L0bS","Rpc3L0bH"});
  addRegions({"Rpc2L1bS","Rpc2L1bH"});
  addRegions({"Rpc3LSS1b"});
  addRegions({"Rpv2L1bH","Rpv2L0b","Rpv2L2bH","Rpv2L2bS","Rpv2L1bS","Rpv2L1bM"});

  //addHistogram("MET",100,0,2000);

}

void SameSignThreeLepton2016::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons = filterCrack(event->getElectrons(10, 2.47, ELooseBLLH|ED0Sigma5)); 
  auto baselineMuons = event->getMuons(10, 2.5, MuMedium); 
  auto jets  = event->getJets(20., 2.8); 
  auto metVec        = event->getMET();

  // Standard SUSY overlap removal      
  auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return 0.1 + 9.6/muon.Pt(); };
  auto radiusCalcEl = [] (const AnalysisObject& electron, const AnalysisObject& ) { return 0.1 + 9.6/electron.Pt(); };

  jets               = overlapRemoval(jets, baselineElectrons, 0.2, NOT(BTag85MV2c10)); 
  jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks); 
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, radiusCalcEl); 
  baselineMuons      = overlapRemoval(baselineMuons, jets, radiusCalcMuon); 

  auto signalElectrons = filterObjects(baselineElectrons, 10, 2.0, EMediumLH|EZ05mm|EIsoFixedCutTight); 
  auto signalMuons     = filterObjects(baselineMuons, 10, 2.5, MuD0Sigma3|MuZ05mm|MuIsoFixedCutTightTrackOnly); 
  
  auto signalLeptons = signalElectrons+signalMuons;

  int nLeptons = signalLeptons.size();
  if (nLeptons<2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];
  
  if (lep0.charge()!=lep1.charge() && nLeptons<3) return;

  // Variables we'll be cutting on
  double met   = metVec.Et();
  int nJets25  = countObjects(jets, 25.);
  int nJets40  = countObjects(jets, 40.);
  int nJets50  = countObjects(jets, 50.);
  int nBJets20 = countObjects(jets, 20., 2.5, BTag70MV2c10); 
  double meff  = sumObjectsPt(jets) + sumObjectsPt(signalLeptons) + met; 
  //Variable needed for the Rpv SRs
  int nNegLep = 0;
  for(unsigned int ilep=0;ilep<signalLeptons.size();ilep++){
    if(signalLeptons[ilep].charge()<0) nNegLep++;
  }
  //Variable needed for the Rpv SRs
  int Zee = 0;
  if (signalElectrons.size()==2){
      if ((signalElectrons[0]+signalElectrons[1]).M() > 81. && (signalElectrons[0]+signalElectrons[1]).M() < 101. && signalElectrons[0].charge()==signalElectrons[1].charge())
        Zee=1;
    }
  if (signalElectrons.size()==3){
    if ((signalElectrons[0]+signalElectrons[1]).M() > 81. && (signalElectrons[0]+signalElectrons[1]).M() < 101. && signalElectrons[0].charge()==signalElectrons[1].charge())
      Zee=1;
    if ((signalElectrons[0]+signalElectrons[2]).M() > 81. && (signalElectrons[0]+signalElectrons[2]).M() < 101. && signalElectrons[0].charge()==signalElectrons[2].charge())
      Zee=1;
    if ((signalElectrons[1]+signalElectrons[2]).M() > 81. && (signalElectrons[1]+signalElectrons[2]).M() < 101. && signalElectrons[1].charge()==signalElectrons[2].charge())
      Zee=1;
  }
  if (signalElectrons.size()==4){
    if ((signalElectrons[0]+signalElectrons[1]).M() > 81. && (signalElectrons[0]+signalElectrons[1]).M() < 101. && signalElectrons[0].charge()==signalElectrons[1].charge())
      Zee=1;
    if ((signalElectrons[0]+signalElectrons[2]).M() > 81. && (signalElectrons[0]+signalElectrons[2]).M() < 101. && signalElectrons[0].charge()==signalElectrons[2].charge())
      Zee=1;
    if ((signalElectrons[1]+signalElectrons[2]).M() > 81. && (signalElectrons[1]+signalElectrons[2]).M() < 101. && signalElectrons[1].charge()==signalElectrons[2].charge())
      Zee=1;
    if ((signalElectrons[0]+signalElectrons[3]).M() > 81. && (signalElectrons[0]+signalElectrons[3]).M() < 101. && signalElectrons[0].charge()==signalElectrons[3].charge())
      Zee=1;
    if ((signalElectrons[1]+signalElectrons[3]).M() > 81. && (signalElectrons[1]+signalElectrons[3]).M() < 101. && signalElectrons[1].charge()==signalElectrons[3].charge())
      Zee=1;
    if ((signalElectrons[2]+signalElectrons[3]).M() > 81. && (signalElectrons[2]+signalElectrons[3]).M() < 101. && signalElectrons[2].charge()==signalElectrons[3].charge())
      Zee=1;
  }

  //Variable needed for the Rpc3LSS1b
  int is3LSS = 0;
  int nPosLep = 0;
  for(unsigned int ilep=0;ilep<signalLeptons.size();ilep++){
    if(signalLeptons[ilep].charge()>0) nPosLep++;
  }
  if (nPosLep>2 || nNegLep>2)
    is3LSS = 1;

  if (nLeptons>=2 && nBJets20>=2 && nJets25>=6 && met>200 && meff>600  && met/meff>0.25 && lep1.Pt()>20              ) accept("Rpc2L2bS");
  if (nLeptons>=2 && nBJets20>=2 && nJets25>=6 &&            meff>1800 && met/meff>0.15 && lep1.Pt()>20              ) accept("Rpc2L2bH");
  if (nLeptons>=2 && nBJets20>=1 && nJets25>=6 && met>100 &&              met/meff>0.3  && lep0.Pt()<100 && lep0.Pt()>20) accept("Rpc2Lsoft1b");
  if (nLeptons>=2 && nBJets20>=2 && nJets25>=6 && met>200 && meff>600  && met/meff>0.25 && lep0.Pt()<100 && lep0.Pt()>20) accept("Rpc2Lsoft2b");
  if (nLeptons>=2 && nBJets20==0 && nJets25>=6 && met>150 &&              met/meff>0.25 && lep1.Pt()>20              ) accept("Rpc2L0bS");
  if (nLeptons>=2 && nBJets20==0 && nJets40>=6 && met>250 && meff>900  &&                  lep1.Pt()>20              ) accept("Rpc2L0bH");
  if (nLeptons>=3 && nBJets20==0 && nJets40>=4 && met>600 &&                               lep1.Pt()>20              ) accept("Rpc3L0bS");
  if (nLeptons>=3 && nBJets20==0 && nJets40>=4 && met>200 && meff>1600 &&                  lep1.Pt()>20              ) accept("Rpc3L0bH");
  if (nLeptons>=2 && nBJets20>=1 && nJets25>=6 && met>150 && meff>600  && met/meff>0.25 && lep1.Pt()>20              ) accept("Rpc2L1bS");
  if (nLeptons>=2 && nBJets20>=1 && nJets25>=6 && met>250 &&              met/meff>0.2  && lep1.Pt()>20              ) accept("Rpc2L1bH");
  if (nLeptons>=3 && nBJets20>=1 &&                                       is3LSS>0      && lep1.Pt()>20     ) accept("Rpc3LSS1b");
  if (nLeptons>=2 && nBJets20>=1 && nJets50>=6 &&            meff>2200 &&                  lep1.Pt()>20              ) accept("Rpv2L1bH");
  if (nLeptons==2 && nBJets20==0 && nJets40>=6 &&            meff>1800 &&                  lep1.Pt()>20 && Zee<1     ) accept("Rpv2L0b");
  if (nLeptons>=2 && nBJets20>=2 && nJets40>=6 &&            meff>2000 &&                  lep1.Pt()>20 && Zee<1     ) accept("Rpv2L2bH");
  if (nLeptons>=2 && nBJets20>=2 && nJets50>=3 &&            meff>1200 &&                  lep1.Pt()>20 && nNegLep>=2) accept("Rpv2L2bS");
  if (nLeptons>=2 && nBJets20>=1 && nJets50>=4 &&            meff>1200 &&                  lep1.Pt()>20 && nNegLep>=2) accept("Rpv2L1bS");
  if (nLeptons>=2 && nBJets20>=1 && nJets50>=4 &&            meff>1800 &&                  lep1.Pt()>20 && nNegLep>=2) accept("Rpv2L1bM");

  return;
}
