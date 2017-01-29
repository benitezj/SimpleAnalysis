#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StopOneLepton2015)

void StopOneLepton2015::Init()
{
  addRegions({"SR1","SR2","SR3"});
  addRegions({"CRT1","CRT2","CRT3"});
  addRegions({"CRW1","CRW2","CRW3"});
  addRegions({"CRST1","CRST2","CRST3"});
}

void StopOneLepton2015::ProcessEvent(AnalysisEvent *event)
{
  auto baseElectrons  = event->getElectrons(7, 2.47, EVeryLooseLH);
  auto baseMuons      = event->getMuons(6,2.7, MuLoose);
  auto baseJets       = event->getJets(20., 4.5); 
  auto baseTaus       = event->getTaus(20., 2.5, TauLoose); 
  auto metVec         = event->getMET();
  float met = metVec.Pt();

  // Reject events with bad jets
  // CHECK: right place?
  if (countObjects(baseJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // overlap removal
  auto radiusCalcLepton = [] (const AnalysisObject& lepton, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/lepton.Pt()); };
  auto muJetSpecial = [] (const AnalysisObject& jet, const AnalysisObject& muon) { 
    if (jet.pass(NOT(BTag77MV2c20)) && (jet.pass(LessThan3Tracks) || muon.Pt()/jet.Pt()>0.7)) return 0.2;
    else return 0.;
  };
  baseMuons     = overlapRemoval(baseMuons, baseElectrons, 0.01, NOT(MuCaloTaggedOnly));
  baseElectrons = overlapRemoval(baseElectrons, baseMuons, 0.01);

  baseJets      = overlapRemoval(baseJets, baseElectrons, 0.2, NOT(BTag77MV2c20));
  baseElectrons = overlapRemoval(baseElectrons, baseJets, 0.2);

  baseJets      = overlapRemoval(baseJets, baseMuons, muJetSpecial, NOT(BTag77MV2c20));
  baseMuons     = overlapRemoval(baseMuons,baseJets, 0.2);

  baseMuons     = overlapRemoval(baseMuons,baseJets, radiusCalcLepton);
  baseElectrons = overlapRemoval(baseElectrons,baseJets, radiusCalcLepton);
  
  baseTaus      = overlapRemoval(baseTaus, baseElectrons, 0.1);

  auto baseLeptons     = baseElectrons+baseMuons;
  auto signalJets      = filterObjects(baseJets, 25, 2.5, JVT50Jet);
  auto signalBJets     = filterObjects(signalJets, 25, 2.5, BTag77MV2c20);
  auto signalElectrons = filterObjects(baseElectrons, 25, 2.47, ELooseLH|ED0Sigma5|EZ05mm|EIsoLooseTrack);
  auto signalMuons     = filterObjects(baseMuons, 25, 2.7, MuD0Sigma3|MuZ05mm|MuIsoLooseTrack);
  auto signalLeptons   = signalElectrons+signalMuons;
  
  auto fatJetsR10 = reclusterJets(signalJets, 1.0, 150, 0.2, 0.05);  
  auto fatJetsR12 = reclusterJets(signalJets, 1.2, 150, 0.2, 0.05);  

  if (met<200 || signalLeptons.size()!=1 || baseLeptons.size()!=1 || signalJets.size()<4) return;

  float mt   = calcMT(signalLeptons[0], metVec);
  float dphiMin = minDphi(metVec, signalJets, 2);
  float mtauT2 = 120;
  if (baseTaus.size()>0&&baseTaus[0].charge()!=signalLeptons[0].charge()) mtauT2 = calcMT2(baseTaus[0],signalLeptons[0],metVec);

  if (mt<30 || dphiMin<0.4 || mtauT2<80) return; //CHECK: Int note says mtauT2<100 while opposite charge is not in paper

  float HTmiss_sig = 25; //FIXME
  float mchitop   = 100; //FIXME
  float topness    = 10.; //FIXME
  
  auto mostBjetLike = signalBJets;
  int idx=0;
  while(mostBjetLike.size()<2) { //FIXME: should be most b-jet like, not highest pt b-jets, but not possible here
    if (signalJets[idx].pass(NOT(BTag77MV2c20))) mostBjetLike.push_back(signalJets[idx]);
    idx++;
  }
  float amt2 = std::min(calcAMT2(mostBjetLike[0]+signalLeptons[0],mostBjetLike[1],metVec,0.,80),
			calcAMT2(mostBjetLike[1]+signalLeptons[0],mostBjetLike[0],metVec,0.,80));

  int nBjets = signalBJets.size();
  float deltaRbl = 0;
  float deltaRbb = 0;
  if (nBjets>0) deltaRbl = signalLeptons[0].DeltaR(signalBJets[0]);
  if (nBjets>1) deltaRbb = signalBJets[1].DeltaR(signalBJets[0]);
  float dphiJet12 = 2.0;
  if (fatJetsR12.size()>=2) dphiJet12 = fabs(metVec.DeltaPhi(fatJetsR12[1]));

  if (signalJets[0].Pt()>80 && signalJets[1].Pt()>50 && signalJets[3].Pt()>40) { //region 1
    if (met>260 && HTmiss_sig>14 && mt>170 && amt2>175 && topness>6.5 && 
        mchitop<270 && deltaRbl<3.0 && nBjets>=1) accept("SR1"); 
    if (met>200 && HTmiss_sig>5 && mt>30 && mt<90 && amt2>100 && amt2<200 && topness>6.5 && 
        mchitop<270 && nBjets>=1) accept("CRT1"); 
    if (met>200 && HTmiss_sig>5 && mt>30 && mt<90 && amt2>100 && topness>6.5 && 
        mchitop<270 && nBjets==0) accept("CRW1"); 
    if (met>200 && HTmiss_sig>5 && mt>30 && mt<120 && amt2>200 && topness>6.5 && 
        mchitop<270 && deltaRbb>1.2 && nBjets>=2) accept("CRST1"); 
  }

  if (signalJets[0].Pt()>120 && signalJets[1].Pt()>80 && 
      signalJets[2].Pt()>50 && signalJets[3].Pt()>25 && fatJetsR12.size()>0) { //region 2
    //CHECK: mt not in line from Javier
    if (met>350 && HTmiss_sig>20 && mt>200 && amt2>175 && deltaRbl<2.5 && nBjets>=1 &&
        fatJetsR12[0].Pt()>200 && fatJetsR12[0].M()>140 && dphiJet12>1) accept("SR2"); 
    if (met>250 && HTmiss_sig>15 && mt>30 && mt<90 && amt2>100 && amt2<200 && nBjets>=1 &&
        fatJetsR12[0].Pt()>200 && fatJetsR12[0].M()>140 && dphiJet12>1) accept("CRT2"); 
    if (met>250 && HTmiss_sig>15 && mt>30 && mt<90 && amt2>100 && nBjets==0 &&
        fatJetsR12[0].Pt()>200 && fatJetsR12[0].M()>140 && dphiJet12>1) accept("CRW2"); 
    if (met>200 && HTmiss_sig>15 && mt>30 && mt<120 && amt2>200 && deltaRbb>1.2 && nBjets>=2 &&
        fatJetsR12[0].Pt()>200 && fatJetsR12[0].M()>0 && dphiJet12>1) accept("CRST2"); 
  }
  if (signalJets[0].Pt()>120 && signalJets[1].Pt()>80 && 
      signalJets[2].Pt()>50 && signalJets[3].Pt()>25 && fatJetsR10.size()>0) { //region 3
    if (met>480 && HTmiss_sig>14 && mt>190 && amt2>175 && topness>9.5 && deltaRbl<2.8 && 
        nBjets>=1 && fatJetsR10[0].Pt()>280 && fatJetsR10[0].M()>70) accept("SR3"); 
    if (met>280 && HTmiss_sig>8 && mt>30 && mt<90 && amt2>100 && amt2<200 && topness>0 && 
        nBjets>=1 && fatJetsR10[0].Pt()>280 && fatJetsR10[0].M()>70) accept("CRT3"); 
    if (met>280 && HTmiss_sig>8 && mt>30 && mt<90 && amt2>100 && topness>0 && 
        nBjets==0 && fatJetsR10[0].Pt()>280 && fatJetsR10[0].M()>70) accept("CRW3"); 
    if (met>200 && HTmiss_sig>5 && mt>30 && mt<120 && amt2>200 && topness>9.5 && deltaRbb>1.2 &&
        nBjets>=2 && fatJetsR10[0].Pt()>280 && fatJetsR10[0].M()>70) accept("CRST3"); 
  }

  return;
}
