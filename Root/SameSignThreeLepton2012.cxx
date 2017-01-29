#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(SameSignThreeLepton2012)

void SameSignThreeLepton2012::Init()
{
  addRegions({"SR0b","SR1b","SR3b","SR3Llow","SR3Lhigh",
		       "SR0b_disc","SR1b_disc","SR3b_disc","SR3Llow_disc","SR3Lhigh_disc"});
}

void SameSignThreeLepton2012::ProcessEvent(AnalysisEvent *event)
{
  auto baselineElectrons  = event->getElectrons(10, 2.47);
  // FIXME: check if it's better to use muons with |eta|<2.4
  auto baselineMuons      = event->getMuons(10, 2.5);
  auto jets_eta28 = event->getJets(20., 2.8);
  auto metVec     = event->getMET();

  // Standard SUSY overlap removal      
  jets_eta28         = overlapRemoval(jets_eta28, baselineElectrons, 0.2);
  baselineElectrons  = overlapRemoval(baselineElectrons, jets_eta28, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, jets_eta28, 0.4);
  baselineElectrons  = overlapRemoval(baselineElectrons, baselineMuons, 0.1);

  if ((baselineElectrons.size() + baselineMuons.size()) < 2) return;

  auto signalElectrons = filterObjects(baselineElectrons, 15.);
  auto signalMuons = filterObjects(baselineMuons, 15.);
  
  AnalysisObjects signalLeptons = signalElectrons+signalMuons;

  int nLep15 = signalLeptons.size();
  if (nLep15 < 2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];

  if ((lep0 + lep1).M() < 12.) return;
  if ((lep0.charge() != lep1.charge()) && nLep15 < 3) return;
  if (lep0.Pt() < 20.) return;

  // Variables we'll be cutting on
  double met = metVec.Et();
  int nJets40 = countObjects(jets_eta28, 40., 2.8);
  int nBJets20 = countObjects(jets_eta28, 20., 2.5, true);
  double mt = calcMT(lep0, metVec);

  double meff = sumObjectsPt(jets_eta28,999,40) + sumObjectsPt(signalLeptons) + met;
  double mllSFOS1 = -1.;
  double mllSFOS2 = -1.;

  bool foundOneSFOS = false;
  if (lep0.type() == lep1.type() && lep0.charge() != lep1.charge()) {
    mllSFOS1 = (lep0 + lep1).M();
    foundOneSFOS = true;
  }
  if (nLep15 > 2) {
    auto lep2 = signalLeptons[2];
    if (lep0.type() == lep2.type() && lep0.charge() != lep2.charge()) {
      if (foundOneSFOS) {
        mllSFOS2 = (lep0 + lep2).M();
      } else {
        mllSFOS1 = (lep0 + lep2).M();
      }
      foundOneSFOS = true;
    }
    if (lep1.type() == lep2.type() && lep1.charge() != lep2.charge()) {
      if (foundOneSFOS) {
        mllSFOS2 = (lep1 + lep2).M();
      } else {
        mllSFOS1 = (lep1 + lep2).M();
      }
      foundOneSFOS = true;
    }
  }

  bool passSR3b     = ( nLep15>=2 && nJets40>= 5 && nBJets20>=3 );
  bool passSR0b     = ( nLep15==2 && met>150 && nJets40>=3 && nBJets20==0 && mt>100 );
  bool passSR1b     = ( nLep15==2 && met>150 && nJets40>=3 && nBJets20>=1 && mt>100 && !passSR3b ); 
  bool passSR3Llow  = ( nLep15>=3 && nJets40>=4 && met > 50 && met<150 && (mllSFOS1 < 84 || mllSFOS1 > 98) && (mllSFOS2 < 84|| mllSFOS2 > 98)  && !passSR3b );
  bool passSR3Lhigh = ( nLep15>=3 && nJets40>=4 && met>150 && !passSR3b );

  bool passSR3b_disc     = passSR3b     && meff > 350;
  bool passSR0b_disc     = passSR0b     && meff > 400;
  bool passSR1b_disc     = passSR1b     && meff > 700;
  bool passSR3Llow_disc  = passSR3Llow  && meff > 400;
  bool passSR3Lhigh_disc = passSR3Lhigh && meff > 400;

  if (passSR0b) accept("SR0b");
  if (passSR1b) accept("SR1b");
  if (passSR3b) accept("SR3b");
  if (passSR3Llow) accept("SR3Llow");
  if (passSR3Lhigh) accept("SR3Lhigh");

  if (passSR0b_disc) accept("SR0b_disc");
  if (passSR1b_disc) accept("SR1b_disc");
  if (passSR3b_disc) accept("SR3b_disc");
  if (passSR3Llow_disc) accept("SR3Llow_disc");
  if (passSR3Lhigh_disc) accept("SR3Lhigh_disc");
}
