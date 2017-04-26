#include "SimpleAnalysis/AnalysisClass.h"
#include <TSystem.h>

DefineAnalysis(Stoph2016)

void Stoph2016::Init()
{
  addRegions({"SR1l4bA", "SR1l4bB", "SR1l4bC"});//define regions here
}

void Stoph2016::ProcessEvent(AnalysisEvent *event)
{

  auto electrons = event->getElectrons(10, 2.47, ELooseLH); //filter on pT, eta, ID
  auto muons     = event->getMuons(10, 2.5, MuMedium);  //filter on pT, eta ID
  auto jets      = event->getJets(25., 2.8);  //filter on pT, eta
  auto metVec    = event->getMET();

  // SUSY overlap removal
  jets               = overlapRemoval(jets, electrons, 0.2);
  jets               = overlapRemoval(jets, muons, 0.4, LessThan3Tracks);
  electrons  = overlapRemoval(electrons, jets, 0.4);
  muons      = overlapRemoval(muons, jets, 0.4);

  // object counting
  auto signalElectrons = filterObjects(electrons, 10, 2.0, EMediumLH | EZ05mm | EIsoGradientLoose);
  auto signalMuons     = filterObjects(muons, 10, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);
  auto signalJets      = filterObjects(jets, 30, 2.5);
  auto bjets           = filterObjects(jets, 30., 2.5, BTag77MV2c10);

  auto nbjet            = countObjects(bjets, 30, 2.5);
  auto njet30           = countObjects(signalJets, 30, 2.5);
  auto njet60           = countObjects(signalJets, 60, 2.5);

  double MET     = metVec.Et();

  AnalysisObjects signalLeptons = signalElectrons + signalMuons;

  int nlep = signalLeptons.size();

  if (nlep > 0 && nlep < 3 && signalJets.size() > 1 && bjets.size()>3 ) {
    auto lep1pT = signalLeptons[0].Pt();
    double bjet1pT = bjets[0].Pt();

    auto mT = sqrt(2. * lep1pT * MET * ( 1. - cos(signalLeptons[0].DeltaPhi(metVec) ) ) );

    double mbb_DR = 0.;
    double ptbb_DR = 0.;
    float DRbb = 999999999.;

    for (UInt_t ijet = 0; ijet < bjets.size(); ijet++) {

      TLorentzVector bjet1;
      bjet1.SetPtEtaPhiM(bjets[ijet].Pt(), bjets[ijet].Eta(), bjets[ijet].Phi(), bjets[ijet].M());

      for (UInt_t ijet2 = 0; ijet2 < bjets.size(); ijet2++) {

        if ( ijet == ijet2 ) continue;

        TLorentzVector bjet2;
        bjet2.SetPtEtaPhiM(bjets[ijet2].Pt(), bjets[ijet2].Eta(), bjets[ijet2].Phi(), bjets[ijet2].M());

        double DRbb_12 = bjet1.DeltaR(bjet2);

        if ( DRbb_12 < DRbb ) {

          mbb_DR = (bjet1 + bjet2).M();
          ptbb_DR = (bjet1 + bjet2).Perp();
          DRbb = DRbb_12;
        }
      }
    }
    auto hadTsumTOT30 = 0.;

    for (UInt_t ijet = 0; ijet < jets.size(); ijet++) {
      if (jets[ijet].Pt() < 30. || fabs(jets[ijet].Eta() > 2.5)) //ADD last condition
        continue;

      TLorentzVector jet1;
      jet1.SetPtEtaPhiM(jets[ijet].Pt(), jets[ijet].Eta(), jets[ijet].Phi(), jets[ijet].M());

      if (jet1.Pt() > 30.) hadTsumTOT30 += jets[ijet].Pt();

    }

    bool check_SR1l4bA = false;
    bool check_SR1l4bB = false;
    bool check_SR1l4bC = false;

    if (nlep < 3 && nbjet > 3 && njet60 > 5 && mbb_DR > 95. && mbb_DR < 155. && MET > 120 && hadTsumTOT30 > 1000 && ptbb_DR > 300.) {
      check_SR1l4bA = true;
    }
    if (nlep < 3 && nbjet > 3 && njet60 > 4 && mT > 150 && MET > 150)  {
      check_SR1l4bB = true;
    }
    if (bjet1pT < 140. && mT > 125  && nlep < 3 && nbjet > 3 && njet30 > 6 && MET > 150) {
      check_SR1l4bC = true;
    }

    if (check_SR1l4bA) accept("SR1l4bA");
    if (check_SR1l4bB) accept("SR1l4bB");
    if (check_SR1l4bC) accept("SR1l4bC");

  }

  return;

}
