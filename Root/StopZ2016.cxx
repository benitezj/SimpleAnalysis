#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StopZ2016)

void StopZ2016::Init()
{
  addRegions({"SR3l1bA", "SR3l1bB", "SR3l1bC"});//define regions here
}

void StopZ2016::ProcessEvent(AnalysisEvent *event)
{
  const float mZ = 91.2;

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

  double MET     = metVec.Et();

  AnalysisObjects signalLeptons = signalElectrons + signalMuons;

  if (signalLeptons.size() > 2 && signalJets.size() > 1 ) {
    auto lep1pT = signalLeptons[0].Pt();

    double bjet1pT;
    if (bjets.size() > 0) bjet1pT = bjets[0].Pt();

    double jet1pT;
    if (signalJets.size() > 0) jet1pT = signalJets[0].Pt();

    //index and flavor of selected two leptons
    auto flavlep1 = -1;
    auto flavlep2 = -1;
    float closest_mll = -1;
    float closest_ptll = 0.;
    float closest_dphill = 0.;

    for (UInt_t ilep = 0; ilep < signalLeptons.size(); ilep++) {

      auto lep1 = signalLeptons[ilep];

      for (UInt_t ilep2 = 0; ilep2 < signalLeptons.size(); ilep2++) {

        if ( ilep == ilep2 ) continue;

        auto lep2 = signalLeptons[ilep2];

        if (std::floor(lep1.M()) == 0)   flavlep1 = 11 * lep1.charge();
        else if (std::floor(lep1.M()) == 105) flavlep1 = 13 * lep1.charge();
        if (std::floor(lep2.M()) == 0)   flavlep2 = 11 * lep2.charge();
        else if (std::floor(lep2.M()) == 105) flavlep2 = 13 * lep2.charge();

        bool isSF = false;
        if (abs(flavlep1 * flavlep2) == 121 || abs(flavlep1 * flavlep2) == 169) isSF = true;

        double mll_12 = (lep1 + lep2).M();
        double isOS_12 = flavlep1 * flavlep2;

        if (lep1.charge()*lep2.charge() < 0 && isSF && (((fabs(mll_12 - mZ) < fabs(closest_mll - mZ)) || (closest_mll == -1) ))) {

          closest_mll = mll_12;
          closest_ptll = (lep1 + lep2).Perp();
          closest_dphill = fabs(lep1.DeltaPhi(lep2));
        }
      }
    }

    bool check_SR3l1bA = false;
    bool check_SR3l1bB = false;
    bool check_SR3l1bC = false;

    if (lep1pT > 40. && jet1pT > 250. && bjet1pT > 40. && closest_mll > 76.2 && closest_mll < 106.2 && nbjet > 0 && MET > 100. && njet30 > 5 && closest_ptll > 150.) {
      check_SR3l1bA = true;
    }
    if (lep1pT > 40. && jet1pT > 80. && bjet1pT > 40. && closest_mll > 76.2 && closest_mll < 106.2 && nbjet > 0 && MET > 180. && njet30 > 5) {
      check_SR3l1bB = true;
    }
    if (lep1pT > 40. && jet1pT > 60. && closest_mll > 76.2 && closest_mll < 106.2 && nbjet > 0 && MET > 140. && njet30 > 4 && closest_ptll < 80.) {
      check_SR3l1bC = true;
    }

    if (check_SR3l1bA) accept("SR3l1bA");
    if (check_SR3l1bB) accept("SR3l1bB");
    if (check_SR3l1bC) accept("SR3l1bC");

  }

  return;
}
