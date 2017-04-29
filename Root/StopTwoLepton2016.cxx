#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StopTwoLepton2016)

void StopTwoLepton2016::Init()
{
  addRegions({"SRASF120","SRADF120","SRASF140","SRADF140","SRASF160","SRADF160","SRASF180","SRADF180","SRBSF120","SRBDF120","SRBSF140","SRBDF140"});

}

void StopTwoLepton2016::ProcessEvent(AnalysisEvent *event)
{

  const float mZ = 91.2;

  auto baselineElectrons = event->getElectrons(10, 2.47, ELooseLH);
  auto baselineMuons = event->getMuons(10, 2.5, MuMedium | MuQoPSignificance);
  auto jets          = event->getJets(20., 2.8, JVT50Jet);
  auto metVec        = event->getMET();

  // SUSY overlap removal
  jets               = overlapRemoval(jets, baselineElectrons, 0.2);
  jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks);
  baselineElectrons  = overlapRemoval(baselineElectrons, jets, 0.4);
  baselineMuons      = overlapRemoval(baselineMuons, jets, 0.4);

  //FIXME: check isolation
  auto signalElectrons = filterObjects(baselineElectrons, 20, 2.47, EMediumLH | EZ05mm | EIsoGradientLoose);
  auto signalMuons     = filterObjects(baselineMuons, 20, 2.4, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);
  auto bjets           = filterObjects(jets, 25., 2.5, BTag77MV2c20);

  auto n_bjets          = countObjects(bjets, 25, 2.5);
  auto n_jets25         = countObjects(jets, 25, 2.5);

  AnalysisObjects signalLeptons = signalElectrons + signalMuons;

  int nLeptons = signalLeptons.size();
  if (nLeptons != 2) return;

  auto lep0 = signalLeptons[0];
  auto lep1 = signalLeptons[1];

  // Variables we'll be cutting on
  double MET                            = metVec.Et();
  double meff                           = sumObjectsPt(jets, 2, 25) + lep0.Pt() + lep1.Pt() + MET;
  auto pbll_TLV                         = metVec + lep0 + lep1;
  auto pll                              = lep0 + lep1;

  double MT2                            = calcMT2(lep0, lep1, metVec);;
  double R1                             = MET / meff;
  double dPhiEtmisspbll                 = fabs(metVec.DeltaPhi(pbll_TLV));

  double mll                            = (lep0 + lep1).M();

  int DX                                = fabs((2 * (lep0.Pz() + lep1.Pz())) / 13000.);

  bool isSF = false;
  if(lep0.type() == lep1.type()) isSF = true;

  //Opposite Sign leptons
  if (lep0.charge() == lep1.charge()) return;

  //Leading lepton with pT>25 GeV
  if (lep0.Pt() < 25.) return;

  //Reject low mass DY etc
  if (mll < 20.) return;

  //bC1 SRs, everything has mT2>100
  if (MT2 < 100.) return;

  if (n_bjets == 0 && isSF && mll> 111.2 && R1 > 0.3 && DX < 0.07 && MT2 > 120. && MT2 < 140.) accept("SRASF120");
  if (n_bjets == 0 && !isSF && DX < 0.07 && MT2 > 120. && MT2 < 140.) accept("SRADF120");
  if (n_bjets == 0 && isSF && mll> 111.2 && R1 > 0.3 && DX < 0.07 && MT2 > 140. && MT2 < 160.) accept("SRASF140");
  if (n_bjets == 0 && !isSF && DX < 0.07 && MT2 > 140. && MT2 < 160.) accept("SRADF140");
  if (n_bjets == 0 && isSF && mll> 111.2 && R1 > 0.3 && DX < 0.07 && MT2 > 160. && MT2 < 180.) accept("SRASF160");
  if (n_bjets == 0 && !isSF && DX < 0.07 && MT2 > 160. && MT2 < 180.) accept("SRADF160");
  if (n_bjets == 0 && isSF && mll> 111.2 && R1 > 0.3 && DX < 0.07 && MT2 > 180.) accept("SRASF180");
  if (n_bjets == 0 && !isSF && DX < 0.07 && MT2 > 180.) accept("SRADF180");

  if (n_bjets > 0 && n_jets25 > 1 && isSF && fabs(mll-mZ)>20. && dPhiEtmisspbll < 1.5 && MT2 > 120. && MT2 < 140.) accept("SRBSF120");
  if (n_bjets > 0 && n_jets25 > 1 && !isSF && dPhiEtmisspbll < 1.5 && MT2 > 120. && MT2 < 140.) accept("SRBDF120");
  if (n_bjets > 0 && n_jets25 > 1 && isSF && fabs(mll-mZ)>20. && dPhiEtmisspbll < 1.5 && MT2 > 140.) accept("SRBSF140");
  if (n_bjets > 0 && n_jets25 > 1 && !isSF && dPhiEtmisspbll < 1.5 && MT2 > 140.) accept("SRBDF140");

  return;
}
