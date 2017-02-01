#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ThreeBjets2016)

void ThreeBjets2016::Init()
{
  addRegions({
              "SR_Gbb_A",
              "SR_Gbb_B",
              "SR_Gbb_C",
              "SR_Gbb_D",
              "SR_Gtt_0l_A",
              "SR_Gtt_0l_B",
              "SR_Gtt_0l_C",
              "SR_Gtt_1l_A",
              "SR_Gtt_1l_B",
              "SR_Gtt_1l_C",
              "SR_0l_Hnj_Hmeff",
              "SR_0l_Hnj_Imeff",
              "SR_0l_Hnj_Lmeff",
              "SR_0l_IStR",
              "SR_0l_Inj_Imeff",
              "SR_0l_Inj_Lmeff",
              "SR_0l_Lnj_Hmeff",
              "SR_0l_Lnj_Imeff",
              "SR_0l_Lnj_Lmeff",
              "SR_1l_Hnj_Hmeff",
              "SR_1l_Hnj_Imeff",
              "SR_1l_Hnj_Lmeff",
              "SR_1l_Inj_Imeff",
              "SR_1l_Inj_Lmeff"
            });

  addRegions({
              "CR_Gbb_A",
              "CR_Gbb_B",
              "CR_Gbb_C",
              "CR_Gbb_D",
              "CR_Gtt_0l_A",
              "CR_Gtt_0l_B",
              "CR_Gtt_0l_C",
              "CR_Gtt_1l_A",
              "CR_Gtt_1l_B",
              "CR_Gtt_1l_C",
              "CR_Hnj_Hmeff",
              "CR_Hnj_Imeff",
              "CR_Hnj_Lmeff",
              "CR_IStR",
              "CR_Inj_Imeff",
              "CR_Inj_Lmeff",
              "CR_Lnj_Hmeff",
              "CR_Lnj_Imeff",
              "CR_Lnj_Lmeff"
            });

  addRegions({
              "VR1_Gtt_0l_A",
              "VR1_Gtt_0l_B",
              "VR1_Gtt_0l_C",
              "VR2_Gtt_0l_A",
              "VR2_Gtt_0l_B",
              "VR2_Gtt_0l_C",
              "VR1_Gtt_1l_A",
              "VR1_Gtt_1l_B",
              "VR1_Gtt_1l_C",
              "VR2_Gtt_1l_A",
              "VR2_Gtt_1l_B",
              "VR2_Gtt_1l_C",
              "VR1_Gbb_A",
              "VR1_Gbb_B",
              "VR1_Gbb_C",
              "VR1_Gbb_D",
              "VR2_Gbb_A",
              "VR2_Gbb_B",
              "VR2_Gbb_C",
              "VR2_Gbb_D",
              "VR0l_Lnj_Lmeff",
              "VR0l_Lnj_Imeff",
              "VR0l_Lnj_Hmeff",
              "VR0l_Inj_Lmeff",
              "VR0l_Inj_Imeff",
              "VR0l_Hnj_Lmeff",
              "VR0l_Hnj_Imeff",
              "VR0l_Hnj_Hmeff",
              "VR0l_IStR",
              "VR1l_Inj_Lmeff",
              "VR1l_Inj_Imeff",
              "VR1l_Hnj_Lmeff",
              "VR1l_Hnj_Imeff",
              "VR1l_Hnj_Hmeff"
            });

  addHistogram("meffi",20,0,2000);
  addHistogram("met",20,0,2000);
  addHistogram("mjsum",20,0,2000);
  addHistogram("mTb_min",20,0,800);
  addHistogram("mT", 20, 0, 500);
  addHistogram("dphiMin4", 50, 0, 5);
  addHistogram("dphi1jet", 50, 0, 5);

  addHistogram("n_jets", 20, -0.5, 19.5);
  addHistogram("n_bjets",20, -0.5, 19.5);
  addHistogram("n_electrons", 20, -0.5, 19.5);
  addHistogram("n_muons", 20, -0.5, 19.5);
  addHistogram("n_leptons", 20, -0.5, 19.5);

  addHistogram("GenMET",40,0,2000);
  addHistogram("GenHT",40,0,2000);
}

void ThreeBjets2016::ProcessEvent(AnalysisEvent *event)
{

  float gen_met      = event->getGenMET();
  float gen_ht       = event->getGenHT();
  int channel_number = event->getMCNumber();

  // handle HT slicing now
  if(channel_number==410000 && gen_ht>600) return;
  if(channel_number==410001 && gen_ht>600) return;
  if(channel_number==410002 && gen_ht>600) return;
  if(channel_number==410003 && gen_ht>600) return;
  if(channel_number==410004 && gen_ht>600) return;

  // baseline electrons are requested to pass the loose likelihood identification criteria
  //    and have pT > 20 GeV and |eta| < 2.47
  auto electrons  = event->getElectrons(20, 2.47, ELooseLH);
  // baseline muons are required to pass the Medium selections and to have pT > 20 GeV, |eta| < 2.5
  auto muons      = event->getMuons(20, 2.5, MuMedium);
  // small-R jets: pT > 20 GeV, |eta| < 2.8
  auto candJets   = event->getJets(20., 2.8);
  // TODO: determine if MET includes TST at all. Doesn't look like it.
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  if(countObjects(candJets, 20, 2.8, NOT(LooseBadJet))!=0) return;
  // No bad muon veto implemented

  // TODO: do we need this? Let's leave it in here as this was here before
  candJets = filterObjects(candJets, 20, 2.8, JVT50Jet);
  // reclusterJets(jets, ReclusterRadius, RCJetPtMin, RCJetSubjetRadius, RCJetPtFrac)
  auto fatJets = reclusterJets(candJets, 0.8, 100, 0.2, 0.10);

  // Overlap removal
  // TODO: double-check that Overlap Removal defined correctly here. I looked at it and it appears to match INT note definitions
  auto radiusCalcJet  = [] (const AnalysisObject& , const AnalysisObject& muon) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  auto radiusCalcElec = [] (const AnalysisObject& elec, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/elec.Pt()); };

  electrons  = overlapRemoval(electrons, muons, 0.01);
  candJets   = overlapRemoval(candJets, electrons, 0.2, NOT(BTag77MV2c10));
  electrons  = overlapRemoval(electrons, candJets, radiusCalcElec);
  candJets   = overlapRemoval(candJets, muons, radiusCalcJet, LessThan3Tracks);
  muons      = overlapRemoval(muons, candJets, radiusCalcMuon);

  // No cosmic muon veto implemented

  // require signal jets to be 30 GeV
  auto signalJets      = filterObjects(candJets, 30);
  // signal electrons are required to pass the Medium likelihood criteria and isolated using LooseTrackOnly
  auto signalElectrons = filterObjects(electrons, 20, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoBoosted);
  // signal muons are required to be isolated using LooseTrackOnly
  auto signalMuons     = filterObjects(muons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoBoosted|MuNotCosmic);
  // combine into signalLeptons for easy counting
  auto signalLeptons   = signalElectrons + signalMuons;

  // get b-jets
  auto bjets = filterObjects(signalJets, 30., 2.5, BTag77MV2c10);

  int n_bjets       = bjets.size();
  int n_jets        = signalJets.size();
  int n_leptons     = signalLeptons.size();

  // require at least 3 b-bjets and 4 signal jets
  if(n_bjets<3 || n_jets<4) return;

  // inclusive - all jets + leptons
  float meffi  = met + sumObjectsPt(signalJets) + sumObjectsPt(signalLeptons);
  // min mT of leading 3 bjets
  float mTbmin    = calcMTmin(bjets, metVec,3);
  // dphimin between leading 4 signal jets and met
  float dphiMin4  = minDphi(metVec, signalJets, 4);
  // leading lepton and met
  float mT = (signalLeptons.size()>0)? calcMT(signalLeptons[0], metVec) : 0.0;
  // sum of leading 4 reclustered jet masses
  float mjsum = sumObjectsM(fatJets, 4);
  // dPhi(j1, MET) for Gbb
  float dphi1jet = minDphi(metVec, signalJets, 1);

  fill("meffi", meffi);
  fill("met", met);
  fill("mjsum", mjsum);
  fill("mTb_min", mTbmin);
  fill("mT", mT);
  fill("dphiMin4", dphiMin4);
  fill("dphi1jet", dphi1jet);

  fill("n_jets", n_jets);
  fill("n_bjets", n_bjets);
  fill("n_electrons", signalElectrons.size());
  fill("n_muons", signalMuons.size());
  fill("n_leptons", n_leptons);

  if(n_leptons == 0 && dphiMin4 > 0.4 && n_jets >= 4         ){
    // handle Gbb 0L SR
    if(                n_bjets >= 3 &&                met > 400 && meffi > 2800               ) accept("SR_Gbb_A");
    if(mTbmin > 155 && n_bjets >= 4 &&                met > 450                               ) accept("SR_Gbb_B");
    if(mTbmin > 100 && n_bjets >= 3 &&                met > 600 &&
       signalJets[0].Pt() > 400 && signalJets[0].pass(NOT(BTag77MV2c10)) && dphi1jet > 2.5    ) accept("SR_Gbb_C");
    if(mTbmin > 90  && n_bjets >= 4 &&                met > 450 && meffi > 1600               ) accept("SR_Gbb_D");
  }

  if(n_leptons == 1 &&                n_jets >= 4 && mT < 150){
    // handle Gbb 0L CR
    if(                n_bjets >= 3 &&                met > 400 && meffi > 2500               ) accept("CR_Gbb_A");
    if(                n_bjets >= 4 &&                met > 375                               ) accept("CR_Gbb_B");
    if(                n_bjets >= 3 &&                met > 600 &&
       signalJets[0].Pt() > 400 && signalJets[0].pass(NOT(BTag77MV2c10)) && dphi1jet > 2.5    ) accept("CR_Gbb_C");
    if(                n_bjets >= 4 &&                met > 300 && meffi > 1600               ) accept("CR_Gbb_D");
  }

  if(n_leptons == 0 && dphiMin4 > 0.4){
    // handle Gbb VR1, VR2
    if(                 n_bjets>=3 && n_jets >=4 && met > 350              && meffi > 1900 && meffi < 2800) accept("VR1_Gbb_A");
    if(mTbmin >= 125 && n_bjets>=4 && n_jets >=4 && met > 350 && met < 450                                ) accept("VR1_Gbb_B");
    if(mTbmin >= 100 && n_bjets>=3 && n_jets >=4 && met > 225 && met < 600
        && signalJets[0].Pt() >= 400 && signalJets[0].pass(NOT(BTag77MV2c10)) && dphi1jet > 2.5           ) accept("VR1_Gbb_C");
    if(mTbmin >=  90 && n_bjets>=4 && n_jets >=4 && met > 250 && met < 450 && meffi > 1600                ) accept("VR1_Gbb_D");
    if(                 n_bjets>=3 && n_jets >=4 && met > 350              && meffi > 1900 && meffi < 2800) accept("VR2_Gbb_A");
    if(mTbmin >= 125 && n_bjets>=4 && n_jets >=4 && met > 350 && met < 450                                ) accept("VR2_Gbb_B");
    if(mTbmin >= 100 && n_bjets>=3 && n_jets >=4 && met > 225 && met < 600
        && signalJets[0].Pt() >= 400 && signalJets[0].pass(NOT(BTag77MV2c10)) && dphi1jet > 2.5           ) accept("VR2_Gbb_C");
    if(mTbmin >=  90 && n_bjets>=4 && n_jets >=4 && met > 250 && met < 450 && meffi > 1600                ) accept("VR2_Gbb_D");
  }

  if(n_leptons == 0 && dphiMin4 > 0.4){
    // handle Gtt 0L SR
    if(mTbmin > 60  && n_bjets >= 3 && n_jets >= 7 && met > 350 && meffi > 2600 && mjsum > 300) accept("SR_Gtt_0l_A");
    if(mTbmin > 120 && n_bjets >= 3 && n_jets >= 7 && met > 500 && meffi > 1800 && mjsum > 200) accept("SR_Gtt_0l_B");
    if(mTbmin > 120 && n_bjets >= 4 && n_jets >= 8 && met > 200 && meffi > 1000 && mjsum > 100) accept("SR_Gtt_0l_C");
  }

  if(n_leptons == 1 && mT < 150){
    // handle Gtt 0L CR
    if(                n_bjets >= 3 && n_jets >= 6 && met > 275 && meffi > 1800 && mjsum > 300) accept("CR_Gtt_0l_A");
    if(                n_bjets >= 3 && n_jets >= 6 && met > 400 && meffi > 1700 && mjsum > 200) accept("CR_Gtt_0l_B");
    if(                n_bjets >= 4 && n_jets >= 7 && met > 250 && meffi > 1000 && mjsum > 100) accept("CR_Gtt_0l_C");
  }

  if(n_leptons == 1 && mT < 150){
    // handle Gtt 0L VR1
    if(mTbmin > 60 && n_bjets >= 3 && n_jets >= 6 && met > 300 && meffi > 1800 && mjsum < 300) accept("VR1_Gtt_0l_A");
    if(mTbmin > 80 && n_bjets >= 3 && n_jets >= 6 && met > 450 && meffi > 1400 && mjsum < 200) accept("VR1_Gtt_0l_B");
    if(mTbmin > 80 && n_bjets >= 4 && n_jets >= 7 && met > 225 && meffi >  850 && mjsum < 100) accept("VR1_Gtt_0l_C");
  }

  if(n_leptons == 0 && dphiMin4 > 0.4){
    // handle Gtt 0L VR2
    if(               n_bjets >= 3 && n_jets >= 6 && met > 250 && meffi > 2000 && mjsum < 300) accept("VR2_Gtt_0l_A");
    if(               n_bjets >= 3 && n_jets >= 6 && met > 450 && meffi > 1400 && mjsum < 200) accept("VR2_Gtt_0l_B");
    if(               n_bjets >= 4 && n_jets >= 7 && met > 250 && meffi > 1000 && mjsum < 100) accept("VR2_Gtt_0l_C");
  }

  if(n_leptons >= 1 && mT > 150){
    // handle Gtt 1L SR
    if(mTbmin > 120 && n_bjets >= 3 && n_jets >= 5 && met > 500 && meffi > 2200 && mjsum > 200) accept("SR_Gtt_1l_A");
    if(mTbmin > 160 && n_bjets >= 3 && n_jets >= 6 && met > 450 && meffi > 1800 && mjsum > 200) accept("SR_Gtt_1l_B");
    if(mTbmin > 160 && n_bjets >= 3 && n_jets >= 7 && met > 350 && meffi > 1000               ) accept("SR_Gtt_1l_C");
  }

  if(n_leptons >= 1 && mT < 150){
    // handle Gtt 1L CR
    if(                n_bjets >= 3 && n_jets == 5 && met > 300 && meffi > 1700 && mjsum > 150) accept("CR_Gtt_1l_A");
    if(                n_bjets >= 3 && n_jets == 6 && met > 400 && meffi > 1500 && mjsum > 100) accept("CR_Gtt_1l_B");
    if(                n_bjets >= 3 && n_jets == 7 && met > 350 && meffi > 1000               ) accept("CR_Gtt_1l_C");
  }

  if(n_leptons >= 1 && mT > 150){
    // handle Gtt 1L VR1
    if(                met > 300 && meffi > 1600 && n_jets >=5 && n_bjets >= 3 && mjsum < 200.0) accept("VR1_Gtt_1l_A");
    if(                met > 250 && meffi > 1200 && n_jets >=6 && n_bjets >= 3 && mjsum < 100.0) accept("VR1_Gtt_1l_B");
    if(mTbmin < 160 && met > 300 && meffi > 1000 && n_jets >=7 && n_bjets >= 3                 ) accept("VR1_Gtt_1l_C");
  }

  if(n_leptons >= 1 && mT < 150){
    // handle Gtt 1L VR2
    if(mTbmin > 120 && met > 400 && meffi > 1400 && n_jets > 5 && n_bjets >= 3 && mjsum > 200) accept("VR2_Gtt_1l_A");
    if(mTbmin > 140 && met > 350 && meffi > 1200 && n_jets > 6 && n_bjets >= 3 && mjsum > 150) accept("VR2_Gtt_1l_B");
    if(mTbmin > 160 && met > 300 && meffi > 1000 && n_jets > 7 && n_bjets >= 3               ) accept("VR2_Gtt_1l_C");
  }

  if(n_leptons == 0 && dphiMin4 > 0.4){
    // handle Multi-Channel 0L SR
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && met > 350 && meffi > 800  && meffi < 1400                && signalJets[3].Pt() > 90) accept("SR_0l_Lnj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && met > 350 && meffi > 1400 && meffi < 2400                && signalJets[3].Pt() > 90) accept("SR_0l_Lnj_Imeff");
    if(                n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && met > 300 && meffi > 2400                                && signalJets[3].Pt() > 90) accept("SR_0l_Lnj_Hmeff");

    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 7 && n_jets <= 8 && met > 300 && meffi > 800  && meffi < 1600                                          ) accept("SR_0l_Inj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 7 && n_jets <= 8 && met > 300 && meffi > 1600 && meffi < 2500 && mjsum > 150                           ) accept("SR_0l_Inj_Imeff");

    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 9                && met > 300 && meffi > 900  && meffi < 1800                                          ) accept("SR_0l_Hnj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 9                && met > 300 && meffi > 1800 && meffi < 2500 && mjsum > 150                           ) accept("SR_0l_Hnj_Imeff");
    if(mTbmin > 100 && n_bjets >= 3 && n_jets >= 7                && met > 400 && meffi > 2500                 && mjsum > 200                           ) accept("SR_0l_Hnj_Hmeff");

    if(mTbmin > 100 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 8 && met > 600                 && meffi < 2200
                    && signalJets[0].pass(BTag77MV2c10)           && dphi1jet > 2.9            && signalJets[0].Pt() > 400                              ) accept("SR_0l_IStR");
  }

  if(n_leptons == 0 && dphiMin4 > 0.4){
    // handle Multi-Channel 0L VR
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && meffi >  800 && meffi < 1250 && signalJets[3].Pt()>90             ) accept("VR0l_Lnj_Lmeff");
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && meffi > 1250 && meffi < 1800 && signalJets[3].Pt()>90             ) accept("VR0l_Lnj_Imeff");
    if(                           met > 200              && n_bjets >= 3 && n_jets >= 4 && n_jets <= 6 && meffi > 2000 && meffi < 2400 && (signalJets[3].Pt()<90 || met<300)) accept("VR0l_Lnj_Hmeff");
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 7 && n_jets <= 8 && meffi >  800 && meffi < 1450                                      ) accept("VR0l_Inj_Lmeff");
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 7 && n_jets <= 8 && meffi > 1450 && meffi < 2000                                      ) accept("VR0l_Inj_Imeff");
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 9                && meffi >  900 && meffi < 1650                                      ) accept("VR0l_Hnj_Lmeff");
    if((mTbmin<140 || met<300) && met > 225              && n_bjets >= 3 && n_jets >= 9                && meffi > 1650 && meffi < 2100                                      ) accept("VR0l_Hnj_Imeff");
    if(                           met > 200 && met < 400 && n_bjets >= 3 && n_jets >= 7                && meffi > 2500                                                      ) accept("VR0l_Hnj_Hmeff");
    if(mTbmin > 100            && met > 250 && met < 600 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 8                 && meffi < 2000
        && dphi1jet>2.9 && signalJets[0].pass(BTag77MV2c10) && signalJets[0].Pt() >= 400                                                                                    ) accept("VR0l_IStR");
  }

  if(n_leptons >= 1 && mT > 150){
    // handle Multi-Channel 1L SR
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && met > 300 && meffi > 800  && meffi < 1600                                          ) accept("SR_1l_Inj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && met > 300 && meffi > 1600 && meffi < 2300 && mjsum > 150                           ) accept("SR_1l_Inj_Imeff");

    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 8                && met > 300 && meffi > 900  && meffi < 1800                                          ) accept("SR_1l_Hnj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 8                && met > 300 && meffi > 1800 && meffi < 2300 && mjsum > 150                           ) accept("SR_1l_Hnj_Imeff");
    if(mTbmin > 120 && n_bjets >= 3 && n_jets >= 6                && met > 500 && meffi > 2300                 && mjsum > 200                           ) accept("SR_1l_Hnj_Hmeff");
  }

  if(n_leptons >= 1 && mT > 150){
    // handle Multi-Channel 1L VR
    if(mTbmin < 140            && met > 300              && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && meffi >  800 && meffi < 1450              ) accept("VR1l_Inj_Lmeff");
    if(mTbmin < 140            && met > 225              && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && meffi > 1450 && meffi < 2000              ) accept("VR1l_Inj_Imeff");
    if(mTbmin < 140            && met > 225              && n_bjets >= 3 && n_jets >= 8                && meffi >  900 && meffi < 1650              ) accept("VR1l_Hnj_Lmeff");
    if((mTbmin<140 || met<300) && met > 200              && n_bjets >= 3 && n_jets >= 8                && meffi > 1600 && meffi < 2100              ) accept("VR1l_Hnj_Imeff");
    if(                           met > 200 && met < 500 && n_bjets >= 3 && n_jets >= 6                && meffi > 2100 && (meffi>2300 || mTbmin<140)) accept("VR1l_Hnj_Hmeff");
  }

  if(n_leptons >= 1 && mT < 150){
    // handle Multi-Channel CRs
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 5 && met > 300 && meffi > 800  && meffi < 1400                && signalJets[3].Pt() > 70) accept("CR_Lnj_Lmeff");
    if(mTbmin > 140 && n_bjets >= 3 && n_jets >= 4 && n_jets <= 5 && met > 300 && meffi > 1400 && meffi < 2000                && signalJets[3].Pt() > 70) accept("CR_Lnj_Imeff");
    if(                n_bjets >= 3 && n_jets >= 4 && n_jets <= 5 && met > 200 && meffi > 2100                                && signalJets[3].Pt() > 30) accept("CR_Lnj_Hmeff");

    if(mTbmin > 130 && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && met > 300 && meffi > 800  && meffi < 1600                                          ) accept("CR_Inj_Lmeff");
    if(mTbmin > 110 && n_bjets >= 3 && n_jets >= 6 && n_jets <= 7 && met > 200 && meffi > 1600 && meffi < 2100 && mjsum > 150                           ) accept("CR_Inj_Imeff");

    if(mTbmin > 130 && n_bjets >= 3 && n_jets >= 8                && met > 250 && meffi > 900  && meffi < 1700                                          ) accept("CR_Hnj_Lmeff");
    if(mTbmin > 60  && n_bjets >= 3 && n_jets >= 8                && met > 200 && meffi > 1700 && meffi < 2100 && mjsum > 150                           ) accept("CR_Hnj_Imeff");
    if(mTbmin > 60  && n_bjets >= 3 && n_jets >= 6                && met > 300 && meffi > 2100                 && mjsum > 150                           ) accept("CR_Hnj_Hmeff");

    if(                n_bjets >= 3 && n_jets >= 4 && n_jets <= 7 && met > 400                 && meffi < 2000
                    && signalJets[0].pass(BTag77MV2c10)           && dphi1jet > 2.9            && signalJets[0].Pt() > 400                              ) accept("CR_IStR");
  }


  ntupVar("gen_met", gen_met);
  ntupVar("gen_ht", gen_ht);
  ntupVar("channel_number", channel_number);

  ntupVar("meffi", meffi);
  ntupVar("met", met);
  ntupVar("mjsum", mjsum);
  ntupVar("mTb_min", mTbmin);
  ntupVar("mT", mT);
  ntupVar("dphiMin4", dphiMin4);
  ntupVar("dphi1jet", dphi1jet);

  ntupVar("n_jets", n_jets);
  ntupVar("n_bjets", n_bjets);
  ntupVar("n_electrons", static_cast<int>(signalElectrons.size()));
  ntupVar("n_muons", static_cast<int>(signalMuons.size()));
  ntupVar("n_leptons", n_leptons);

  ntupVar("signalJets0_pass_BTag77MV2c10", signalJets[0].pass(BTag77MV2c10));
  ntupVar("signalJets0_Pt", signalJets[0].Pt());
  ntupVar("signalJets3_Pt", signalJets[3].Pt());

  return;
}
