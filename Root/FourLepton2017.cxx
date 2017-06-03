#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(FourLepton2017)

void FourLepton2017::Init() {
  addRegions({"SR0A","SR0B","SR0C","SR0D"});

}

void FourLepton2017::ProcessEvent(AnalysisEvent *event) {
    auto hardElectrons = event->getElectrons(7, 2.47, EMediumLH);
    auto softElectrons = event->getElectrons(7, 2.47, EVeryLooseLH);
    auto hardMuons = event->getMuons(5, 2.5, MuMedium);
    auto softMuons = event->getMuons(5, 2.5, MuMedium);
    auto preJets = event->getJets(20., 4.5);
    auto metVec = event->getMET();
    float met = metVec.Pt();

    // Reject events with bad jets
    if (countObjects(preJets, 20, 4.5, NOT(LooseBadJet)) != 0) return;

    // hard lepton overlap removal and signal lepton/jet definitions
    auto hardJets = overlapRemoval(preJets, hardElectrons, 0.2);
    hardElectrons = overlapRemoval(hardElectrons, hardJets, 0.4);
    hardElectrons = overlapRemoval(hardElectrons, hardMuons, 0.01);
    hardJets = overlapRemoval(hardJets, hardMuons, 0.4, LessThan3Tracks);
    hardMuons = overlapRemoval(hardMuons, hardJets, 0.4);
    hardJets = filterObjects(hardJets, 30, 2.8, JVT50Jet);
    auto hardBJets = filterObjects(hardJets, 30, 2.5, BTag77MV2c10);
    auto hardSignalElectrons = filterObjects(hardElectrons, 35, 2.47, ETightLH | ED0Sigma5 | EZ05mm | EIsoGradientLoose);
    auto hardSignalMuons = filterObjects(hardMuons, 35, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);


    // Soft lepton overlap removal and signal lepton/jet definitions
    auto softJets = overlapRemoval(preJets, softElectrons, 0.2, NOT(BTag77MV2c10));
    softElectrons = overlapRemoval(softElectrons, softJets, 0.4);
    softElectrons = overlapRemoval(softElectrons, softMuons, 0.01);
    softJets = overlapRemoval(softJets, softMuons, 0.4, LessThan3Tracks);
    softMuons = overlapRemoval(softMuons, softJets, 0.4);
    softJets = filterObjects(softJets, 30, 2.8, JVT50Jet);
    auto softBJets = filterObjects(softJets, 30, 2.5, BTag77MV2c10);
    auto softSignalElectrons = filterObjects(softElectrons, 7, 2.47, ETightLH | ED0Sigma5 | EZ05mm | EIsoGradientLoose);
    auto softSignalMuons = filterObjects(softMuons, 6, 2.5, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);

    // Soft Lepton signal and control regions
    auto softLeptons = softElectrons + softMuons;
    auto softSignalLeptons = softSignalElectrons + softSignalMuons;
    if (softLeptons.size() == 1 && softSignalLeptons.size() == 1 && softSignalLeptons[0].Pt() < (35 > (5 * softJets.size()) ? 5 * (softJets.size()) : 35)) {
        float mt = calcMT(softSignalLeptons[0], metVec);
        float ht = sumObjectsPt(softJets, 999, 30);
        float meff = ht + met + softSignalLeptons[0].Pt();
        //float Ap   = aplanarity(softJets);

        // b-Tag signal region
        if (softBJets.size() > 0) {
            if (softJets.size() >= 2 && met > 430 && mt > 100 && met / meff > 0.25 && meff > 700) {
                if (meff < 1100) accept("SRs2j1bT");
                else if (meff < 1500) accept("SRs2j2bT");
                else if (meff < 1900) accept("SRs2j3bT");
                else accept("SRs2j4bT");
            }
        }
        // b-Veto signal region
        else {
            if (softJets.size() >= 2 && met > 430 && mt > 100 && met / meff > 0.25 && meff > 700) {
                if (meff < 1100) accept("SRs2j1bV");
                else if (meff < 1500) accept("SRs2j2bV");
                else if (meff < 1900) accept("SRs2j3bV");
                else accept("SRs2j4bV");
            }
        }
        // b-inclusive discovery SR
        if (softJets.size() >= 2 && met > 430 && mt > 100 && met / meff > 0.25 && meff > 700) {
            if (meff > 1100) accept("SRs2j");
        }
    }

    return;
}
