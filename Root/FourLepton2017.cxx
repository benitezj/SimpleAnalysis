#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(FourLepton2017)

const float Z_Mass = 91.2;


void FourLepton2017::Init() {
    addRegions( { "SR0A", "SR0B", "SR0C", "SR0D", "SR1", "SR2" });

}

void FourLepton2017::ProcessEvent(AnalysisEvent *event) {
    auto softElectrons = event->getElectrons(7, 2.47, EVeryLooseLH);
    auto softMuons = event->getMuons(5, 2.5, MuMedium);
    auto softTaus = event->getTaus(20, 2.47, TauOneProng | TauThreeProng);
    auto preJets = event->getJets(20., 4.5);
    auto metVec = event->getMET();
    float met = metVec.Pt();

    // Reject events with bad jets
    if (countObjects(preJets, 20, 4.5, NOT(LooseBadJet)) != 0) return;
    auto baselineTaus = overlapRemoval(softTaus, softElectrons, 0.2);
    // Impose the pt<50 GeV and not combined muon at some point
    baselineTaus = overlapRemoval(baselineTaus, softMuons,[](const AnalysisObject &t, const AnalysisObject & m){
       if (t.Pt() > 50 && ( m.pass(MuCaloTaggedOnly) || fabs(m.Eta()) > 2.5) )return 1.e5;
        return t.DeltaR(m);
    },0.2);
    auto baselineMuons = overlapRemoval(softMuons, softElectrons, 0.01, MuCaloTaggedOnly);
    auto baselineElectrons = overlapRemoval(softElectrons, baselineMuons, 0.01);
    auto baselineJets = overlapRemoval(preJets, baselineElectrons, 0.2);
    baselineElectrons = overlapRemoval(baselineElectrons, baselineJets, 0.4);

    baselineJets = overlapRemoval(baselineJets, baselineMuons, 0.4, LessThan3Tracks);
    baselineMuons = overlapRemoval(baselineMuons, baselineJets, 0.4);

    auto signalTaus = filterObjects(baselineTaus, 20, 2.47, TauMedium);
    baselineJets = overlapRemoval(baselineJets, signalTaus, 0.4);
    //Overlap removal done
    //now start the low-mass removal
    auto passLMR_electrons = lowMassRemoval(baselineElectrons + baselineMuons, [](const AnalysisObject& l, const AnalysisObject& l1) {return l.charge() != l1.charge();}, 0, 4., ELECTRON);
    auto passLMR_muons = lowMassRemoval(baselineElectrons + baselineMuons, [](const AnalysisObject& l, const AnalysisObject& l1) {return l.charge() != l1.charge();}, 0, 4., MUON);
    passLMR_electrons = lowMassRemoval(passLMR_electrons, IsSFOS, 8.4, 10.4);
    passLMR_muons = lowMassRemoval(passLMR_muons, IsSFOS, 8.4, 10.4);
    //final signal selections
    auto signalJets = filterObjects(baselineJets, 20., 2.8, JVT50Jet);
    auto signalElectrons = filterObjects(passLMR_electrons, 7, 2.47, EMediumLH | ED0Sigma5 | EZ05mm | EIsoGradientLoose);
    auto signalMuons = filterObjects(passLMR_muons, 5, 2.7, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);
    //Count the number of signal leptons
    unsigned int N_SignalEle = signalElectrons.size();
    unsigned int N_SignalMuo = signalMuons.size();
    unsigned int N_SignalTau = signalTaus.size();
    unsigned int N_SignalLep = N_SignalEle + N_SignalMuo;
    float Ht_Lep = sumObjectsPt(signalElectrons, -1, 999) + sumObjectsPt(signalMuons, -1, 999) + sumObjectsPt(signalTaus, -1, 999);
    float Ht_Jet = sumObjectsPt(signalJets, 40, 999);
    float Meff = met + Ht_Lep + Ht_Jet;
    bool ZVeto = PassZVeto(signalElectrons, signalMuons);
    std::pair<float, float> DiZ = DiZSelection(signalElectrons, signalMuons);
    if (N_SignalLep >= 4 && N_SignalTau == 0) {
        if (ZVeto && Meff > 600) accept("SR0A");
        if (ZVeto && Meff > 1100) accept("SR0B");
        if (fabs(DiZ.first - Z_Mass) < 10. && 61.2 < DiZ.second && DiZ.second < 101.2) {
            if (met > 50) accept("SR0C");
            if (met > 100) accept("SR0D");
        }
    } else if (N_SignalLep == 3 && N_SignalTau >= 1 && Meff > 700) accept("SR1");
    else if (N_SignalLep == 2 && N_SignalTau >= 2 && Meff > 650) accept("SR2");
    return;
}
