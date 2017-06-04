#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(FourLepton2017)

void FourLepton2017::Init() {
    addRegions( { "SR0A", "SR0B", "SR0C", "SR0D" });

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
    baselineTaus = overlapRemoval(baselineTaus, softMuons, 0.2);
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
    auto signalJets = filterObjects(baselineJets, 20.,2.8,JVT50Jet);
    auto signalElectrons = filterObjects(passLMR_electrons,7,2.47, EMediumLH | ED0Sigma5 | EZ05mm | EIsoGradientLoose );
    auto signalMuons = filterObjects(passLMR_muons,5,2.7, MuD0Sigma3 | MuZ05mm | MuIsoGradientLoose);

    //Count the number of signal leptons
    unsigned int N_SignalEle = signalElectrons.size();
    unsigned int N_SignalMuo = signalMuons.size();
    unsigned int N_SignalTau  = signalTaus.size();
    unsigned int N_SignalLep = N_SignalEle + N_SignalMuo;
    float Ht_Lep = sumObjectsPt(signalElectrons,-1,999) + sumObjectsPt(signalMuons,-1,999) + sumObjectsPt(signalTaus,-1,999);
    float Ht_Jet = sumObjectsPt(signalJets,40, 999);
    float Meff = met + Ht_Lep + Ht_Jet;
    bool ZVeto = PassZVeto(signalElectrons,signalMuons);

    if (ZVeto && N_SignalLep >= 4 && N_SignalTau == 0) {
        if (Meff > 600) accept("SROA");
        if (Meff > 1100) accept("SR0B");
    }


    //                else accept("SRs2j4bT");
//            }
//        }
//        // b-Veto signal region
//        else {
//            if (softJets.size() >= 2 && met > 430 && mt > 100 && met / meff > 0.25 && meff > 700) {
//                if (meff < 1100) accept("SRs2j1bV");
//                else if (meff < 1500) accept("SRs2j2bV");
//                else if (meff < 1900) accept("SRs2j3bV");
//                else accept("SRs2j4bV");
//            }
//        }
//        // b-inclusive discovery SR
//        if (softJets.size() >= 2 && met > 430 && mt > 100 && met / meff > 0.25 && meff > 700) {
//            if (meff > 1100) accept("SRs2j");
//        }
//    }

    return;
}
