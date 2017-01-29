#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ThreeBjets2012)

void ThreeBjets2012::Init()
{
  addRegions({"SR_0l_4j_A","SR_0l_4j_B","SR_0l_4j_C",
	               "SR_0l_7j_A","SR_0l_7j_B","SR_0l_7j_C",
	               "SR_1l_6j_A","SR_1l_6j_B","SR_1l_6j_C"});
}

void ThreeBjets2012::ProcessEvent(AnalysisEvent *event)
{
    
    auto electrons  = event->getElectrons(20,2.47);
    auto muons      = event->getMuons(10,2.5);
    auto jets       = event->getJets(20.,2.8);
    auto metVec     = event->getMET();
    double met      = metVec.Et();
    
    // Standard SUSY overlap removal
    jets       = overlapRemoval(jets,electrons,0.2);
    electrons  = overlapRemoval(electrons,jets,0.4);
    muons      = overlapRemoval(muons,jets,0.4);

    // get b-jets
    auto bjets = filterObjects(jets,30.,2.5,true);
    
    // Count electrons and muons
    int nb_baseline_muons     = muons.size();
    int nb_baseline_electrons = electrons.size();
    int nb_signal_muons       = countObjects(muons,25.);
    int nb_signal_electrons   = countObjects(electrons,25.);
    
    // count jets
    int nb_jets_30  = countObjects(jets,30.);
    int nb_jets_50  = countObjects(jets,50.);
    int nb_bjets_30 = countObjects(bjets,30.);
    int nb_bjets_50 = countObjects(bjets,50.);

    //baseline selection
    if(nb_jets_30 < 4) return;
    else if(nb_bjets_30 < 3) return;
    else if(jets.at(0).Pt() < 90.) return;
    
    // Define Ht, Meff and MetSig
    float Ht_4j    = sumObjectsPt(jets,4,30);
    float Ht_incl  = sumObjectsPt(jets,999,30);
    float dphi_min = minDphi(metVec,jets,4);
    
    AnalysisObjects leptons = electrons+muons;
    Ht_incl += sumObjectsPt(leptons,999,20);
    
    float meff_4j    = met + Ht_4j;
    float meff_incl  = met + Ht_incl;
    float met_sig_4j = 0;
    if (Ht_4j>0.) met_sig_4j = met/sqrt(Ht_4j);
    
    /// mT for the leading signal lepton
    float mT = 0.;
    if(leptons.size()>0) mT=calcMT(leptons[0],metVec);
    
    // Apply lepton veto
    bool passLeptonVeto = ((nb_baseline_muons + nb_baseline_electrons)==0);
    // Apply leading jet is anti b-tagged
    bool passAntiBtag = (jets.at(0).Pt()>bjets.at(0).Pt());
    
    // 0-lepton signal regions
    if(passLeptonVeto && dphi_min>0.5 && met/meff_4j>0.2) {        
        if (nb_jets_50>=4 && nb_bjets_50>=3 && met>250.  && meff_4j>1300.) accept("SR_0l_4j_A");
        if (nb_jets_50>=4 && nb_bjets_50>=3 && met>350.  && meff_4j>1100.) accept("SR_0l_4j_B");
        if (passAntiBtag && met>400.  && meff_4j>1000. && met_sig_4j>16)   accept("SR_0l_4j_C");

        if (nb_jets_30>=7 && met>200.  && meff_incl>1000.) accept("SR_0l_7j_A");
        if (nb_jets_30>=7 && met>350.  && meff_incl>1000.) accept("SR_0l_7j_B");
        if (nb_jets_30>=7 && met>250.  && meff_incl>1500.) accept("SR_0l_7j_C");
    }
    // 1-lepton signal regions
    if((nb_signal_electrons + nb_signal_muons) >= 1) {
        if (nb_jets_30>=6 && mT> 140. && met>175.  && meff_incl>700.) accept("SR_1l_6j_A");
        if (nb_jets_30>=6 && mT> 140. && met>225.  && meff_incl>800.) accept("SR_1l_6j_B");
        if (nb_jets_30>=6 && mT> 160. && met>275.  && meff_incl>900.) accept("SR_1l_6j_C");
    }    
    return;
}
