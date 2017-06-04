#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/mt2_bisect.h"
#include "SimpleAnalysis/MT2.h"

#include "TMatrixDSym.h"
#include "TVectorD.h"

#include <algorithm>

#define FASTJET
#ifdef FASTJET
// Jets
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#endif

const float Z_Mass = 91.2;

AnalysisObject operator+(const AnalysisObject& lhs, const AnalysisObject& rhs) {
    const TLorentzVector &tlhs = lhs;
    return AnalysisObject(tlhs + rhs, lhs.charge() + rhs.charge(), 0, COMBINED, 0);
}

AnalysisObjects operator+(const AnalysisObjects& lhs, const AnalysisObjects& rhs) {
    AnalysisObjects combined;
    for (const auto &cand : lhs)
        combined.push_back(cand);
    for (const auto &cand : rhs)
        combined.push_back(cand);
    AnalysisClass::sortObjectsByPt(combined);
    return combined;
}

std::vector<AnalysisClass*> *getAnalysisList() {
    static std::vector<AnalysisClass*> *list = new std::vector<AnalysisClass*>;
    return list;
}

void AnalysisClass::ntupVar(const std::string &label, AnalysisObject &object, bool saveMass, bool saveType) {
    _output->ntupVar(label + "_pt", object.Pt());
    _output->ntupVar(label + "_eta", object.Eta());
    _output->ntupVar(label + "_phi", object.Phi());
    _output->ntupVar(label + "_charge", object.charge());
    if (saveMass) _output->ntupVar(label + "_m", object.M());
    if (saveType) _output->ntupVar(label + "_type", object.type());
    _output->ntupVar(label + "_id", object.id());
}

void AnalysisClass::ntupVar(const std::string &label, AnalysisObjects &objects, bool saveMass, bool saveType) {
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<int> charge;
    std::vector<int> type;
    std::vector<int> id;
    std::vector<float> mass;
    for (const auto& object : objects) {
        pt.push_back(object.Pt());
        eta.push_back(object.Eta());
        phi.push_back(object.Phi());
        charge.push_back(object.charge());
        mass.push_back(object.M());
        type.push_back(object.type());
        id.push_back(object.id());
    }
    _output->ntupVar(label + "_pt", pt);
    _output->ntupVar(label + "_eta", eta);
    _output->ntupVar(label + "_phi", phi);
    _output->ntupVar(label + "_charge", charge);
    if (saveMass) _output->ntupVar(label + "_m", mass);
    if (saveType) _output->ntupVar(label + "_type", type);
    _output->ntupVar(label + "_id", id);

}

AnalysisObjects AnalysisClass::overlapRemoval(const AnalysisObjects &cands, const AnalysisObjects &others, std::function<float(const AnalysisObject&, const AnalysisObject&)> radiusFunc, int passId) {
    AnalysisObjects reducedList;

    for (const auto& cand : cands) {
        bool overlap = false;
        for (const auto& other : others) {
            if (cand.DeltaR(other) < radiusFunc(cand, other) && cand != other && cand.pass(passId)) {
                overlap = true;
                break;
            }
        }
        if (!overlap) reducedList.push_back(cand);
    }
    return reducedList;
}
AnalysisObjects AnalysisClass::overlapRemoval(const AnalysisObjects &cands, const AnalysisObjects &others, float deltaR, int passId) {
    return overlapRemoval(cands, others, [deltaR] (const AnalysisObject& , const AnalysisObject&) {return deltaR;}, passId);
}

AnalysisObjects AnalysisClass::lowMassRemoval(const AnalysisObjects & cand, std::function<bool(const AnalysisObject&, const AnalysisObject&)> RemoveSwitch, float MinMass, float MaxMass, int type) {
    AnalysisObjects reducedList;
    unsigned int NObj = cand.size();
    AnalysisObjects::const_iterator begin_cand = cand.begin();
    std::vector<bool> Good(NObj, true);
    for (unsigned int obj = 0; obj != NObj; ++obj) {
        //Check if the object is already rejected
        if (!Good.at(obj)) continue;
        const AnalysisObject& Object = *(begin_cand + obj);
        for (unsigned int obj1 = 0; obj1 != obj; ++obj) {
            if (!Good.at(obj1)) continue;
            //Remove the pair
            const AnalysisObject& Object1 = *(begin_cand + obj1);
            float InvMass = (Object + Object1).M();
            if (MinMass < InvMass && InvMass < MaxMass && RemoveSwitch(Object, Object1)) {
                Good.at(obj) = false;
                Good.at(obj1) = false;
                break;
            }
        }
    }
    for (unsigned int o = 0; o < NObj; ++o) {
        if ((type == -1  || cand.at(o).type() == type) && Good.at(o)) reducedList.push_back(cand.at(o));
    }
    return reducedList;
}
AnalysisObjects AnalysisClass::filterObjects(const AnalysisObjects& cands, float ptCut, float etaCut, int id) {
    AnalysisObjects reducedList;
    for (const auto& cand : cands) {
        if ((cand.Pt() >= ptCut) && (fabs(cand.Eta()) < etaCut) && (cand.pass(id))) reducedList.push_back(cand);
    }
    return reducedList;
}

AnalysisObjects AnalysisClass::filterCrack(const AnalysisObjects& cands, float minEta, float maxEta) {
    AnalysisObjects reducedList;
    for (const auto& cand : cands) {
        if (fabs(cand.Eta()) < minEta || fabs(cand.Eta()) > maxEta) reducedList.push_back(cand);
    }
    return reducedList;
}

int AnalysisClass::countObjects(const AnalysisObjects& cands, float ptCut, float etaCut, int id) {
    int count = 0;
    for (const auto& cand : cands) {
        if (cand.Pt() >= ptCut && fabs(cand.Eta()) < etaCut && cand.pass(id)) count++;
    }
    return count;
}

float AnalysisClass::sumObjectsPt(const AnalysisObjects& cands, unsigned int maxNum, float ptCut) {
    float sum = 0;
    for (int ii = 0; ii < std::min((int) maxNum, (int) cands.size()); ii++)
        if (cands[ii].Pt() > ptCut) sum += cands[ii].Pt();
    return sum;
}

float AnalysisClass::sumObjectsM(const AnalysisObjects& cands, unsigned int maxNum, float mCut) {
    float sum = 0;
    for (int ii = 0; ii < std::min((int) maxNum, (int) cands.size()); ii++)
        if (cands[ii].M() > mCut) sum += cands[ii].M();
    return sum;
}

struct pt_sort {
        bool operator()(const TLorentzVector& v1, const TLorentzVector& v2) const {
            return v1.Pt() > v2.Pt();
        }
};
void AnalysisClass::sortObjectsByPt(AnalysisObjects& cands) {
    std::sort(cands.begin(), cands.end(), pt_sort());
}

float AnalysisClass::calcMCT(const AnalysisObject& o1, const AnalysisObject& o2) {
    float mCT = pow(o1.Et() + o2.Et(), 2) - pow(o1.Px() - o2.Px(), 2) - pow(o1.Py() - o2.Py(), 2);
    mCT = (mCT >= 0.) ? sqrt(mCT) : sqrt(-mCT);
    return mCT;
}

float AnalysisClass::calcMT(const AnalysisObject& lepton, const AnalysisObject& met) {
    float mT = 2 * lepton.Pt() * met.Et() * (1 - TMath::Cos(lepton.Phi() - met.Phi()));
    mT = (mT >= 0.) ? sqrt(mT) : sqrt(-mT);
    return mT;
}

float AnalysisClass::calcMTmin(const AnalysisObjects& cands, const AnalysisObject& met, int maxNum) {
    float mtmin = 1e9;

    for (int ii = 0; ii < std::min((int) maxNum, (int) cands.size()); ii++)
        mtmin = std::min(mtmin, calcMT(cands[ii], met));
    return mtmin;
}

float AnalysisClass::calcMT2(const AnalysisObject &o1, const AnalysisObject &o2, const AnalysisObject &met) {
    double patr[3] = { o1.M(), o1.Px(), o1.Py() };
    double pbtr[3] = { o2.M(), o2.Px(), o2.Py() };
    double pmisstr[3] = { 0, met.Px(), met.Py() };

    mt2_bisect::mt2 mt2calculator;
    mt2calculator.set_momenta(patr, pbtr, pmisstr);
    mt2calculator.set_mn(0);
    return mt2calculator.get_mt2();
}

float AnalysisClass::calcAMT2(const AnalysisObject &o1, const AnalysisObject &o2, const AnalysisObject &met, float m1, float m2) {
    return asymm_mt2_lester_bisect::get_mT2(o1.M(), o1.Px(), o1.Py(), o2.M(), o2.Px(), o2.Py(), met.Px(), met.Py(), m1, m2);
}

float AnalysisClass::minDphi(const AnalysisObject &met, const AnalysisObjects& cands, unsigned int maxNum, float ptCut) {
    float dphi_min = 999.;
    for (int ii = 0; ii < std::min((int) maxNum, (int) cands.size()); ii++) {
        float dphi = fabs(met.DeltaPhi(cands[ii]));
        if (dphi < dphi_min && cands[ii].Pt() > ptCut) dphi_min = dphi;
    }
    return dphi_min;
}

float AnalysisClass::minDR(const AnalysisObjects& cands, unsigned int maxNum, float ptCut) {
    float dr_min = 999.;
    for (int ii = 0; ii < std::min((int) maxNum, (int) cands.size()); ii++) {
        for (int ij = ii + 1; ij < std::min((int) maxNum, (int) cands.size()); ij++) {
            float dr = fabs(cands[ii].DeltaR(cands[ij]));
            if (dr < dr_min && cands[ii].Pt() > ptCut) dr_min = dr;
        }
    }
    return dr_min;
}

AnalysisObjects AnalysisClass::reclusterJets(const AnalysisObjects &jets, float radius, float ptmin, float rclus, float ptfrac) {
    std::vector<fastjet::PseudoJet> JetVector;
    for (const auto& jet : jets) {
        fastjet::PseudoJet Pjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
        JetVector.push_back(Pjet);
    }
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius, recomb_scheme, strategy);

    std::vector<fastjet::PseudoJet> fat_jets;
    fastjet::ClusterSequence clust_seq(JetVector, jet_def);
    if (!JetVector.empty()) {
        fat_jets = clust_seq.inclusive_jets(ptmin);
    }
    AnalysisObjects fatJets;
    fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, rclus), fastjet::SelectorPtFractionMin(ptfrac));
    for (const auto& fat_jet : fat_jets) {
        if (ptfrac >= 0) {
            fastjet::PseudoJet trimmed_jet = trimmer(fat_jet);
            if (trimmed_jet.pt() > ptmin) fatJets.push_back(AnalysisObject(trimmed_jet.px(), trimmed_jet.py(), trimmed_jet.pz(), trimmed_jet.E(), 0, 0, COMBINED, 0));
        } else {
            fatJets.push_back(AnalysisObject(fat_jet.px(), fat_jet.py(), fat_jet.pz(), fat_jet.E(), 0, 0, COMBINED, 0));
        }
    }
    sortObjectsByPt(fatJets);
    return fatJets;
}

float AnalysisClass::aplanarity(const AnalysisObjects& jets) {
    TMatrixDSym momTensor(3);
    TVectorD eigenValues(3);
    momTensor.Zero();
    if (jets.size() < 2) return 0;
    double norm = 0;
    for (const auto& jet : jets) {
        momTensor(0, 0) += jet.Px() * jet.Px();
        momTensor(0, 1) += jet.Px() * jet.Py();
        momTensor(0, 2) += jet.Px() * jet.Pz();
        momTensor(1, 0) += jet.Py() * jet.Px();
        momTensor(1, 1) += jet.Py() * jet.Py();
        momTensor(1, 2) += jet.Py() * jet.Pz();
        momTensor(2, 0) += jet.Pz() * jet.Px();
        momTensor(2, 1) += jet.Pz() * jet.Py();
        momTensor(2, 2) += jet.Pz() * jet.Pz();
        norm += jet.Vect().Mag2();
    }
    momTensor *= 1. / norm;
    momTensor.EigenVectors(eigenValues);
    return 1.5 * eigenValues(2);
}
bool AnalysisClass::IsSFOS(const AnalysisObject& L, const AnalysisObject& L1) {
    return L1.type() == L.type() && L1.charge() != L.charge();
}
std::pair<float, float> AnalysisClass::DiZSelection(const AnalysisObjects& electrons, const AnalysisObjects& muons) {
    std::pair<float, float> ZMasses(-1, -1);
    float BestDZ = 1.e25;
    auto Leptons = electrons + muons;
    AnalysisObjects::const_iterator Lep_begin = Leptons.begin();
    AnalysisObjects::const_iterator Lep_end = Leptons.end();
    unsigned int NumLep = Leptons.size();

    for (AnalysisObjects::const_iterator L = Lep_begin + 1; L != Lep_end; ++L) {
        for (AnalysisObjects::const_iterator L1 = Lep_begin; L1 != L; ++L1) {
            //Loop over the first two leptons
            const AnalysisObject& FirstLep = (*L);
            const AnalysisObject& SecondLep = (*L1);
            //No SFOS pair
            if (!IsSFOS(FirstLep, SecondLep)) continue;
            float M1 = (FirstLep + SecondLep).M();
            float dM1 = fabs(M1 - Z_Mass);
            //Find the second pair in the row
            for (AnalysisObjects::const_iterator L2 = Lep_begin + 1; L2 != Lep_end; ++L2) {
                //Do not use the same object twice
                if (L2 == L1 || L2 == L) continue;
                for (AnalysisObjects::const_iterator L3 = Lep_begin; L3 != L2; ++L3) {
                    if (L3 == L1 || L3 == L) continue;
                    const AnalysisObject& ThirdLep = (*L2);
                    const AnalysisObject& FourthLep = (*L3);
                    if (!IsSFOS(ThirdLep, FourthLep)) continue;
                    float M2 = (ThirdLep + FourthLep).M();
                    float dM2 = fabs(M2 - Z_Mass);
                    float dZTest = dM1 + dM2;
                    //The sum of the distances to M_Z is smaller than anything known
                    if (dZTest < BestDZ) {
                        //Order them by distance closest -> leading
                        if (dM1 < dM2) ZMasses = std::pair<float, float>(M1, M2);
                        else ZMasses = std::pair<float, float>(M2, M1);
                        BestDZ = dZTest;
                    }
                } //End of the second pair finder

                //No way to get two Zs of the event. Find the best one nervertheless
                if (NumLep < 4) {
                    if (dM1 < BestDZ) {
                        BestDZ = dM1;
                        ZMasses.first = M1;
                    }
                }
            }
        }
    }
    return ZMasses;
}
bool AnalysisClass::PassZVeto(const AnalysisObjects& electrons, const AnalysisObjects &muons, float Window) {
    auto Leptons = electrons + muons;
    AnalysisObjects::const_iterator Lep_begin = Leptons.begin();
    AnalysisObjects::const_iterator Lep_end = Leptons.end();
    for (AnalysisObjects::const_iterator L = Lep_begin + 1; L != Lep_end; ++L) {
        for (AnalysisObjects::const_iterator L1 = Lep_begin; L1 != L; ++L1) {
            //Loop over the first two leptons
            const AnalysisObject& FirstLep = (*L);
            const AnalysisObject& SecondLep = (*L1);
            AnalysisObject DiLep = (FirstLep + SecondLep);
            if (IsSFOS(FirstLep, SecondLep) && fabs(DiLep.M() - Z_Mass) < Window) return false;
            for (AnalysisObjects::const_iterator L2 = Lep_begin; L2 != L1; ++L2) {
                const AnalysisObject& ThirdLep = (*L2);
                AnalysisObject TriLep = DiLep + ThirdLep;
                //Radiated photons might be converted into leptons
                bool TriLepSFOS = IsSFOS(FirstLep, SecondLep) || IsSFOS(FirstLep, ThirdLep) || IsSFOS(SecondLep, ThirdLep);
                if (TriLepSFOS && fabs(ThirdLep.M() - Z_Mass) < Window) return false;
                //Check forZ ->llll
                for (AnalysisObjects::const_iterator L3 = Lep_begin; L3 != L2; ++L3) {
                    const AnalysisObject& FourthLep = (*L3);
                    bool FourLepSFOS = (IsSFOS(FirstLep, SecondLep) && IsSFOS(ThirdLep, FourthLep)) //eemumu
                    || (IsSFOS(FirstLep, ThirdLep) && IsSFOS(SecondLep, FourthLep)) //e mu e mu
                            || (IsSFOS(FirstLep, FourthLep) && IsSFOS(ThirdLep, SecondLep)); //e mu mu e
                    if (FourLepSFOS && fabs((TriLep + FourthLep).M() - Z_Mass) < Window) return false;
                }
            }
        }
    }
    return true;
}

