#ifndef TRUTHEVENT_H
#define TRUTHEVENT_H

#include "SimpleAnalysis/AnalysisClass.h"
#include <cmath>

class TruthEvent : public AnalysisEvent
{
 public:
 TruthEvent(float sumet, double Ex, double Ey) : _sumet(sumet), _met(Ex,Ey,0,sqrt(Ex*Ex+Ey*Ey),0,0,MET,-1,0) {};
  virtual AnalysisObjects getElectrons(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _baseElectrons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getMuons(float ptCut,float etaCut,int isolation)  { return AnalysisClass::filterObjects( _baseMuons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getTaus(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _baseTaus, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getPhotons(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _basePhotons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getJets(float ptCut,float etaCut,int btag) { return AnalysisClass::filterObjects( _baseJets, ptCut, etaCut, btag); };
  virtual AnalysisObjects getFatJets(float ptCut,float etaCut,int btag) { return AnalysisClass::filterObjects( _baseFatJets, ptCut, etaCut, btag); };
  virtual AnalysisObject  getMET() { return _met; };
  virtual float  getSumET() { return _sumet; };

  virtual float  getGenMET() { return _genmet; };
  virtual float  getGenHT() { return _genht; };

  virtual int    getMCNumber() { return _mcChannel; }; //Temporary until better solution found
  virtual int    getSUSYChannel() { return _susyChannel; }; //Temporary until better solution found
  virtual std::vector<float>  getMCWeights() { return _mcWeights; }; //Temporary until better solution found

  virtual void   setChannelInfo(int mcChannel, int susyChannel) {  _mcChannel=mcChannel; _susyChannel=susyChannel; };

  void addElectron(double Px, double Py, double Pz, double E, int charge, int iso, int idx) {
    _baseElectrons.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,ELECTRON,-1,idx));
  };
  void addElectron(TLorentzVector tlv, int charge,int iso, int idx) {
    _baseElectrons.push_back(AnalysisObject(tlv,charge,iso,ELECTRON,-1,idx));
  };
  void addMuon(double Px, double Py, double Pz, double E, int charge,int iso, int idx) {
    _baseMuons.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,MUON,-1,idx));
  };
  void addMuon(TLorentzVector tlv, int charge, int iso, int idx) {
    _baseMuons.push_back(AnalysisObject(tlv,charge,iso,MUON,-1,idx));
  };
  void addTau(double Px, double Py, double Pz, double E, int charge,int iso, int idx) {
    _baseTaus.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,TAU,-1,idx));
  };
  void addTau(TLorentzVector tlv, int charge, int iso, int idx) {
    _baseTaus.push_back(AnalysisObject(tlv,charge,iso,TAU,-1,idx));
  };
  void addPhoton(double Px, double Py, double Pz, double E,int iso, int idx) {
    _basePhotons.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,PHOTON,-1,idx));
  };
  void addPhoton(TLorentzVector tlv, int iso, int idx) {
    _basePhotons.push_back(AnalysisObject(tlv,0,iso,PHOTON,-1,idx));
  };
  void addJet(double Px, double Py, double Pz, double E, int iso, int truth_id, int idx) {
    _baseJets.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,JET,truth_id,idx));
  };
  void addJet(TLorentzVector tlv, int iso, int truth_id, int idx) {
    _baseJets.push_back(AnalysisObject(tlv,0,iso,JET,truth_id,idx));
  };
  void addFatJet(double Px, double Py, double Pz, double E, int iso, int idx) {
    _baseFatJets.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,FATJET,-1,idx));
  };
  void addFatJet(TLorentzVector tlv, int iso, int idx) {
    _baseFatJets.push_back(AnalysisObject(tlv,0,iso,FATJET,-1,idx));
  };

  virtual void setGenMET(float genMET=0.){
    _genmet = genMET;
  };

  virtual void setGenHT(float genHT=0.){
    _genht = genHT;
  };

  virtual void setMCWeights(std::vector<float> ws){
    _mcWeights = ws;
  }

 private:
  AnalysisObjects _baseElectrons;
  AnalysisObjects _baseMuons;
  AnalysisObjects _baseTaus;
  AnalysisObjects _basePhotons;
  AnalysisObjects _baseJets;
  AnalysisObjects _baseFatJets;
  float           _sumet;
  AnalysisObject  _met;
  float           _genmet;
  float           _genht;

  int             _mcChannel;
  int             _susyChannel;

  std::vector<float> _mcWeights;

 public:
  void sortObjects() {
    AnalysisClass::sortObjectsByPt(_baseElectrons);
    AnalysisClass::sortObjectsByPt(_baseMuons);
    AnalysisClass::sortObjectsByPt(_baseTaus);
    AnalysisClass::sortObjectsByPt(_basePhotons);
    AnalysisClass::sortObjectsByPt(_baseJets);
    AnalysisClass::sortObjectsByPt(_baseFatJets);
  }

};

#endif
