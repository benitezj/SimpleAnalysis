#ifndef TRUTHEVENT_H
#define TRUTHEVENT_H

#include "SimpleAnalysis/AnalysisClass.h"
#include <cmath>
#include "xAODTruth/TruthParticleContainer.h"

class TruthEvent : public AnalysisEvent
{
 public:
 TruthEvent(float sumet, double Ex, double Ey) : _sumet(sumet), _met(Ex,Ey,0,sqrt(Ex*Ex+Ey*Ey),0,0,MET,0) {};
  virtual AnalysisObjects getElectrons(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _baseElectrons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getMuons(float ptCut,float etaCut,int isolation)  { return AnalysisClass::filterObjects( _baseMuons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getTaus(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _baseTaus, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getPhotons(float ptCut,float etaCut,int isolation) { return AnalysisClass::filterObjects( _basePhotons, ptCut, etaCut, isolation); };
  virtual AnalysisObjects getJets(float ptCut,float etaCut,int btag) { return AnalysisClass::filterObjects( _baseJets, ptCut, etaCut, btag); };
  virtual AnalysisObjects getFatJets(float ptCut,float etaCut,int btag) { return AnalysisClass::filterObjects( _baseFatJets, ptCut, etaCut, btag); };
  virtual AnalysisObject  getMET() { return _met; };
  virtual const xAOD::TruthParticleContainer * getTruthParticles(){return _truthParticles;}
  void setTruthParticles(const xAOD::TruthParticleContainer * truthparticles){_truthParticles=truthparticles;}

  virtual float  getSumET() { return _sumet; };

  virtual float  getGenMET() { return _genmet; };
  virtual float  getGenHT() { return _genht; };

  virtual int    getMCNumber() { return _mcChannel; }; //Temporary until better solution found
  virtual int    getSUSYChannel() { return _susyChannel; }; //Temporary until better solution found
  virtual std::vector<float>  getMCWeights() { return _mcWeights; }; //Temporary until better solution found

  virtual void   setChannelInfo(int mcChannel, int susyChannel) {  _mcChannel=mcChannel; _susyChannel=susyChannel; };
  virtual int   getPDF_id1()  { return _id1; };
  virtual float getPDF_x1()   { return _x1; };
  virtual float getPDF_pdf1() { return _pdf1; };
  virtual int   getPDF_id2()  { return _id2; };
  virtual float getPDF_x2()   { return _x2; };
  virtual float getPDF_pdf2() { return _pdf2; };
  virtual float getPDF_scale() { return _scale; };
  virtual void setPDFInfo(int id1, float x1, float pdf1,
			  int id2, float x2, float pdf2, float scale) {
    _id1=id1; _x1=x1; _pdf1=pdf1;
    _id2=id2; _x2=x2; _pdf2=pdf2; _scale=scale;
  };

  void addElectron(double Px, double Py, double Pz, double E, int charge, int iso, int idx) {
    _baseElectrons.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,ELECTRON,idx));
  };
  void addElectron(TLorentzVector tlv, int charge,int iso, int idx) {
    _baseElectrons.push_back(AnalysisObject(tlv,charge,iso,ELECTRON,idx));
  };
  void addMuon(double Px, double Py, double Pz, double E, int charge,int iso, int idx) {
    _baseMuons.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,MUON,idx));
  };
  void addMuon(TLorentzVector tlv, int charge, int iso, int idx) {
    _baseMuons.push_back(AnalysisObject(tlv,charge,iso,MUON,idx));
  };
  void addTau(double Px, double Py, double Pz, double E, int charge,int iso, int idx) {
    _baseTaus.push_back(AnalysisObject(Px,Py,Pz,E,charge,iso,TAU,idx));
  };
  void addTau(TLorentzVector tlv, int charge, int iso, int idx) {
    _baseTaus.push_back(AnalysisObject(tlv,charge,iso,TAU,idx));
  };
  void addPhoton(double Px, double Py, double Pz, double E,int iso, int idx) {
    _basePhotons.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,PHOTON,idx));
  };
  void addPhoton(TLorentzVector tlv, int iso, int idx) {
    _basePhotons.push_back(AnalysisObject(tlv,0,iso,PHOTON,idx));
  };
  void addJet(double Px, double Py, double Pz, double E, int iso, int idx) {
    _baseJets.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,JET,idx));
  };
  void addJet(TLorentzVector tlv, int iso, int idx) {
    _baseJets.push_back(AnalysisObject(tlv,0,iso,JET,idx));
  };
  void addFatJet(double Px, double Py, double Pz, double E, int iso, int idx) {
    _baseFatJets.push_back(AnalysisObject(Px,Py,Pz,E,0,iso,FATJET,idx));
  };
  void addFatJet(TLorentzVector tlv, int iso, int idx) {
    _baseFatJets.push_back(AnalysisObject(tlv,0,iso,FATJET,idx));
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
  const xAOD::TruthParticleContainer * _truthParticles=NULL;

  int             _mcChannel;
  int             _susyChannel;
  int             _id1;
  float           _x1;
  float           _pdf1;
  int             _id2;
  float           _x2;
  float           _pdf2;
  float           _scale;

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
