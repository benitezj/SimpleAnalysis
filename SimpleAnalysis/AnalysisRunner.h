#ifndef ANALYSISRUNNER_H
#define ANALYSISRUNNER_H

#include <vector>
#include <string>

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/TruthSmear.h"
#include "SimpleAnalysis/PDFReweight.h"


class AnalysisRunner
{
 public:
 AnalysisRunner(std::vector<AnalysisClass*>& analysisList) : _analysisList(analysisList), _reweighter(0) {};
  void init() { for(const auto& analysis : _analysisList) analysis->Init(); };
  void final() { for(const auto& analysis : _analysisList) analysis->Final(); };
  void SetSmearing(TruthSmear *smear) { _smear=smear; };
  void SetReweighting(PDFReweighter *reweighter) { _reweighter=reweighter; };
  void SetMCWeightIndex(int mcwidx) { _mcwindex=mcwidx; };

  int getMCWeightIndex(){ return _mcwindex; };

  void processEvent(TruthEvent *event,int eventNumber) { 
    if (_smear) event = _smear->smearEvent(event);
    double weight = 1.;
    if (_mcwindex >= int(event->getMCWeights().size())){
      throw std::runtime_error("The specified MC weight index is out of range! ");
    }
    if (_mcwindex>=0) weight = event->getMCWeights()[_mcwindex];
    event->sortObjects();
    if (_reweighter) weight *= _reweighter->reweightEvent(event);
    for(const auto& analysis : _analysisList) {
      analysis->getOutput()->setEventWeight(weight);
      analysis->getOutput()->ntupVar("Event", eventNumber);
      analysis->ProcessEvent(event);
      analysis->getOutput()->ntupFill();
    }
    if (_smear) delete event;
  };

 private:
  std::vector<AnalysisClass*>& _analysisList;
  TruthSmear *_smear;
  PDFReweighter *_reweighter;
  int _mcwindex;
};

class Reader
{
 public:
  Reader(std::vector<AnalysisClass*>& analysisList) {_analysisRunner=new AnalysisRunner(analysisList); }
  virtual ~Reader() {};
  virtual void SetSmearing(TruthSmear *smear) { _analysisRunner->SetSmearing(smear); };
  virtual void SetReweighting(PDFReweighter *reweighter) { _analysisRunner->SetReweighting(reweighter); };
  virtual void SetMCWeightIndex(int mcwidx) { _analysisRunner->SetMCWeightIndex(mcwidx); };

  virtual int getMCWeightIndex(){ return _analysisRunner->getMCWeightIndex(); };

  virtual void processFiles(const std::vector<std::string>& inputNames) {
    _analysisRunner->init();
    processFilesInternal(inputNames);
    _analysisRunner->final();
  }
 protected:
  virtual void processFilesInternal(const std::vector<std::string>& inputNames) = 0;
  AnalysisRunner *_analysisRunner;
 private:
		    
};


#endif
