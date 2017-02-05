#ifndef ANALYSISRUNNER_H
#define ANALYSISRUNNER_H

#include <vector>
#include <string>

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/TruthSmear.h"


class AnalysisRunner
{
 public:
 AnalysisRunner(std::vector<AnalysisClass*>& analysisList) : _analysisList(analysisList) {};
  void init() { for(const auto& analysis : _analysisList) analysis->Init(); };
  void final() { for(const auto& analysis : _analysisList) analysis->Final(); };
  void SetSmearing(TruthSmear *smear) { _smear=smear; };
  void processEvent(TruthEvent *event,double weight,int eventNumber) { 
    if (_smear) event = _smear->smearEvent(event);
    event->sortObjects();
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
};

class Reader
{
 public:
  Reader(std::vector<AnalysisClass*>& analysisList) {_analysisRunner=new AnalysisRunner(analysisList); }
  virtual ~Reader() {};
  virtual void SetSmearing(TruthSmear *smear) { _analysisRunner->SetSmearing(smear); };
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
