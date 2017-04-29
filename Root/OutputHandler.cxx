#include <stdlib.h>
#include <math.h>

#include "SimpleAnalysis/OutputHandler.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

int OutputHandler::addEntry(const std::string &name) {
  int num=eventCounts.size();
  if (label2idx.find(name)!=label2idx.end()) {
    std::cerr<<"Duplicate signal region label: "<<name<<std::endl;
    exit(1);
  }
  if (name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890_") != std::string::npos) {
    std::cerr << "Illegal signal region name: "<<name<<std::endl;
    std::cerr << "Signal region names should only have alphanumeric characters and '_'\n"<<std::endl;
    exit(1);
  }
  label2idx[name]=num;
  idx2label[num]=name;
  eventCounts.push_back(0);
  weightedSums.push_back(0);
  weightedSquaredSums.push_back(0);

  return num;
}

void OutputHandler::addEntries(const std::vector<std::string> &labels) {
  for(const auto& label: labels ) addEntry(label);
}

void OutputHandler::addHistogram(const std::string &label,int bins,float min,float max) {
  TH1 *hist=new TH1D((_title+label).c_str(),label.c_str(),bins,min,max);
  hist->SetDirectory(_oFile);
  histograms[label]=hist; 
}

void OutputHandler::addHistogram(const std::string &label,int bins,float *edges) {
  TH1 *hist=new TH1D((_title+label).c_str(),label.c_str(),bins,edges);
  hist->SetDirectory(_oFile);
  histograms[label]=hist;
}

void OutputHandler::addHistogram(const std::string &label,std::vector<float> &edges) {
  float* my_edges = new float[ edges.size() ];
  for (unsigned int e=0;e<edges.size();++e) my_edges[e] = edges[e];
  TH1 *hist=new TH1D((_title+label).c_str(),label.c_str(),edges.size(),my_edges);
  hist->SetDirectory(_oFile);
  histograms[label]=hist;
}

void OutputHandler::addHistogram(const std::string &label,
				 int binsX,float minX,float maxX,
				 int binsY,float minY,float maxY) {
  TH1 *hist=new TH2D((_title+label).c_str(),label.c_str(),binsX,minX,maxX,binsY,minY,maxY);
  hist->SetDirectory(_oFile);
  histograms[label]=hist; 
}

void OutputHandler::pass(int num,double weight) {
  weight*=_eventWeight;
  if (idx2label.find(num)==idx2label.end()) {
    std::cerr<<"Unknown signal region number: "<<num<<std::endl;
    exit(1);
  }
  if (num==0) {
    ntupVar("eventWeight",weight);
  } else {
    ntupVar(idx2label[num],weight);
  }
  eventCounts[num]++;
  weightedSums[num]+=weight;
  weightedSquaredSums[num]+=weight*weight;
}

void OutputHandler::pass(const std::string & name,double weight) {
  if (label2idx.find(name)==label2idx.end()) {
    std::cerr<<"Unknown signal region label: "<<name<<std::endl;
    exit(1);
  }
  pass(label2idx[name],weight);
}

void OutputHandler::fillHistogram(const std::string &label,
				 double x) {
  if (histograms.find(label)==histograms.end()) {
    std::cerr<<"Unknown histogram label: "<<label<<std::endl;
    exit(1);
  }
  histograms[label]->Fill(x,_eventWeight);
}

void OutputHandler::fillHistogram(const std::string &label,
				  double x,double y) {
  if (histograms.find(label)==histograms.end()) {
    std::cerr<<"Unknown histogram label: "<<label<<std::endl;
    exit(1);
  }
  ((TH2* ) histograms[label])->Fill(x,y,_eventWeight);
}

void OutputHandler::setEventWeight(double weight) {
  for(auto & ntup : ntupInts) *ntup=0;
  for(auto & ntup : ntupFloats) *ntup=0;
  for(auto & ntupVec: ntupVectorFloats) ntupVec->clear();
  for(auto & ntupVec: ntupVectorInts) ntupVec->clear();

  _eventWeight = weight;
  pass(0); //count all events
}

void OutputHandler::createNtuple() {
  _ntuple=new TTree((_title+"ntuple").c_str(),"Simple Analysis ntuple");
  _ntuple->SetDirectory(_oFile);   
}

void OutputHandler::ntupVar(const std::string &label,int value) {
  if (!_doNtuple) return;
  if (!_ntuple) createNtuple();
  if (ntupInt.find(label)==ntupInt.end()) {
    ntupInt[label]=0;
    TBranch *branch=_ntuple->Branch(label.c_str(),&(ntupInt[label]),(label+"/I").c_str());
    ntupInts.push_back(&(ntupInt[label]));
    for(int ii=0;ii<_ntuple->GetEntries();++ii) branch->Fill(); //backfill
  }
  ntupInt[label]=value;
}

void OutputHandler::ntupVar(const std::string &label,float value) {
  if (!_doNtuple) return;
  if (!_ntuple) createNtuple();
  if (ntupFloat.find(label)==ntupFloat.end()) {
    ntupFloat[label]=0;
    TBranch *branch=_ntuple->Branch(label.c_str(),&(ntupFloat[label]),(label+"/F").c_str());
    ntupFloats.push_back(&(ntupFloat[label]));
    for(int ii=0;ii<_ntuple->GetEntries();++ii) branch->Fill(); //backfill
  }
  ntupFloat[label]=value;
}

void OutputHandler::ntupVar(const std::string &label,std::vector<float> &values) {
  if (!_doNtuple) return;
  if (!_ntuple) createNtuple();
  if (ntupVectorFloat.find(label)==ntupVectorFloat.end()) {
    ntupVectorFloat[label]=new std::vector<float>;
    TBranch *branch=_ntuple->Branch(label.c_str(),ntupVectorFloat[label]);
    ntupVectorFloats.push_back(ntupVectorFloat[label]);
    for(int ii=0;ii<_ntuple->GetEntries();++ii) branch->Fill(); //backfill
  }
  std::vector<float> *dest=ntupVectorFloat[label];
  for(float value : values)
    dest->push_back(value);
}

void OutputHandler::ntupVar(const std::string &label,std::vector<int> &values) {
  if (!_doNtuple) return;
  if (!_ntuple) createNtuple();
  if (ntupVectorInt.find(label)==ntupVectorInt.end()) {
    ntupVectorInt[label]=new std::vector<int>;
    TBranch *branch=_ntuple->Branch(label.c_str(),ntupVectorInt[label]);
    ntupVectorInts.push_back(ntupVectorInt[label]);
    for(int ii=0;ii<_ntuple->GetEntries();++ii) branch->Fill(); //backfill
  }
  std::vector<int> *dest=ntupVectorInt[label];
  for(int value : values)
    dest->push_back(value);
}


void OutputHandler::saveRegions(std::ostream& filehandle,bool header) {
  //  csv_ofile << "SR,Acceptance,Sum,WeightedSum,WeightedSquareSum"<<std::endl;
  if (header) {
    filehandle<<"SR,events,acceptance,err"<<std::endl;
    // Not so obvious, but we only want to save to ROOT file once
    _oFile->Write();
    _oFile->Close();
  }
  for(const auto& idx: idx2label) {
    filehandle<< _title+idx.second<<","<<eventCounts[idx.first];
    if (idx.first==0) {
      filehandle<<","<<weightedSums[0]<<","<<weightedSquaredSums[0];
    } else {
      double acc=weightedSums[idx.first]/weightedSums[0];
      double err=sqrt(weightedSquaredSums[idx.first])/weightedSums[0];
      filehandle<<","<<acc<<","<<err;
    }
    filehandle<<std::endl;
  }
}

