#include <iostream>
#include <sstream>
#include <algorithm>

#include <TFile.h>
#include <TChain.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "SimpleAnalysis/SlimReader.h"
#include "SimpleAnalysis/OldSlimReader.h"
#include "SimpleAnalysis/xAODTruthReader.h"
#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/OutputHandler.h"
#include "SimpleAnalysis/AnalysisRunner.h"

static void splitCommaString(const std::string& names,std::vector<std::string>& result) {
  std::stringstream ss(names);
  while( ss.good() ) {
    std::string substr;
    getline( ss, substr, ',' );
    result.push_back( substr );
  }
}

int main(int argc, char **argv) {
  if ( ! xAOD::Init().isSuccess() ) {
    throw std::runtime_error("Cannot initialise xAOD access !");
  }

  po::options_description desc("Run one or more truth-level analyses");
  std::string outputName;
  desc.add_options()
    ("help,h", "print usage and exit")
    ("output,o",    po::value<std::string>(&outputName), "Output name - if not supplied use analyses names")
    ("analyses,a",  po::value<std::string>(), "Comma-separated list of analyses to run")
    ("listanalyses,l", "List available analyses and exit")
    ("input-files", po::value< vector<std::string> >(), "Comma-separated list of input files")
    ("ntuple,n", "Fill ntuple")
    ;
  po::positional_options_description p;
  p.add("input-files", -1);
  
  po::variables_map vm;
  po::store( po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
  
  if (vm.count("listanalyses")) {
    std::cout<<"Select analyses among:"<<std::endl;
    for(auto *analysis : *getAnalysisList()) {
      std::cout<<" "<<analysis->name()<<std::endl;
    }
    return 0;
  }
  
  if (vm.count("help") || (vm.count("input-files")==0)) {
    std::cout << desc << std::endl;
    return 1;
  }
  bool doNtuple = vm.count("ntuple");
  bool mergedOutput = vm.count("output");

  std::vector<std::string> inputFileNames;
  std::vector<std::string> analysisNames;
  
  if (vm.count("analyses")) splitCommaString(vm["analyses"].as<std::string>(),analysisNames);
  for(const auto& fileNames: vm["input-files"].as<vector<std::string> >()) {
    splitCommaString(fileNames,inputFileNames);
  }
  
  std::cout<<"Analyses to run: ";
  for(const auto& name: analysisNames) std::cout<<name<<" ";
  if (analysisNames.size()==0)  std::cout<<"all";
  std::cout<<std::endl;
  std::cout<<"Files to analyze: ";
  for(const auto& fileName: inputFileNames) std::cout<<fileName<<" ";
  std::cout<<std::endl;
  if (mergedOutput) 
    std::cout<<"Output merged into: "<<outputName<<".[txt|root]"<<std::endl;
  else
    std::cout<<"Output split per analysis"<<std::endl;

  TFile *oRoot = 0;
  if (mergedOutput) oRoot=new TFile((outputName+".root").c_str(),"RECREATE");

  std::vector<AnalysisClass*> analysisList;
  bool selectAnalysis=analysisNames.size()>0;
  for(auto *analysis : *getAnalysisList()) {
    if (selectAnalysis) {
      const auto namePtr=std::find(analysisNames.begin(),analysisNames.end(),analysis->name());
      if (namePtr==analysisNames.end()) continue;
      analysisNames.erase(namePtr);
    }
    if (!mergedOutput) {
      oRoot = new TFile((analysis->name()+".root").c_str(),"RECREATE");
    }
    OutputHandler* output=new OutputHandler(oRoot,doNtuple);
    if (outputName.size()) output->title(analysis->name());
    analysis->setOutput(output);
    analysisList.push_back(analysis);
  }
  if (analysisNames.size()!=0) {
    std::cerr<<"Unknown analysis requested: ";
    for(const auto& name: analysisNames) std::cerr<<name<<" ";
    std::cerr<<std::endl;
    return 1;
  }

  TFile *fh=TFile::Open(inputFileNames[0].c_str());
  if (fh==0) {
    std::cerr<<"Failed to open the first file: "<<inputFileNames[0]<<std::endl;
    return 2;
  }
  Reader *reader=0;
  if (fh->FindKey("truth")) {
    std::cout<<"Reading (old-style) slimmed input"<<std::endl;
    reader=new OldSlimReader(analysisList);
  } else if (fh->FindKey("ntuple")) {
    std::cout<<"Reading slimmed input"<<std::endl;
    reader=new SlimReader(analysisList);
  } else if (fh->FindKey("CollectionTree")) {
    std::cout<<"Reading xAOD input"<<std::endl;
    reader=new xAODTruthReader(analysisList);
  } else {
    std::cerr<<"Unknown input format in: "<<inputFileNames[0]<<std::endl;
    return 2;
  }
  reader->processFiles(inputFileNames);
  delete reader;

  std::ofstream oFile;
  if (mergedOutput) oFile.open(outputName+".txt");
  bool first=true;
  for(const auto& analysis : analysisList) {
    if (!mergedOutput) {
      oFile.close();
      oFile.open(analysis->name()+".txt");
    }
    analysis->getOutput()->saveRegions(oFile,first);
    if (mergedOutput) first=false;
  }
  
  return 0;
}
