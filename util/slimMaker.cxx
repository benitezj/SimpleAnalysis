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
#include "SimpleAnalysis/NtupleMaker.h"

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

  po::options_description desc("Slim xAOD truth to small ntuple");
  std::string outputName="slimNtuple.root";
  double minElecPt=3;
  double minMuonPt=3;
  double minTauPt=15;
  double minPhotonPt=15;
  double minJetPt=15;
  double minFatJetPt=100;
  desc.add_options()
    ("help,h", "print usage and exit")
    ("output,o",    po::value<std::string>(&outputName), "Output name [default: slimNtuple.root]")
    ("minElectronPt",  po::value<double>(&minElecPt), "Minimum electron pt [default: 3 GeV]")
    ("minMuonPt",  po::value<double>(&minMuonPt), "Minimum muon pt [default: 3 GeV]")
    ("minTauPt",  po::value<double>(&minTauPt), "Minimum tau pt [default: 15 GeV]")
    ("minPhotonPt",  po::value<double>(&minPhotonPt), "Minimum photon pt [default: 15 GeV]")
    ("minJetPt",  po::value<double>(&minJetPt), "Minimum jet pt [default: 15 GeV]")
    ("minFatJetPt",  po::value<double>(&minFatJetPt), "Minimum fat jet pt [default: 100 GeV]")
    ("input-files", po::value< vector<std::string> >(), "Comma-separated list of input files")
    ;
  po::positional_options_description p;
  p.add("input-files", -1);
  
  po::variables_map vm;
  po::store( po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
   
  if (vm.count("help") || (vm.count("input-files")==0)) {
    std::cout << desc << std::endl;
    return 1;
  }
  std::vector<std::string> inputFileNames;

  for(const auto& fileNames: vm["input-files"].as<vector<std::string> >()) {
    splitCommaString(fileNames,inputFileNames);
  }
  
  std::cout<<"Files to slim: ";
  for(const auto& fileName: inputFileNames) std::cout<<fileName<<" ";
  std::cout<<std::endl;

  TFile *oRoot = new TFile((outputName).c_str(),"RECREATE");

  AnalysisClass* writer= new NtupleMaker(minElecPt,minMuonPt,minTauPt,minPhotonPt,minJetPt,minFatJetPt);
  OutputHandler* output=new OutputHandler(oRoot,true);
  std::vector<AnalysisClass*> analysisList;
  writer->setOutput(output);
  analysisList.push_back(writer);

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

  std::stringstream oFile;
  writer->getOutput()->saveRegions(oFile,true);
  
  return 0;
}
