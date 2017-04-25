#ifndef PDFREWEIGHT_H
#define PDFREWEIGHT_H

#include <vector>
#include <string>

class AnalysisEvent;
namespace LHAPDF {
  class PDF;
}

class PDFReweighter
{
 public:
  PDFReweighter(std::vector<std::string>& options);
  double reweightEvent(AnalysisEvent *event);

 private:
  double inputEnergy;
  double outputEnergy;
  LHAPDF::PDF* pdf;
};


#endif
