#include <iostream>
#include <fstream>
#include <algorithm>    // std::max
#include <stdlib.h>
#include <iomanip>
#include "TSystemDirectory.h"

// Average neutrino decay position in beam coordinates as a function of off axis x: Will be used to set the decay position event-by-event.
double OffAxisPoints[14] = {-2,      0.5,     3,      5.5,     8,       10.5,    13,      15.5,    18,      20.5,   23,      25.5,    28,      30.5};
double meanPDPZ[14]      = {93.6072, 93.362,  90.346, 85.6266, 81.1443, 76.6664, 73.0865, 69.8348, 67.5822, 65.005, 62.4821, 60.8336, 59.1433, 57.7352};

using namespace std;

double My_meanPDPZ(double x) {

  double result = 0.0;

  // Linear interpolation
  if ( x >= -2                && x < OffAxisPoints[1] ) result = meanPDPZ[0] + ( x - OffAxisPoints[0] )*( meanPDPZ[1]  - meanPDPZ[0] )/( OffAxisPoints[1]  - OffAxisPoints[0] );
  if ( x >= OffAxisPoints[1]  && x < OffAxisPoints[2] ) result = meanPDPZ[1] + ( x - OffAxisPoints[1] )*( meanPDPZ[2]  - meanPDPZ[1] )/( OffAxisPoints[2]  - OffAxisPoints[1] );
  if ( x >= OffAxisPoints[2]  && x < OffAxisPoints[3] ) result = meanPDPZ[2] + ( x - OffAxisPoints[2] )*( meanPDPZ[3]  - meanPDPZ[2] )/( OffAxisPoints[3]  - OffAxisPoints[2] );
  if ( x >= OffAxisPoints[3]  && x < OffAxisPoints[4] ) result = meanPDPZ[3] + ( x - OffAxisPoints[3] )*( meanPDPZ[4]  - meanPDPZ[3] )/( OffAxisPoints[4]  - OffAxisPoints[3] );
  if ( x >= OffAxisPoints[4]  && x < OffAxisPoints[5] ) result = meanPDPZ[4] + ( x - OffAxisPoints[4] )*( meanPDPZ[5]  - meanPDPZ[4] )/( OffAxisPoints[5]  - OffAxisPoints[4] );
  if ( x >= OffAxisPoints[5]  && x < OffAxisPoints[6] ) result = meanPDPZ[5] + ( x - OffAxisPoints[5] )*( meanPDPZ[6]  - meanPDPZ[5] )/( OffAxisPoints[6]  - OffAxisPoints[5] );
  if ( x >= OffAxisPoints[6]  && x < OffAxisPoints[7] ) result = meanPDPZ[6] + ( x - OffAxisPoints[6] )*( meanPDPZ[7]  - meanPDPZ[6] )/( OffAxisPoints[7]  - OffAxisPoints[6] );
  if ( x >= OffAxisPoints[7]  && x < OffAxisPoints[8] ) result = meanPDPZ[7] + ( x - OffAxisPoints[7] )*( meanPDPZ[8]  - meanPDPZ[7] )/( OffAxisPoints[8]  - OffAxisPoints[7] );
  if ( x >= OffAxisPoints[8]  && x < OffAxisPoints[9] ) result = meanPDPZ[8] + ( x - OffAxisPoints[8] )*( meanPDPZ[9]  - meanPDPZ[8] )/( OffAxisPoints[9]  - OffAxisPoints[8] );
  if ( x >= OffAxisPoints[9]  && x < 30.5  )            result = meanPDPZ[9] + ( x - OffAxisPoints[9] )*( meanPDPZ[10] - meanPDPZ[9] )/( OffAxisPoints[10] - OffAxisPoints[9] );

  return result;
}

void addfilesMany(TChain *ch, const std::vector<string>& v, const TString ext=".root") {
  bool verbose(true);
  int nfiles=0;
  for (std::string dirname : v) {
    TSystemDirectory dir(dirname.c_str(), dirname.c_str());
    TList *files = dir.GetListOfFiles();
    if (files) {
      if (verbose) std::cout << " Found files" << std::endl;
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
        fname = file->GetName();
        if (verbose) std::cout << " file name: " << dirname + fname << std::endl;
        if (!file->IsDirectory() && fname.BeginsWith(ext)) {
          nfiles++;
          if (verbose) std::cout << " adding #" << nfiles << ": " << dirname + fname << std::endl;
          ch->Add(dirname + fname);
        }
      }
    }
  }
}
