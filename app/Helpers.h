#include <iostream>
#include <fstream>
#include <algorithm>    // std::max
#include <stdlib.h>
#include <iomanip>
#include "TSystemDirectory.h"

// Average neutrino decay position in beam coordinates as a function of off axis x: Will be used to set the decay position event-by-event.
double OffAxisPoints[] = {-2,      0.5,     3,      5.5,     8,       10.5,    13,      15.5,    18,      20.5,   23,      25.5,    28,      30.5}; // unit: meters
double meanPDPZ[]      = {93.6072, 93.362,  90.346, 85.6266, 81.1443, 76.6664, 73.0865, 69.8348, 67.5822, 65.005, 62.4821, 60.8336, 59.1433, 57.7352};

float offset[] = { 0., 305., 5. };
float collarLo[] = {-320., -120., 30. };
float collarHi[] = { 320.,  120., 470.};

using namespace std;

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
