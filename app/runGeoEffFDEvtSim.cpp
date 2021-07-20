#include "geoEff.h"
#include <iostream>
#include <iomanip>
using namespace std;
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TFitResultPtr.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TEfficiency.h>
#include <TMath.h>
#include "TLorentzVector.h"

#include "Helpers.h"   // Linear interpolate average neutrino decay position as a function of off axis x

using namespace std;

int main(){

  Int_t Run;
  Int_t SubRun;
  Int_t Event;

  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  t->Add("/dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root"); // ntuple path on FNAL dunegpvm machine

  t->SetBranchAddress("Run",    &Run);
  t->SetBranchAddress("SubRun", &SubRun);
	t->SetBranchAddress("Event",  &Event);

  geoEff * eff = new geoEff(314, true);
  eff->setNthrows(1024);

  // Loop over FD events
  int nentries;
  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;

  for ( int i = 0; i < nentries; i++ ) {
    t->GetEntry(i);
    if ( (i % 3) == 0  ) std::cout << "Looking at Events " << i << std::endl;
    std::cout << "run: " << Run << ", subrun: " << SubRun << ", event: " << Event << std::endl;
  } // end loop over events

}
