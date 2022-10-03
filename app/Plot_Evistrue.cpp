// ROOT includes
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
#include <TRandom3.h>
#include "TSystem.h"
#include "TROOT.h"
#include <TGraph2D.h>
#include <TRandom.h>
#include <TF2.h>

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector> // Need this for generate dictionary for nested vectors
using namespace std;

// Include customized functions and constants
#include "Helpers.h"

void Plot_Evistrue() // /pnfs/dune/persistent/users/flynnguo/myFDntuples/myntuple_61454381_*.root
{
  gROOT->Reset();
  // Input FDroot file
  TString FileIn = "/pnfs/dune/persistent/users/flynnguo/myFDntuples/myntuple_61454381_991.root";
  //
  // Read branch from input trees
  //
  TChain *t_E = new TChain("MyEnergyAnalysis/MyTree");
  t_E->Add(FileIn.Data());

  double E_vis_true;                 // True vis energy
  t_E->SetBranchAddress("E_vis_true",                      &E_vis_true);

  // Create TCanvas
  TCanvas *c1 = new TCanvas("E_vis_true","E_vis_true",700,500);
  c1->cd(1);
  TH1D *hist_E_vis_true = new TH1D("hist_E_vis_true","hist_E_vis_true",100,0,100);

  // Loop over all events
  int nentries = 0; // Total input events
  nentries = t_E->GetEntries();
  for ( int ientry = 0; ientry < nentries; ientry++ )
  {
    hist_E_vis_true->Fill(E_vis_true);
  }

  // Add labels;
  hist_E_vis_true->GetYaxis()->SetTitle("# of events");
  hist_E_vis_true->GetXaxis()->SetTitle("E_vis_true [MeV]");
  hist_E_vis_true->Draw();


  //Save into root file
  c1->SaveAs("E_vis_true.pdf");


} // end Plot_Evisture
