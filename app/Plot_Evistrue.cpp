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
  //
  // Read branch from input trees
  //
  TChain *t_effValues = new TChain("effTreeND");
  t_effValues->Add("/pnfs/dune/persistent/users/flynnguo/FDGeoEffinND/FDGeoEff_39254532_99.root");

  // double E_vis_true;                 // True vis energy
  double ND_Gen_numu_E;

  // t_E->SetBranchAddress("E_vis_true",                      &E_vis_true);
  t_effValues->SetBranchAddress("ND_Gen_numu_E",                      &ND_Gen_numu_E);


  // TH1D *hist_E_vis_true = new TH1D("hist_E_vis_true","FD_hist_E_vis_true",100,0,100);
  TH1D *hist_ND_Gen_numu_E = new TH1D("hist_ND_Gen_numu_E","hist_ND_Gen_numu_E",100,0,100);

  // Loop over all events
  int nentries = 0; // Total input events
  int ientry = 0;
  nentries = t_effValues->GetEntries();
  cout<< "nentries:" << nentries<<endl;
  for ( int i = 0; i < (nentries/330); i++ )
  {
    ientry = i*330;
    t_effValues->GetEntry(ientry);
    // hist_E_vis_true->Fill(E_vis_true);
    hist_ND_Gen_numu_E->Fill(ND_Gen_numu_E);
    cout << "ientry: " <<ientry<< ", ND_Gen_numu_E: " << ND_Gen_numu_E << endl;
  }

  // Add labels;
  // Create TCanvas
  TCanvas *c1 = new TCanvas("E_true","E_true",700,500);
  c1->cd(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  // hist_E_vis_true->GetYaxis()->SetTitle("# of events");
  // hist_E_vis_true->GetXaxis()->SetTitle("E_vis_true [MeV]");
  // hist_E_vis_true->Draw();
  hist_ND_Gen_numu_E->GetYaxis()->SetTitle("# of events");
  hist_ND_Gen_numu_E->GetXaxis()->SetTitle("E_ND_Gen_numu [MeV]");
  hist_ND_Gen_numu_E->Draw();


  //Save into root file

  c1->SaveAs("E_true.pdf");


} // end Plot_Evisture
