// Run this macro:
// root -l -b -q FDEffCalc.C
//
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

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>

// Include customized functions and constants
#include "Helpers.h"

void FDEffCalc() {

  bool debug           = false; // Print out for debug purpose
  int nentries         = 0;     // Total input events
  TString TitleX       = "Vertex x [cm]";
  TString TitleY       = "Hadron contain efficiency";
  TString TitleZ       = "Events";
  TString DrawOption   = "COLZ TEXT";
  Float_t Offset       = 0.7;   // 2D histogram Z axis title offset
  Double_t eff_binning = 0.05;  // Efficiency binning

  //
  // Read hadron throw result branch from input tree
  //

  vector<double> *ND_off_axis_pos_vec = 0; // unit: m
  vector<double> *Sim_mu_start_vx     = 0; // unit: cm
  // Nested vector: ND off axis position x, evt vtx x, vetoSize, vetoEnergy, many 64-throw-result chunks
  vector<vector<vector<vector<vector<uint64_t> > > > > *Sim_hadron_throw_result = 0; // Need initialize 0 here to avoid error

  // Read output after transformations are done on FD evts
  TChain *mychain = new TChain("effTreeFD");
  // Output throw result path on FNAL dunegpvm machine
  mychain->Add("/dune/app/users/weishi/NDEff/DUNE_ND_GeoEff/bin/Output_FDGeoEff.root");
  mychain->SetBranchAddress("ND_off_axis_pos_vec",     &ND_off_axis_pos_vec); // is a vector for a event: entries = written evts * ND_off_axis_pos_steps
  mychain->SetBranchAddress("Sim_mu_start_vx",         &Sim_mu_start_vx);     // is a vector for a event: entries = written evts * vtx_vx_steps, equivalent to b_vtx_vx_vec
  mychain->SetBranchAddress("Sim_hadron_throw_result", &Sim_hadron_throw_result);

  //
  // 2D histograms for events for each ND off axis pos
  //

  // Set the number of histograms
  int nplot = ( OffAxisPoints[13] - OffAxisPoints[0] ) / ND_off_axis_pos_stepsize + 1; // Use variables defined in Helpers.h
  // Create a dynamic array of pointers
  TH2F** NdOffAxisPos = new TH2F*[nplot];
  for ( int ihist = 0; ihist < nplot; ihist++ ) {
    // Calculate the off axis pos for each histogram
    NdOffAxisPos[ihist] = new TH2F(TString::Format( "NdOffAxisPos_%.2f_m", ihist*ND_off_axis_pos_stepsize + OffAxisPoints[0] ), TString::Format( "ND off-axis position: %.2f m", ihist*ND_off_axis_pos_stepsize + OffAxisPoints[0] ), ( ND_local_x_max - ND_local_x_min ) / ND_local_x_stepsize + 1, ND_local_x_min, ND_local_x_max + ND_local_x_stepsize, int(1/eff_binning)+1, 0., 1 + eff_binning);
  }

  //
  // Loop over written events
  //

  nentries = mychain->GetEntries();

  std::cout << "Tot written evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < nentries; ientry++ ) {

    mychain->GetEntry(ientry);

    if ( ientry%10000 == 0 ) std::cout << "Looking at entry " << ientry << std::endl;

    // First loop over ND off_axis vector
    //   then loop over evt vtx x
    //     then [0][0][k-th chunk]
    //     or [i][j][k-th chunk], i for vetoSize, j for vetoEnergy

    if ( debug ) std::cout << "ND off axis size: " << Sim_hadron_throw_result->size() << std::endl;

    int counter1 = 0;

    for ( vector<vector<vector<vector<vector<uint64_t> > > > >::iterator it_nd_off_axis_pos = Sim_hadron_throw_result->begin(); it_nd_off_axis_pos != Sim_hadron_throw_result->end(); ++it_nd_off_axis_pos ) {

      counter1++;
      if ( debug ) std::cout << "  #" << counter1 << ", evt vtx x size: " << it_nd_off_axis_pos->size() << std::endl;

      int counter2 = 0;

      for ( vector<vector<vector<vector<uint64_t> > > >::iterator it_evt_vtx_x = it_nd_off_axis_pos->begin(); it_evt_vtx_x != it_nd_off_axis_pos->end(); ++it_evt_vtx_x ) {

        counter2++;
        if ( debug ) std::cout << "    #" << counter2 << ", vetoSize size: " << it_evt_vtx_x->size() << std::endl;

        int counter3 = 0;

        for ( vector<vector<vector<uint64_t> > >::iterator it_veto_size = it_evt_vtx_x->begin(); it_veto_size != it_evt_vtx_x->end(); ++it_veto_size ) {

          counter3++;
          if ( debug ) std::cout << "      #" << counter3 << ", vetoEnergy size: " << it_veto_size->size() << std::endl;

          int counter4 = 0;

          for ( vector<vector<uint64_t> >::iterator it_veto_energy = it_veto_size->begin(); it_veto_energy != it_veto_size->end(); ++it_veto_energy ) {

            // Every 64 throw result is a chunk
            // current test case: each evt has 128 throws, so 2 chunks
            counter4++;
            if ( debug ) std::cout << "        #" << counter4 << ", 64-throw size (chunks): " << it_veto_energy->size() << std::endl;

            int counter5 = 0;
            int pass = 0; // count no. of throws that pass containment requirement

            for ( vector<uint64_t>::iterator it_chunk = it_veto_energy->begin(); it_chunk != it_veto_energy->end(); ++it_chunk ) {

              counter5++;
              if ( debug ) std::cout << "          chunk #" << counter5 << ": " << *it_chunk << std::endl;

              for ( unsigned int ithrow = 0; ithrow < 64; ithrow++ ) {

                // Access per throw result for the evt
                // Example (in ROOT):
                //   In:  uint64_t number = 18446744073709551615 (this is 2^64-1, the max value for uint64_t)
                //   In:  number & ((uint64_t)1) << 0
                //   Out: (unsigned long long) 1
                //   In:  number & ((uint64_t)1) << 1
                //   Out: (unsigned long long) 2
                //   ...
                //   In:  number & ((uint64_t)1) << 63
                //   Out: (unsigned long long) 9223372036854775808
                // A non-zero result indicates the throw result is true, zero indicates the throw result is false

                uint64_t throw_result = (*it_chunk) & ( ((uint64_t)1)<<(ithrow%64) );

                if ( debug ) std::cout << "                  throw #" << ithrow+1 << ": " << throw_result << std::endl;

                // Count no. of throws that pass containment requirement
                if ( throw_result != 0 ) pass++;

              } // loop over 64 throw results in a chunk

            }   // end loop over 64-bit chunks

            // Calculate per-event hadron containment efficiency: Sum(N throw booleans)/N
            double hadron_contain_eff = pass*1.0/(64*counter5); // Total throws N = 64 * counter5
            if ( debug ) std::cout << "        Passed throws: " << pass << ", tot. throws: " << 64*counter5 << ", eff: " << hadron_contain_eff << ", evt vtx x [cm]: " << Sim_mu_start_vx->at( counter2 -1 ) << ", nd off-axis pos [m]: " << ND_off_axis_pos_vec->at( counter1 -1 ) << std::endl;

            //
            // 2D histogram: efficiency - evt vtx x, for ND-off-axis[1]-[14]
            //

            // Skip first element of evt vtx x as it is random
            if ( counter2 - 1 > 0 ) {

              // Fill histograms
              for ( int ifill = 0; ifill < nplot; ifill++ ) {

                // Also skip first element of ND-off-axis as it is random
                if ( counter1 - 1 == ifill + 1 ) NdOffAxisPos[ifill]->Fill(Sim_mu_start_vx->at( counter2 -1 ), hadron_contain_eff);

              } // end for fill

            } // end skip first element of evt vtx x

            //
            // 2D histogram: efficiency - hadron E, for one ND-off-axis and one evt vtx x
            //

          }     // end loop over veto energy

        }       // end loop over veto size

      }         // end loop over evt vtx x

    }           // end loop over ND off axis pos

  }             // end loop over events


  // As a function of Sim nu E

  // Write trees
  TFile * outFile = new TFile("FDEffCalc.root", "RECREATE");

  // Create a dynamic array of pointers
  TCanvas** c_NdOffAxisPos = new TCanvas*[nplot];

  for ( int ic = 0; ic < nplot; ic++ ) {

    // Draw the histogram
    c_NdOffAxisPos[ic] = new TCanvas(TString::Format( "c_NdOffAxisPos_%.2f_m", ic*ND_off_axis_pos_stepsize + OffAxisPoints[0] ), TString::Format( "ND off-axis position: %.2f m", ic*ND_off_axis_pos_stepsize + OffAxisPoints[0] ), 700, 500);
    c_NdOffAxisPos[ic]->cd();
    NdOffAxisPos[ic]->GetXaxis()->SetTitle( TitleX.Data() );
    NdOffAxisPos[ic]->GetYaxis()->SetTitle( TitleY.Data() );
    NdOffAxisPos[ic]->GetZaxis()->SetTitle( TitleZ.Data() );
    NdOffAxisPos[ic]->GetZaxis()->SetTitleOffset(Offset);
    NdOffAxisPos[ic]->SetStats(0);
    NdOffAxisPos[ic]->Draw( DrawOption.Data() );
    c_NdOffAxisPos[ic]->Write();
  }

  outFile->Close();

} // end macro
