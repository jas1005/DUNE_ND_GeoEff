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

void FDEffCalc() {

  bool debug = true;

  int nentries = 0;                   // Total input events

  //
  // Read hadron throw result branch from input tree
  //

  // Result of hadron containment is stored in nested vectors, need to generate dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*uint64_t*"); // Remove old dictionary if exists
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*vector*uint64_t*");
  // Generate new dictionary
  gInterpreter->GenerateDictionary("vector<vector<uint64_t> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<uint64_t> > > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<vector<uint64_t> > > > >", "vector");
  // Nested vector: ND off axis position x, evt vtx x, vetoSize, vetoEnergy, many 64-throw-result chunks
  vector<vector<vector<vector<vector<uint64_t> > > > > *Sim_hadron_throw_result = 0; // Need initialize 0 here to avoid error

  // Read output after transformations are done on FD evts
  TChain *mychain = new TChain("effTreeFD");
  // Output throw result path on FNAL dunegpvm machine
  mychain->Add("/dune/app/users/weishi/NDEff/DUNE_ND_GeoEff/bin/Output_FDGeoEff.root");

  mychain->SetBranchAddress("Sim_hadron_throw_result",     &Sim_hadron_throw_result);

  //
  // Test with 128 throws/evt, so 2 chunks, only one vetoSize and one vetoEnergy: [0][0][0] and [0][0][1]
  //

  // First loop over ND off_axis vector
  //   then loop over evt vtx x
  //     then [0][0][k-th chunk]
  //     or [i][j][k-th chunk], i for vetoSize, j for vetoEnergy

  // Per-evt hadron containment efficiency: Sum(128 booleans)/128
  // Example in ROOT:
  //   In:  uint64_t number = 18446744073709551615 (this is 2^64-1, the max value for uint64_t)
  //   In:  number & ((uint64_t)1) << 0
  //   Out: (unsigned long long) 1
  //   In:  number & ((uint64_t)1) << 1
  //   Out: (unsigned long long) 2
  //   ...
  //   In:  number & ((uint64_t)1) << 63
  //   Out: (unsigned long long) 9223372036854775808
  // A non-zero result indicates the throw result is true, zero indicates the throw result is false

  // 2D histogram: efficiency - evt vtx x, for ND-off-axis[1]-[6], [0] is random, skip for now

  nentries = mychain->GetEntries();

  if (debug) std::cout << "Tot written evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < 3; ientry++ ) {

    mychain->GetEntry(ientry);

    if (debug) std::cout << "Looking at entry " << ientry << std::endl;

    if (debug) std::cout << "ND off axis size: " << Sim_hadron_throw_result->size() << std::endl;

    int counter1 = 0;

    for ( vector<vector<vector<vector<vector<uint64_t> > > > >::iterator it_nd_off_axis_pos = Sim_hadron_throw_result->begin(); it_nd_off_axis_pos != Sim_hadron_throw_result->end(); ++it_nd_off_axis_pos ) {

      counter1++;
      if (debug) std::cout << "  #" << counter1 << ", evt vtx x size: " << it_nd_off_axis_pos->size() << std::endl;

      int counter2 = 0;

      for ( vector<vector<vector<vector<uint64_t> > > >::iterator it_evt_vtx_x = it_nd_off_axis_pos->begin(); it_evt_vtx_x != it_nd_off_axis_pos->end(); ++it_evt_vtx_x ) {

        counter2++;
        if (debug) std::cout << "    #" << counter2 << ", vetoSize size: " << it_evt_vtx_x->size() << std::endl;

        int counter3 = 0;

        for ( vector<vector<vector<uint64_t> > >::iterator it_veto_size = it_evt_vtx_x->begin(); it_veto_size != it_evt_vtx_x->end(); ++it_veto_size ) {

          counter3++;
          if (debug) std::cout << "      #" << counter3 << ", vetoEnergy size: " << it_veto_size->size() << std::endl;

          int counter4 = 0;

          for ( vector<vector<uint64_t> >::iterator it_veto_energy = it_veto_size->begin(); it_veto_energy != it_veto_size->end(); ++it_veto_energy ) {

            // Every 64 throw result is a chunk, current test case: each evt has 128 throws, 2 chunks
            counter4++;
            if (debug) std::cout << "        #" << counter4 << ", 64-throw size (chunks): " << it_veto_energy->size() << std::endl;

            int counter5 = 0;

            for ( vector<uint64_t>::iterator it_chunk = it_veto_energy->begin(); it_chunk != it_veto_energy->end(); ++it_chunk ) {

              counter5++;
              if (debug) std::cout << "          chunk #" << counter5 << ": " << *it_chunk << std::endl;

              for ( unsigned int ithrow = 0; ithrow < 64; ithrow++ ) {

                uint64_t throw_result = (*it_chunk) & ( ((uint64_t)1)<<(ithrow%64) );

                if (debug) std::cout << "                  throw #" << ithrow+1 << ": " << throw_result << std::endl;

              } // loop over 64 throw results in a chunk

            } // end loop over 64-bit chunks

          } // end loop over veto energy

        } // end loop over veto size

      } // end loop over evt vtx x

    } // end loop over ND off axis pos

  } // end loop over events


  // As a function of Sim nu E

} // end macro
