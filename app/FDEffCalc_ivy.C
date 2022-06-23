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

// Library methods
#include "geoEff.h"

// Include customized functions and constants
#include "Helpers.h"

void FDEffCalc_ivy() {

  bool debug             = false; // Print out for debug purpose
  // Inout file on FNAL dunegpvm machine
  TString FileIn         = "/home/fyguo/NDEff/DUNE_ND_GeoEff/bin/Output_FDGeoEff_ivy.root";
  double effthreshold    = 0.;    // b/t 0 and 1: only use evts with hadron contain eff above the threshold in some plots
  Double_t eff_binning   = 0.05;  // A positive number
  Double_t Vtx_x_binning = 1.;    // cm, a positive number
  Double_t Vtx_y_binning = 1.;
  Double_t Vtx_z_binning = 1.;
  Double_t PlotMinY      = 0.5;   // events
  Double_t PlotMaxY      = 2000.; // events
  double LepE_low        = 0.;
  double LepE_up         = 30.;   // GeV
  double LepE_binning    = 1;     // GeV, a positive number
  double HadE_low        = 0.;
  double HadE_up         = 2000.; // MeV
  double HadE_binning    = 100;   // MeV, a positive number
  int nentries           = 0;     // Total input events

  //
  // Read branch from input trees
  //

  // Initialize
  double ND_Gen_numu_E;
  vector<double> *ND_vtx_vx_vec=0;
  vector<double> *ND_off_axis_pos_vec=0;
  vector<vector<vector<vector<vector<uint64_t> > > > > *ND_Sim_hadron_throw_result=0; // ...... ....... .. ... .... ........ .. ... ... .. ......... .........., 64-throw-result chunks
  vector<vector<vector<uint64_t> > > *hadron_throw_result=0;
  // array: (x,y,z)<->dim=(0,1,2)
  double ND_OffAxis_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_OnAxis_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_OffAxis_Sim_mu_start_p[3]; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double ND_OffAxis_Sim_mu_start_E;
  double ND_Sim_hadronic_Edep_a2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]

  // Random throws
  vector<float> *throwVtxY = 0; // unit: cm
  vector<float> *throwVtxZ = 0;

  // Read ND info
  TChain *myhadronchain = new TChain("effTreeFD");
  myhadronchain->Add( FileIn.Data() );
  myhadronchain->SetBranchAddress("ND_Gen_numu_E",                          &ND_Gen_numu_E);
  myhadronchain->SetBranchAddress("ND_off_axis_pos_vec",                    &ND_off_axis_pos_vec);
  myhadronchain->SetBranchAddress("ND_vtx_vx_vec",                          &ND_vtx_vx_vec);
  // myhadronchain->SetBranchAddress("ND_Sim_hadron_throw_result",             &ND_Sim_hadron_throw_result);
  myhadronchain->SetBranchAddress("ND_OffAxis_Sim_mu_start_v",              &ND_OffAxis_Sim_mu_start_v);
  myhadronchain->SetBranchAddress("ND_OnAxis_Sim_mu_start_v",               &ND_OnAxis_Sim_mu_start_v);
  myhadronchain->SetBranchAddress("ND_OffAxis_Sim_mu_start_p",              &ND_OffAxis_Sim_mu_start_p);
  myhadronchain->SetBranchAddress("ND_OffAxis_Sim_mu_start_E",              &ND_OffAxis_Sim_mu_start_E);
  myhadronchain->SetBranchAddress("ND_Sim_hadronic_Edep_a2",                &ND_Sim_hadronic_Edep_a2);    // tot hadron deposited E in the evt
  myhadronchain->SetBranchAddress("hadron_throw_result",                     &hadron_throw_result);
  // Read hadron throw result
  TChain *myhadronthrowchain = new TChain("effThrow");
  myhadronthrowchain->Add( FileIn.Data() );
  myhadronthrowchain->SetBranchAddress("ND_Sim_hadron_throw_result",             &ND_Sim_hadron_throw_result);

  // Read random throws each FD evt
  TChain *mythrowchain = new TChain("ThrowsFD");
  mythrowchain->Add( FileIn.Data() );
  mythrowchain->SetBranchAddress("throwVtxY", &throwVtxY); // vector<float>: entries = [ (int)(written evts / 100) + 1 ] * N_throws
  mythrowchain->SetBranchAddress("throwVtxZ", &throwVtxZ);
  //
  // geoEff * eff = new geoEff(314, false); // set verbose to true for debug
  // Add output txt file
  ofstream myfile;
   myfile.open ("Output_GeoEff_ivy.txt");
   bool verbose = true;
  //
  // Loop over written events
  //
  nentries = myhadronchain->GetEntries();

  if ( debug ) std::cout << "Tot written evts: " << nentries << std::endl;
  myfile << "Tot written evts: " << nentries << '\n';
  for ( int iwritten = 0; iwritten < nentries; iwritten++ ) {


    myhadronchain->GetEntry(iwritten);
    if ( verbose ) myfile << "Looking at iwritten " << iwritten << '\n';
    cout << " Looking at iwritten " << iwritten<< endl;
    int ND_off_axis_pos_vec_counter = iwritten%(ND_off_axis_pos_vec->size());
    int ND_vtx_vx_vec_counter = iwritten%(ND_vtx_vx_vec->size());

    if ( debug ) std::cout << "N_throws set #" << iwritten/100 << std::endl;

    mythrowchain->GetEntry(iwritten/100);
    if ( debug ) std::cout << "throwVtxY size: " << throwVtxY->size() << std::endl; // this should be N_throws set by user
    // First loop over ND off_axis vector
    //   then loop over evt vtx x
    //     then [0][0][k-th chunk]
    //     or [i][j][k-th chunk], i for vetoSize, j for vetoEnergy

    myhadronthrowchain->GetEntry(iwritten/9);
    std::cout << "ND off axis size: " << ND_Sim_hadron_throw_result->size() << std::endl;
    if ( verbose ) myfile << "ND off axis size: " << ND_Sim_hadron_throw_result->size() << '\n';
    // cout << "hadron_throw_result->size(): " << hadron_throw_result->size() << endl;
    int counter1 = 0;

        int counter3 = 0;

        for ( vector<vector<vector<uint64_t> > >::iterator it_veto_size = hadron_throw_result->begin(); it_veto_size != hadron_throw_result->end(); ++it_veto_size ) {

          counter3++;
          if ( debug ) std::cout << "      #" << counter3 << ", vetoEnergy size: " << it_veto_size->size() << std::endl;

          int counter4 = 0;

          for ( vector<vector<uint64_t> >::iterator it_veto_energy = it_veto_size->begin(); it_veto_energy != it_veto_size->end(); ++it_veto_energy ) {

            counter4++;

            // Every 64 throw result is a chunk
            // current test case: each evt has 128 throws, so 2 chunks
            if ( debug ) std::cout << "        #" << counter4 << ", 64-throw size (chunks): " << it_veto_energy->size() << std::endl;

            int counter5    = 0;
            int validthrows = 0; // count no. of throws that meet ND FV cut for this evt
            int hadronpass  = 0; // count no. of throws that meet hadron containment cut
            double hadron_contain_eff = 0.; // initialize eff to zero

            for ( vector<uint64_t>::iterator it_chunk = it_veto_energy->begin(); it_chunk != it_veto_energy->end(); ++it_chunk ) {

              counter5++;
              if ( debug ) std::cout << "          chunk #" << counter5 << ": " << *it_chunk << std::endl;
              if ( debug ) std::cout << "                    throws passed ND FV cut" << std::endl;

              for ( unsigned int ithrow = 0; ithrow < 64; ithrow++ ) {
                // myfile << "IsInNDFV: " << ND_vtx_vx_vec->at( counter2 - 1 ) << ", "<< throwVtxY->at( (counter5-1)*64 + ithrow )-5.5 << ", "<< throwVtxZ->at( (counter5-1)*64 + ithrow )-411 << '\n';

                // For the numerator, only consider throws where throwed FD evt vtx x/y/z is in ND FV, same as what is done for ND evts
                // For now, we use mu start pos as evt vtx pos, random throws for y/z are stored in the ThrowsFD tree
                if ( FDEffCalc_cfg::IsInNDFV(ND_vtx_vx_vec->at(ND_vtx_vx_vec_counter), throwVtxY->at( (counter5-1)*64 + ithrow )-geoEff::getCurrentOffset(1), throwVtxZ->at( (counter5-1)*64 + ithrow )-geoEff::getCurrentOffset(2)) ) {
                  // if ( FDEffCalc_cfg::IsInNDFV(ND_vtx_vx_vec->at(ND_vtx_vx_vec_counter), throwVtxY->at( (counter5-1)*64 + ithrow )-5.5, throwVtxZ->at( (counter5-1)*64 + ithrow )-411. )) {
                    validthrows++;

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

                  if ( debug ) std::cout << "                    throw #" << ithrow+1 << ": " << throw_result << std::endl;

                  // Count no. of throws passed hadron containment requirement
                  if ( throw_result != 0 ) hadronpass++;

                } // end if FD event passed ND FV cut

              }   // end loop over 64 throws in a chunk

            }     // end loop over 64-throw chunks
            // cout << "counter1: " << counter1 << endl;
            // cout << "counter2: " << counter2 << endl;
            // cout << "counter3: " << counter3 << endl;
            // cout << "counter4: " << counter4 << endl;
            // cout << "counter5: " << counter5 << endl;
            //
            // Calculate per-event hadron containment efficiency from throws
            //

            // If a throwed evt vtx is in ND dead region, it is not going to be detected by ND, but probably will be detected in FD (assume)
            // In principle, we should NOT exclude such throws in the denominator

            //hadron_contain_eff = hadronpass*1.0/N_throws; // N_throws also equals 64 * counter5
            //if ( debug ) std::cout << "        Passed throws: " << hadronpass << ", eff: " << hadron_contain_eff << ", evt vtx x [cm]: " << ND_OffAxis_Sim_mu_start_v[0]->at( counter2 - 1 ) << ", nd off-axis pos [m]: " << ND_off_axis_pos_vec->at( counter1 -1 ) << std::endl;

            // But for the ND FV cut, we could probably factorize it analytically instead of using these random throws to evaluate
            // therefore we exclude such throws from the denominator as well
            if ( validthrows > 0 ) hadron_contain_eff = hadronpass*1.0/validthrows;
            // myfile << "iwritten: " << iwritten << '\n';
            // myfile << "hadron_contain_eff: " << hadron_contain_eff << '\n';
            if ( debug ) std::cout << "        Passed throws: " << hadronpass << ", tot. valid throws: " << validthrows << ", eff: " << hadron_contain_eff << ", nd off-axis pos (evt vtx x) [cm]: " << ND_OffAxis_Sim_mu_start_v[0] << std::endl;
            myfile << "        Passed throws: " << hadronpass << ", tot. valid throws: " << validthrows << ", eff: " << hadron_contain_eff << ", evt vtx x [cm]: " << ND_vtx_vx_vec->at(ND_vtx_vx_vec_counter)<< ", nd off-axis pos [m]: " << ND_off_axis_pos_vec->at(ND_off_axis_pos_vec_counter)<< "\n";

          }     // end loop over veto energy
        }       // end loop over veto size

  }             // end loop over events


  myfile.close();

} // end macro
