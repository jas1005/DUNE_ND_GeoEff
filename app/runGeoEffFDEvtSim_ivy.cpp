// Library methods
#include "geoEff.h"

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
#include <vector> // Need this for generate dictionary for nested vectors

// Include customized functions and constants
#include "Helpers.h"




int main(){

  //
  // Branches to be read from n-tuple produced from FD MC
  //
  int FD_Run; // # of the run being processed
  int FD_SubRun; // # of the sub-run being processed
  int FD_Event; // # of the event being processed
  int FD_Sim_nNumu; // # of Sim muon neutrinos (numu and numubar)
  double FD_Gen_numu_E; // Energy of generator level neutrino [GeV]
  int FD_Sim_nMu; // # of Sim muons (mu+/mu-)
  int FD_CCNC_truth; // 0 =CC 1 =NC
  int FD_neuPDG; // Generator level neutrino PDG
  double FD_Sim_mu_start_vx; // Position of the muon trajectory at start point on the x-axis [cm]
  double FD_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double FD_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  double FD_Sim_mu_end_vx; // Position of the muon trajectory at end point on the x-axis [cm]
  double FD_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double FD_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double FD_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double FD_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double FD_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  double FD_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  double FD_Sim_mu_end_px; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double FD_Sim_mu_end_py; // Momentum of the muon trajectory at end point on the y-axis [GeV]
  double FD_Sim_mu_end_pz; // Momentum of the muon trajectory at end point on the z-axis [GeV]
  double FD_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  double FD_Sim_hadronic_Edep_a2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  int FD_Sim_n_hadronic_Edep_a; // # of hadronic energy deposits
  //int FD_Sim_n_hadronic_Edep_a_vector;  // # of add up hadronic energy deposits
  vector<float> *FD_Sim_hadronic_hit_Edep_a2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_a = 0; // Position of each energy deposit on the x-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_y_a = 0; // Position of each energy deposit on the y-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_z_a = 0; // Position of each energy deposit on the z-axis [cm]

  // Read ntuple from FD MC
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // Ntuple path on FNAL dunegpvm machine
  // For Ivy machine:
  t->Add("/home/fyguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple_ivy.root");
  // For FNAL machine:
  // t->Add("/dune/app/users/flynnguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root");

  t->SetBranchAddress("Run",                      &FD_Run);
  t->SetBranchAddress("SubRun",                   &FD_SubRun);
  t->SetBranchAddress("Event",                    &FD_Event);
  t->SetBranchAddress("Sim_nNumu",                &FD_Sim_nNumu);
  t->SetBranchAddress("Gen_numu_E",               &FD_Gen_numu_E);
  t->SetBranchAddress("Sim_nMu",                  &FD_Sim_nMu);
  t->SetBranchAddress("CCNC_truth",               &FD_CCNC_truth);
  t->SetBranchAddress("neuPDG",                   &FD_neuPDG);
  t->SetBranchAddress("Sim_mu_start_vx",          &FD_Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",          &FD_Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",          &FD_Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",            &FD_Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",            &FD_Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",            &FD_Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",          &FD_Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",          &FD_Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",          &FD_Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_start_E",           &FD_Sim_mu_start_E);
  t->SetBranchAddress("Sim_mu_end_px",            &FD_Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",            &FD_Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",            &FD_Sim_mu_end_pz);
  t->SetBranchAddress("Sim_mu_end_E",             &FD_Sim_mu_end_E);
  t->SetBranchAddress("Sim_hadronic_Edep_a2",     &FD_Sim_hadronic_Edep_a2);
  t->SetBranchAddress("Sim_n_hadronic_Edep_a",    &FD_Sim_n_hadronic_Edep_a);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_a2", &FD_Sim_hadronic_hit_Edep_a2);
  t->SetBranchAddress("Sim_hadronic_hit_x_a",     &FD_Sim_hadronic_hit_x_a);
  t->SetBranchAddress("Sim_hadronic_hit_y_a",     &FD_Sim_hadronic_hit_y_a);
  t->SetBranchAddress("Sim_hadronic_hit_z_a",     &FD_Sim_hadronic_hit_z_a);

  //
  // Declare variables used in this program
  //
  int nentries = 0; // Total input events
  float vetoEnergyFD; // Total hadron deposited energy in FD veto region
  float vetoEnergyFD_new; // Add new restriction to vetoEnergy
  int iwritten = 0; // Output event counter
  float decayZbeamCoord; // Decay point (neutrino production point) in beam coordinate [cm]
  float decayXdetCoord; // Decay point (neutrino production point) in detector coordinate on the x-axis [cm]
  float decayYdetCoord; // Decay point (neutrino production point) in detector coordinate on the y-axis [cm]
  float decayZdetCoord; // Decay point (neutrino production point) in detector coordinate on the z-axis [cm]
  vector<float> HadronHitEdeps; // Hadron hit segment energy deposits [MeV]
  vector<float> HadronHitPoss;  // Hadron hit segment energy deposits position [cm]

  vector<double> ND_off_axis_pos_vec = {0,-7*100,-30*100}; // ND off-axis positions' choices for each FD evt [cm]
  vector<double> ND_vtx_vx_vec={-200,0,200};          // Vtx x choices for each FD evt in ND volume [cm]
  //vector<double> ND_off_axis_pos_vec; // unit: meters, ND off-axis choices for each FD evt: 1st element is randomized for each evt
  //vector<double> ND_vtx_vx_vec;          // unit: cm, vtx x choices for each FD evt in ND volume: 1st element is randomized for each evt

  // int ND_off_axis_pos_steps = 0;
  // int vtx_vx_steps = 0;
  //
  // // Initialize first element as -999, to be replaced by a random off-axis nd pos in each evt below
  // ND_off_axis_pos_vec.emplace_back(-999.);
  //
  // if ( ND_off_axis_pos_stepsize > 0 && ND_off_axis_pos_stepsize <= OffAxisPoints[13] ) {
  //   ND_off_axis_pos_steps = ( OffAxisPoints[13] - OffAxisPoints[0] ) / ND_off_axis_pos_stepsize;
  // }
  // else std::cout << "Error: please set the ND_off_axis_pos_stepsize above 0 and below max element of OffAxisPoints." << std::endl;
  //
  // if (verbose) std::cout << "ND_off_axis_pos_steps: " << ND_off_axis_pos_steps << std::endl;
  //
  // // The rest elements follow fixed increments from min ND local x
  // for ( int i_ND_off_axis_pos_step = 0; i_ND_off_axis_pos_step < ND_off_axis_pos_steps + 1; i_ND_off_axis_pos_step++ ){
  //   ND_off_axis_pos_vec.emplace_back( i_ND_off_axis_pos_step*ND_off_axis_pos_stepsize + OffAxisPoints[0] );
  // }
  //
  // if (verbose) std::cout << "ND_off_axis_pos_vec size: "<< ND_off_axis_pos_vec.size() << std::endl;
  //
  // // Initialize first element as -999, to be replaced by a random vtx x in each evt below
  // ND_vtx_vx_vec.emplace_back(-999.);
  //
  // if ( ND_local_x_stepsize > 0 && ND_local_x_stepsize <= ND_local_x_max ) {
  //   vtx_vx_steps = ( ND_local_x_max - ND_local_x_min ) / ND_local_x_stepsize;
  // }
  // else std::cout << "Error: please set the ND_local_x_stepsize above 0 and below ND_local_x_max." << std::endl;
  //
  // if (verbose) std::cout << "vtx_vx_steps: " << vtx_vx_steps << std::endl;
  //
  // // The rest elements follow fixed increments from min ND local x
  // for ( int i_vtx_vx_step = 0; i_vtx_vx_step < vtx_vx_steps + 1; i_vtx_vx_step++ ){
  //   ND_vtx_vx_vec.emplace_back( i_vtx_vx_step*ND_local_x_stepsize + ND_local_x_min );
  // }
  //
  // if (verbose) std::cout << "ND_vtx_vx_vec size: "<< ND_vtx_vx_vec.size() << std::endl;
  //


  //
  // Move all info from FD to ND
  //
  // Lepton info: expressed in ND coordinate sys, do not confuse with branches read above in FD coordinate sys
  double ND_Gen_numu_E; // Energy of generator level neutrino [GeV]
  // Add veto E with two different restrictions
  double ND_vetoEnergyFD; // Total hadron deposited energy in FD veto region
  double ND_vetoEnergyFD_new; // Add new restriction to vetoEnergy
  // Muon info
  vector<double> ND_Sim_mu_start_vx; // Vector corresponds to randomized/stepwise vtx x, position of the muon trajectory at start point on the x-axis [cm]
  vector<double> ND_Sim_mu_end_vx; // Vector corresponds to randomized/stepwise vtx x, position of the muon trajectory at end point on the x-axis [cm]
  double ND_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double ND_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  double ND_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double ND_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double ND_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double ND_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double ND_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  double ND_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  double ND_Sim_mu_end_px; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double ND_Sim_mu_end_py; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double ND_Sim_mu_end_pz; // Momentum of the muon trajectory at end point on the x-axis [GeV]
  double ND_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  // Hadron info
  double ND_Sim_hadronic_Edep_a2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  // Result of hadron containment is stored in nested vectors, need to generate dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*bool*"); // Remove old dictionary if exists
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*bool*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*bool*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*vector*uint64_t*");
  // Generate new dictionary
  gInterpreter->GenerateDictionary("vector<vector<bool> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<bool> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<bool> > > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<uint64_t> > > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<vector<uint64_t> > > > >", "vector");
  // Nested vector (vetoSize > vetoEnergy) returned from function getHadronContainmentOrigin
  vector<vector<bool> > hadron_contain_result_before_throw;
  vector<vector<vector<bool> > > hadron_contain_result_before_throw_vec_for_vtx_vx; // One more vector: randomized/stepwise evt vtx x
  vector<vector<vector<vector<bool> > > > ND_Sim_hadron_contain_result_before_throw; // ................ ................... ND off axis position x
  // Nested vector (vetoSize > vetoEnergy > 64_bit_throw_result) returned from function getHadronContainmentThrows
  vector<vector<vector<uint64_t> > > hadron_throw_result;
  vector<vector<vector<vector<uint64_t> > > > hadron_throw_result_vec_for_vtx_vx;   // One more vector: randomized/stepwise evt vtx x
  vector<vector<vector<vector<vector<uint64_t> > > > > ND_Sim_hadron_throw_result;   // ................ ................... ND off axis position x
  vector<double> ND_Sim_hadronic_hit_x; // Position of each energy deposit on the x-axis [cm]
  vector<double> ND_Sim_hadronic_hit_y; // Position of each energy deposit on the y-axis [cm]
  vector<double> ND_Sim_hadronic_hit_z; // Position of each energy deposit on the z-axis [cm]
  // //
  // Put event back to the origin (the on-axis position where beam center goes through) [cm]
  //
  vector<double> ND_OffAxis_Sim_mu_start_vx; // Position of the muon trajectory at start point on the x-axis [cm]
  double ND_OffAxis_Sim_mu_start_vy; // Position of the muon trajectory at start point on the y-axis [cm]
  double ND_OffAxis_Sim_mu_start_vz; // Position of the muon trajectory at start point on the z-axis [cm]
  vector<double> ND_OffAxis_Sim_mu_end_vx; // Position of the muon trajectory at end point on the x-axis [cm]
  double ND_OffAxis_Sim_mu_end_vy; // Position of the muon trajectory at end point on the y-axis [cm]
  double ND_OffAxis_Sim_mu_end_vz; // Position of the muon trajectory at end point on the z-axis [cm]
  double ND_OnAxis_Sim_mu_start_px; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  double ND_OnAxis_Sim_mu_start_py; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  double ND_OnAxis_Sim_mu_start_pz; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  vector<float> ND_OffAxis_Sim_hadronic_hit_x; // Position of each energy deposit on the x-axis [cm]
  vector<float> ND_OffAxis_Sim_hadronic_hit_y; // Position of each energy deposit on the y-axis [cm]
  vector<float> ND_OffAxis_Sim_hadronic_hit_z; // Position of each energy deposit on the z-axis [cm]

  // double ND_OnAxis_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  // double ND_OnAxis_Sim_mu_end_E; // Energy of leading mu at end point [GeV]

  //
  // Event after moving from ND(On-Axis) to ND(Off-Axis)
  //
  float ND_OffAxis_Sim_mu_end_vx_rotated; // Position of the muon trajectory at end point on the x-axis [cm]
  float ND_OffAxis_Sim_mu_end_vy_rotated; // Position of the muon trajectory at end point on the y-axis [cm]
  float ND_OffAxis_Sim_mu_end_vz_rotated; // Position of the muon trajectory at end point on the z-axis [cm]
  float ND_OffAxis_Sim_mu_start_px_rotated; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  float ND_OffAxis_Sim_mu_start_py_rotated; // Momentum of the muon trajectory at start point on the y-axis [GeV]
  float ND_OffAxis_Sim_mu_start_pz_rotated; // Momentum of the muon trajectory at start point on the z-axis [GeV]
  vector<float>  ND_OffAxis_Sim_hadronic_hit_x_rotated; // Position of each energy deposit on the x-axis [cm]
  vector<float>  ND_OffAxis_Sim_hadronic_hit_y_rotated; // Position of each energy deposit on the y-axis [cm]
  vector<float>  ND_OffAxis_Sim_hadronic_hit_z_rotated; // Position of each energy deposit on the z-axis [cm]

  // double ND_OffAxis_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  // double ND_OffAxis_Sim_mu_end_E; // Energy of leading mu at end point [GeV]

  // Check the momentum conservation
  float ND_OnAxis_Sim_mu_start_p; // Total momentum of the muon trajectory at start point [GeV]
  float ND_OffAxis_Sim_mu_start_p; // Total momentum of the muon trajectory at start point [GeV]
  // Add output txt file
  ofstream myfile;
   myfile.open ("OutputDataCheck.txt");



  //
  // A tree to store lepton info (for NN training)
  // and result of hadron containment after applying transformations
  //
  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  effTreeFD->Branch("ND_Gen_numu_E",                             &ND_Gen_numu_E,            "ND_Gen_numu_E/D");
  // Add two vetoE branches
  effTreeFD->Branch("ND_vetoEnergyFD",                           &ND_vetoEnergyFD,          "ND_vetoEnergyFD/D");
  effTreeFD->Branch("ND_vetoEnergyFD_new",                       &ND_vetoEnergyFD_new,      "ND_vetoEnergyFD_new/D");
  effTreeFD->Branch("ND_off_axis_pos_vec",                       &ND_off_axis_pos_vec);                             // vector<double>: entries = written evts * ND_off_axis_pos_steps
  effTreeFD->Branch("ND_Sim_mu_start_vx",                        &ND_Sim_mu_start_vx);                              // ............... entries = written evts * vtx_vx_steps, equivalent to ND_vtx_vx_vec
  effTreeFD->Branch("ND_Sim_mu_start_vy",                        &ND_Sim_mu_start_vy,       "ND_Sim_mu_start_vy/D");   // entries = written evts
  effTreeFD->Branch("ND_Sim_mu_start_vz",                        &ND_Sim_mu_start_vz,       "ND_Sim_mu_start_vz/D");
  effTreeFD->Branch("ND_Sim_mu_end_vx",                          &ND_Sim_mu_end_vx);
  effTreeFD->Branch("ND_Sim_mu_end_vy",                          &ND_Sim_mu_end_vy,         "ND_Sim_mu_end_vy/D");
  effTreeFD->Branch("ND_Sim_mu_end_vz",                          &ND_Sim_mu_end_vz,         "ND_Sim_mu_end_vz/D");
  effTreeFD->Branch("ND_Sim_mu_start_px",                        &ND_Sim_mu_start_px,       "ND_Sim_mu_start_px/D");
  effTreeFD->Branch("ND_Sim_mu_start_py",                        &ND_Sim_mu_start_py,       "ND_Sim_mu_start_py/D");
  effTreeFD->Branch("ND_Sim_mu_start_pz",                        &ND_Sim_mu_start_pz,       "ND_Sim_mu_start_pz/D");
  effTreeFD->Branch("ND_Sim_mu_start_E",                         &ND_Sim_mu_start_E,        "ND_Sim_mu_start_E/D");
  effTreeFD->Branch("ND_Sim_mu_end_px",                          &ND_Sim_mu_end_px,         "ND_Sim_mu_end_px/D");
  effTreeFD->Branch("ND_Sim_mu_end_py",                          &ND_Sim_mu_end_py,         "ND_Sim_mu_end_py/D");
  effTreeFD->Branch("ND_Sim_mu_end_pz",                          &ND_Sim_mu_end_pz,         "ND_Sim_mu_end_pz/D");
  effTreeFD->Branch("ND_Sim_mu_end_E",                           &ND_Sim_mu_end_E,          "ND_Sim_mu_end_E/D");
  effTreeFD->Branch("ND_Sim_hadron_contain_result_before_throw", &ND_Sim_hadron_contain_result_before_throw);          // nested vector
  effTreeFD->Branch("ND_Sim_hadron_throw_result",                &ND_Sim_hadron_throw_result);
  effTreeFD->Branch("ND_Sim_hadronic_Edep_a2",                   &ND_Sim_hadronic_Edep_a2,  "ND_Sim_hadronic_Edep_a2/D"); // entries = written evts
  effTreeFD->Branch("FD_Sim_mu_start_vy",                        &FD_Sim_mu_start_vy,       "FD_Sim_mu_start_vy/D");
  effTreeFD->Branch("FD_Sim_mu_start_vz",                        &FD_Sim_mu_start_vz,       "FD_Sim_mu_start_vz/D");
  effTreeFD->Branch("ND_Sim_hadronic_hit_x",                     &ND_Sim_hadronic_hit_x);
  effTreeFD->Branch("ND_Sim_hadronic_hit_y",                     &ND_Sim_hadronic_hit_y);
  effTreeFD->Branch("ND_Sim_hadronic_hit_z",                     &ND_Sim_hadronic_hit_z);
  // On-Axis events
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_vx",                 &ND_OffAxis_Sim_mu_start_vx, "ND_OffAxis_Sim_mu_start_vx/D");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_vy",                 &ND_OffAxis_Sim_mu_start_vy, "ND_OffAxis_Sim_mu_start_vy/D");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_vz",                 &ND_OffAxis_Sim_mu_start_vz, "ND_OffAxis_Sim_mu_start_vz/D");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vx",                   &ND_OffAxis_Sim_mu_end_vx,  "ND_OffAxis_Sim_mu_end_vx/D");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vy",                   &ND_OffAxis_Sim_mu_end_vy,  "ND_OffAxis_Sim_mu_end_vy/D");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vz",                   &ND_OffAxis_Sim_mu_end_vz,  "ND_OffAxis_Sim_mu_end_vz/D");
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_x",              &ND_OffAxis_Sim_hadronic_hit_x);
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_y",              &ND_OffAxis_Sim_hadronic_hit_y);
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_z",              &ND_OffAxis_Sim_hadronic_hit_z);
  // Off-Axis events
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vx_rotated",                       &ND_OffAxis_Sim_mu_end_vx_rotated,      "ND_OffAxis_Sim_mu_end_vx_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vy_rotated",                       &ND_OffAxis_Sim_mu_end_vy_rotated,      "ND_OffAxis_Sim_mu_end_vy_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_vz_rotated",                       &ND_OffAxis_Sim_mu_end_vz_rotated,      "ND_OffAxis_Sim_mu_end_vz_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_px_rotated",                     &ND_OffAxis_Sim_mu_start_px_rotated,    "ND_OffAxis_Sim_mu_start_px_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_py_rotated",                     &ND_OffAxis_Sim_mu_start_py_rotated,    "ND_OffAxis_Sim_mu_start_py_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_pz_rotated",                     &ND_OffAxis_Sim_mu_start_pz_rotated,    "ND_OffAxis_Sim_mu_start_pz_rotated/F");
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_x_rotated",              &ND_OffAxis_Sim_hadronic_hit_x_rotated);
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_y_rotated",              &ND_OffAxis_Sim_hadronic_hit_y_rotated);
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_z_rotated",              &ND_OffAxis_Sim_hadronic_hit_z_rotated);




  //
  // A separate tree to store translations and rotations of throws
  // which will be applied to leptons before NN training
  //

  vector<float> throwVtxY;
  vector<float> throwVtxZ;
  vector<float> throwRot;

  TTree * ThrowsFD = new TTree("ThrowsFD", "FD Throws");
  ThrowsFD->Branch("throwVtxY", &throwVtxY); // vector<float>: entries = [ (int)(written evts / 100) + 1 ] * N_throws
  ThrowsFD->Branch("throwVtxZ", &throwVtxZ);
  ThrowsFD->Branch("throwRot",  &throwRot);

  // Mean neutrino production point (beam coordinate) on z axis as a function of ND off-axis position
  TGraph* gDecayZ = new TGraph(14, OffAxisPoints, meanPDPZ);

  //
  // Get beam parameters: eventually need to read from XML file: which XML for FD? same as ND?
  //

  // Clockwise rotate beamline around ND local x axis
  double beamLineRotation = -0.101;           // unit: rad

  // Put back into beam center(0.0, 0.05387, 6.66) tanslation 1 [cm]
  float ND_OnAxis_Sim_mu_start_vx = 0.0;
  float ND_OnAxis_Sim_mu_start_vy = 0.05387*100;
  float ND_OnAxis_Sim_mu_start_vz = 6.6*100;

  // Coordinate transformation, units: meters
  double beamRefDetCoord[3] = {ND_OnAxis_Sim_mu_start_vx/100., ND_OnAxis_Sim_mu_start_vy/100., ND_OnAxis_Sim_mu_start_vz/100.}; // (0, 0, 0) is ND detector origin
  double detRefBeamCoord[3] = {0., 0., 574.}; // (0, 0, 0) is beam origin
  // Calculate neutrino production x in detector coordinate, y/z later as they depend on ND off-axis position
  decayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0];

  //
  // Initialize geometric efficiency module
  //

  geoEff * eff = new geoEff(314, false); // set verbose to true for debug
  eff->setNthrows(N_throws);
  // Rotate w.r.t. neutrino direction, rather than fixed beam direction
  eff->setUseFixedBeamDir(false);

  // 30 cm veto
  eff->setVetoSizes(vector<float>(1, 30.));
  // 30 MeV
  eff->setVetoEnergyThresholds(vector<float>(1, 30.));

  // Active detector dimensions for ND
  eff->setActiveX(NDActiveVol_min[0], NDActiveVol_max[0]);
  eff->setActiveY(NDActiveVol_min[1], NDActiveVol_max[1]);
  eff->setActiveZ(NDActiveVol_min[2], NDActiveVol_max[2]);

  // Range for translation throws. Use full active volume but fix X.
  eff->setRangeX(-1, -1);
  eff->setRandomizeX(false);
  eff->setRangeY(NDActiveVol_min[1], NDActiveVol_max[1]);
  eff->setRangeZ(NDActiveVol_min[2], NDActiveVol_max[2]);

  // Set offset between MC coordinate system and det volumes
  eff->setOffsetX(offset[0]);
  eff->setOffsetY(offset[1]);
  eff->setOffsetZ(offset[2]);

  // Set On-Axis vertex where center beam cross
  eff->setOnAxisVertex(ND_OnAxis_Sim_mu_start_vx,ND_OnAxis_Sim_mu_start_vy,ND_OnAxis_Sim_mu_start_vz);

  // Add hist of veto E
  TH1F *hist_vetoEnergyFD = new TH1F("hist_vetoEnergyFD", "hist_vetoEnergyFD", 1500, 0, 1500);
  TH1F *hist_vetoEnergyFD_new = new TH1F("hist_vetoEnergyFD_new", "hist_vetoEnergyFD_new", 1500, 0, 1500);

  //
  // Loop over FD events
  //

  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;
  myfile << "Tot evts: " << nentries << "\n";
  for ( int ientry = 0; ientry < 3; ientry++ ) {

    t->GetEntry(ientry);
    if ( ientry%10000 == 0 )
    {
      std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;
      myfile << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << "\n\n\n";
    }
    myfile << "ientry: " << ientry <<"\n";


    //
    // Skip events without muon/hadronic deposits
    //

    if ( FD_Sim_nMu == 0 || FD_Sim_n_hadronic_Edep_a == 0 ) continue;
    myfile << "FD_Sim_n_hadronic_Edep_a: " << FD_Sim_n_hadronic_Edep_a <<"\n";
    if ( FD_CCNC_truth == 1) continue;   // only use CC events
    if ( abs(FD_neuPDG) != 14 ) continue;       // only use muon neu

    //
    // Calculate total hadron E in FD veto region
    //
    // Old veto region restriction
    vetoEnergyFD = 0.;
    // Loop over hadron E deposits
    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){

      // Veto region size: 30 cm from the active volume
      if ( ( FD_Sim_hadronic_hit_x_a->at(ihadronhit) > FDActiveVol_min[0] && FD_Sim_hadronic_hit_x_a->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
           ( FD_Sim_hadronic_hit_y_a->at(ihadronhit) > FDActiveVol_min[1] && FD_Sim_hadronic_hit_y_a->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
           ( FD_Sim_hadronic_hit_z_a->at(ihadronhit) > FDActiveVol_min[2] && FD_Sim_hadronic_hit_z_a->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
           ( FD_Sim_hadronic_hit_x_a->at(ihadronhit) > FDActiveVol_max[0] - 30 && FD_Sim_hadronic_hit_x_a->at(ihadronhit) < FDActiveVol_max[0] ) ||
           ( FD_Sim_hadronic_hit_y_a->at(ihadronhit) > FDActiveVol_max[1] - 30 && FD_Sim_hadronic_hit_y_a->at(ihadronhit) < FDActiveVol_max[1] ) ||
           ( FD_Sim_hadronic_hit_z_a->at(ihadronhit) > FDActiveVol_max[2] - 30 && FD_Sim_hadronic_hit_z_a->at(ihadronhit) < FDActiveVol_max[2] )
         ){
           vetoEnergyFD += FD_Sim_hadronic_hit_Edep_a2->at(ihadronhit);
      } // end if hadron deposit in FD veto region

    } // end loop over hadron E deposits
    ND_vetoEnergyFD = vetoEnergyFD;
    //add a vetoEnergyFD histogram in the root file

    hist_vetoEnergyFD->Fill(vetoEnergyFD);



    // New veto region restriction
    vetoEnergyFD_new = 0.;
    // Loop over hadron E deposits
    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){

      // Veto region size: 30 cm from the active volume
      if ( ( FD_Sim_hadronic_hit_x_a->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
           ( FD_Sim_hadronic_hit_y_a->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
           ( FD_Sim_hadronic_hit_z_a->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
           ( FD_Sim_hadronic_hit_x_a->at(ihadronhit) > FDActiveVol_max[0] - 30 ) ||
           ( FD_Sim_hadronic_hit_y_a->at(ihadronhit) > FDActiveVol_max[1] - 30 ) ||
           ( FD_Sim_hadronic_hit_z_a->at(ihadronhit) > FDActiveVol_max[2] - 30 )
         ){
           vetoEnergyFD_new += FD_Sim_hadronic_hit_Edep_a2->at(ihadronhit);
      } // end if hadron deposit in FD veto region

    } // end loop over hadron E deposits
    ND_vetoEnergyFD_new = vetoEnergyFD_new;
    //add a vetoEnergyFD histogram in the root file
    hist_vetoEnergyFD_new->Fill(vetoEnergyFD_new);

    /*// Add vetoEnergyFD ntuple before threshold
    ND_vetoEnergyFD = vetoEnergyFD;
    vetoE->Fill();*/

    //
    // Skip FD event if the total hadron E in veto region exceeds vetoEnergy [MeV]
    //
    if ( vetoEnergyFD > 30 ) continue; // 30 MeV

    // Renew throws every 100th written event to save file size, i.e., if N = 128,
    // for written evt 0-99:    same 128 transformations for each event,
    // for written evt 100-199: same but renewed 128 transformations for each evt
    // so on so forth...
    // These transformations will be applied to leptons in the event, so need to keep track of iwritten
    if ( iwritten % 100 == 0 ) {

      // Produce N throws defined at setNthrows(N)
      // Same throws applied for hadron below
      eff->throwTransforms(); // Does not depend on evt vtx
      throwVtxY.clear();
      throwVtxZ.clear();
      throwRot.clear();
      throwVtxY = eff->getCurrentThrowTranslationsY();
      throwVtxZ = eff->getCurrentThrowTranslationsZ();
      throwRot  = eff->getCurrentThrowRotations();
      ThrowsFD->Fill();
    }

    // Local y-z axes in FD and ND are rotated due to Earth curvature, x direction is not change
    // FD event coordinates, if unchanged, would represent different event in ND coordinate sys.
    // Apply an active transformation matrix R_x(theta): rotate each point counterclockwise by theta around x-axis in ND
    // theta is 2*|beamLineRotation|
    // Transform FD relative coordinate, x coordinate unchanged
    //
    //              [ 1          0             0
    // R_x(theta) =   0      cos(theta)   -sin(theta)
    //                0      sin(theta)    cos(theta) ]
    //

    // Rotation affects mu start/end position and hadron positions below
    ND_Sim_mu_start_vx.emplace_back(FD_Sim_mu_start_vx);
    ND_Sim_mu_start_vy = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vy - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vz;
    ND_Sim_mu_start_vz = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vy + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vz;

    ND_Sim_mu_end_vx.emplace_back(FD_Sim_mu_end_vx);
    ND_Sim_mu_end_vy = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vy - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vz;
    ND_Sim_mu_end_vz = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vy + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vz;

    // Initialize for the event
    ND_Sim_hadron_contain_result_before_throw.clear();
    ND_Sim_hadron_throw_result.clear();

    // Branches that are not affected by ND off axis position and vtx x (loops below)
    ND_Gen_numu_E      = FD_Gen_numu_E;

    // X momentum is not affected by coordinate rotation
    ND_Sim_mu_start_px = FD_Sim_mu_start_px;
    ND_Sim_mu_end_px   = FD_Sim_mu_end_px;

    // Rotation affects mu start/end momentum vector in Y and Z axes, using the same rotation matrix below.
    ND_Sim_mu_start_py = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_py - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_pz;
    ND_Sim_mu_start_pz = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_py + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_pz;
    ND_Sim_mu_end_py   = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_py - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_pz;
    ND_Sim_mu_end_pz   = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_py + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_pz;

    // Energy is a scalar, doesn't change
    ND_Sim_mu_start_E  = FD_Sim_mu_start_E;
    ND_Sim_mu_end_E    = FD_Sim_mu_end_E;
    ND_Sim_hadronic_Edep_a2 = FD_Sim_hadronic_Edep_a2;

    // Hardon hit positions
    HadronHitEdeps.clear();
    HadronHitPoss.clear();
    HadronHitEdeps.reserve(FD_Sim_n_hadronic_Edep_a);
    HadronHitPoss.reserve(FD_Sim_n_hadronic_Edep_a*3);

    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
      ND_Sim_hadronic_hit_x.emplace_back(FD_Sim_hadronic_hit_x_a->at(ihadronhit));
      cout<< "FD hadron hit_x" << FD_Sim_hadronic_hit_x_a->at(ihadronhit)<<endl;
      cout<< "ND hadron hit_x" <<ND_Sim_hadronic_hit_x.at(ihadronhit)<<endl;
      ND_Sim_hadronic_hit_y.emplace_back(cos( 2*abs(beamLineRotation) )*FD_Sim_hadronic_hit_y_a->at(ihadronhit) - sin( 2*abs(beamLineRotation) )*FD_Sim_hadronic_hit_z_a->at(ihadronhit));
      ND_Sim_hadronic_hit_z.emplace_back(sin( 2*abs(beamLineRotation) )*FD_Sim_hadronic_hit_y_a->at(ihadronhit) + cos( 2*abs(beamLineRotation) )*FD_Sim_hadronic_hit_z_a->at(ihadronhit));
    }

    //
    // Find the new coordinates of each vector after putting events back to the OnAxis_vector
    //
    // ND_Sim_mu_start_vector
    ND_OffAxis_Sim_mu_start_vx.clear();
    ND_OffAxis_Sim_mu_start_vy = ND_OnAxis_Sim_mu_start_vy;
    ND_OffAxis_Sim_mu_start_vz = ND_OnAxis_Sim_mu_start_vz;
    // ND_Sim_mu_end_vector
    ND_OffAxis_Sim_mu_end_vx.clear();
    ND_OffAxis_Sim_mu_end_vy = ND_OffAxis_Sim_mu_start_vy - ( ND_Sim_mu_start_vy - ND_Sim_mu_end_vy);
    ND_OffAxis_Sim_mu_end_vz = ND_OffAxis_Sim_mu_start_vz - ( ND_Sim_mu_start_vz - ND_Sim_mu_end_vz);
    // ND_OnAxis_Sim_mu_start_momentum
    ND_OnAxis_Sim_mu_start_px = ND_Sim_mu_start_px;
    ND_OnAxis_Sim_mu_start_py = ND_Sim_mu_start_py;
    ND_OnAxis_Sim_mu_start_pz = ND_Sim_mu_start_pz;
    // Hadron hit energy is a scale, doesn't change
    // ND_OnAxis_Sim_mu_start_E = ND_Sim_mu_start_E;
    // ND_OnAxis_Sim_mu_end_E   = ND_Sim_mu_end_E;
    // Hardon hit positions
    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ )
    {
      ND_OffAxis_Sim_hadronic_hit_x.emplace_back(ND_Sim_hadronic_hit_x.at(ihadronhit) - ( FD_Sim_mu_start_vx - ND_OnAxis_Sim_mu_start_vx ));
      HadronHitPoss.emplace_back(ND_OffAxis_Sim_hadronic_hit_x.at(ihadronhit));
      ND_OffAxis_Sim_hadronic_hit_y.emplace_back(ND_Sim_hadronic_hit_y.at(ihadronhit) - ( ND_Sim_mu_start_vy - ND_OffAxis_Sim_mu_start_vy ));
      HadronHitPoss.emplace_back(ND_OffAxis_Sim_hadronic_hit_y.at(ihadronhit));
      ND_OffAxis_Sim_hadronic_hit_z.emplace_back(ND_Sim_hadronic_hit_z.at(ihadronhit) - ( ND_Sim_mu_start_vz - ND_OffAxis_Sim_mu_start_vz ));
      HadronHitPoss.emplace_back(ND_OffAxis_Sim_hadronic_hit_z.at(ihadronhit));
      HadronHitEdeps.emplace_back( FD_Sim_hadronic_hit_Edep_a2->at(ihadronhit) );
    }



    //
    // Loop over ND_off_axis_pos_vec: random off_axis_pos or every ND_off_axis_pos_stepsize
    // Don't put it outside event loop to avoid looping over all events multiple times
    //

    int ND_off_axis_pos_counter = 0;
    for ( double i_ND_off_axis_pos : ND_off_axis_pos_vec ) {

      ND_off_axis_pos_counter++;
      //
      // Loop over vtx x: random x or stepwise increased x
      // Don't put it outside event loop to avoid looping over all events multiple times
      //
      int vtx_vx_counter = 0;
      hadron_throw_result_vec_for_vtx_vx.clear(); // Need to initialize before loop over vtx x vec
      hadron_contain_result_before_throw_vec_for_vtx_vx.clear();

      for ( double i_vtx_vx : ND_vtx_vx_vec ) {

        vtx_vx_counter++;

        // Interpolate event neutrino production point (beam coordinate)
        decayZbeamCoord = gDecayZ->Eval( i_ND_off_axis_pos + i_vtx_vx - detRefBeamCoord[0] );

        // Calculate neutrino production point in detector coordinate
        decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation);
        decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation);
        // Set production point in unit: cm
        eff->setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.);



        //if (verbose) std::cout << "nd off_axis x #" << ND_off_axis_pos_counter << ": " << i_ND_off_axis_pos << " m" << std::endl;
        // std::cout << "Nd Off-Axis x #" << ND_off_axis_pos_counter << ": " << i_ND_off_axis_pos << " cm" << std::endl;
        // std::cout << "Nd Off-Axis x in detector #" << vtx_vx_counter << ": " << i_vtx_vx << " cm" << std::endl;
        // std::cout << "Event position #" << vtx_vx_counter << ": " << i_ND_off_axis_pos + i_vtx_vx << " cm" << std::endl;
       myfile << "Nd Off-Axis x #" << ND_off_axis_pos_counter << "[cm]: " << i_ND_off_axis_pos << "\n";
       myfile << "Nd Off-Axis x in detector #" << vtx_vx_counter << "[cm]: " << i_vtx_vx << "\n";
       myfile << "Event position #" << vtx_vx_counter << "[cm]: " << i_ND_off_axis_pos + i_vtx_vx << "\n\n";


        // ND off-axis position does not affect evt vx, so only fill branches below once when loop over ND off-axis vec
        // if ( ND_off_axis_pos_counter == 1 ) {
          ND_OffAxis_Sim_mu_start_vx.emplace_back( i_ND_off_axis_pos + i_vtx_vx );
          ND_OffAxis_Sim_mu_end_vx.emplace_back( FD_Sim_mu_end_vx - FD_Sim_mu_start_vx + i_ND_off_axis_pos + i_vtx_vx ); // w.r.t. mu start x
        // }
        // Mu_start_vy and vx do not change
        ND_OffAxis_Sim_mu_start_vy = ND_OffAxis_Sim_mu_start_vy;
        ND_OffAxis_Sim_mu_start_vz = ND_OffAxis_Sim_mu_start_vz;

        // Evt vtx pos in unit: cm
        eff->setVertex( i_ND_off_axis_pos + i_vtx_vx, ND_OffAxis_Sim_mu_start_vy, ND_OffAxis_Sim_mu_start_vz);
        eff->setOffAxisVertex(i_ND_off_axis_pos + i_vtx_vx, ND_OffAxis_Sim_mu_start_vy, ND_OffAxis_Sim_mu_start_vz);

        // Find ND_Sim_mu_end_v_af
        eff->setMuEndV( ND_OnAxis_Sim_mu_start_vx - ( FD_Sim_mu_start_vx - FD_Sim_mu_end_vx ) + i_ND_off_axis_pos + i_vtx_vx, ND_OffAxis_Sim_mu_end_vy, ND_OffAxis_Sim_mu_end_vz);
        if (i_ND_off_axis_pos+i_vtx_vx==0)
        {
          ND_OffAxis_Sim_mu_end_vx_rotated = ND_OnAxis_Sim_mu_start_vx - ( FD_Sim_mu_start_vx - FD_Sim_mu_end_vx ) + i_ND_off_axis_pos + i_vtx_vx;
          ND_OffAxis_Sim_mu_end_vy_rotated = ND_OffAxis_Sim_mu_end_vy;
          ND_OffAxis_Sim_mu_end_vz_rotated = ND_OffAxis_Sim_mu_end_vz;
        }
        else{
          ND_OffAxis_Sim_mu_end_vx_rotated = eff->getOffAxisMuEndV(0);
          ND_OffAxis_Sim_mu_end_vy_rotated = eff->getOffAxisMuEndV(1);
          ND_OffAxis_Sim_mu_end_vz_rotated = eff->getOffAxisMuEndV(2);}
        myfile << "ND_OffAxis_Sim_mu_start_vx[cm]: " << i_ND_off_axis_pos + i_vtx_vx << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_vy[cm]: " << ND_OffAxis_Sim_mu_start_vy << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_vz[cm]: " << ND_OffAxis_Sim_mu_start_vz << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vx[cm]: " << ND_OnAxis_Sim_mu_start_vx - ( FD_Sim_mu_start_vx - FD_Sim_mu_end_vx ) + i_ND_off_axis_pos + i_vtx_vx << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vy[cm]: " << ND_OffAxis_Sim_mu_end_vy << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vz[cm]: " << ND_OffAxis_Sim_mu_end_vz << "\n";
        myfile << "Distance between OffAxis_Sim_mu_v end-start[cm]: " << sqrt(pow(ND_OnAxis_Sim_mu_start_vx - ( FD_Sim_mu_start_vx - FD_Sim_mu_end_vx ),2)+pow(ND_OffAxis_Sim_mu_end_vy -  ND_OffAxis_Sim_mu_start_vy,2)+pow(ND_OffAxis_Sim_mu_end_vz -  ND_OffAxis_Sim_mu_start_vz,2)) << "\n\n";
        myfile << "ND_OffAxis_Sim_mu_start_vx_w/rotation[cm]: " << i_ND_off_axis_pos + i_vtx_vx << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_vy_w/rotation[cm]: " << ND_OffAxis_Sim_mu_start_vy << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_vz_w/rotation[cm]: " << ND_OffAxis_Sim_mu_start_vz << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vx_w/rotation[cm]: " << ND_OffAxis_Sim_mu_end_vx_rotated << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vy_w/rotation[cm]: " << ND_OffAxis_Sim_mu_end_vy_rotated << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_vz_w/rotation[cm]: " << ND_OffAxis_Sim_mu_end_vz_rotated << "\n";
        myfile << "Distance between OffAxis_Sim_mu_v end-start_w/rotation[cm]: " << sqrt(pow(ND_OffAxis_Sim_mu_end_vx_rotated - (i_ND_off_axis_pos + i_vtx_vx),2)+pow(ND_OffAxis_Sim_mu_end_vy_rotated -  ND_OffAxis_Sim_mu_start_vy,2)+pow(ND_OffAxis_Sim_mu_end_vz_rotated -  ND_OffAxis_Sim_mu_start_vz,2)) << "\n\n";



        // Find ND_Sim_mu_start_momentum
        ND_OnAxis_Sim_mu_start_p = sqrt(pow(ND_Sim_mu_start_px,2)+pow(ND_Sim_mu_start_py,2)+pow(ND_Sim_mu_start_pz,2));
        myfile << "ND_OffAxis_Sim_mu_start_px[GeV]: " << ND_Sim_mu_start_px << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_py[GeV]: " << ND_Sim_mu_start_py << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_pz[GeV]: " << ND_Sim_mu_start_pz << "\n";
        myfile << "Total ND_OnAxis_Sim_mu_start_p[GeV]: " << ND_OnAxis_Sim_mu_start_p << "\n";

        // Find ND_Sim_mu_end_p_af
        eff->setMuStartP(ND_OnAxis_Sim_mu_start_px, ND_OnAxis_Sim_mu_start_py, ND_OnAxis_Sim_mu_start_pz);
        if (i_ND_off_axis_pos+i_vtx_vx==0)
        {ND_OffAxis_Sim_mu_start_px_rotated = ND_OnAxis_Sim_mu_start_px;
        ND_OffAxis_Sim_mu_start_py_rotated = ND_OnAxis_Sim_mu_start_py;
        ND_OffAxis_Sim_mu_start_pz_rotated = ND_OnAxis_Sim_mu_start_pz;}
        else{
        ND_OffAxis_Sim_mu_start_px_rotated = eff->getOffAxisMuStartP(0);
        ND_OffAxis_Sim_mu_start_py_rotated = eff->getOffAxisMuStartP(1);
        ND_OffAxis_Sim_mu_start_pz_rotated = eff->getOffAxisMuStartP(2);}
        ND_OffAxis_Sim_mu_start_p = sqrt(pow(ND_OffAxis_Sim_mu_start_px_rotated,2)+pow(ND_OffAxis_Sim_mu_start_py_rotated,2)+pow(ND_OffAxis_Sim_mu_start_pz_rotated,2));
        myfile << "ND_OffAxis_Sim_mu_start_px_w/rotation[GeV]: " << ND_OffAxis_Sim_mu_start_px_rotated << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_py_w/rotation[GeV]: " << ND_OffAxis_Sim_mu_start_py_rotated << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_pz_w/rotation[GeV]: " << ND_OffAxis_Sim_mu_start_pz_rotated << "\n";
        myfile << "Total ND_OffAxis_Sim_mu_start_p_w/rotation[GeV]: " << ND_OffAxis_Sim_mu_start_p << "\n\n\n";



        // Find ND_Sim_mu_hadronic_hit
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ )
        {

          eff->setHadronHitV(ND_OffAxis_Sim_hadronic_hit_x.at(ihadronhit),ND_OffAxis_Sim_hadronic_hit_y.at(ihadronhit),ND_OffAxis_Sim_hadronic_hit_z.at(ihadronhit));
          if (i_ND_off_axis_pos+i_vtx_vx==0)
          {
            ND_OffAxis_Sim_hadronic_hit_x_rotated.emplace_back(ND_OffAxis_Sim_hadronic_hit_x.at(ihadronhit));
            ND_OffAxis_Sim_hadronic_hit_y_rotated.emplace_back(ND_OffAxis_Sim_hadronic_hit_y.at(ihadronhit));
            ND_OffAxis_Sim_hadronic_hit_z_rotated.emplace_back(ND_OffAxis_Sim_hadronic_hit_z.at(ihadronhit));
          }
          else{
            ND_OffAxis_Sim_hadronic_hit_x_rotated.emplace_back(eff->getOffAxisHadronHitV(0));
            ND_OffAxis_Sim_hadronic_hit_y_rotated.emplace_back(eff->getOffAxisHadronHitV(1));
            ND_OffAxis_Sim_hadronic_hit_z_rotated.emplace_back(eff->getOffAxisHadronHitV(2));
          }
          myfile << "Hit #" << ihadronhit <<"/" << FD_Sim_n_hadronic_Edep_a << "\n";
          myfile << "Distance between OffAxis_Sim_hadronic_hit-OffAxisVertex[cm]: " << sqrt(pow(ND_OffAxis_Sim_hadronic_hit_x.at(ihadronhit)-(i_ND_off_axis_pos + i_vtx_vx),2)+pow(ND_OffAxis_Sim_hadronic_hit_y.at(ihadronhit)-ND_OffAxis_Sim_mu_start_vy,2)+pow(ND_OffAxis_Sim_hadronic_hit_z.at(ihadronhit)-ND_OffAxis_Sim_mu_start_vz,2)) << "\n";
          myfile << "Distance between OffAxis_Sim_hadronic_hit-OffAxisVertex_w/rotation[cm]: " << sqrt(pow(ND_OffAxis_Sim_hadronic_hit_x_rotated.at(ihadronhit)-(i_ND_off_axis_pos + i_vtx_vx),2)+pow(ND_OffAxis_Sim_hadronic_hit_y_rotated.at(ihadronhit)-ND_OffAxis_Sim_mu_start_vy,2)+pow(ND_OffAxis_Sim_hadronic_hit_z_rotated.at(ihadronhit)-ND_OffAxis_Sim_mu_start_vz,2)) << "\n\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_x[cm]: " << ND_OffAxis_Sim_hadronic_hit_x.at(ihadronhit) << "\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_y[cm]: " << ND_OffAxis_Sim_hadronic_hit_y.at(ihadronhit) << "\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_z[cm]: " << ND_OffAxis_Sim_hadronic_hit_z.at(ihadronhit) << "\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_x_rotated[cm]: " << ND_OffAxis_Sim_hadronic_hit_x_rotated.at(ihadronhit) << "\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_y_rotated[cm]: " << ND_OffAxis_Sim_hadronic_hit_y_rotated.at(ihadronhit) << "\n";
        myfile << "ND_OffAxis_Sim_hadronic_hit_z_rotated[cm]: " << ND_OffAxis_Sim_hadronic_hit_z_rotated.at(ihadronhit) << "\n\n";
        }





        // // Energy doesn't change
        // ND_OffAxis_Sim_mu_start_E = ND_OnAxis_Sim_mu_start_E;
        // ND_OffAxis_Sim_mu_end_E   = ND_OnAxis_Sim_mu_end_E;


        // Need to loop deposits many times as the hadron containment below need vtx x first
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){

          // Relative to muon start pos in ND coordinate sys: (i_vtx_vx, ND_Sim_mu_start_vy, ND_Sim_mu_start_vz)
          HadronHitPoss.emplace_back( ND_OffAxis_Sim_hadronic_hit_x_rotated.at(ihadronhit)); // w.r.t. mu start x
          // Again, need to apply R_x(theta) for hadron y/z, do not affect x
          HadronHitPoss.emplace_back( ND_OffAxis_Sim_hadronic_hit_y_rotated.at(ihadronhit) );
          HadronHitPoss.emplace_back( ND_OffAxis_Sim_hadronic_hit_z_rotated.at(ihadronhit) );
          HadronHitEdeps.emplace_back( FD_Sim_hadronic_hit_Edep_a2->at(ihadronhit) );
        } // end for loop

        eff->setHitSegEdeps(HadronHitEdeps);
        eff->setHitSegPoss(HadronHitPoss); // this is converted hadrom deposit pos in ND coordinate sys.

        //
        // Get hadron containment result after everything is set to ND coordinate sys
        //

        // Store hadron containment result for the FD evt before doing random throws
        hadron_contain_result_before_throw = eff->getHadronContainmentOrigin();
        hadron_contain_result_before_throw_vec_for_vtx_vx.emplace_back(hadron_contain_result_before_throw);

        // Do random throws regardless whether FD evt is contained in ND volume by setting a false flag
        hadron_throw_result = eff->getHadronContainmentThrows(false); // Every 64 throw results stored into a 64 bit unsigned int: 0101101...
        hadron_throw_result_vec_for_vtx_vx.emplace_back(hadron_throw_result);

        if (verbose) std::cout << "vtx x #" << vtx_vx_counter << ": " << i_vtx_vx << " cm, throw result[0][0][0]: " << hadron_throw_result[0][0][0] << std::endl;
        effTreeFD->Fill();
        ND_OffAxis_Sim_hadronic_hit_x_rotated.clear();
        ND_OffAxis_Sim_hadronic_hit_y_rotated.clear();
        ND_OffAxis_Sim_hadronic_hit_z_rotated.clear();

      } // end Loop over ND_vtx_vx_vec


      ND_Sim_hadron_contain_result_before_throw.emplace_back(hadron_contain_result_before_throw_vec_for_vtx_vx);
      ND_Sim_hadron_throw_result.emplace_back(hadron_throw_result_vec_for_vtx_vx);

    }   // end Loop over ND_off_axis_pos_vec


    iwritten++;

    ND_Sim_hadronic_hit_x.clear();
    ND_Sim_hadronic_hit_y.clear();
    ND_Sim_hadronic_hit_z.clear();
    ND_OffAxis_Sim_hadronic_hit_x.clear();
    ND_OffAxis_Sim_hadronic_hit_y.clear();
    ND_OffAxis_Sim_hadronic_hit_z.clear();



  }     // end loop over events

  std::cout << "Written evts: " << iwritten << std::endl;
  myfile << "Written evts: " << iwritten << "\n";
  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff_ivy.root", "RECREATE");
  ThrowsFD->Write();
  effTreeFD->Write();
  hist_vetoEnergyFD->Write();
  hist_vetoEnergyFD_new->Write();

  myfile.close();
  outFile->Close();

} // end main
