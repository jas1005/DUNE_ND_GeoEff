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

  int FD_Run;
  int FD_SubRun;
  int FD_Event;
  int FD_Sim_nNumu;
  double FD_Gen_numu_E;
  int FD_Sim_nMu;
  int FD_CCNC_truth; // Choose CCNC
  int FD_neuPDG; // neutrino PDG
  double FD_Sim_mu_start_vx; // unit: cm?
  double FD_Sim_mu_start_vy;
  double FD_Sim_mu_start_vz;
  double FD_Sim_mu_end_vx;
  double FD_Sim_mu_end_vy;
  double FD_Sim_mu_end_vz;
  double FD_Sim_mu_start_px;
  double FD_Sim_mu_start_py;
  double FD_Sim_mu_start_pz;
  double FD_Sim_mu_start_E;
  double FD_Sim_mu_end_px;
  double FD_Sim_mu_end_py;
  double FD_Sim_mu_end_pz;
  double FD_Sim_mu_end_E;
  double FD_Sim_hadronic_Edep_a2;
  int FD_Sim_n_hadronic_Edep_a;
  vector<float> *FD_Sim_hadronic_hit_Edep_a2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_a     = 0; // Same as mu pos unit: cm?
  vector<float> *FD_Sim_hadronic_hit_y_a     = 0;
  vector<float> *FD_Sim_hadronic_hit_z_a     = 0;

  // Read ntuple from FD MC
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // ntuple path on FNAL dunegpvm machine
  // For Ivy machine:
  t->Add("/home/fyguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root");
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

  int nentries = 0;        // Total input events
  float vetoEnergyFD;      // Total hadron deposited energy in FD veto region
  float vetoEnergyFD_new;  // Add vetoEnergy for new restriction
  int iwritten = 0;        // Output event counter
  float decayZbeamCoord;
  float decayXdetCoord;
  float decayYdetCoord;
  float decayZdetCoord;
  vector<float> HadronHitEdeps;
  vector<float> HadronHitPoss;

  vector<double> ND_off_axis_pos_vec = {0,7,30}; // unit: cm, ND off-axis choices for each FD evt: 1st element is randomized for each evt
  vector<double> ND_vtx_vx_vec={-0.02,0,0.02};          // unit: m, vtx x choices for each FD evt in ND volume: 1st element is randomized for each evt
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
  // Lepton info: expressed in ND coordinate sys, do not confuse with branches read above in FD coordinate sys
  double ND_Gen_numu_E;
  // Add veto E with two different restrictions
  double ND_vetoEnergyFD;
  double ND_vetoEnergyFD_new;

  vector<double> ND_Sim_mu_start_vx; // Vector corresponds to randomized/stepwise vtx x
  vector<double> ND_Sim_mu_end_vx;
  double ND_Sim_mu_start_vy;         // Do not use float!
  double ND_Sim_mu_start_vz;
  double ND_Sim_mu_end_vy;
  double ND_Sim_mu_end_vz;
  double ND_Sim_mu_start_px;
  double ND_Sim_mu_start_py;
  double ND_Sim_mu_start_pz;
  double ND_Sim_mu_start_E;
  double ND_Sim_mu_end_px;
  double ND_Sim_mu_end_py;
  double ND_Sim_mu_end_pz;
  double ND_Sim_mu_end_E;

  // Tot hadron E dep
  double ND_Sim_hadronic_Edep_a2;
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
  //
  vector<float> ND_Sim_mu_end_vx_af;
  vector<float> ND_Sim_mu_end_vy_af;
  vector<float> ND_Sim_mu_end_vz_af;
  vector<float> ND_Sim_mu_start_px_af;
  vector<float> ND_Sim_mu_start_py_af;
  vector<float> ND_Sim_mu_start_pz_af;
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
  effTreeFD->Branch("ND_Sim_mu_end_vx_af",                       &ND_Sim_mu_end_vx_af,      "ND_Sim_mu_end_vx_af/F");
  effTreeFD->Branch("ND_Sim_mu_end_vy_af",                       &ND_Sim_mu_end_vy_af,      "ND_Sim_mu_end_vy_af/F");
  effTreeFD->Branch("ND_Sim_mu_end_vz_af",                       &ND_Sim_mu_end_vz_af,      "ND_Sim_mu_end_vz_af/F");
  effTreeFD->Branch("ND_Sim_mu_start_px_af",                     &ND_Sim_mu_start_px_af,    "ND_Sim_mu_start_px_af/F");
  effTreeFD->Branch("ND_Sim_mu_start_py_af",                     &ND_Sim_mu_start_py_af,    "ND_Sim_mu_start_py_af/F");
  effTreeFD->Branch("ND_Sim_mu_start_pz_af",                     &ND_Sim_mu_start_pz_af,    "ND_Sim_mu_start_pz_af/F");
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
  // Coordinate transformation, units: meters
  double beamRefDetCoord[3] = {0.0, 0.05387, 6.66}; // (0, 0, 0) is ND detector origin
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

  // Add hist of veto E
  TH1F *hist_vetoEnergyFD = new TH1F("hist_vetoEnergyFD", "hist_vetoEnergyFD", 1500, 0, 1500);
  TH1F *hist_vetoEnergyFD_new = new TH1F("hist_vetoEnergyFD_new", "hist_vetoEnergyFD_new", 1500, 0, 1500);



  //
  // Loop over FD events
  //

  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < nentries; ientry++ ) {

    t->GetEntry(ientry);
    if ( ientry%10000 == 0 ) std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;

    //
    // Skip events without muon/hadronic deposits
    //

    if ( FD_Sim_nMu == 0 || FD_Sim_n_hadronic_Edep_a == 0 ) continue;
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
    ND_Sim_mu_start_vy = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vy - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vz;
    ND_Sim_mu_start_vz = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vy + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_start_vz;

    ND_Sim_mu_end_vy = cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vy - sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vz;
    ND_Sim_mu_end_vz = sin( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vy + cos( 2*abs(beamLineRotation) )*FD_Sim_mu_end_vz;

    // Initialize for the event
    ND_Sim_mu_start_vx.clear();
    ND_Sim_mu_end_vx.clear();
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

    // Put back into beam center(0.0, 0.05387, 6.66) tanslation 1
    float ND_Sim_mu_start_OnAxis_vx = 0.0;
    float ND_Sim_mu_start_OnAxis_vy = 0.05387*100;
    float ND_Sim_mu_start_OnAxis_vz = 6.6*100;





    eff->setOnAxisVertex(ND_Sim_mu_start_OnAxis_vx,ND_Sim_mu_start_OnAxis_vy,ND_Sim_mu_start_OnAxis_vz);
    eff->setMuStartP(ND_Sim_mu_start_px,ND_Sim_mu_start_py,ND_Sim_mu_start_pz);

    //
    // Two options for setting ND off-axis position
    //
    // Option 1: randomize ND off-axis position once for each event
    // Option 2: stepwise increments along ND detector hall off axis range
    //
    // If only want option 1, set random_ND_off_axis_pos to true in Helpers.h; default is false (use both options)
    //

    //Initialize random number generator
    //This needs to be inside the event loop to make sure each event has a different random number
    TRandom3 *r3_OffAxisPoint = new TRandom3();
    // Set the seed (required to avoid repeated random numbers in each sequence)
    r3_OffAxisPoint->SetSeed(0);
    ND_off_axis_pos_vec.at(0) = r3_OffAxisPoint->Uniform(OffAxisPoints[0], OffAxisPoints[13]);

    if (verbose) std::cout << "random OffAxisPoint [meters]: " << ND_off_axis_pos_vec.at(0) << std::endl;

    //
    // Similarly, two options for setting event vtx x position
    //
    // Option 1: randomize x once for each event
    // Option 2: stepwise increments along ND local width in x: -2m to 2m
    //
    // If only want option 1, set random_vtx_vx to true in Helpers.h; default is false (use both options)
    //

    TRandom3 *r3_vtx_x = new TRandom3();
    r3_vtx_x->SetSeed(0);
    ND_vtx_vx_vec.at(0) = r3_vtx_x->Uniform(ND_local_x_min, ND_local_x_max);

    if (verbose) std::cout << "random vtx_x [cm]: " << ND_vtx_vx_vec.at(0) << std::endl;

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
        std::cout << "nd off_axis x #" << ND_off_axis_pos_counter << ": " << i_ND_off_axis_pos << " m" << std::endl;
        std::cout << "nd off_axis x in detector #" << vtx_vx_counter << ": " << i_vtx_vx << " m" << std::endl;
        std::cout << "Event position # : " << i_ND_off_axis_pos + i_vtx_vx << " m" << std::endl;


        // ND off-axis position does not affect evt vx, so only fill branches below once when loop over ND off-axis vec
        // if ( ND_off_axis_pos_counter == 1 ) {
          ND_Sim_mu_start_vx.emplace_back( i_vtx_vx+i_ND_off_axis_pos );
          ND_Sim_mu_end_vx.emplace_back( FD_Sim_mu_end_vx - FD_Sim_mu_start_vx + i_vtx_vx+i_ND_off_axis_pos ); // w.r.t. mu start x
        // }

        // Evt vtx pos in unit: cm
        eff->setVertex( i_vtx_vx+i_ND_off_axis_pos, ND_Sim_mu_start_vy, ND_Sim_mu_start_vz );
        eff->setNewVertexBF(i_vtx_vx+i_ND_off_axis_pos, ND_Sim_mu_start_OnAxis_vy, ND_Sim_mu_start_OnAxis_vz);
        eff->setMuEndV(FD_Sim_mu_end_vx - FD_Sim_mu_start_vx + i_vtx_vx+i_ND_off_axis_pos,ND_Sim_mu_end_vy,ND_Sim_mu_end_vz);
        ND_Sim_mu_end_vx_af = eff->getRotMuEndV_AF_X();
        ND_Sim_mu_end_vy_af = eff->getRotMuEndV_AF_Y();
        ND_Sim_mu_end_vz_af = eff->getRotMuEndV_AF_Z();
        ND_Sim_mu_start_px_af = eff->getRotMuStartP_AF_X();
        ND_Sim_mu_start_py_af = eff->getRotMuStartP_AF_Y();
        ND_Sim_mu_start_pz_af = eff->getRotMuStartP_AF_Z();


        HadronHitEdeps.clear();
        HadronHitPoss.clear();
        HadronHitEdeps.reserve(FD_Sim_n_hadronic_Edep_a);
        HadronHitPoss.reserve(FD_Sim_n_hadronic_Edep_a*3);

        // Need to loop deposits many times as the hadron containment below need vtx x first
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){

          // Relative to muon start pos in ND coordinate sys: (i_vtx_vx, ND_Sim_mu_start_vy, ND_Sim_mu_start_vz)
          HadronHitPoss.emplace_back( FD_Sim_hadronic_hit_x_a->at(ihadronhit) - FD_Sim_mu_start_vx + i_vtx_vx + +i_ND_off_axis_pos); // w.r.t. mu start x
          // Again, need to apply R_x(theta) for hadron y/z, do not affect x
          HadronHitPoss.emplace_back( cos( 2*abs(beamLineRotation) )*( FD_Sim_hadronic_hit_y_a->at(ihadronhit) - FD_Sim_mu_start_vy + ND_Sim_mu_start_vy ) - sin( 2*abs(beamLineRotation) )*( FD_Sim_hadronic_hit_z_a->at(ihadronhit) - FD_Sim_mu_start_vz + ND_Sim_mu_start_vz ) );
          HadronHitPoss.emplace_back( sin( 2*abs(beamLineRotation) )*( FD_Sim_hadronic_hit_y_a->at(ihadronhit) - FD_Sim_mu_start_vy + ND_Sim_mu_start_vy ) + cos( 2*abs(beamLineRotation) )*( FD_Sim_hadronic_hit_z_a->at(ihadronhit) - FD_Sim_mu_start_vz + ND_Sim_mu_start_vz ) );

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

      } // end Loop over ND_vtx_vx_vec

      ND_Sim_hadron_contain_result_before_throw.emplace_back(hadron_contain_result_before_throw_vec_for_vtx_vx);
      ND_Sim_hadron_throw_result.emplace_back(hadron_throw_result_vec_for_vtx_vx);

    }   // end Loop over ND_off_axis_pos_vec

    effTreeFD->Fill();
    iwritten++;

  }     // end loop over events

  std::cout << "Written evts: " << iwritten << std::endl;

  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff_test.root", "RECREATE");
  ThrowsFD->Write();
  effTreeFD->Write();
  hist_vetoEnergyFD->Write();
  hist_vetoEnergyFD_new->Write();

  outFile->Close();

} // end main
