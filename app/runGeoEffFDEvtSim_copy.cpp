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

  Int_t Run;
  Int_t SubRun;
  Int_t Event;
  Int_t Sim_nNumu;
  double Gen_numu_E;
  Int_t Sim_nMu;
  double Sim_mu_start_vx; // unit: cm?
  double Sim_mu_start_vy;
  double Sim_mu_start_vz;
  double Sim_mu_end_vx;
  double Sim_mu_end_vy;
  double Sim_mu_end_vz;
  double Sim_mu_start_px;
  double Sim_mu_start_py;
  double Sim_mu_start_pz;
  double Sim_mu_start_E;
  double Sim_mu_end_px;
  double Sim_mu_end_py;
  double Sim_mu_end_pz;
  double Sim_mu_end_E;
  double Sim_hadronic_Edep_a2;
  Int_t Sim_n_hadronic_Edep_a;
  vector<float> *Sim_hadronic_hit_Edep_a2 = 0; // Need initialize 0 here to avoid error
  vector<float> *Sim_hadronic_hit_x_a     = 0; // Same as mu pos unit: cm?
  vector<float> *Sim_hadronic_hit_y_a     = 0;
  vector<float> *Sim_hadronic_hit_z_a     = 0;

  // Read ntuple from FD MC
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // ntuple path on FNAL dunegpvm machine
  t->Add("/home/fyguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root");

  t->SetBranchAddress("Run",                      &Run);
  t->SetBranchAddress("SubRun",                   &SubRun);
  t->SetBranchAddress("Event",                    &Event);
  t->SetBranchAddress("Sim_nNumu",                &Sim_nNumu);
  t->SetBranchAddress("Gen_numu_E",               &Gen_numu_E);
  t->SetBranchAddress("Sim_nMu",                  &Sim_nMu);
  t->SetBranchAddress("Sim_mu_start_vx",          &Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",          &Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",          &Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",            &Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",            &Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",            &Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",          &Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",          &Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",          &Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_start_E",           &Sim_mu_start_E);
  t->SetBranchAddress("Sim_mu_end_px",            &Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",            &Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",            &Sim_mu_end_pz);
  t->SetBranchAddress("Sim_mu_end_E",             &Sim_mu_end_E);
  t->SetBranchAddress("Sim_hadronic_Edep_a2",     &Sim_hadronic_Edep_a2);
  t->SetBranchAddress("Sim_n_hadronic_Edep_a",    &Sim_n_hadronic_Edep_a);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_a2", &Sim_hadronic_hit_Edep_a2);
  t->SetBranchAddress("Sim_hadronic_hit_x_a",     &Sim_hadronic_hit_x_a);
  t->SetBranchAddress("Sim_hadronic_hit_y_a",     &Sim_hadronic_hit_y_a);
  t->SetBranchAddress("Sim_hadronic_hit_z_a",     &Sim_hadronic_hit_z_a);

  //
  // Declare variables used in this program
  //

  int nentries = 0;        // Total input events
  float vetoEnergyFD; // Total hadron deposited energy in FD veto region
  float b_vetoEnergyFD;   // Add vetoEnergy ntuple
  int iwritten = 0;        // Output event counter
  float decayZbeamCoord;
  float decayXdetCoord;
  float decayYdetCoord;
  float decayZdetCoord;
  vector<float> HadronHitEdeps;
  vector<float> HadronHitPoss;

  vector<double> b_ND_off_axis_pos_vec; // unit: meters, ND off-axis choices for each FD evt: 1st element is randomized for each evt
  vector<double> b_vtx_vx_vec;          // unit: cm, vtx x choices for each FD evt in ND volume: 1st element is randomized for each evt
  int ND_off_axis_pos_steps = 0;
  int vtx_vx_steps = 0;

  // Initialize first element as -999, to be replaced by a random off-axis nd pos in each evt below
  b_ND_off_axis_pos_vec.emplace_back(-999.);

  if ( ND_off_axis_pos_stepsize > 0 && ND_off_axis_pos_stepsize <= OffAxisPoints[13] ) {
    ND_off_axis_pos_steps = ( OffAxisPoints[13] - OffAxisPoints[0] ) / ND_off_axis_pos_stepsize;
  }
  else std::cout << "Error: please set the ND_off_axis_pos_stepsize above 0 and below max element of OffAxisPoints." << std::endl;

  if (verbose) std::cout << "ND_off_axis_pos_steps: " << ND_off_axis_pos_steps << std::endl;

  // The rest elements follow fixed increments from min ND local x
  for ( int i_ND_off_axis_pos_step = 0; i_ND_off_axis_pos_step < ND_off_axis_pos_steps + 1; i_ND_off_axis_pos_step++ ){
    b_ND_off_axis_pos_vec.emplace_back( i_ND_off_axis_pos_step*ND_off_axis_pos_stepsize + OffAxisPoints[0] );
  }

  if (verbose) std::cout << "b_ND_off_axis_pos_vec size: "<< b_ND_off_axis_pos_vec.size() << std::endl;

  // Initialize first element as -999, to be replaced by a random vtx x in each evt below
  b_vtx_vx_vec.emplace_back(-999.);

  if ( ND_local_x_stepsize > 0 && ND_local_x_stepsize <= ND_local_x_max ) {
    vtx_vx_steps = ( ND_local_x_max - ND_local_x_min ) / ND_local_x_stepsize;
  }
  else std::cout << "Error: please set the ND_local_x_stepsize above 0 and below ND_local_x_max." << std::endl;

  if (verbose) std::cout << "vtx_vx_steps: " << vtx_vx_steps << std::endl;

  // The rest elements follow fixed increments from min ND local x
  for ( int i_vtx_vx_step = 0; i_vtx_vx_step < vtx_vx_steps + 1; i_vtx_vx_step++ ){
    b_vtx_vx_vec.emplace_back( i_vtx_vx_step*ND_local_x_stepsize + ND_local_x_min );
  }

  if (verbose) std::cout << "b_vtx_vx_vec size: "<< b_vtx_vx_vec.size() << std::endl;

  // Lepton info: expressed in ND coordinate sys, do not confuse with branches read above in FD coordinate sys
  double b_Gen_numu_E;
  vector<double> b_Sim_mu_start_vx; // Vector corresponds to randomized/stepwise vtx x
  vector<double> b_Sim_mu_end_vx;
  double b_Sim_mu_start_vy;         // Do not use float!
  double b_Sim_mu_start_vz;
  double b_Sim_mu_end_vy;
  double b_Sim_mu_end_vz;
  double b_Sim_mu_start_px;
  double b_Sim_mu_start_py;
  double b_Sim_mu_start_pz;
  double b_Sim_mu_start_E;
  double b_Sim_mu_end_px;
  double b_Sim_mu_end_py;
  double b_Sim_mu_end_pz;
  double b_Sim_mu_end_E;
  // Tot hadron E dep
  double b_Sim_hadronic_Edep_a2;
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
  vector<vector<vector<vector<bool> > > > b_Sim_hadron_contain_result_before_throw; // ................ ................... ND off axis position x
  // Nested vector (vetoSize > vetoEnergy > 64_bit_throw_result) returned from function getHadronContainmentThrows
  vector<vector<vector<uint64_t> > > hadron_throw_result;
  vector<vector<vector<vector<uint64_t> > > > hadron_throw_result_vec_for_vtx_vx;   // One more vector: randomized/stepwise evt vtx x
  vector<vector<vector<vector<vector<uint64_t> > > > > b_Sim_hadron_throw_result;   // ................ ................... ND off axis position x

  //
  // A tree to store lepton info (for NN training)
  // and result of hadron containment after applying transformations
  //

  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  effTreeFD->Branch("Gen_numu_E",                             &b_Gen_numu_E,            "Gen_numu_E/D");
  effTreeFD->Branch("ND_off_axis_pos_vec",                    &b_ND_off_axis_pos_vec);                             // vector<double>: entries = written evts * ND_off_axis_pos_steps
  effTreeFD->Branch("Sim_mu_start_vx",                        &b_Sim_mu_start_vx);                                 // ............... entries = written evts * vtx_vx_steps, equivalent to b_vtx_vx_vec
  effTreeFD->Branch("Sim_mu_start_vy",                        &b_Sim_mu_start_vy,       "Sim_mu_start_vy/D");      // entries = written evts
  effTreeFD->Branch("Sim_mu_start_vz",                        &b_Sim_mu_start_vz,       "Sim_mu_start_vz/D");
  effTreeFD->Branch("Sim_mu_end_vx",                          &b_Sim_mu_end_vx);
  effTreeFD->Branch("Sim_mu_end_vy",                          &b_Sim_mu_end_vy,         "Sim_mu_end_vy/D");
  effTreeFD->Branch("Sim_mu_end_vz",                          &b_Sim_mu_end_vz,         "Sim_mu_end_vz/D");
  effTreeFD->Branch("Sim_mu_start_px",                        &b_Sim_mu_start_px,       "Sim_mu_start_px/D");
  effTreeFD->Branch("Sim_mu_start_py",                        &b_Sim_mu_start_py,       "Sim_mu_start_py/D");
  effTreeFD->Branch("Sim_mu_start_pz",                        &b_Sim_mu_start_pz,       "Sim_mu_start_pz/D");
  effTreeFD->Branch("Sim_mu_start_E",                         &b_Sim_mu_start_E,        "Sim_mu_start_E/D");
  effTreeFD->Branch("Sim_mu_end_px",                          &b_Sim_mu_end_px,         "Sim_mu_end_px/D");
  effTreeFD->Branch("Sim_mu_end_py",                          &b_Sim_mu_end_py,         "Sim_mu_end_py/D");
  effTreeFD->Branch("Sim_mu_end_pz",                          &b_Sim_mu_end_pz,         "Sim_mu_end_pz/D");
  effTreeFD->Branch("Sim_mu_end_E",                           &b_Sim_mu_end_E,          "Sim_mu_end_E/D");
  effTreeFD->Branch("Sim_hadron_contain_result_before_throw", &b_Sim_hadron_contain_result_before_throw);          // nested vector
  effTreeFD->Branch("Sim_hadron_throw_result",                &b_Sim_hadron_throw_result);
  effTreeFD->Branch("Sim_hadronic_Edep_a2",                   &b_Sim_hadronic_Edep_a2,  "Sim_hadronic_Edep_a2/D"); // entries = written evts
  // Add VetoE_FD branch
  TBranch *vetoE = effTreeFD->Branch("VetoEnergyFD",                           &b_vetoEnergyFD,            "vetoEnergyFD/F");
  
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
  double beamRefDetCoord[3] = {0., 0., 15.5}; // (0, 0, 0) is ND detector origin
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

  //
  // Loop over FD events
  //

  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < nentries; ientry++ ) {

    t->GetEntry(ientry);
    if ( ientry%10000 == 0 ) std::cout << "Looking at entry " << ientry << ", run: " << Run << ", subrun: " << SubRun << ", event: " << Event << std::endl;

    //
    // Skip events without muon/hadronic deposits
    //

    if ( Sim_nMu == 0 || Sim_n_hadronic_Edep_a == 0 ) continue;

    //
    // Calculate total hadron E in FD veto region
    //

    vetoEnergyFD = 0.;
    // Loop over hadron E deposits
    for ( int ihadronhit = 0; ihadronhit < Sim_n_hadronic_Edep_a; ihadronhit++ ){

      // Veto region size: 30 cm from the active volume
      if ( ( Sim_hadronic_hit_x_a->at(ihadronhit) < FDActiveVol_min[0] + 30 ) ||
           ( Sim_hadronic_hit_y_a->at(ihadronhit) < FDActiveVol_min[1] + 30 ) ||
           ( Sim_hadronic_hit_z_a->at(ihadronhit) < FDActiveVol_min[2] + 30 ) ||
           ( Sim_hadronic_hit_x_a->at(ihadronhit) > FDActiveVol_max[0] - 30 ) ||
           ( Sim_hadronic_hit_y_a->at(ihadronhit) > FDActiveVol_max[1] - 30 ) ||
           ( Sim_hadronic_hit_z_a->at(ihadronhit) > FDActiveVol_max[2] - 30 )
         ){
           vetoEnergyFD += Sim_hadronic_hit_Edep_a2->at(ihadronhit);
      } // end if hadron deposit in FD veto region

    } // end loop over hadron E deposits

    // Add vetoEnergyFD ntuple before threshold 
    b_vetoEnergyFD = vetoEnergyFD; 
    vetoE->Fill();
    
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

    //
    // Use relative coordinate for FD event (relative to muon start position)
    // Set muon start pos in ND to (random x/stepwise increased x, random y, random z),
    // eventually should set this for actual event vtx, not mu start pos
    //

    // Branches that are not affected by ND off axis position and vtx x (loops below)
    b_Gen_numu_E      = Gen_numu_E;

    // Use random y/z for each FD evt in ND volume, it will be randomly translated later anyway
    TRandom3 *r3_mu_start_vy_nd = new TRandom3(); // Initialize random number generator, put inside the event loop so each event is different
    r3_mu_start_vy_nd->SetSeed(0);                // Set the seed (required to avoid repeated random numbers in each sequence)
    b_Sim_mu_start_vy = r3_mu_start_vy_nd->Uniform( NDActiveVol_min[1], NDActiveVol_max[1] );

    TRandom3 *r3_mu_start_vz_nd = new TRandom3();
    r3_mu_start_vz_nd->SetSeed(0);
    b_Sim_mu_start_vz = r3_mu_start_vz_nd->Uniform( NDActiveVol_min[2], NDActiveVol_max[2] );

    b_Sim_mu_end_vy   = Sim_mu_end_vy - Sim_mu_start_vy + b_Sim_mu_start_vy; // w.r.t. mu start random y in ND
    b_Sim_mu_end_vz   = Sim_mu_end_vz - Sim_mu_start_vz + b_Sim_mu_start_vz;
    b_Sim_mu_start_px = Sim_mu_start_px;                                     // momentum is not affected by coordinate
    b_Sim_mu_start_py = Sim_mu_start_py;
    b_Sim_mu_start_pz = Sim_mu_start_pz;
    b_Sim_mu_start_E  = Sim_mu_start_E;
    b_Sim_mu_end_px   = Sim_mu_end_px;
    b_Sim_mu_end_py   = Sim_mu_end_py;
    b_Sim_mu_end_pz   = Sim_mu_end_pz;
    b_Sim_mu_end_E    = Sim_mu_end_E;

    b_Sim_hadronic_Edep_a2 = Sim_hadronic_Edep_a2;

    // Local y-z axes in FD and ND are rotated due to Earth curvature
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
    b_Sim_mu_start_vy = cos( 2*abs(beamLineRotation) )*b_Sim_mu_start_vy - sin( 2*abs(beamLineRotation) )*b_Sim_mu_start_vz;
    b_Sim_mu_start_vz = sin( 2*abs(beamLineRotation) )*b_Sim_mu_start_vy + cos( 2*abs(beamLineRotation) )*b_Sim_mu_start_vz;

    b_Sim_mu_end_vy = cos( 2*abs(beamLineRotation) )*b_Sim_mu_end_vy - sin( 2*abs(beamLineRotation) )*b_Sim_mu_end_vz;
    b_Sim_mu_end_vz = sin( 2*abs(beamLineRotation) )*b_Sim_mu_end_vy + cos( 2*abs(beamLineRotation) )*b_Sim_mu_end_vz;

    // Initialize for the event
    b_Sim_mu_start_vx.clear();
    b_Sim_mu_end_vx.clear();
    b_Sim_hadron_contain_result_before_throw.clear();
    b_Sim_hadron_throw_result.clear();

    //
    // Two options for setting ND off-axis position
    //
    // Option 1: randomize ND off-axis position once for each event
    // Option 2: stepwise increments along ND detector hall off axis range
    //
    // If only want option 1, set random_ND_off_axis_pos to true in Helpers.h; default is false (use both options)
    //

    // Initialize random number generator
    // This needs to be inside the event loop to make sure each event has a different random number
    TRandom3 *r3_OffAxisPoint = new TRandom3();
    // Set the seed (required to avoid repeated random numbers in each sequence)
    r3_OffAxisPoint->SetSeed(0);
    b_ND_off_axis_pos_vec.at(0) = r3_OffAxisPoint->Uniform(OffAxisPoints[0], OffAxisPoints[13]);

    if (verbose) std::cout << "random OffAxisPoint [meters]: " << b_ND_off_axis_pos_vec.at(0) << std::endl;

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
    b_vtx_vx_vec.at(0) = r3_vtx_x->Uniform(ND_local_x_min, ND_local_x_max);

    if (verbose) std::cout << "random vtx_x [cm]: " << b_vtx_vx_vec.at(0) << std::endl;

    //
    // Loop over b_ND_off_axis_pos_vec: random off_axis_pos or every ND_off_axis_pos_stepsize
    // Don't put it outside event loop to avoid looping over all events multiple times
    //

    int ND_off_axis_pos_counter = 0;
    for ( double i_ND_off_axis_pos : b_ND_off_axis_pos_vec ) {

      // Skip the stepwise increased option if only want a random off axis ND position per event to save file size
      if ( random_ND_off_axis_pos && ND_off_axis_pos_counter != 0 ) continue;

      ND_off_axis_pos_counter++;

      // Interpolate event neutrino production point (beam coordinate)
      decayZbeamCoord = gDecayZ->Eval( i_ND_off_axis_pos - detRefBeamCoord[0] );

      // Calculate neutrino production point in detector coordinate
      decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation);
      decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation);
      // Set production point in unit: cm
      eff->setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.);

      if (verbose) std::cout << "nd off_axis x #" << ND_off_axis_pos_counter << ": " << i_ND_off_axis_pos << " m" << std::endl;

      //
      // Loop over vtx x: random x or stepwise increased x
      // Don't put it outside event loop to avoid looping over all events multiple times
      //

      int vtx_vx_counter = 0;
      hadron_throw_result_vec_for_vtx_vx.clear(); // Need to initialize before loop over vtx x vec
      hadron_contain_result_before_throw_vec_for_vtx_vx.clear();

      for ( double i_vtx_vx : b_vtx_vx_vec ) {

        // Skip the stepwise increased option if only want a random evt vtx x to save file size
        if ( random_vtx_vx && vtx_vx_counter != 0 ) continue;

        vtx_vx_counter++;

        // ND off-axis position does not affect evt vx, so only fill branches below once when loop over ND off-axis vec
        if ( ND_off_axis_pos_counter == 1 ) {
          b_Sim_mu_start_vx.emplace_back( i_vtx_vx );
          b_Sim_mu_end_vx.emplace_back( Sim_mu_end_vx - Sim_mu_start_vx + i_vtx_vx ); // w.r.t. mu start x
        }

        // Evt vtx pos in unit: cm
        eff->setVertex( i_vtx_vx, b_Sim_mu_start_vy, b_Sim_mu_start_vz );

        HadronHitEdeps.clear();
        HadronHitPoss.clear();
        HadronHitEdeps.reserve(Sim_n_hadronic_Edep_a);
        HadronHitPoss.reserve(Sim_n_hadronic_Edep_a*3);

        // Need to loop deposits many times as the hadron containment below need vtx x first
        for ( int ihadronhit = 0; ihadronhit < Sim_n_hadronic_Edep_a; ihadronhit++ ){

          // Relative to muon start pos in ND coordinate sys: (i_vtx_vx, b_Sim_mu_start_vy, b_Sim_mu_start_vz)
          HadronHitPoss.emplace_back( Sim_hadronic_hit_x_a->at(ihadronhit) - Sim_mu_start_vx + i_vtx_vx ); // w.r.t. mu start x
          // Again, need to apply R_x(theta) for hadron y/z, do not affect x
          HadronHitPoss.emplace_back( cos( 2*abs(beamLineRotation) )*( Sim_hadronic_hit_y_a->at(ihadronhit) - Sim_mu_start_vy + b_Sim_mu_start_vy ) - sin( 2*abs(beamLineRotation) )*( Sim_hadronic_hit_z_a->at(ihadronhit) - Sim_mu_start_vz + b_Sim_mu_start_vz ) );
          HadronHitPoss.emplace_back( sin( 2*abs(beamLineRotation) )*( Sim_hadronic_hit_y_a->at(ihadronhit) - Sim_mu_start_vy + b_Sim_mu_start_vy ) + cos( 2*abs(beamLineRotation) )*( Sim_hadronic_hit_z_a->at(ihadronhit) - Sim_mu_start_vz + b_Sim_mu_start_vz ) );

          HadronHitEdeps.emplace_back( Sim_hadronic_hit_Edep_a2->at(ihadronhit) );

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

        if (verbose) std::cout << "                           vtx x #" << vtx_vx_counter << ": " << i_vtx_vx << " cm, throw result[0][0][0]: " << hadron_throw_result[0][0][0] << std::endl;

      } // end Loop over b_vtx_vx_vec

      b_Sim_hadron_contain_result_before_throw.emplace_back(hadron_contain_result_before_throw_vec_for_vtx_vx);
      b_Sim_hadron_throw_result.emplace_back(hadron_throw_result_vec_for_vtx_vx);

    }   // end Loop over b_ND_off_axis_pos_vec

    effTreeFD->Fill();
    iwritten++;

  }     // end loop over events

  std::cout << "Written evts: " << iwritten << std::endl;

  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff.root", "RECREATE");
  ThrowsFD->Write();
  effTreeFD->Write();

  outFile->Close();

} // end main
