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

// The program which do all translations and rotations when moving event from FD to ND:
// 0. FD: read event from FD MC ntuple: before earth curvature rotation
// 1. FD to ND: after earth curvature rotation
// 2. ND: move back to the beam center
// 3. ND to ND: translate from OnAxis to OffAxis
// 4. ND: get after eigen rotated vectors for step 3
// 5. ND: generate random throws
int main()
{
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 0. FD: read event from FD MC ntuple: before earth curvature rotation
  //
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  // Ntuple path on FNAL dunegpvm machine
  // For Ivy machine:
  t->Add("/home/fyguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple_ivy.root");
  // For FNAL machine:
  // t->Add("/dune/app/users/flynnguo/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root");

  // Define variables for FD event
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
  vector<float> *FD_Sim_hadronic_hit_Edep_a2 = 0; // Need initialize 0 here to avoid error
  vector<float> *FD_Sim_hadronic_hit_x_a = 0; // Position of each energy deposit on the x-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_y_a = 0; // Position of each energy deposit on the y-axis [cm]
  vector<float> *FD_Sim_hadronic_hit_z_a = 0; // Position of each energy deposit on the z-axis [cm]

  // Extract event info from ntuple
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
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Declare variables used in this program
  //
  int nentries = 0; // Total input events
  float vetoEnergyFD; // Total hadron deposited energy in FD veto region
  int iwritten = 0; // Output event counter
  float decayZbeamCoord; // Decay point (neutrino production point) in beam coordinate [cm]
  float decayXdetCoord; // Decay point (neutrino production point) in detector coordinate on the x-axis [cm]
  float decayYdetCoord; // Decay point (neutrino production point) in detector coordinate on the y-axis [cm]
  float decayZdetCoord; // Decay point (neutrino production point) in detector coordinate on the z-axis [cm]
  vector<float> HadronHitEdeps; // Hadron hit segment energy deposits [MeV]
  vector<float> HadronHitPoss;  // Hadron hit segment energy deposits position [cm]

  vector<double> ND_off_axis_pos_vec = {0,-7*100,-30*100}; // ND off-axis positions' choices for each FD evt [cm]
  vector<double> ND_vtx_vx_vec={-200,0,200};          // Vtx x choices for each FD evt in ND volume [cm]
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Result of hadron containment is stored in nested vectors, need to generate dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*float*"); // Remove old dictionary if exists
  gSystem->Exec("rm -f AutoDict*vector*vector*bool*"); // Remove old dictionary if exists
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*bool*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*bool*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*uint64_t*");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*vector*uint64_t*");
  // Generate new dictionary
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<bool> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<bool> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<bool> > > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<uint64_t> > > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<vector<uint64_t> > > > >", "vector");
  // Nested vector (vetoSize > vetoEnergy > 64_bit_throw_result) returned from function getHadronContainmentThrows
  vector<vector<vector<uint64_t> > > hadron_throw_result;
  vector<vector<vector<vector<uint64_t> > > > hadron_throw_result_vec_for_vtx_vx;   // One more vector: randomized/stepwise evt vtx x
  vector<vector<vector<vector<vector<uint64_t> > > > > ND_Sim_hadron_throw_result;   // ................ ................... ND off axis position x
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 1. FD to ND: after earth curvature rotation
  //
  // Lepton info: expressed in ND coordinate sys, do not confuse with branches read above in FD coordinate sys
  double ND_Gen_numu_E; // Energy of generator level neutrino [GeV]
  // Muon info
  // array: (x,y,z)<->dim=(0,1,2)
  double ND_RandomVtx_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_RandomVtx_Sim_mu_end_v[3]; // Position of the muon trajectory at end point [cm]
  double ND_RandomVtx_Sim_mu_start_p[3]; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  // double ND_RandomVtx_Sim_mu_start_E; // Energy of leading mu at start point [GeV]
  // double ND_RandomVtx_Sim_mu_end_E; // Energy of leading mu at end point [GeV]
  // Hadron info
  double ND_Sim_hadronic_Edep_a2; // Total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
  vector<vector<float>> ND_RandomVtx_Sim_hadronic_hit; // Position of each energy deposit [cm]
  vector<float> ND_RandomVtx_Sim_hadronic_hit_xyz; // Position of each energy deposit [cm]
  // Add veto E with two different restrictions
  double ND_vetoEnergyFD; // Total hadron deposited energy in FD veto region
  // Momentum conservation check
  float ND_RandomVtx_Sim_mu_start_p_total; // Total momentum of ND_OnAxis_Sim_mu_start_p_total

  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 2. ND: move back to the beam center
  //
  // Muon info
  // array: (x,y,z)<->dim=(0,1,2)
  double ND_OnAxis_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_OnAxis_Sim_mu_end_v[3]; // Position of the muon trajectory at end point [cm]
  double ND_OnAxis_Sim_mu_start_p[3]; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  // Hadron info
  vector<vector<float>> ND_OnAxis_Sim_hadronic_hit; // Position of each energy deposit [cm]
  vector<float> ND_OnAxis_Sim_hadronic_hit_xyz;

  // Momentum conservation check
  float ND_OnAxis_Sim_mu_start_p_total; // Total momentum of ND_OnAxis_Sim_mu_start_p_total
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 3. ND to ND: translate from OnAxis to OffAxis
  //
  // Muon info
  // array: (x,y,z)<->dim=(0,1,2)
  double ND_OffAxis_Unrotated_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_OffAxis_Unrotated_Sim_mu_end_v[3]; // Position of the muon trajectory at end point [cm]
  double ND_OffAxis_Unrotated_Sim_mu_start_p[3]; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  // Hadron info
  vector<vector<float>> ND_OffAxis_Unrotated_Sim_hadronic_hit; // Position of each energy deposit [cm]
  vector<float> ND_OffAxis_Unrotated_Sim_hadronic_hit_xyz; // Position of each energy deposit [cm]
  // Momentum conservation check
  float ND_OffAxis_Unrotated_Sim_mu_start_p_total; // Total momentum of ND_OffAxis_Unrotated_Sim_mu_start_p
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // 4. ND: get after eigen rotated vectors for step 4
  //
  // Muon info
  // array: (x,y,z)<->dim=(0,1,2)
  double ND_OffAxis_Sim_mu_start_v[3]; // Position of the muon trajectory at start point [cm]
  double ND_OffAxis_Sim_mu_end_v[3]; // Position of the muon trajectory at end point [cm]
  double ND_OffAxis_Sim_mu_start_p[3]; // Momentum of the muon trajectory at start point on the x-axis [GeV]
  // Hadron info
  vector<vector<float>> ND_OffAxis_Sim_hadronic_hit; // Position of each energy deposit [cm]
  vector<float> ND_OffAxis_Sim_hadronic_hit_xyz;
  // Momentum conservation check
  float ND_OffAxis_Sim_mu_start_p_total; // Total momentum of ND_OffAxis_Sim_mu_start_p
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Store variables into a tree
  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  effTreeFD->Branch("ND_Gen_numu_E",                             &ND_Gen_numu_E,            "ND_Gen_numu_E/D");
  effTreeFD->Branch("ND_vetoEnergyFD",                           &ND_vetoEnergyFD,          "ND_vetoEnergyFD/D");
  effTreeFD->Branch("ND_Sim_n_hadronic_Edep_a",                  &FD_Sim_n_hadronic_Edep_a,            "FD_Sim_n_hadronic_Edep_a/I");
  // 1. FD to ND: after earth curvature rotation
  effTreeFD->Branch("ND_RandomVtx_Sim_mu_start_v",               ND_RandomVtx_Sim_mu_start_v,       "ND_RandomVtx_Sim_mu_start_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_RandomVtx_Sim_mu_end_v",                 ND_RandomVtx_Sim_mu_end_v,         "ND_RandomVtx_Sim_mu_end_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_RandomVtx_Sim_mu_start_p",               ND_RandomVtx_Sim_mu_start_p,       "ND_RandomVtx_Sim_mu_start_p[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_Sim_hadronic_Edep_a2",                   &ND_Sim_hadronic_Edep_a2,          "ND_Sim_hadronic_Edep_a2/D"); // entries = written evts
  effTreeFD->Branch("ND_RandomVtx_Sim_hadronic_hit_xyz",         &ND_RandomVtx_Sim_hadronic_hit);

  // 2. ND: move back to the beam center
  effTreeFD->Branch("ND_OnAxis_Sim_mu_start_v",               ND_OnAxis_Sim_mu_start_v,       "ND_OnAxis_Sim_mu_start_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OnAxis_Sim_mu_end_v",                 ND_OnAxis_Sim_mu_end_v,         "ND_OnAxis_Sim_mu_end_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OnAxis_Sim_mu_start_p",               ND_OnAxis_Sim_mu_start_p,       "ND_OnAxis_Sim_mu_start_p[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OnAxis_Sim_hadronic_hit",             &ND_OnAxis_Sim_hadronic_hit);
  // 3. ND to ND: translate from OnAxis to OffAxis
  effTreeFD->Branch("ND_off_axis_pos_vec",                       &ND_off_axis_pos_vec);                             // vector<double>: entries = written evts * ND_off_axis_pos_steps
  effTreeFD->Branch("ND_vtx_vx_vec",                             &ND_vtx_vx_vec);
  effTreeFD->Branch("ND_OffAxis_Unrotated_Sim_mu_start_v",               ND_OffAxis_Unrotated_Sim_mu_start_v,       "ND_OffAxis_Unrotated_Sim_mu_start_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Unrotated_Sim_mu_end_v",                 ND_OffAxis_Unrotated_Sim_mu_end_v,         "ND_OffAxis_Unrotated_Sim_mu_end_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Unrotated_Sim_mu_start_p",               ND_OffAxis_Unrotated_Sim_mu_start_p,       "ND_OffAxis_Unrotated_Sim_mu_start_p[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Unrotated_Sim_hadronic_hit_xyz",         &ND_OffAxis_Unrotated_Sim_hadronic_hit);
  // 4. ND: get after eigen rotated vectors for step 4
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_v",               ND_OffAxis_Sim_mu_start_v,       "ND_OffAxis_Sim_mu_start_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Sim_mu_end_v",                 ND_OffAxis_Sim_mu_end_v,         "ND_OffAxis_Sim_mu_end_v[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Sim_mu_start_p",               ND_OffAxis_Sim_mu_start_p,       "ND_OffAxis_Sim_mu_start_p[3]/D");   // entries = written evts*3
  effTreeFD->Branch("ND_OffAxis_Sim_hadronic_hit_xyz",         &ND_OffAxis_Sim_hadronic_hit);
  // 5. ND: generate random throws
  effTreeFD->Branch("ND_Sim_hadron_throw_result",              &ND_Sim_hadron_throw_result);
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
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
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Get beam parameters: eventually need to read from XML file: which XML for FD? same as ND?
  //
  // Clockwise rotate beamline around ND local x axis
  double beamLineRotation = -0.101;           // unit: rad

  // Put back into beam center(0.0, 0.05387, 6.66) tanslation 1 [cm]
  ND_OnAxis_Sim_mu_start_v[0] = 0.0*100.;
  ND_OnAxis_Sim_mu_start_v[1] = 0.05387*100.;
  ND_OnAxis_Sim_mu_start_v[2] = 6.6*100.;
  // Coordinate transformation, units: meters
  double beamRefDetCoord[3] = {ND_OnAxis_Sim_mu_start_v[0]/100., ND_OnAxis_Sim_mu_start_v[1]/100., ND_OnAxis_Sim_mu_start_v[2]/100.}; // (0, 0, 0) is ND detector origin
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
  eff->setOnAxisVertex(ND_OnAxis_Sim_mu_start_v[0],ND_OnAxis_Sim_mu_start_v[1],ND_OnAxis_Sim_mu_start_v[2]);

  // Add hist of veto E
  TH1F *hist_vetoEnergyFD = new TH1F("hist_vetoEnergyFD", "hist_vetoEnergyFD", 1500, 0, 1500);
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Add output txt file
  ofstream myfile;
   myfile.open ("Output_FDGeoEff_DataCheck_ivy.txt");
  //
  //
  // Loop over FD events
  //
  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;
  myfile << "Tot evts: " << nentries << "\n";
  for ( int ientry = 0; ientry < nentries; ientry++ )
  {
    t->GetEntry(ientry);
    if ( ientry%10000 == 0 )
    {
      std::cout << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << std::endl;
      myfile << "Looking at entry " << ientry << ", FD_run: " << FD_Run << ", FD_subrun: " << FD_SubRun << ", FD_event: " << FD_Event << "\n\n\n";
    }
    myfile << "\n ientry: " << ientry <<"\n\n";
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
    //
    // Skip FD event if the total hadron E in veto region exceeds vetoEnergy [MeV]
    //
    if ( vetoEnergyFD > 30 ) continue; // 30 MeV
    //
    // Renew throws every 100th written event to save file size, i.e., if N = 128,
    // for written evt 0-99:    same 128 transformations for each event,
    // for written evt 100-199: same but renewed 128 transformations for each evt
    // so on so forth...
    // These transformations will be applied to leptons in the event, so need to keep track of iwritten
    //
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
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //
    // 1. FD to ND: after earth curvature rotation
    // Rotation affects mu start/end position and hadron positions below
    // Mu_start_v
    double FD_Sim_mu_start_v[3] = {FD_Sim_mu_start_vx, FD_Sim_mu_start_vy, FD_Sim_mu_start_vz};
    for(int i=0; i<3; i++)
    {
      ND_RandomVtx_Sim_mu_start_v[i] = eff->getEarthCurvature(FD_Sim_mu_start_v, beamLineRotation, i);
    }
    myfile << "FD_Sim_mu_start_vx[cm]: " << FD_Sim_mu_start_vx << "\n";
    myfile << "FD_Sim_mu_start_vy[cm]: " << FD_Sim_mu_start_vy << "\n";
    myfile << "FD_Sim_mu_start_vz[cm]: " << FD_Sim_mu_start_vz << "\n\n";
    myfile << "ND_RandomVtx_Sim_mu_start_vx[cm]: " << ND_RandomVtx_Sim_mu_start_v[0] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_start_vy[cm]: " << ND_RandomVtx_Sim_mu_start_v[1] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_start_vz[cm]: " << ND_RandomVtx_Sim_mu_start_v[2] << "\n\n";
    // Mu_end_v
    double FD_Sim_mu_end_v[3] = {FD_Sim_mu_end_vx, FD_Sim_mu_end_vy, FD_Sim_mu_end_vz};
    for(int i=0; i<3; i++)
    {
      ND_RandomVtx_Sim_mu_end_v[i] = eff->getEarthCurvature(FD_Sim_mu_end_v, beamLineRotation, i);
    }
    myfile << "FD_Sim_mu_end_vx[cm]: " << FD_Sim_mu_end_vx << "\n";
    myfile << "FD_Sim_mu_end_vy[cm]: " << FD_Sim_mu_end_vy << "\n";
    myfile << "FD_Sim_mu_end_vz[cm]: " << FD_Sim_mu_end_vz << "\n";
    double FD_Sim_mu_v_end_start = eff->getDistance(FD_Sim_mu_start_v,FD_Sim_mu_end_v);
    myfile << "Distance between FD_Sim_mu_v_end_start[cm]: " << FD_Sim_mu_v_end_start << "\n\n";
    myfile << "ND_RandomVtx_Sim_mu_end_vx[cm]: " << ND_RandomVtx_Sim_mu_end_v[0] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_end_vy[cm]: " << ND_RandomVtx_Sim_mu_end_v[1] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_end_vz[cm]: " << ND_RandomVtx_Sim_mu_end_v[2] << "\n";
    double ND_RandomVtx_Sim_mu_v_end_start = eff->getDistance(ND_RandomVtx_Sim_mu_start_v,ND_RandomVtx_Sim_mu_end_v);
    myfile << "Distance between ND_RandomVtx_Sim_mu_v_end_start[cm]: " << ND_RandomVtx_Sim_mu_v_end_start << "\n\n";

    // Mu_start_p
    double FD_Sim_mu_start_p[3] = {FD_Sim_mu_start_px, FD_Sim_mu_start_py, FD_Sim_mu_start_pz};
    for(int i=0; i<3; i++)
    {
      ND_RandomVtx_Sim_mu_start_p[i] = eff->getEarthCurvature(FD_Sim_mu_start_p, beamLineRotation, i);
    }
    ND_RandomVtx_Sim_mu_start_p_total = eff->getTotalMomentum(ND_RandomVtx_Sim_mu_start_p);
    double FD_Sim_mu_start_p_total = eff->getTotalMomentum(FD_Sim_mu_start_p);
    myfile << "FD_Sim_mu_start_px[cm]: " << FD_Sim_mu_start_px << "\n";
    myfile << "FD_Sim_mu_start_py[cm]: " << FD_Sim_mu_start_py << "\n";
    myfile << "FD_Sim_mu_start_pz[cm]: " << FD_Sim_mu_start_pz << "\n";
    myfile << "FD_Sim_mu_start_p_total[GeV]: " << FD_Sim_mu_start_p_total << "\n\n";
    myfile << "ND_RandomVtx_Sim_mu_start_px[GeV]: " << ND_RandomVtx_Sim_mu_start_p[0] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_start_py[GeV]: " << ND_RandomVtx_Sim_mu_start_p[1] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_start_pz[GeV]: " << ND_RandomVtx_Sim_mu_start_p[2] << "\n";
    myfile << "ND_RandomVtx_Sim_mu_start_p_total[GeV]: " << ND_RandomVtx_Sim_mu_start_p_total << "\n\n";

    // Initialize for the event
    ND_Sim_hadron_throw_result.clear();
    // Branches that are not affected by ND off axis position and vtx x (loops below)
    ND_Gen_numu_E      = FD_Gen_numu_E;
    ND_Sim_hadronic_Edep_a2 = FD_Sim_hadronic_Edep_a2;


    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
      double FD_Sim_hadronic_hit[3] = {FD_Sim_hadronic_hit_x_a->at(ihadronhit),FD_Sim_hadronic_hit_y_a->at(ihadronhit), FD_Sim_hadronic_hit_z_a->at(ihadronhit)};
      for (int i =0; i<3;i++)
      {
        ND_RandomVtx_Sim_hadronic_hit_xyz.emplace_back((float)eff->getEarthCurvature(FD_Sim_hadronic_hit, beamLineRotation, i));
      }
      ND_RandomVtx_Sim_hadronic_hit.emplace_back(ND_RandomVtx_Sim_hadronic_hit_xyz);
      ND_RandomVtx_Sim_hadronic_hit_xyz.clear();
    }
    //
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //
    // 2. ND: move back to the beam center
    // ND_OnAxis_Sim_mu_end_v
    // double geoEff::getTranslations(double v_bf[3], double vtx_bf[3], double vtx_af[3], int dim)
    for(int i=0; i<3; i++)
    {
      ND_OnAxis_Sim_mu_end_v[i] = eff->getTranslations(ND_RandomVtx_Sim_mu_end_v, ND_RandomVtx_Sim_mu_start_v, ND_OnAxis_Sim_mu_start_v, i);
    }
    myfile << "ND_OnAxis_Sim_mu_start_vx[cm]: " << ND_OnAxis_Sim_mu_start_v[0] << "\n";
    myfile << "ND_OnAxis_Sim_mu_start_vy[cm]: " << ND_OnAxis_Sim_mu_start_v[1] << "\n";
    myfile << "ND_OnAxis_Sim_mu_start_vz[cm]: " << ND_OnAxis_Sim_mu_start_v[2] << "\n\n";
    myfile << "ND_OnAxis_Sim_mu_end_vx[cm]: " << ND_OnAxis_Sim_mu_end_v[0] << "\n";
    myfile << "ND_OnAxis_Sim_mu_end_vy[cm]: " << ND_OnAxis_Sim_mu_end_v[1] << "\n";
    myfile << "ND_OnAxis_Sim_mu_end_vz[cm]: " << ND_OnAxis_Sim_mu_end_v[2] << "\n";
    double ND_OnAxis_Sim_mu_v_end_start = eff->getDistance(ND_OnAxis_Sim_mu_start_v,ND_OnAxis_Sim_mu_end_v);
    myfile << "Distance between ND_OnAxis_Sim_mu_v_end_start[cm]: " << ND_OnAxis_Sim_mu_v_end_start << "\n\n";
    // ND_OnAxis_Sim_mu_start_p
    for(int i=0; i<3; i++)
    {
      ND_OnAxis_Sim_mu_start_p[i] = eff->RemainUnchanged(ND_RandomVtx_Sim_mu_start_p[i]);
    }
    ND_OnAxis_Sim_mu_start_p_total = eff->getTotalMomentum(ND_OnAxis_Sim_mu_start_p);
    myfile << "ND_OnAxis_Sim_mu_start_px[cm]: " << ND_OnAxis_Sim_mu_start_p[0] << "\n";
    myfile << "ND_OnAxis_Sim_mu_start_py[cm]: " << ND_OnAxis_Sim_mu_start_p[1] << "\n";
    myfile << "ND_OnAxis_Sim_mu_start_pz[cm]: " << ND_OnAxis_Sim_mu_start_p[2] << "\n";
    myfile << "ND_OnAxis_Sim_mu_start_p_total[GeV]: " << ND_OnAxis_Sim_mu_start_p_total << "\n\n";
    // ND_OnAxis_Sim_hadronic_hit
    for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
      double ND_RandomVtx_Sim_hadronic_hit_array[3] = {ND_RandomVtx_Sim_hadronic_hit[ihadronhit][0],ND_RandomVtx_Sim_hadronic_hit[ihadronhit][1],ND_RandomVtx_Sim_hadronic_hit[ihadronhit][2]};
      for (int i =0; i<3;i++)
      {
        ND_OnAxis_Sim_hadronic_hit_xyz.emplace_back((float)eff->getTranslations(ND_RandomVtx_Sim_hadronic_hit_array, ND_RandomVtx_Sim_mu_start_v, ND_OnAxis_Sim_mu_start_v, i));
      }
      ND_OnAxis_Sim_hadronic_hit.emplace_back(ND_OnAxis_Sim_hadronic_hit_xyz);
      ND_OnAxis_Sim_hadronic_hit_xyz.clear();
    }
    //
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //
    // Loop over ND_off_axis_pos_vec: random off_axis_pos or every ND_off_axis_pos_stepsize
    // Don't put it outside event loop to avoid looping over all events multiple times
    //

    int ND_off_axis_pos_counter = 0;
    for ( double i_ND_off_axis_pos : ND_off_axis_pos_vec )
    {
      ND_off_axis_pos_counter++;
      //
      // Loop over vtx x: random x or stepwise increased x
      // Don't put it outside event loop to avoid looping over all events multiple times
      //
      int vtx_vx_counter = 0;
      hadron_throw_result_vec_for_vtx_vx.clear(); // Need to initialize before loop over vtx x vec
      //
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //------------------------------------------------------------------------------
      //
      for ( double i_vtx_vx : ND_vtx_vx_vec )
      {

        vtx_vx_counter++;

        // Interpolate event neutrino production point (beam coordinate)
        decayZbeamCoord = gDecayZ->Eval( i_ND_off_axis_pos + i_vtx_vx - detRefBeamCoord[0] );

        // Calculate neutrino production point in detector coordinate
        decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation);
        decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation);
        // Set production point in unit: cm
        eff->setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.);

        if (verbose) std::cout << "nd off_axis x #" << ND_off_axis_pos_counter << ": " << i_ND_off_axis_pos << " m" << std::endl;
        myfile << "Nd Off-Axis x #" << ND_off_axis_pos_counter << "[cm]: " << i_ND_off_axis_pos << "\n";
        myfile << "Nd Off-Axis x in detector #" << vtx_vx_counter << "[cm]: " << i_vtx_vx << "\n";
        myfile << "Event position #" << vtx_vx_counter << "[cm]: " << i_ND_off_axis_pos + i_vtx_vx << "\n\n";

        //
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //
        // 3. ND to ND: translate from OnAxis to OffAxis
        //
        // ND_OffAxis_Unrotated_Sim_mu_start_v
        ND_OffAxis_Unrotated_Sim_mu_start_v[0] = ND_OnAxis_Sim_mu_start_v[0] + i_ND_off_axis_pos + i_vtx_vx;
        ND_OffAxis_Unrotated_Sim_mu_start_v[1] = eff->RemainUnchanged(ND_OnAxis_Sim_mu_start_v[1]);
        ND_OffAxis_Unrotated_Sim_mu_start_v[2] = eff->RemainUnchanged(ND_OnAxis_Sim_mu_start_v[2]);
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_vx[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_v[0] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_vy[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_v[1] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_vz[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_v[2] << "\n\n";
        // ND_OffAxis_Unrotated_Sim_mu_end_v
        // double geoEff::getTranslations(double v_bf[3], double vtx_bf[3], double vtx_af[3], int dim)
        for(int i=0; i<3; i++)
        {
          ND_OffAxis_Unrotated_Sim_mu_end_v[i] = eff->getTranslations(ND_OnAxis_Sim_mu_end_v, ND_OnAxis_Sim_mu_start_v, ND_OffAxis_Unrotated_Sim_mu_start_v, i);
        }
        myfile << "ND_OffAxis_Unrotated_Sim_mu_end_vx[cm]: " << ND_OffAxis_Unrotated_Sim_mu_end_v[0] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_end_vy[cm]: " << ND_OffAxis_Unrotated_Sim_mu_end_v[1] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_end_vz[cm]: " << ND_OffAxis_Unrotated_Sim_mu_end_v[2] << "\n";
        double ND_OffAxis_Unrotated_Sim_mu_v_end_start = eff->getDistance(ND_OffAxis_Unrotated_Sim_mu_end_v,ND_OffAxis_Unrotated_Sim_mu_start_v);
        myfile << "Distance between ND_OffAxis_Unrotated_Sim_mu_v_end_start[cm]: " << ND_OffAxis_Unrotated_Sim_mu_v_end_start << "\n\n";
        // ND_OffAxis_Unrotated_Sim_mu_start_p
        for(int i=0; i<3; i++)
        {
          ND_OffAxis_Unrotated_Sim_mu_start_p[i] = eff->RemainUnchanged(ND_OnAxis_Sim_mu_start_p[i]);
        }
        ND_OffAxis_Unrotated_Sim_mu_start_p_total = eff->getTotalMomentum(ND_OffAxis_Unrotated_Sim_mu_start_p);
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_px[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_p[0] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_py[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_p[1] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_pz[cm]: " << ND_OffAxis_Unrotated_Sim_mu_start_p[2] << "\n";
        myfile << "ND_OffAxis_Unrotated_Sim_mu_start_p_total[GeV]: " << ND_OffAxis_Unrotated_Sim_mu_start_p_total << "\n\n";
        // ND_OffAxis_Unrotated_Sim_hadronic_hit
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
          double ND_OnAxis_Sim_hadronic_hit_array[3] = {ND_OnAxis_Sim_hadronic_hit[ihadronhit][0],ND_OnAxis_Sim_hadronic_hit[ihadronhit][1],ND_OnAxis_Sim_hadronic_hit[ihadronhit][2]};
          for (int i =0; i<3;i++)
          {
            ND_OffAxis_Unrotated_Sim_hadronic_hit_xyz.emplace_back((float)eff->getTranslations(ND_OnAxis_Sim_hadronic_hit_array, ND_OnAxis_Sim_mu_start_v, ND_OffAxis_Unrotated_Sim_mu_start_v, i));
          }
          ND_OffAxis_Unrotated_Sim_hadronic_hit.emplace_back(ND_OffAxis_Unrotated_Sim_hadronic_hit_xyz);
          ND_OffAxis_Unrotated_Sim_hadronic_hit_xyz.clear();
        }
        //
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //
        // 4. ND: get after eigen rotated vectors for step 3
        // ND_OffAxis_Sim_mu_start_v
        ND_OffAxis_Sim_mu_start_v[0] = eff->RemainUnchanged(ND_OffAxis_Unrotated_Sim_mu_start_v[0]);
        ND_OffAxis_Sim_mu_start_v[1] = eff->RemainUnchanged(ND_OffAxis_Unrotated_Sim_mu_start_v[1]);
        ND_OffAxis_Sim_mu_start_v[2] = eff->RemainUnchanged(ND_OffAxis_Unrotated_Sim_mu_start_v[2]);
        myfile << "ND_OffAxis_Sim_mu_start_v[cm]: " << ND_OffAxis_Sim_mu_start_v[0] << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_v[cm]: " << ND_OffAxis_Sim_mu_start_v[1] << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_v[cm]: " << ND_OffAxis_Sim_mu_start_v[2] << "\n\n";
        // setVertex
        eff->setOffAxisVertex(ND_OffAxis_Sim_mu_start_v[0], ND_OffAxis_Sim_mu_start_v[1], ND_OffAxis_Sim_mu_start_v[2]);
        // ND_OffAxis_Sim_mu_end_v
        eff->setMuEndV(ND_OffAxis_Unrotated_Sim_mu_end_v[0],ND_OffAxis_Unrotated_Sim_mu_end_v[1],ND_OffAxis_Unrotated_Sim_mu_end_v[2]);
        for(int i=0; i<3; i++)
        {
          ND_OffAxis_Sim_mu_end_v[i] = eff->getOffAxisMuEndV(i);
        }
        myfile << "ND_OffAxis_Sim_mu_end_v[cm]: " << ND_OffAxis_Sim_mu_end_v[0] << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_v[cm]: " << ND_OffAxis_Sim_mu_end_v[1] << "\n";
        myfile << "ND_OffAxis_Sim_mu_end_v[cm]: " << ND_OffAxis_Sim_mu_end_v[2] << "\n";
        double ND_OffAxis_Sim_mu_v_end_start = eff->getDistance(ND_OffAxis_Sim_mu_end_v,ND_OffAxis_Sim_mu_start_v);
        myfile << "Distance between ND_OffAxis_Sim_mu_v_end_start[cm]: " << ND_OffAxis_Sim_mu_v_end_start << "\n\n";
        // ND_OffAxis_Sim_mu_start_p
        eff->setMuStartP(ND_OffAxis_Unrotated_Sim_mu_start_p[0],ND_OffAxis_Unrotated_Sim_mu_start_p[1],ND_OffAxis_Unrotated_Sim_mu_start_p[2]);
        for(int i=0; i<3; i++)
        {
          ND_OffAxis_Sim_mu_start_p[i] = eff->getOffAxisMuStartP(i);
        }
        ND_OffAxis_Sim_mu_start_p_total = eff->getTotalMomentum(ND_OffAxis_Sim_mu_start_p);
        myfile << "ND_OffAxis_Sim_mu_start_px[cm]: " << ND_OffAxis_Sim_mu_start_p[0] << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_py[cm]: " << ND_OffAxis_Sim_mu_start_p[1] << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_pz[cm]: " << ND_OffAxis_Sim_mu_start_p[2] << "\n";
        myfile << "ND_OffAxis_Sim_mu_start_p_total[GeV]: " << ND_OffAxis_Sim_mu_start_p_total << "\n\n";
        // ND_OffAxis_Sim_hadronic_hit
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
          eff->setHadronHitV(ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][0],ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][1],ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][2]);
          for (int i =0; i<3;i++)
          {
            ND_OffAxis_Sim_hadronic_hit_xyz.emplace_back(eff->getOffAxisHadronHitV(i));
          }
          ND_OffAxis_Sim_hadronic_hit.emplace_back(ND_OffAxis_Sim_hadronic_hit_xyz);
          ND_OffAxis_Sim_hadronic_hit_xyz.clear();
        }

        // Add output
        for (int ihadronhit = FD_Sim_n_hadronic_Edep_a-2; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++) {
          myfile << "Hit #" << ihadronhit <<"/" << FD_Sim_n_hadronic_Edep_a << "\n";
          myfile<<"FD_Sim_hadronic_hit_x[cm]: "<<FD_Sim_hadronic_hit_x_a->at(ihadronhit)<<"\n";
          myfile<<"FD_Sim_hadronic_hit_y[cm]: "<<FD_Sim_hadronic_hit_y_a->at(ihadronhit)<<"\n";
          myfile<<"FD_Sim_hadronic_hit_z[cm]: "<<FD_Sim_hadronic_hit_z_a->at(ihadronhit)<<"\n";
          double FD_Sim_hadronic_hit[3] = {FD_Sim_hadronic_hit_x_a->at(ihadronhit),FD_Sim_hadronic_hit_y_a->at(ihadronhit), FD_Sim_hadronic_hit_z_a->at(ihadronhit)};
          double FD_Sim_hadronic_hit_end_start = eff->getDistance(FD_Sim_hadronic_hit,FD_Sim_mu_start_v);
          myfile << "Distance between FD_Sim_hadronic_hit_end_start[cm]:" << FD_Sim_hadronic_hit_end_start <<"\n\n";

          for (int j = 0; j < 3; j++)
          {
            myfile << "ND_RandomVtx_Sim_hadronic_hit_["<<j<<"][cm]: "<< ND_RandomVtx_Sim_hadronic_hit[ihadronhit][j] << "\n";
          }
          double ND_RandomVtx_Sim_hadronic_hit_array[3] = {ND_RandomVtx_Sim_hadronic_hit[ihadronhit][0],ND_RandomVtx_Sim_hadronic_hit[ihadronhit][1],ND_RandomVtx_Sim_hadronic_hit[ihadronhit][2]};
          double ND_RandomVtx_Sim_hadronic_hit_end_start = eff->getDistance(ND_RandomVtx_Sim_hadronic_hit_array,ND_RandomVtx_Sim_mu_start_v);
          myfile << "Distance between ND_RandomVtx_Sim_hadronic_hit_end_start[cm]:" << ND_RandomVtx_Sim_hadronic_hit_end_start <<"\n\n";

          for (int j = 0; j < 3; j++)
          {
            myfile << "ND_OnAxis_Sim_hadronic_hit_["<<j<<"][cm]: "<< ND_OnAxis_Sim_hadronic_hit[ihadronhit][j] << "\n";
          }
          double ND_OnAxis_Sim_hadronic_hit_array[3] = {ND_OnAxis_Sim_hadronic_hit[ihadronhit][0],ND_OnAxis_Sim_hadronic_hit[ihadronhit][1],ND_OnAxis_Sim_hadronic_hit[ihadronhit][2]};
          double ND_OnAxis_Sim_hadronic_hit_end_start = eff->getDistance(ND_OnAxis_Sim_hadronic_hit_array,ND_OnAxis_Sim_mu_start_v);
          myfile << "Distance between ND_OnAxis_Sim_hadronic_hit_end_start[cm]:" << ND_OnAxis_Sim_hadronic_hit_end_start <<"\n\n";

          for (int j = 0; j < 3; j++)
          {
            myfile << "ND_OffAxis_Unrotated_Sim_hadronic_hit_["<<j<<"][cm]: "<< ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][j] << "\n";
          }
          double ND_OffAxis_Unrotated_Sim_hadronic_hit_array[3] = {ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][0],ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][1],ND_OffAxis_Unrotated_Sim_hadronic_hit[ihadronhit][2]};
          double ND_OffAxis_Unrotated_Sim_hadronic_hit_end_start = eff->getDistance(ND_OffAxis_Unrotated_Sim_hadronic_hit_array,ND_OffAxis_Unrotated_Sim_mu_start_v);
          myfile << "Distance between ND_OffAxis_Unrotated_Sim_hadronic_hit_end_start[cm]:" << ND_OffAxis_Unrotated_Sim_hadronic_hit_end_start <<"\n\n";

          for (int j = 0; j < 3; j++)
          {
            myfile << "ND_OffAxis_Sim_hadronic_hit_["<<j<<"][cm]: "<< ND_OffAxis_Sim_hadronic_hit[ihadronhit][j] << "\n";
          }
          double ND_OffAxis_Sim_hadronic_hit_array[3] = {ND_OffAxis_Sim_hadronic_hit[ihadronhit][0],ND_OffAxis_Sim_hadronic_hit[ihadronhit][1],ND_OffAxis_Sim_hadronic_hit[ihadronhit][2]};
          double ND_OffAxis_Sim_hadronic_hit_end_start = eff->getDistance(ND_OffAxis_Sim_hadronic_hit_array,ND_OffAxis_Sim_mu_start_v);
          myfile << "Distance between ND_OffAxis_Sim_hadronic_hit_end_start[cm]:" << ND_OffAxis_Sim_hadronic_hit_end_start <<"\n\n";
        }
        //
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //
        // 5. ND: generate random throws
        // Hardon hit positions
        HadronHitEdeps.clear();
        HadronHitPoss.clear();
        HadronHitEdeps.reserve(FD_Sim_n_hadronic_Edep_a);
        HadronHitPoss.reserve(FD_Sim_n_hadronic_Edep_a*3);
        // Set HadronHitPoss
        for ( int ihadronhit = 0; ihadronhit < FD_Sim_n_hadronic_Edep_a; ihadronhit++ ){
          for (int i =0; i<3;i++)
          {
            HadronHitPoss.emplace_back(ND_OffAxis_Sim_hadronic_hit[ihadronhit][i]);
          }
          HadronHitEdeps.emplace_back( FD_Sim_hadronic_hit_Edep_a2->at(ihadronhit) );
        }

        eff->setVertex(ND_OffAxis_Sim_mu_start_v[0], ND_OffAxis_Sim_mu_start_v[1], ND_OffAxis_Sim_mu_start_v[2]);
        eff->setHitSegEdeps(HadronHitEdeps);
        eff->setHitSegPoss(HadronHitPoss); // this is converted hadrom deposit pos in ND coordinate sys.
        //
        // Get hadron containment result after everything is set to ND coordinate sys
        // Do random throws regardless whether FD evt is contained in ND volume by setting a false flag
        hadron_throw_result = eff->getHadronContainmentThrows(false); // Every 64 throw results stored into a 64 bit unsigned int: 0101101...
        hadron_throw_result_vec_for_vtx_vx.emplace_back(hadron_throw_result);

        if (verbose) std::cout << "vtx x #" << vtx_vx_counter << ": " << i_vtx_vx << " cm, throw result[0][0][0]: " << hadron_throw_result[0][0][0] << std::endl;
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        //
        effTreeFD->Fill();
        ND_OffAxis_Unrotated_Sim_hadronic_hit.clear();
        ND_OffAxis_Sim_hadronic_hit.clear();
     } // end Loop over ND_vtx_vx_vec

      ND_Sim_hadron_throw_result.emplace_back(hadron_throw_result_vec_for_vtx_vx);

    }   // end Loop over ND_off_axis_pos_vec
    iwritten++;
    ND_RandomVtx_Sim_hadronic_hit.clear();
    ND_OnAxis_Sim_hadronic_hit.clear();

  } // end loop over events entries
  std::cout << "Written evts: " << iwritten << std::endl;
  myfile << "Written evts: " << iwritten << "\n";
  //
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //------------------------------------------------------------------------------
  //
  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff_ivy.root", "RECREATE");
  ThrowsFD->Write();
  effTreeFD->Write();
  hist_vetoEnergyFD->Write();

  myfile.close();
  outFile->Close();
} // end main
