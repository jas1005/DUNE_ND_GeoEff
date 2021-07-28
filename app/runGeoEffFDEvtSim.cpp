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
  // Declare variables used in this program
  //

  int nentries = 0;                 // Total input events
  int iwritten = 0;                 // Output event counter
  vector<float> HadronHitEdeps;
  vector<float> HadronHitPoss;
  double Random_OffAxisPoint;       // Random off-axis ND x position for each event
  float decayZbeamCoord;
  float decayXdetCoord;
  float decayYdetCoord;
  float decayZdetCoord;

  //
  // Branches to be read from n-tuple produced from FD MC
  //

  Int_t Run;
  Int_t SubRun;
  Int_t Event;
  Int_t Sim_nMu;
  double Sim_mu_start_vx; // cm?
  double Sim_mu_start_vy;
  double Sim_mu_start_vz;
  double Sim_mu_end_vx;
  double Sim_mu_end_vy;
  double Sim_mu_end_vz;
  double Sim_mu_start_px;
  double Sim_mu_start_py;
  double Sim_mu_start_pz;
  double Sim_mu_end_px;
  double Sim_mu_end_py;
  double Sim_mu_end_pz;
  Int_t Sim_n_hadronic_Edep_a;
  vector<float> *Sim_hadronic_hit_Edep_a2 = 0; // Need initialize 0 here to avoid error
  vector<float> *Sim_hadronic_hit_x_a     = 0; // Same as mu pos unit: cm?
  vector<float> *Sim_hadronic_hit_y_a     = 0;
  vector<float> *Sim_hadronic_hit_z_a     = 0;

  // Read ntuple from FD MC
  TChain *t = new TChain("MyEnergyAnalysis/MyTree");
  t->Add("/dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis/myntuple.root"); // ntuple path on FNAL dunegpvm machine

  t->SetBranchAddress("Run",    &Run);
  t->SetBranchAddress("SubRun", &SubRun);
  t->SetBranchAddress("Event",  &Event);
  t->SetBranchAddress("Sim_nMu",          &Sim_nMu);
  t->SetBranchAddress("Sim_mu_start_vx",  &Sim_mu_start_vx);
  t->SetBranchAddress("Sim_mu_start_vy",  &Sim_mu_start_vy);
  t->SetBranchAddress("Sim_mu_start_vz",  &Sim_mu_start_vz);
  t->SetBranchAddress("Sim_mu_end_vx",    &Sim_mu_end_vx);
  t->SetBranchAddress("Sim_mu_end_vy",    &Sim_mu_end_vy);
  t->SetBranchAddress("Sim_mu_end_vz",    &Sim_mu_end_vz);
  t->SetBranchAddress("Sim_mu_start_px",  &Sim_mu_start_px);
  t->SetBranchAddress("Sim_mu_start_py",  &Sim_mu_start_py);
  t->SetBranchAddress("Sim_mu_start_pz",  &Sim_mu_start_pz);
  t->SetBranchAddress("Sim_mu_end_px",    &Sim_mu_end_px);
  t->SetBranchAddress("Sim_mu_end_py",    &Sim_mu_end_py);
  t->SetBranchAddress("Sim_mu_end_pz",    &Sim_mu_end_pz);
  t->SetBranchAddress("Sim_n_hadronic_Edep_a",    &Sim_n_hadronic_Edep_a);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_a2", &Sim_hadronic_hit_Edep_a2);
  t->SetBranchAddress("Sim_hadronic_hit_x_a",     &Sim_hadronic_hit_x_a);
  t->SetBranchAddress("Sim_hadronic_hit_y_a",     &Sim_hadronic_hit_y_a);
  t->SetBranchAddress("Sim_hadronic_hit_z_a",     &Sim_hadronic_hit_z_a);

  // A tree to store lepton info (for NN training)
  // and result of hadron containment after applying transformations

  // Lepton info
  vector<double> b_Sim_mu_start_vx; // Vector corresponds to randomized/stepwise vtx x
  vector<double> b_Sim_mu_end_vx;
  double b_Sim_mu_start_vy;         // Do not use float!
  double b_Sim_mu_start_vz;
  double b_Sim_mu_end_vy;
  double b_Sim_mu_end_vz;
  double b_Sim_mu_start_px;
  double b_Sim_mu_start_py;
  double b_Sim_mu_start_pz;
  double b_Sim_mu_end_px;
  double b_Sim_mu_end_py;
  double b_Sim_mu_end_pz;

  // Result of hadron containment is stored in nested vectors, need to generate dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*uint64_t*"); // Remove old dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*vector*uint64_t*");
  // Generate new dictionary
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<vector<uint64_t> > > >", "vector");
  vector<vector<vector<uint64_t> > > hadron_throw_result;
  vector<vector<vector<vector<uint64_t> > > > b_Sim_hadron_throw_result; // One more vector corresponds to randomized/stepwise vtx x

  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  effTreeFD->Branch("Sim_hadron_throw_result", &b_Sim_hadron_throw_result);
  effTreeFD->Branch("Sim_mu_start_vx",         &b_Sim_mu_start_vx); // vector
  effTreeFD->Branch("Sim_mu_start_vy",         &b_Sim_mu_start_vy,       "Sim_mu_start_vy/D");
  effTreeFD->Branch("Sim_mu_start_vz",         &b_Sim_mu_start_vz,       "Sim_mu_start_vz/D");
  effTreeFD->Branch("Sim_mu_end_vx",           &b_Sim_mu_end_vx);   // vector
  effTreeFD->Branch("Sim_mu_end_vy",           &b_Sim_mu_end_vy,         "Sim_mu_end_vy/D");
  effTreeFD->Branch("Sim_mu_end_vz",           &b_Sim_mu_end_vz,         "Sim_mu_end_vz/D");
  effTreeFD->Branch("Sim_mu_start_px",         &b_Sim_mu_start_px,       "Sim_mu_start_px/D");
  effTreeFD->Branch("Sim_mu_start_py",         &b_Sim_mu_start_py,       "Sim_mu_start_py/D");
  effTreeFD->Branch("Sim_mu_start_pz",         &b_Sim_mu_start_pz,       "Sim_mu_start_pz/D");
  effTreeFD->Branch("Sim_mu_end_px",           &b_Sim_mu_end_px,         "Sim_mu_end_px/D");
  effTreeFD->Branch("Sim_mu_end_py",           &b_Sim_mu_end_py,         "Sim_mu_end_py/D");
  effTreeFD->Branch("Sim_mu_end_pz",           &b_Sim_mu_end_pz,         "Sim_mu_end_pz/D");

  // A separate tree to store translations and rotations of throws
  // which will be applied to leptons before NN training

  vector<float> throwVtxY;
  vector<float> throwVtxZ;
  vector<float> throwRot;

  TTree * ThrowsFD = new TTree("ThrowsFD", "FD Throws");
  ThrowsFD->Branch("throwVtxY", &throwVtxY);
  ThrowsFD->Branch("throwVtxZ", &throwVtxZ);
  ThrowsFD->Branch("throwRot",  &throwRot);

  //
  // Get beam parameters: eventually need to read it from XML file: which XML for FD? same?
  //

  // What is this? rotate in ND y-z plane?
  double beamLineRotation = -0.101;           // unit: rad
  // Coordinate transformation, units: meters
  double beamRefDetCoord[3] = {0., 0., 15.5}; // (0, 0, 0) is ND detector origin
  double detRefBeamCoord[3] = {0., 0., 574.}; // (0, 0, 0) is beam origin

  //
  // Initialize geometric efficiency module
  //

  geoEff * eff = new geoEff(314, false); // set verbose to true for debug
  eff->setNthrows(128);
  // Rotate w.r.t. neutrino direction, rather than fixed beam direction
  eff->setUseFixedBeamDir(false);

  // 30 cm veto
  eff->setVetoSizes(vector<float>(1, 30.));
  // 30 MeV
  eff->setVetoEnergyThresholds(vector<float>(1, 30.));

  // Active detector dimensions for ND
  eff->setActiveX(collarLo[0]-30, collarHi[0]+30);
  eff->setActiveY(collarLo[1]-30, collarHi[1]+30);
  eff->setActiveZ(collarLo[2]-30, collarHi[2]+30);

  // Range for translation throws. Use full active volume but fix X.
  eff->setRangeX(-1, -1);
  eff->setRandomizeX(false);
  eff->setRangeY(collarLo[1]-30, collarHi[1]+30);
  eff->setRangeZ(collarLo[2]-30, collarHi[2]+30);

  // Set offset between MC coordinate system and volumes defined above: FD need changes?
  eff->setOffsetX(offset[0]);
  eff->setOffsetY(offset[1]);
  eff->setOffsetZ(offset[2]);

  //
  // Loop over FD events
  //

  nentries = t->GetEntries();
  if (verbose) std::cout << "Tot evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < nentries; ientry++ ) {

    t->GetEntry(ientry);
    if (verbose) std::cout << "Looking at entry " << ientry << ", run: " << Run << ", subrun: " << SubRun << ", event: " << Event << std::endl;

    // Skip events without hadronic deposits or muon
    if ( Sim_nMu == 0 || Sim_n_hadronic_Edep_a == 0 ) continue;

    // Renew throws every 100th written event to save file size, i.e., if N = 128,
    // for written evt 0-99:   same 128 transformations for each event,
    // for written evt 99-200: same but renewed 128 transformations for each event
    // These transformations will be applied to leptons in the event, so need to keep track of iwritten
    if ( iwritten % 100 == 0 ) {

      // Produce N throws defined at setNthrows(N)
      eff->throwTransforms(); // Doesn't depend on event vtx
      throwVtxY.clear();
      throwVtxZ.clear();
      throwRot.clear();
      throwVtxY = eff->getCurrentThrowTranslationsY();
      throwVtxZ = eff->getCurrentThrowTranslationsZ();
      throwRot  = eff->getCurrentThrowRotations();
      ThrowsFD->Fill();
    }

    //
    // Rotate event around x-axis +/-0.101*2? to bring it into ND system before doing ND constraint
    //

    //
    // The virtual ND can be at a random off-axis position for each event
    //

    // Initialize random number generator
    TRandom3 *r3_OffAxisPoint = new TRandom3();
    // Set the seed (required to avoid repeated random numbers in each sequence)
    r3_OffAxisPoint->SetSeed(0);
    Random_OffAxisPoint = r3_OffAxisPoint->Uniform(OffAxisPoints[0], OffAxisPoints[13]);
    if (verbose) std::cout << "random OffAxisPoint [meters]: " << Random_OffAxisPoint << std::endl;

    // Mean neutrino production point (beam coordinate) on z axis as a function of ND off-axis position
    TGraph* gDecayZ = new TGraph(14, OffAxisPoints, meanPDPZ);
    // Interpolate event neutrino production point (beam coordinate)
    decayZbeamCoord = gDecayZ->Eval( Random_OffAxisPoint - detRefBeamCoord[0] );

    // Calculate neutrino production point in detector coordinate
    decayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0];
    decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation);
    decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation);
    // Set production point in unit: cm
    eff->setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.);

    //
    // Use relative coordinate for FD event (relative to muon start)
    // Set muon start to (x, 0, 0), eventually should use (x, 0, 0) for event vtx
    //

    // Branches that are not affected by vtx x
    b_Sim_mu_start_vy = 0.; // y and z will be randomly translated later, so just use 0
    b_Sim_mu_start_vz = 0.;
    b_Sim_mu_end_vy   = Sim_mu_end_vy - Sim_mu_start_vy; // + 0 omitted
    b_Sim_mu_end_vz   = Sim_mu_end_vz - Sim_mu_start_vz; // + 0 omitted
    b_Sim_mu_start_px = Sim_mu_start_px; // momentum is not affected by x
    b_Sim_mu_start_py = Sim_mu_start_py;
    b_Sim_mu_start_pz = Sim_mu_start_pz;
    b_Sim_mu_end_px   = Sim_mu_end_px;
    b_Sim_mu_end_py   = Sim_mu_end_py;
    b_Sim_mu_end_pz   = Sim_mu_end_pz;

    // Two options to decide x, user can switch via random_vtx_vx
    //
    // Option 1: randomize x once for each event
    // Option 2: 10cm steps along ND local width in x: -2m to 2m, 40 steps (larger output file)

    TRandom3 *r3_vtx_x = new TRandom3();
    r3_vtx_x->SetSeed(0);
    vtx_vx_vec.at(0) = r3_vtx_x->Uniform(ND_local_x_min, ND_local_x_max);
    if (verbose) std::cout << "vtx_vx_vec size: "<< vtx_vx_vec.size() << ", random vtx_x [cm]: " << vtx_vx_vec.at(0) << std::endl;

    b_Sim_mu_start_vx.clear();
    b_Sim_mu_end_vx.clear();

    // Loop over vtx x: random x (ix=0) or stepwise increased x (ix>0)
    // Can't put it outside event loop as that will loop over all events more than once
    for ( double i_vtx_vx : vtx_vx_vec ) {

      // Here in b_Sim_mu_start_vx and b_Sim_mu_end_vx, every (ix_max-ix_min) chunk is for one event
      // Output branch entry will be (ix_max-ix_min) * written evts
      b_Sim_mu_start_vx.emplace_back( i_vtx_vx );
      b_Sim_mu_end_vx.emplace_back( Sim_mu_end_vx - Sim_mu_start_vx + i_vtx_vx ); // w.r.t. mu start x
      // Event vertex pos in unit: cm
      eff->setVertex( i_vtx_vx, b_Sim_mu_start_vy, b_Sim_mu_start_vz );

      HadronHitEdeps.clear();
      HadronHitPoss.clear();
      HadronHitEdeps.reserve(Sim_n_hadronic_Edep_a);
      HadronHitPoss.reserve(Sim_n_hadronic_Edep_a*3);

      // Need to loop deposits many times as the hadron containment below need vtx x first
      for ( int ihadronhit = 0; ihadronhit < Sim_n_hadronic_Edep_a; ihadronhit++ ){
        // Relative to muon start position
        HadronHitPoss.emplace_back( Sim_hadronic_hit_x_a->at(ihadronhit) - Sim_mu_start_vx + i_vtx_vx );
        HadronHitPoss.emplace_back( Sim_hadronic_hit_y_a->at(ihadronhit) - Sim_mu_start_vy ); // + 0 omitted
        HadronHitPoss.emplace_back( Sim_hadronic_hit_z_a->at(ihadronhit) - Sim_mu_start_vz ); // + 0 omitted
        HadronHitEdeps.emplace_back( Sim_hadronic_hit_Edep_a2->at(ihadronhit) );
      }

      eff->setHitSegEdeps(HadronHitEdeps);
      eff->setHitSegPoss(HadronHitPoss);

      // Get hadron containment result after setVertex and setDecayPos
      hadron_throw_result = eff->getHadronContainmentThrows();     // Every 64 throw results is stored into a 64 bit unsigned int: 0101101...
      b_Sim_hadron_throw_result.emplace_back(hadron_throw_result); // This is now a class of vector<vector<vector<vector<uint64_t> > > >

      if (verbose) std::cout << "vtx x: " << i_vtx_vx << std::endl;
      if (verbose) std::cout << "Throw result, 0,0,0: " << hadron_throw_result[0][0][0] << std::endl;

    } // end Loop over vtx x

    effTreeFD->Fill();
    iwritten++;

  } // end loop over events

  if (verbose) std::cout << "Written evts: " << iwritten << std::endl;

  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff.root", "RECREATE");
  ThrowsFD->Write();
  effTreeFD->Write();

  outFile->Close();

}   // end main
