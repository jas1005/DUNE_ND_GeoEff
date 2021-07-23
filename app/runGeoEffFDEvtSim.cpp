// Library methods
#include "geoEff.h"

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>

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

// Include customized functions and constants
#include "Helpers.h"
#include "TROOT.h"
#include <vector> // Need this for generate dictionary for nested vectors

int main(){

  //
  // Declare variables used in this program
  //

  int nentries = 0;                 // Total input events
  int iwritten = 0;                 // Output event counter
  vector<float> HadronHitEdeps;
  vector<float> HadronHitPoss;
  double RnOffAxisPoint;            // Random off-axis ND x position for each event
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
  double Sim_mu_start_vx;
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
  vector<float> *Sim_hadronic_hit_x_a     = 0;
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

  //
  // A tree to store lepton info (for NN training)
  // and result of hadron containment after applying transformations
  //

  // Lepton info
  double b_Sim_mu_start_vx;   // Do not use float!
  double b_Sim_mu_start_vy;
  double b_Sim_mu_start_vz;
  double b_Sim_mu_end_vx;
  double b_Sim_mu_end_vy;
  double b_Sim_mu_end_vz;
  double b_Sim_mu_start_px;
  double b_Sim_mu_start_py;
  double b_Sim_mu_start_pz;
  double b_Sim_mu_end_px;
  double b_Sim_mu_end_py;
  double b_Sim_mu_end_pz;

  // Result of hadron containment is stored in nested vectors, need to generate dictionary
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*uint64_t*"); // First remove old dictionary
  // Generate new dictionary
  gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t> > >", "vector");
  vector<vector<vector<uint64_t> > > b_Sim_hadron_throw_result;

  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  effTreeFD->Branch("Sim_hadron_throw_result", &b_Sim_hadron_throw_result);
  effTreeFD->Branch("Sim_mu_start_vx",         &b_Sim_mu_start_vx, "Sim_mu_start_vx/D");
  effTreeFD->Branch("Sim_mu_start_vy",         &b_Sim_mu_start_vy, "Sim_mu_start_vy/D");
  effTreeFD->Branch("Sim_mu_start_vz",         &b_Sim_mu_start_vz, "Sim_mu_start_vz/D");
  effTreeFD->Branch("Sim_mu_end_vx",           &b_Sim_mu_end_vx,   "Sim_mu_end_vx/D");
  effTreeFD->Branch("Sim_mu_end_vy",           &b_Sim_mu_end_vy,   "Sim_mu_end_vy/D");
  effTreeFD->Branch("Sim_mu_end_vz",           &b_Sim_mu_end_vz,   "Sim_mu_end_vz/D");
  effTreeFD->Branch("Sim_mu_start_px",         &b_Sim_mu_start_px, "Sim_mu_start_px/D");
  effTreeFD->Branch("Sim_mu_start_py",         &b_Sim_mu_start_py, "Sim_mu_start_py/D");
  effTreeFD->Branch("Sim_mu_start_pz",         &b_Sim_mu_start_pz, "Sim_mu_start_pz/D");
  effTreeFD->Branch("Sim_mu_end_px",           &b_Sim_mu_end_px,   "Sim_mu_end_px/D");
  effTreeFD->Branch("Sim_mu_end_py",           &b_Sim_mu_end_py,   "Sim_mu_end_py/D");
  effTreeFD->Branch("Sim_mu_end_pz",           &b_Sim_mu_end_pz,   "Sim_mu_end_pz/D");


  //
  // A separate tree to store translations and rotations of throws
  // which will be applied to leptons before NN training
  //

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
  double beamLineRotation = -0.101;           // unit: rad,
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
  std::cout << "Tot evts: " << nentries << std::endl;
  for ( int ientry = 0; ientry < nentries; ientry++ ) {

    t->GetEntry(ientry);
    std::cout << "Looking at entry " << ientry << ", run: " << Run << ", subrun: " << SubRun << ", event: " << Event << std::endl;

    // Skip events without hadronic deposits or muon
    if ( Sim_nMu == 0 || Sim_n_hadronic_Edep_a == 0 ) continue;

    HadronHitEdeps.clear();
    HadronHitPoss.clear();
    HadronHitEdeps.reserve(Sim_n_hadronic_Edep_a);
    HadronHitPoss.reserve(Sim_n_hadronic_Edep_a*3);

    //
    // We want to put the FD event as if it's in ND, which can be off-axis
    // Need a random generated off-axis position for each event
    //

    // Initialize random number generator
    TRandom3 *r3 = new TRandom3();
    // Set the seed (required to avoid repeated random numbers in each sequence)
    r3->SetSeed(0);
    RnOffAxisPoint = r3->Uniform(OffAxisPoints[0], OffAxisPoints[13]);

    // Mean neutrino production point (beam coordinate) on z axis as a function of ND off-axis position
    TGraph* gDecayZ = new TGraph(14, OffAxisPoints, meanPDPZ);
    // Interpolate event neutrino production point (beam coordinate)
    decayZbeamCoord = gDecayZ->Eval( RnOffAxisPoint - detRefBeamCoord[0] );

    // Calculate neutrino production point in detector coordinate
    decayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0];
    decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation);
    decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation);
    // Set production point using cm unit
    eff->setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.);

    //
    // Event vertex pos in unit: cm
    // Here we use lepton track start position as vertex position, eventually will need to change to true evt vtx for NC events
    // These vertex vx, vy, vz are in FD coordinate sys?, need to translate to ND box?
    // vx can be randomized in ND? Every 5cm
    //
    
    eff->setVertex(Sim_mu_start_vx/10., Sim_mu_start_vy/10., Sim_mu_start_vz/10.);

    //
    // Renew throws every 100th written event to save file size, i.e., if N = 128,
    // for written evt 0-99:   same 128 transformations for each event,
    // for written evt 99-200: same but renewed 128 transformations for each event
    // These transformations will be applied to leptons in the event, so need to keep track of iwritten
    //

    if ( iwritten % 100 == 0 ) {
      // Produce N throws defined at setNthrows(N)
      eff->throwTransforms();
      throwVtxY.clear();
      throwVtxZ.clear();
      throwRot.clear();
      throwVtxY = eff->getCurrentThrowTranslationsY();
      throwVtxZ = eff->getCurrentThrowTranslationsZ();
      throwRot  = eff->getCurrentThrowRotations();
      ThrowsFD->Fill();
    }

    // Read E and position info for each sim hadron hit
    for ( int ihadronhit = 0; ihadronhit < Sim_n_hadronic_Edep_a; ihadronhit++ ){

      HadronHitEdeps.emplace_back(Sim_hadronic_hit_Edep_a2->at(ihadronhit));
      HadronHitPoss.emplace_back(Sim_hadronic_hit_x_a->at(ihadronhit)/10.);
      HadronHitPoss.emplace_back(Sim_hadronic_hit_y_a->at(ihadronhit)/10.);
      HadronHitPoss.emplace_back(Sim_hadronic_hit_z_a->at(ihadronhit)/10.);

      //
      // In the future consider the ND size for FD and where to put ND volume in FD
      //

      // A rough active detector dimensions for ND:
      // eff->setActiveX(collarLo[0]-30, collarHi[0]+30);
      // eff->setActiveY(collarLo[1]-30, collarHi[1]+30);
      // eff->setActiveZ(collarLo[2]-30, collarHi[2]+30);
    }

    eff->setHitSegEdeps(HadronHitEdeps);
    eff->setHitSegPoss(HadronHitPoss);

    b_Sim_hadron_throw_result = eff->getHadronContainmentThrows();
    std::cout << "Throw result, 0,0,1: " << b_Sim_hadron_throw_result[0][0][1] << std::endl;
    b_Sim_mu_start_vx = Sim_mu_start_vx;
    b_Sim_mu_start_vy = Sim_mu_start_vy;
    b_Sim_mu_start_vz = Sim_mu_start_vz;
    b_Sim_mu_end_vx   = Sim_mu_end_vx;
    b_Sim_mu_end_vy   = Sim_mu_end_vy;
    b_Sim_mu_end_vz   = Sim_mu_end_vz;
    b_Sim_mu_start_px = Sim_mu_start_px;
    b_Sim_mu_start_py = Sim_mu_start_py;
    b_Sim_mu_start_pz = Sim_mu_start_pz;
    b_Sim_mu_end_px   = Sim_mu_end_px;
    b_Sim_mu_end_py   = Sim_mu_end_py;
    b_Sim_mu_end_pz   = Sim_mu_end_pz;

    effTreeFD->Fill();
    iwritten++;

  } // end loop over events

  std::cout << "Written evts: " << iwritten << std::endl;

  // Write trees
  TFile * outFile = new TFile("Output_FDGeoEff.root", "RECREATE");
  effTreeFD->Write();
  ThrowsFD->Write();

  outFile->Close();

}   // end main
