#include "geoEff.h" // All methods are in DUNE_ND_GeoEff/src/geoEff.cpp
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
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

#include "Helpers.h"   // Linear interpolate average neutrino decay position as a function of off axis x
#include "TROOT.h"
#include <vector>
//#ifdef __MAKECINT__
//#pragma link C++ class std::vector<std::vector<std::vector<uint64_t>>>+;
//#endif

int main(){

  //
  // Declare variables used in this program
  //

  unsigned long int nentries;
  vector<float> HadronHitEdeps;
  vector<float> HadronHitPoss;
  float decayZbeamCoord; // in cm
  float decayXdetCoord;
  float decayYdetCoord;
  float decayZdetCoord;

  //
  // Branches to be read from Ntuple produced from FD MC
  //
  Int_t Run;
  Int_t SubRun;
  Int_t Event;
  Int_t Sim_nMu;
  float Sim_mu_start_vx;
  float Sim_mu_start_vy;
  float Sim_mu_start_vz;
  Int_t Sim_n_hadronic_Edep_a;
  vector<float> * Sim_hadronic_hit_Edep_a2;
  vector<float> * Sim_hadronic_hit_x_a;
  vector<float> * Sim_hadronic_hit_y_a;
  vector<float> * Sim_hadronic_hit_z_a;

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
  t->SetBranchAddress("Sim_n_hadronic_Edep_a",    &Sim_n_hadronic_Edep_a);
  t->SetBranchAddress("Sim_hadronic_hit_Edep_a2", &Sim_hadronic_hit_Edep_a2);
  t->SetBranchAddress("Sim_hadronic_hit_x_a",     &Sim_hadronic_hit_x_a);
	t->SetBranchAddress("Sim_hadronic_hit_y_a",     &Sim_hadronic_hit_y_a);
  t->SetBranchAddress("Sim_hadronic_hit_z_a",     &Sim_hadronic_hit_z_a);

  //
  // Branches to be created and write to new tree effTreeFD below
  //

  //vector<vector<vector<uint64_t>>> HadronContainThrowResult;
  TTree * effTreeFD = new TTree("effTreeFD", "FD eff Tree");
  //effTreeFD->Branch("HadronContainThrowResult", &HadronContainThrowResult);

  //
  // Branches to be created and write to new tree ThrowsFD below
  //

  vector<float> throwVtxY;
  vector<float> throwVtxZ;
  vector<float> throwRot;

  // Separate tree to store translations and rotations of throws: for muon NN training?
  TTree * ThrowsFD = new TTree("ThrowsFD", "FD Throws");
  ThrowsFD->Branch("throwVtxY", &throwVtxY);
  ThrowsFD->Branch("throwVtxZ", &throwVtxZ);
  ThrowsFD->Branch("throwRot",  &throwRot);

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

  // Active detector dimensions ND?
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
  // Get beam parameters: eventually need to read it from XML file
  //

  double beamLineRotation = -0.101; // unit: rad
  double beamRefDetCoord[3] = {0., 0.05387, 6.66};
  double detRefBeamCoord[3] = {0., 0.,      562.1179};

  nentries = t->GetEntries();
  std::cout << "Tot evts: " << nentries << std::endl;
  // Loop over FD events
  for ( unsigned long int ientry = 0; ientry < nentries; ientry++ ) {
    t->GetEntry(ientry);
    std::cout << "Looking at entry " << ientry << ", run: " << Run << ", subrun: " << SubRun << ", event: " << Event << std::endl;

    // Skip events without hadronic deposits or muon
    if ( Sim_nMu == 0 || Sim_n_hadronic_Edep_a == 0 ) continue;

    HadronHitEdeps.clear();
    HadronHitPoss.clear();

    HadronHitEdeps.reserve(Sim_n_hadronic_Edep_a);
    HadronHitPoss.reserve(Sim_n_hadronic_Edep_a*3);

    //
    // Use vertex to determine mean decay point
    //

    // Here we use lepton track start position as vertex position, eventually will need to change to true evt vtx for NC events
    TGraph* gDecayZ = new TGraph(14, OffAxisPoints, meanPDPZ);
    decayZbeamCoord = gDecayZ->Eval(Sim_mu_start_vx/1000. - detRefBeamCoord[0])*100.; // in cm
    decayXdetCoord = -1*detRefBeamCoord[0]*100 + beamRefDetCoord[0]*100;
    decayYdetCoord = (-1*detRefBeamCoord[1]*100.*cos(beamLineRotation)) + (-1*detRefBeamCoord[2]*100. + decayZbeamCoord)*sin(beamLineRotation) + beamRefDetCoord[1]*100.;
    decayZdetCoord = (-1*detRefBeamCoord[2]*100. + decayZbeamCoord)*cos(beamLineRotation) - (-1*detRefBeamCoord[1]*100.*sin(beamLineRotation)) + beamRefDetCoord[2]*100.;

    eff->setDecayPos(decayXdetCoord, decayYdetCoord, decayZdetCoord);
    eff->setVertex(Sim_mu_start_vx/10., Sim_mu_start_vy/10., Sim_mu_start_vz/10.);

    if ( ientry%100 == 0 ) {
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
    }

    eff->setHitSegEdeps(HadronHitEdeps);
    eff->setHitSegPoss(HadronHitPoss);

    vector<vector<vector<uint64_t>>> HadronContainThrowResultList = eff->getHadronContainmentThrows();
    std::cout << "Throw result, 0,0,1: " << HadronContainThrowResultList[0][0][1] << std::endl;

    effTreeFD->Fill();

  } // end loop over events

  // Write trees
  TFile * outFile = new TFile("./Output_FDGeoEff.root", "RECREATE");
  effTreeFD->Write();
  ThrowsFD->Write();

  outFile->Close();

}   // end main
