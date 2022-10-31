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
#include <TGraph2D.h>
#include <TRandom.h>
#include <TF2.h>

// C++ includes
#include <iostream>
#include <iomanip>
using namespace std;
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector> // Need this for generate dictionary for nested vectors
using namespace std;

// Include customized functions and constants
#include "Helpers.h"

float vSize = 30.;


void ReadHadronHitNtuple()
{
  gROOT->Reset();


  // Input FDroot file
  TString FileIn = "/home/fyguo/NDEff/DUNE_ND_GeoEff/bin/Test_Output_FDGeoEff_ivy.root";
  //
  // Read branch from input trees
  //
  // Read effValues
  TChain *t_effValues = new TChain("effValues");
  t_effValues->Add(FileIn.Data());
  Int_t iwritten;
  Double_t ND_OffAxis_pos;
  Double_t ND_LAr_pos;
  Double_t ND_OffAxis_eff;
  t_effValues->SetBranchAddress("iwritten",         &iwritten);
  t_effValues->SetBranchAddress("ND_OffAxis_pos",   &ND_OffAxis_pos);
  t_effValues->SetBranchAddress("ND_LAr_pos",       &ND_LAr_pos);
  t_effValues->SetBranchAddress("ND_OffAxis_eff",   &ND_OffAxis_eff);
  //hardon info
  TChain *t_effTree = new TChain("effTreeND");
  t_effTree->Add(FileIn.Data());
  int ND_Sim_n_hadronic_Edep_a;
  vector<float> *HadronHitEdeps =0; // Hadron hit segment energy deposits [MeV]
  vector<vector<float>> *CurrentThrowDepsX = 0; // Coordinates of hadron hits X after random throws
  vector<vector<float>> *CurrentThrowDepsY =0; // Coordinates of hadron hits Y after random throws
  vector<vector<float>> *CurrentThrowDepsZ = 0; // Coordinates of hadron hits Z after random throws
  vector<float> *CurrentThrowVetoE = 0;
  vector<float> *CurrentThrowTotE = 0;
  vector<vector<float>> *ND_OffAxis_Sim_hadronic_hit_xyz=0; // coordinates of hadron hits before random throws

  t_effTree->SetBranchAddress("ND_Sim_n_hadronic_Edep_a",         &ND_Sim_n_hadronic_Edep_a);
  t_effTree->SetBranchAddress("CurrentThrowDepsX",         &CurrentThrowDepsX);
  t_effTree->SetBranchAddress("CurrentThrowDepsY",         &CurrentThrowDepsY);
  t_effTree->SetBranchAddress("CurrentThrowDepsZ",         &CurrentThrowDepsZ);
  t_effTree->SetBranchAddress("CurrentThrowVetoE",         &CurrentThrowVetoE);
  t_effTree->SetBranchAddress("CurrentThrowTotE",          &CurrentThrowTotE);
  t_effTree->SetBranchAddress("HadronHitEdeps",            &HadronHitEdeps);
  t_effTree->SetBranchAddress("ND_OffAxis_Sim_hadronic_hit_xyz",            &ND_OffAxis_Sim_hadronic_hit_xyz);

  // Read PosVec
  TChain *t_PosVec = new TChain("PosVec");
  t_PosVec->Add(FileIn.Data());
  vector<Double_t> *ND_OffAxis_pos_vec = 0; // unit: cm, ND off-axis choices for each FD evt
  vector<Double_t> *ND_LAr_pos_vec = 0;       // unit: cm, vtx x choices for each FD evt in ND volume
  vector<Int_t> *iwritten_vec = 0;
  TBranch *b_ND_OffAxis_pos_vec = 0;
  TBranch *b_ND_LAr_pos_vec = 0;
  TBranch *b_iwritten_vec = 0;
  t_PosVec->SetBranchAddress("ND_OffAxis_pos_vec",         &ND_OffAxis_pos_vec , &b_ND_OffAxis_pos_vec);
  t_PosVec->SetBranchAddress("ND_LAr_pos_vec",             &ND_LAr_pos_vec,      &b_ND_LAr_pos_vec);
  t_PosVec->SetBranchAddress("iwritten_vec",               &iwritten_vec,        &b_iwritten_vec);



  Long64_t tentry = t_PosVec->LoadTree(0);
  b_ND_OffAxis_pos_vec->GetEntry(tentry);
  b_ND_LAr_pos_vec->GetEntry(tentry);
  b_iwritten_vec->GetEntry(tentry);

  Int_t ND_OffAxis_pos_vec_size = ND_OffAxis_pos_vec->size();
  Int_t ND_LAr_pos_vec_size = ND_LAr_pos_vec->size();
  Int_t iwritten_vec_size = iwritten_vec->size();
  Int_t tot_size = ND_OffAxis_pos_vec_size*ND_LAr_pos_vec_size;
  Int_t hadronhit_n_plots = tot_size * N_throws;

  Int_t nentries = t_effValues->GetEntries();

  // Output
  TFile * outFile = new TFile("HadronHitPlots_ivy.root", "RECREATE");
  TDirectory *IPhadronhit =(TDirectory*)outFile->mkdir("hadron hits");//create a new folder in the root file
  TDirectory *IPhadronhit_offaxis =(TDirectory*)outFile->mkdir("OffAxis hadron hits");//create a new folder in the root file


  // Canvas

  TCanvas** c_hadronhit = new TCanvas*[hadronhit_n_plots];
  TH2F** h_hadronhit_xy = new TH2F*[hadronhit_n_plots];
  TH2F** h_hadronhit_zx = new TH2F*[hadronhit_n_plots];
  TH2F** h_hadronhit_zy = new TH2F*[hadronhit_n_plots];


  TCanvas** c_offaxis_hadronhit = new TCanvas*[tot_size];
  TH2F** h_offaxis_hadronhit_xy = new TH2F*[tot_size];
  TH2F** h_offaxis_hadronhit_zx = new TH2F*[tot_size];
  TH2F** h_offaxis_hadronhit_zy = new TH2F*[tot_size];

  // Set Palette
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.85);		//Stat box x position (top right hand corner)
  gStyle->SetStatY(0.9); 		//Stat box y position
  gStyle->SetStatW(0.2);	 		//Stat box width as fraction of pad size
  gStyle->SetStatH(0.15);	 		//Size of each line in stat box
  // gStyle->SetStatFontSize(0.02);	 		//Size of each line in stat box
  vector<Int_t> a_ND_off_axis_pos_vec = {200,-50,-300, -2550,-3050};
  vector<Int_t> a_ND_LAr_pos_vec = {-264, -216, -168, 24, 168, 216, 264};
  // Store info
  ofstream myfile;
   myfile.open ("Output_HadronhitCheck_ivy.txt");

  // Loop all events
  for (Int_t i_iwritten : *iwritten_vec)
  {
    if(myfileVerbose) myfile<< "i_iwritten: " << i_iwritten << "\n\n";
    cout << "i_iwritten: " << i_iwritten << "\n";
    Int_t canvas_counter = 0;
    for (Int_t i_ND_OffAxis_pos: a_ND_off_axis_pos_vec)
    {
      for (Int_t i_ND_LAr_pos: a_ND_LAr_pos_vec)
      {
        Int_t i_entry = tot_size * i_iwritten;
        for (i_entry ; i_entry < tot_size * (i_iwritten+1); i_entry++ )
        {
          t_effTree->GetEntry(i_entry);
          t_effValues->GetEntry(i_entry);
          if (ND_OffAxis_pos == i_ND_OffAxis_pos && ND_LAr_pos == i_ND_LAr_pos)
          {
            Int_t offset_X = i_ND_OffAxis_pos;
            vector<float> ND_OffAxis_hadronic_hit_xyz;
            vector<vector<float>> ND_OffAxis_hadronic_hit;

            ND_OffAxis_hadronic_hit.clear();
            int it_hadronhit_counter =0;
            for (vector<vector<float>>::iterator it_hadronhit = ND_OffAxis_Sim_hadronic_hit_xyz->begin(); it_hadronhit!=ND_OffAxis_Sim_hadronic_hit_xyz->end(); ++it_hadronhit)
            {
              for (Int_t i = 0; i < 3; i++)
              {
                if(verbose) cout << "it_hadronhit_counter: " << it_hadronhit_counter << ", it_xyz_counter: " << i << ", hit_x,y,z: " << it_hadronhit->at(i) << endl;
                ND_OffAxis_hadronic_hit_xyz.emplace_back(it_hadronhit->at(i));
              }
              ND_OffAxis_hadronic_hit.emplace_back(ND_OffAxis_hadronic_hit_xyz);
              ND_OffAxis_hadronic_hit_xyz.clear();
              it_hadronhit_counter++;
            }

            TString h_hadronhit_xy_name = Form("hadronhitXY_event_%d_OffAxis_%d_cm_LAr_%d_cm", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos);
            h_offaxis_hadronhit_xy[canvas_counter] = new TH2F(h_hadronhit_xy_name,h_hadronhit_xy_name,200,-500,500,200,-500,500);
            h_offaxis_hadronhit_xy[canvas_counter]->GetXaxis()->SetTitle("X [cm]");
            h_offaxis_hadronhit_xy[canvas_counter]->GetYaxis()->SetTitle("Y [cm]");
            h_offaxis_hadronhit_xy[canvas_counter]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");

            TString h_hadronhit_zx_name = Form("hadronhitZX_event_%d_OffAxis_%d_cm_LAr_%d_cm", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos);
            h_offaxis_hadronhit_zx[canvas_counter] = new TH2F(h_hadronhit_zx_name,h_hadronhit_zx_name,200,-100,600,200,-500,500);
            h_offaxis_hadronhit_zx[canvas_counter]->GetXaxis()->SetTitle("Z [cm]");
            h_offaxis_hadronhit_zx[canvas_counter]->GetYaxis()->SetTitle("X [cm]");
            h_offaxis_hadronhit_zx[canvas_counter]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");

            TString h_hadronhit_zy_name = Form("hadronhitZY_event_%d_OffAxis_%d_cm_LAr_%d_cm", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos);
            h_offaxis_hadronhit_zy[canvas_counter] = new TH2F(h_hadronhit_zy_name,h_hadronhit_zy_name,200,-100,600,200,-500,500);
            h_offaxis_hadronhit_zy[canvas_counter]->GetXaxis()->SetTitle("Z [cm]");
            h_offaxis_hadronhit_zy[canvas_counter]->GetYaxis()->SetTitle("Y [cm]");
            h_offaxis_hadronhit_zy[canvas_counter]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");


            Double_t vetoEnergyND = 0.;
            Double_t totEnergyND = 0.;

            for(Int_t ihadronhit = 0; ihadronhit < ND_OffAxis_hadronic_hit.size(); ihadronhit++)
            {
              h_offaxis_hadronhit_xy[canvas_counter]->Fill(ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X,ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1],HadronHitEdeps->at(ihadronhit));
              h_offaxis_hadronhit_zx[canvas_counter]->Fill(ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2],ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X,HadronHitEdeps->at(ihadronhit));
              h_offaxis_hadronhit_zy[canvas_counter]->Fill(ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2],ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1],HadronHitEdeps->at(ihadronhit));

              totEnergyND += HadronHitEdeps->at(ihadronhit);;
              if ( ( ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X                      > NDActiveVol_min[0]         && ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X                      < NDActiveVol_min[0] + vSize ) ||
                   ( ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1]        > NDActiveVol_min[1]         && ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1]        < NDActiveVol_min[1] + vSize ) ||
                   ( ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2]        > NDActiveVol_min[2]         && ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2]        < NDActiveVol_min[2] + vSize ) ||
                   ( ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X                      > NDActiveVol_max[0] - vSize && ND_OffAxis_hadronic_hit[ihadronhit][0]-offset_X                      < NDActiveVol_max[0] ) ||
                   ( ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1]        > NDActiveVol_max[1] - vSize && ND_OffAxis_hadronic_hit[ihadronhit][1]-NDLAr_OnAxis_offset[1]        < NDActiveVol_max[1] ) ||
                   ( ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2]        > NDActiveVol_max[2] - vSize && ND_OffAxis_hadronic_hit[ihadronhit][2]-NDLAr_OnAxis_offset[2]        < NDActiveVol_max[2] )
                 ){
                   vetoEnergyND += HadronHitEdeps->at(ihadronhit);
              } // end if hadron deposit in FD veto region
              // add outputs
              if(myfileVerbose)
              {
                myfile << "ND_OffAxis_pos: " << ND_OffAxis_pos << ", ND_LAr_pos: " << ND_LAr_pos << ", ihadronhit: " << ihadronhit << "\n";
                myfile << "ND_OffAxis_hadronic_hit_X: " << ND_OffAxis_hadronic_hit[ihadronhit][0] - offset_X << "\n";
                myfile << "ND_OffAxis_hadronic_hit_Y: " << ND_OffAxis_hadronic_hit[ihadronhit][1] - NDLAr_OnAxis_offset[1] << "\n";
                myfile << "ND_OffAxis_hadronic_hit_Z: " << ND_OffAxis_hadronic_hit[ihadronhit][2] - NDLAr_OnAxis_offset[2] << "\n\n";
              }
            }
            TString energy_name = Form("VetoE_%.2f_MeV, TotE_%.2f_MeV", vetoEnergyND, totEnergyND);

            // create canvas
            TString c_hadronhit_name = Form("c_hadronhit_event_%d_OffAxis_%d_cm_LAr_%d_cm", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos);
            TString c_hadronhit_title = Form("hadronhit event_%d_OffAxis_%d_cm_LAr_%d_cm", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos);
            c_offaxis_hadronhit[canvas_counter] = new TCanvas(c_hadronhit_name, c_hadronhit_title, 0,53,995,597);
            c_offaxis_hadronhit[canvas_counter]->Clear();
            c_offaxis_hadronhit[canvas_counter]->SetLeftMargin(0.10);
            c_offaxis_hadronhit[canvas_counter]->SetRightMargin(0.10);
            c_offaxis_hadronhit[canvas_counter]->Divide(2,2);

            c_offaxis_hadronhit[canvas_counter]->cd(1);
            c_offaxis_hadronhit[canvas_counter]->GetPad(1)->SetRightMargin(.15);
            h_offaxis_hadronhit_xy[canvas_counter]->Draw("COLZ");
            auto *xy_box = new TBox(NDActiveVol_min[0],NDActiveVol_min[1],NDActiveVol_max[0],NDActiveVol_max[1]);
            xy_box->SetLineColor(kBlack);
            xy_box->SetLineWidth(2);
            xy_box->SetFillStyle(0);
            xy_box->Draw();
            auto *xy_box1 = new TBox(NDActiveVol_min[0]+30,NDActiveVol_min[1]+30,NDActiveVol_max[0]-30,NDActiveVol_max[1]-30);
            xy_box1->SetLineColor(kBlue);
            xy_box1->SetLineWidth(2);
            xy_box1->SetFillStyle(0);
            xy_box1->Draw();
            auto *xy_box2 = new TBox(ND_FV_min[0],ND_FV_min[1],ND_FV_max[0],ND_FV_max[1]);
            xy_box2->SetLineColor(kRed);
            xy_box2->SetLineWidth(2);
            xy_box2->SetFillStyle(0);
            xy_box2->Draw();
            TLatex xy_text(-400,400,energy_name);
            xy_text.DrawClone();

            c_offaxis_hadronhit[canvas_counter]->cd(2);
            c_offaxis_hadronhit[canvas_counter]->GetPad(2)->SetRightMargin(.15);
            h_offaxis_hadronhit_zx[canvas_counter]->Draw("COLZ");
            // ND active volume
            auto *zx_box = new TBox(NDActiveVol_min[2],NDActiveVol_min[0],NDActiveVol_max[2],NDActiveVol_max[0]);
            zx_box->SetLineColor(kBlack);
            zx_box->SetLineWidth(2);
            zx_box->SetFillStyle(0);
            zx_box->Draw();
            // ND veto volume
            auto *zx_box1 = new TBox(NDActiveVol_min[2]+30,NDActiveVol_min[0]+30,NDActiveVol_max[2]-30,NDActiveVol_max[0]-30);
            zx_box1->SetLineColor(kBlue);
            zx_box1->SetLineWidth(2);
            zx_box1->SetFillStyle(0);
            zx_box1->Draw();
            // ND fiducial volume
            auto *zx_box2 = new TBox(ND_FV_min[2],ND_FV_min[0],ND_FV_max[2],ND_FV_max[0]);
            zx_box2->SetLineColor(kRed);
            zx_box2->SetLineWidth(2);
            zx_box2->SetFillStyle(0);
            zx_box2->Draw();
            TLatex xz_text(-50,400,energy_name);
            xz_text.DrawClone();

            c_offaxis_hadronhit[canvas_counter]->cd(3);
            c_offaxis_hadronhit[canvas_counter]->GetPad(3)->SetRightMargin(.15);
            h_offaxis_hadronhit_zy[canvas_counter]->Draw("COLZ");
            auto *yz_box = new TBox(NDActiveVol_min[2],NDActiveVol_min[1],NDActiveVol_max[2],NDActiveVol_max[1]);
            yz_box->SetLineColor(kBlack);
            yz_box->SetLineWidth(2);
            yz_box->SetFillStyle(0);
            yz_box->Draw();
            auto *yz_box1 = new TBox(NDActiveVol_min[2]+30,NDActiveVol_min[1]+30,NDActiveVol_max[2]-30,NDActiveVol_max[1]-30);
            yz_box1->SetLineColor(kBlue);
            yz_box1->SetLineWidth(2);
            yz_box1->SetFillStyle(0);
            yz_box1->Draw();
            auto *yz_box2 = new TBox(ND_FV_min[2],ND_FV_min[1],ND_FV_max[2],ND_FV_max[1]);
            yz_box2->SetLineColor(kRed);
            yz_box2->SetLineWidth(2);
            yz_box2->SetFillStyle(0);
            yz_box2->Draw();
            TLatex yz_text(-50,400,energy_name);
            yz_text.DrawClone();


            gPad->Update();
            gPad->Modified();
            gSystem->ProcessEvents();


            outFile->cd("OffAxis hadron hits");
            c_offaxis_hadronhit[canvas_counter]->Write();
            c_offaxis_hadronhit[canvas_counter]->Close();

            if(verbose) cout << "canvas_counter: " << canvas_counter << ", offset_X: " << offset_X << endl;
            canvas_counter++;
          }
        }
      }
    }
    //
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //
    // Create canvas for hadron hit
    Int_t n_plot = 0;
    Int_t i_n_plot = 0;
    for (Int_t i_ND_OffAxis_pos: a_ND_off_axis_pos_vec)
    {
      for (Int_t i_ND_LAr_pos: a_ND_LAr_pos_vec)
      {
        Int_t i_entry = tot_size * i_iwritten;
        for (i_entry ; i_entry < tot_size * (i_iwritten+1); i_entry++ )
        {
          t_effTree->GetEntry(i_entry);
          t_effValues->GetEntry(i_entry);

          if (ND_OffAxis_pos == i_ND_OffAxis_pos && ND_LAr_pos == i_ND_LAr_pos)
          {
            Int_t offset_X = i_ND_OffAxis_pos;
            // Store hadron hit info into vector
            vector<float> ThrowDepsX_hit;
            vector<vector<float>> ThrowDepsX;
            ThrowDepsX.clear();
            int it_throw_x_counter =0;
            for (vector<vector<float>>::iterator it_throw = CurrentThrowDepsX->begin(); it_throw!=CurrentThrowDepsX->end(); ++it_throw)
            {
              for (Int_t ihadronhit =0; ihadronhit < ND_Sim_n_hadronic_Edep_a; ihadronhit++)
              {
                if(verbose) cout << "it_throw_x_counter: " << it_throw_x_counter << ", ihadronhit: " << ihadronhit << ", hit_x: " << it_throw->at(ihadronhit) << endl;
                ThrowDepsX_hit.emplace_back(it_throw->at(ihadronhit));
              }
              ThrowDepsX.emplace_back(ThrowDepsX_hit);
              ThrowDepsX_hit.clear();
              it_throw_x_counter++;
            }
            vector<float> ThrowDepsY_hit;
            vector<vector<float>> ThrowDepsY;
            ThrowDepsY.clear();
            int it_throw_y_counter =0;
            for (vector<vector<float>>::iterator it_throw = CurrentThrowDepsY->begin(); it_throw!=CurrentThrowDepsY->end(); ++it_throw)
            {
              for (Int_t ihadronhit =0; ihadronhit < ND_Sim_n_hadronic_Edep_a; ihadronhit++)
              {
                if(verbose) cout << "it_throw_x_counter: " << it_throw_x_counter << ", ihadronhit: " << ihadronhit << ", hit_y: " << it_throw->at(ihadronhit) << endl;
                ThrowDepsY_hit.emplace_back(it_throw->at(ihadronhit));
              }
              ThrowDepsY.emplace_back(ThrowDepsY_hit);
              ThrowDepsY_hit.clear();
              it_throw_y_counter++;
            }
            vector<float> ThrowDepsZ_hit;
            vector<vector<float>> ThrowDepsZ;
            ThrowDepsZ.clear();
            int it_throw_z_counter =0;
            for (vector<vector<float>>::iterator it_throw = CurrentThrowDepsZ->begin(); it_throw!=CurrentThrowDepsZ->end(); ++it_throw)
            {
              for (Int_t ihadronhit =0; ihadronhit < ND_Sim_n_hadronic_Edep_a; ihadronhit++)
              {
                if(verbose) cout<< "ND_OffAxis_pos: " << ND_OffAxis_pos << ", ND_LAr_pos: " << ND_LAr_pos << ", ithrow: " << it_throw_z_counter << ", ihadronhit: " << ihadronhit << ", hit_z: " << it_throw->at(ihadronhit) << endl;
                ThrowDepsZ_hit.emplace_back(it_throw->at(ihadronhit));
              }
              ThrowDepsZ.emplace_back(ThrowDepsZ_hit);
              ThrowDepsZ_hit.clear();
              it_throw_z_counter++;
            }

            // Draw plots
            // for(Int_t ithrow = 0; ithrow < N_throws; ithrow++)
            for(Int_t ithrow = 0; ithrow < 20; ithrow++)
            {
              cout << "ithrow: " << ithrow <<endl;
              n_plot = i_n_plot*N_throws + ithrow;

              TString h_hadronhit_xy_name = Form("hadronhitXY_event_%d_OffAxis_%d_cm_LAr_%d_cm_throw_%d", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos,ithrow);
              h_hadronhit_xy[n_plot] = new TH2F(h_hadronhit_xy_name,h_hadronhit_xy_name,200,-500,500,200,-500,500);
              h_hadronhit_xy[n_plot]->GetXaxis()->SetTitle("X [cm]");
              h_hadronhit_xy[n_plot]->GetYaxis()->SetTitle("Y [cm]");
              h_hadronhit_xy[n_plot]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");

              TString h_hadronhit_zx_name = Form("hadronhitZX_event_%d_OffAxis_%d_cm_LAr_%d_cm_throw_%d", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos,ithrow);
              h_hadronhit_zx[n_plot] = new TH2F(h_hadronhit_zx_name,h_hadronhit_zx_name,200,-100,600,200,-500,500);
              h_hadronhit_zx[n_plot]->GetXaxis()->SetTitle("Z [cm]");
              h_hadronhit_zx[n_plot]->GetYaxis()->SetTitle("X [cm]");
              h_hadronhit_zx[n_plot]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");

              TString h_hadronhit_zy_name = Form("hadronhitZY_event_%d_OffAxis_%d_cm_LAr_%d_cm_throw_%d", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos,ithrow);
              h_hadronhit_zy[n_plot] = new TH2F(h_hadronhit_zy_name,h_hadronhit_zy_name,200,-100,600,200,-500,500);
              h_hadronhit_zy[n_plot]->GetXaxis()->SetTitle("Z [cm]");
              h_hadronhit_zy[n_plot]->GetYaxis()->SetTitle("Y [cm]");
              h_hadronhit_zy[n_plot]->GetZaxis()->SetTitle("HadronHitEdeps [MeV]");

              for(Int_t ihadronhit =0; ihadronhit < ND_Sim_n_hadronic_Edep_a; ihadronhit++)
              {
                h_hadronhit_xy[n_plot]->Fill(ThrowDepsX[ithrow][ihadronhit] - offset_X,ThrowDepsY[ithrow][ihadronhit]-NDLAr_OnAxis_offset[1],HadronHitEdeps->at(ihadronhit));
                h_hadronhit_zx[n_plot]->Fill(ThrowDepsZ[ithrow][ihadronhit] - NDLAr_OnAxis_offset[2],ThrowDepsX[ithrow][ihadronhit] - offset_X,HadronHitEdeps->at(ihadronhit));
                h_hadronhit_zy[n_plot]->Fill(ThrowDepsZ[ithrow][ihadronhit] - NDLAr_OnAxis_offset[2],ThrowDepsY[ithrow][ihadronhit] - NDLAr_OnAxis_offset[1],HadronHitEdeps->at(ihadronhit));
                if(myfileVerbose)
                {
                  myfile << "ND_OffAxis_pos: " << ND_OffAxis_pos << ", ND_LAr_pos: " << ND_LAr_pos << ", ThrowDepsX at ithrow" << ithrow << ", at ihadronhit: " << ihadronhit << ", is: " << ThrowDepsX[ithrow][ihadronhit] - offset_X<< "\n";
                  myfile << "ND_OffAxis_pos: " << ND_OffAxis_pos << ", ND_LAr_pos: " << ND_LAr_pos << ", ThrowDepsY at ithrow" << ithrow << ", at ihadronhit: " << ihadronhit << ", is: " << ThrowDepsY[ithrow][ihadronhit] - NDLAr_OnAxis_offset[1] << "\n";
                  myfile << "ND_OffAxis_pos: " << ND_OffAxis_pos << ", ND_LAr_pos: " << ND_LAr_pos << ", ThrowDepsZ at ithrow" << ithrow << ", at ihadronhit: " << ihadronhit << ", is: " << ThrowDepsZ[ithrow][ihadronhit] - NDLAr_OnAxis_offset[2] << "\n\n";
                }
              }
              TString energy_name = Form("VetoE_%.2f_MeV, TotE_%.2f_MeV", CurrentThrowVetoE->at(ithrow), CurrentThrowTotE->at(ithrow));

              TString c_hadronhit_name = Form("c_hadronhit_event_%d_OffAxis_%d_cm_LAr_%d_cm_throw_%d", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos,ithrow);
              TString c_hadronhit_title = Form("hadronhit event_%d_OffAxis_%d_cm_LAr_%d_cm_throw_%d", i_iwritten,i_ND_OffAxis_pos,i_ND_LAr_pos,ithrow);
              c_hadronhit[n_plot] = new TCanvas(c_hadronhit_name, c_hadronhit_title, 0,53,995,597);
              c_hadronhit[n_plot]->Clear();
              c_hadronhit[n_plot]->Range(-666.6667,-625,1000,625);
              c_hadronhit[n_plot]->SetLeftMargin(0.10);
              c_hadronhit[n_plot]->SetRightMargin(0.10);
              c_hadronhit[n_plot]->Divide(2,2);

              c_hadronhit[n_plot]->cd(1);
              c_hadronhit[n_plot]->GetPad(1)->SetRightMargin(.15);
              h_hadronhit_xy[n_plot]->Draw("COLZ");
              auto *xy_box = new TBox(NDActiveVol_min[0],NDActiveVol_min[1],NDActiveVol_max[0],NDActiveVol_max[1]);
              xy_box->SetLineColor(kBlack);
              xy_box->SetLineWidth(2);
              xy_box->SetFillStyle(0);
              xy_box->Draw();
              auto *xy_box1 = new TBox(NDActiveVol_min[0]+30,NDActiveVol_min[1]+30,NDActiveVol_max[0]-30,NDActiveVol_max[1]-30);
              xy_box1->SetLineColor(kBlue);
              xy_box1->SetLineWidth(2);
              xy_box1->SetFillStyle(0);
              xy_box1->Draw();
              auto *xy_box2 = new TBox(ND_FV_min[0],ND_FV_min[1],ND_FV_max[0],ND_FV_max[1]);
              xy_box2->SetLineColor(kRed);
              xy_box2->SetLineWidth(2);
              xy_box2->SetFillStyle(0);
              xy_box2->Draw();
              TLatex xy_text(-400,400,energy_name);
              xy_text.DrawClone();

              c_hadronhit[n_plot]->cd(2);
              c_hadronhit[n_plot]->GetPad(2)->SetRightMargin(.15);
              h_hadronhit_zx[n_plot]->Draw("COLZ");
              // ND active volume
              auto *zx_box = new TBox(NDActiveVol_min[2],NDActiveVol_min[0],NDActiveVol_max[2],NDActiveVol_max[0]);
              zx_box->SetLineColor(kBlack);
              zx_box->SetLineWidth(2);
              zx_box->SetFillStyle(0);
              zx_box->Draw();
              // ND veto volume
              auto *zx_box1 = new TBox(NDActiveVol_min[2]+30,NDActiveVol_min[0]+30,NDActiveVol_max[2]-30,NDActiveVol_max[0]-30);
              zx_box1->SetLineColor(kBlue);
              zx_box1->SetLineWidth(2);
              zx_box1->SetFillStyle(0);
              zx_box1->Draw();
              // ND fiducial volume
              auto *zx_box2 = new TBox(ND_FV_min[2],ND_FV_min[0],ND_FV_max[2],ND_FV_max[0]);
              zx_box2->SetLineColor(kRed);
              zx_box2->SetLineWidth(2);
              zx_box2->SetFillStyle(0);
              zx_box2->Draw();
              TLatex xz_text(-50,400,energy_name);
              xz_text.DrawClone();

              c_hadronhit[n_plot]->cd(3);
              c_hadronhit[n_plot]->GetPad(3)->SetRightMargin(.15);
              h_hadronhit_zy[n_plot]->Draw("COLZ");
              auto *yz_box = new TBox(NDActiveVol_min[2],NDActiveVol_min[1],NDActiveVol_max[2],NDActiveVol_max[1]);
              yz_box->SetLineColor(kBlack);
              yz_box->SetLineWidth(2);
              yz_box->SetFillStyle(0);
              yz_box->Draw();
              auto *yz_box1 = new TBox(NDActiveVol_min[2]+30,NDActiveVol_min[1]+30,NDActiveVol_max[2]-30,NDActiveVol_max[1]-30);
              yz_box1->SetLineColor(kBlue);
              yz_box1->SetLineWidth(2);
              yz_box1->SetFillStyle(0);
              yz_box1->Draw();
              auto *yz_box2 = new TBox(ND_FV_min[2],ND_FV_min[1],ND_FV_max[2],ND_FV_max[1]);
              yz_box2->SetLineColor(kRed);
              yz_box2->SetLineWidth(2);
              yz_box2->SetFillStyle(0);
              yz_box2->Draw();
              TLatex yz_text(-50,400,energy_name);
              yz_text.DrawClone();



              gPad->Update();
              gPad->Modified();
              gSystem->ProcessEvents();

              outFile->cd("hadron hits");
              c_hadronhit[n_plot]->Write();
              c_hadronhit[n_plot]->Close();

            }
            i_n_plot++;
          }
        }
      }
    }


  } // end iwritten_vec



  delete[] h_hadronhit_xy;
  delete[] h_hadronhit_zx;
  delete[] h_hadronhit_zy;
  delete[] c_hadronhit;
  delete[] h_offaxis_hadronhit_zy;
  delete[] h_offaxis_hadronhit_zx;
  delete[] h_offaxis_hadronhit_xy;
  delete[] c_offaxis_hadronhit;

  myfile.close();
  outFile->Close();
} // end ReadNtuple
