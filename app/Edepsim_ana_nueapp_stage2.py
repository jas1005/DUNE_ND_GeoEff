#! /usr/bin/env python
"""
Runs over edepsim files and produce paired events in near and far detector.
"""
import time
start_time = time.time()

import numpy as np
import os
import os.path
import sys
import torch
from scipy.spatial.transform import Rotation as R
from muonEffModel import muonEffModel
import ROOT
from ROOT import TG4Event, TFile, TTree, TGraph
from ROOT import gROOT # for creating the output file
from array import array
from math import cos, sin
import random

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("stage1out")
parser.add_argument("--config", type=str, default="UserConfig.py")
parser.add_argument(
    "--out_dir", type=str, default="/dune/app/users/weishi/testn2fd/DUNE_ND_GeoEff/app/output"
)
args = parser.parse_args()
stage1file, config_file, out_path = args.stage1out, args.config, args.out_dir

with open(config_file) as infile:
    exec(infile.read())
import pyGeoEff

ROOT.gROOT.ProcessLine('#include<vector>') # for outputting data to ROOT file

# Load muon neural network
net = muonEffModel()
net.load_state_dict(torch.load("./muonEff30.nn", map_location=torch.device('cpu')))
net.eval()

#########
# I/O
#########
# Default directory: should exist from stage 1
rootpath = out_path + "/nue_output"
if not os.path.exists( out_path):
    print("out_path '" + out_path + "' should exist from stage 1, but it's not!")
    sys.exit()
if not os.path.exists( rootpath):
    print("rootpath '" + rootpath + "' should exist from stage 1, but it's not!")
    sys.exit()
print(" \n") # separate output

##################################
# Open the stage 1 file
##################################
# Here it assumes to look in the default rootpath folder
numu_stage1_file = TFile('{}/{}'.format(rootpath, stage1file), "OPEN")
if not numu_stage1_file:
    print ("Error: could not open file {}/{}".format(rootpath, stage1file))
    sys.exit()
inputstage1Tree = numu_stage1_file.Get("myStage1Events")
stage1entries = inputstage1Tree.GetEntriesFast()

# Output
outfilename = "n2fd_nueapp_paired_stage2"
# Create the ROOT file that will hold the output of this script
f_out = TFile('{}/{}.root'.format(rootpath, outfilename), 'RECREATE')

# Copy stage 1 tree and append stage 2 output
myStage2Events = inputstage1Tree.CloneTree(0)
#myStage2Events = TTree('myStage2Events', 'myStage2Events')
maxEdeps = 100000 # max number of edeps for initialize arrays

##################################
# electron LArBath evt info
##################################
nElectronEdeps = array('i', [0])
myStage2Events.Branch('nElectronEdeps', nElectronEdeps, 'nElectronEdeps/I')
electrondeps_trkID = np.zeros((maxEdeps,), dtype=np.int32)
myStage2Events.Branch('electrondeps_trkID', electrondeps_trkID, 'electrondeps_trkID[nElectronEdeps]/I')
electrondeps_parentID = np.zeros((maxEdeps,), dtype=np.int32)
myStage2Events.Branch('electrondeps_parentID', electrondeps_parentID, 'electrondeps_parentID[nElectronEdeps]/I')
electrondeps_pdg = np.zeros((maxEdeps,), dtype=np.int32)
myStage2Events.Branch('electrondeps_pdg', electrondeps_pdg, 'electrondeps_pdg[nElectronEdeps]/I')
electrondeps_E_MeV = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('electrondeps_E_MeV', electrondeps_E_MeV, 'electrondeps_E_MeV[nElectronEdeps]/F')
electrondeps_start_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('electrondeps_start_t_us', electrondeps_start_t_us, 'electrondeps_start_t_us[nElectronEdeps]/F')
electrondeps_stop_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('electrondeps_stop_t_us', electrondeps_stop_t_us, 'electrondeps_stop_t_us[nElectronEdeps]/F')
# edeps generated in LArBath (start point)
larbath_electrondeps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_start_x_cm', larbath_electrondeps_start_x_cm, 'larbath_electrondeps_start_x_cm[nElectronEdeps]/F') # larbath edeps x
larbath_electrondeps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_start_y_cm', larbath_electrondeps_start_y_cm, 'larbath_electrondeps_start_y_cm[nElectronEdeps]/F')
larbath_electrondeps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_start_z_cm', larbath_electrondeps_start_z_cm, 'larbath_electrondeps_start_z_cm[nElectronEdeps]/F')
# edeps generated in LArBath (stop point)
larbath_electrondeps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_stop_x_cm', larbath_electrondeps_stop_x_cm, 'larbath_electrondeps_stop_x_cm[nElectronEdeps]/F')
larbath_electrondeps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_stop_y_cm', larbath_electrondeps_stop_y_cm, 'larbath_electrondeps_stop_y_cm[nElectronEdeps]/F')
larbath_electrondeps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('larbath_electrondeps_stop_z_cm', larbath_electrondeps_stop_z_cm, 'larbath_electrondeps_stop_z_cm[nElectronEdeps]/F')

##################################
# FD paired evt with nd non-ecc
##################################
# Random thrown vertex in FD
nFDnueEdeps = array('i', [0])
myStage2Events.Branch('nFDnueEdeps', nFDnueEdeps, 'nFDnueEdeps/I')
fd_nuevtx_cm_pair_nd_nonecc = array('f', 3*[0.0])
myStage2Events.Branch('fd_nuevtx_cm_pair_nd_nonecc', fd_nuevtx_cm_pair_nd_nonecc, 'fd_nuevtx_cm_pair_nd_nonecc[3]/F')
# Edeps in FD random throw (start points)
fd_nuedeps_start_x_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_start_x_cm_pair_nd_nonecc', fd_nuedeps_start_x_cm_pair_nd_nonecc, 'fd_nuedeps_start_x_cm_pair_nd_nonecc[nFDnueEdeps]/F')
fd_nuedeps_start_y_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_start_y_cm_pair_nd_nonecc', fd_nuedeps_start_y_cm_pair_nd_nonecc, 'fd_nuedeps_start_y_cm_pair_nd_nonecc[nFDnueEdeps]/F')
fd_nuedeps_start_z_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_start_z_cm_pair_nd_nonecc', fd_nuedeps_start_z_cm_pair_nd_nonecc, 'fd_nuedeps_start_z_cm_pair_nd_nonecc[nFDnueEdeps]/F')
# Edeps in FD random throw (stop points)
fd_nuedeps_stop_x_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_stop_x_cm_pair_nd_nonecc', fd_nuedeps_stop_x_cm_pair_nd_nonecc, 'fd_nuedeps_stop_x_cm_pair_nd_nonecc[nFDnueEdeps]/F')
fd_nuedeps_stop_y_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_stop_y_cm_pair_nd_nonecc', fd_nuedeps_stop_y_cm_pair_nd_nonecc, 'fd_nuedeps_stop_y_cm_pair_nd_nonecc[nFDnueEdeps]/F')
fd_nuedeps_stop_z_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myStage2Events.Branch('fd_nuedeps_stop_z_cm_pair_nd_nonecc', fd_nuedeps_stop_z_cm_pair_nd_nonecc, 'fd_nuedeps_stop_z_cm_pair_nd_nonecc[nFDnueEdeps]/F')

fd_nue_throws_passed = array('i', [0])
myStage2Events.Branch('fd_nue_throws_passed', fd_nue_throws_passed, 'fd_nue_throws_passed/I')

##########################################################
## Loop over stage 1 events and collate electron deposits
##########################################################
for ientry in range(stage1entries):

    ####################################################
    # For event ientry, read the corresponding electron edepsim file
    ####################################################
    electron_edep_file = TFile('{}/edep_LArBath_electron_evtnb_{}.root'.format(rootpath, ientry), "OPEN")
    if not electron_edep_file:
        print ("Error: could not open file {}/edep_LArBath_electron_evtnb_{}.root".format(rootpath, ientry))
        sys.exit()

    inputTree = electron_edep_file.Get("EDepSimEvents")

    event = TG4Event()
    inputTree.SetBranchAddress("Event", event)

    entries = inputTree.GetEntriesFast()
    if entries != 1:
        print("Edep-sim tree has more than one entries!")
        sys.exit()

    # Only has one event
    inputTree.GetEntry(0)  # Only has one event
    # initialize
    # Define here but do not write out
    electron_dep_startpos_list = list()
    electron_dep_stoppos_list = list()
    electron_dep_trkID_list  = list()
    electron_dep_parentID_list = list()
    electron_dep_pdg_list = list()
    electron_edep_list   = list()
    electron_dep_starttime_list = list()
    electron_dep_stoptime_list = list()
    electron_dep_pos_list = list()
    nElectronEdeps[0] = 0
    fd_nuevtx_cm_pair_nd_nonecc[0] = 0; fd_nuevtx_cm_pair_nd_nonecc[1] = 0; fd_nuevtx_cm_pair_nd_nonecc[2] = 0;
    fd_nue_throws_passed[0] = 0
    trajectories_parentid = np.empty(len(event.Trajectories), dtype=np.int32)
    trajectories_pdg = np.empty(len(event.Trajectories), dtype=np.int32)

    for iTraj, trajectory in enumerate(event.Trajectories):
        # iTraj is same as TrackId, starts from #0
        trajectories_parentid[iTraj] = trajectory.GetParentId()
        trajectories_pdg[iTraj] = trajectory.GetPDGCode()
        #print("iTraj #", iTraj, " parentID: ", trajectories_parentid[iTraj], " pdg: ", trajectories_pdg[iTraj] )

    for containerName, hitSegments in event.SegmentDetectors:
        # iHit is the index, hitSEgment is the data stored at the index in the second item in event.SegementDectectors
        for iHit, hitSegment in enumerate(hitSegments):

            # Energy deposit from primary particles
            edep_start_x = hitSegment.GetStart().X() * edep2cm
            edep_start_y = hitSegment.GetStart().Y() * edep2cm
            edep_start_z = hitSegment.GetStart().Z() * edep2cm
            edep_start_t = hitSegment.GetStart().T() * edep2us
            edep_stop_x = hitSegment.GetStop().X() * edep2cm
            edep_stop_y = hitSegment.GetStop().Y() * edep2cm
            edep_stop_z = hitSegment.GetStop().Z() * edep2cm
            edep_stop_t = hitSegment.GetStop().T() * edep2us
            edep_x = (edep_start_x + edep_stop_x)/2 # use this for hadronic veto
            edep_y = (edep_start_y + edep_stop_y)/2
            edep_z = (edep_start_z + edep_stop_z)/2
            edep_trkID = hitSegment.Contrib[0] # what is different if use GetPrimaryId()?
            edep_parentID = trajectories_parentid[edep_trkID]
            edep_pdg = trajectories_pdg[edep_trkID]
            edep = hitSegment.GetEnergyDeposit()

            electron_dep_startpos_list.append(edep_start_x)
            electron_dep_startpos_list.append(edep_start_y)
            electron_dep_startpos_list.append(edep_start_z)
            electron_dep_starttime_list.append(edep_start_t)
            electron_dep_stoppos_list.append(edep_stop_x)
            electron_dep_stoppos_list.append(edep_stop_y)
            electron_dep_stoppos_list.append(edep_stop_z)
            electron_dep_stoptime_list.append(edep_stop_t)
            electron_dep_trkID_list.append(edep_trkID)
            electron_dep_parentID_list.append(edep_parentID)
            electron_dep_pdg_list.append(edep_pdg)
            electron_edep_list.append(edep)

    nElectronEdeps[0] = len(electron_edep_list)

    ####################################################
    # End read the electron edepsim file
    ####################################################

    # NOW get all ecc non-leptonic/hadronic deposits from stage 1
    inputstage1Tree.GetEntry(ientry)
    # initialize
    had_dep_pos_list = list()
    had_edep_list = list()
    all_dep_startpos_list = list()
    all_dep_stoppos_list = list()
    all_edep_list = list()
    nFDnueEdeps[0] = 0
    # Convert stage 1 ecc nonlepton edeps np.array to list
    had_dep_startpos_x_list = list(inputstage1Tree.nd_nonlepdeps_start_x_cm_ecc)
    had_dep_startpos_y_list = list(inputstage1Tree.nd_nonlepdeps_start_y_cm_ecc)
    had_dep_startpos_z_list = list(inputstage1Tree.nd_nonlepdeps_start_z_cm_ecc)
    had_dep_stoppos_x_list = list(inputstage1Tree.nd_nonlepdeps_stop_x_cm_ecc)
    had_dep_stoppos_y_list = list(inputstage1Tree.nd_nonlepdeps_stop_y_cm_ecc)
    had_dep_stoppos_z_list = list(inputstage1Tree.nd_nonlepdeps_stop_z_cm_ecc)
    had_edep_list = list(inputstage1Tree.nonlepdeps_E_MeV)
    all_edep_list = list(inputstage1Tree.nonlepdeps_E_MeV)

    # Add hadronic/non-leptonic edeps from stage 1
    for ihaddep in range(len(had_edep_list)):
        all_dep_startpos_list.append(had_dep_startpos_x_list[ihaddep])
        all_dep_startpos_list.append(had_dep_startpos_y_list[ihaddep])
        all_dep_startpos_list.append(had_dep_startpos_z_list[ihaddep])
        all_dep_stoppos_list.append(had_dep_stoppos_x_list[ihaddep])
        all_dep_stoppos_list.append(had_dep_stoppos_y_list[ihaddep])
        all_dep_stoppos_list.append(had_dep_stoppos_z_list[ihaddep])
        had_dep_pos_list.append((had_dep_startpos_x_list[ihaddep] + had_dep_stoppos_x_list[ihaddep])/2)
        had_dep_pos_list.append((had_dep_startpos_y_list[ihaddep] + had_dep_stoppos_y_list[ihaddep])/2)
        had_dep_pos_list.append((had_dep_startpos_z_list[ihaddep] + had_dep_stoppos_z_list[ihaddep])/2)

    # Append/concatenate electron edeps
    all_dep_startpos_list = all_dep_startpos_list + electron_dep_startpos_list
    all_dep_stoppos_list = all_dep_stoppos_list + electron_dep_stoppos_list
    all_edep_list = all_edep_list + electron_edep_list
    nFDnueEdeps[0] = len(all_edep_list)

    ################
    # Paired FD evt
    ################
    # By this point we should have the stitched event after earth curvature correction (ecc)
    # (i.e. electron generated with same muon kinematics + non-muon part deposits from stage 1)
    # Only random translate the event in fd to obtain paired FD event

    ##################################################################
    # Initialize geometric efficiency module to manipulate energy deps
    ##################################################################
    seed = random.randrange(1024)
    print ("-- seed:", seed)

    geoEff = pyGeoEff.geoEff(seed, False)

    # Use neutrino decay position, rather than fixed neutrino direction as symmetry axis
    geoEff.setUseFixedBeamDir(False)

    # 30 cm veto region for hadronic veto
    geoEff.setVetoSizes([30])

    # 30 MeV E threshold for hadronic veto
    geoEff.setVetoEnergyThresholds([30])

    # Far detector active dimensions for hadronic veto
    geoEff.setFDActiveX(FDActiveVol_min[0], FDActiveVol_max[0])
    geoEff.setFDActiveY(FDActiveVol_min[1], FDActiveVol_max[1])
    geoEff.setFDActiveZ(FDActiveVol_min[2], FDActiveVol_max[2])

    # Range for random translation throws in FD fiducial volume
    geoEff.setRangeXFD(FD_FV_min[0], FD_FV_max[0])
    geoEff.setRangeYFD(FD_FV_min[1], FD_FV_max[1])
    geoEff.setRangeZFD(FD_FV_min[2], FD_FV_max[2])

    # Further conversion of list to Matrix3Xf class
    all_startposdep_fdorig_matrix = geoEff.Vecf2EigenMat3Xf(all_dep_startpos_list)
    all_stopposdep_fdorig_matrix = geoEff.Vecf2EigenMat3Xf(all_dep_stoppos_list)
    had_posdep_fdorig_matrix = geoEff.Vecf2EigenMat3Xf(had_dep_pos_list)

    # Tell the module where the vertex is in FD
    geoEff.setVertexFD(0, 0, 0) # it's at FD origin because we moved it to origin in stage 1 and then just rotated at there

    # Randomly throw the stitched event in FD
    tot_fd_throw_pair_nd_nonecc = 0

    # If after max_fd_throws still don't pass at nd, stop and move to next event (otherwise too much computing resources)
    while tot_fd_throw_pair_nd_nonecc < max_fd_throws:
        print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc)
        ##########################################################################################
        # Below do random throw (translate only) in FD similar to ND: only one throw in FD at a time
        ##########################################################################################
        geoEff.setNthrowsFD(1)
        geoEff.throwTransformsFD() # this randomly generates new vtx position in FD FV

        fd_nuevtx_x_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsX()
        fd_nuevtx_y_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsY()
        fd_nuevtx_z_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsZ()

        # Check if it passes FD hadronic veto
        geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list
        fdthrowresulthad_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(had_posdep_fdorig_matrix)

        if (fdthrowresulthad_pair_nd_nonecc.containresult[0][0][0] != 0):
            print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc, "passed fd had veto")
            print ("---- throw x: ", fd_nuevtx_x_cm_pair_nd_nonecc[0], "y: ", fd_nuevtx_y_cm_pair_nd_nonecc[0], ", z: ", fd_nuevtx_z_cm_pair_nd_nonecc[0])
            ###########################################################
            # FD rand throw passes veto, write paired evt info
            ###########################################################

            # Now change to the full list of edeps
            # the random thrown x/y/z should reamin the same because throw is done above already
            geoEff.setHitSegEdeps(all_edep_list)
            fdthrowresultall_start_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(all_startposdep_fdorig_matrix)
            # Repeat for edepsim stop points !!!
            geoEff.setHitSegEdeps(all_edep_list)
            fdthrowresultall_stop_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(all_stopposdep_fdorig_matrix)

            print ("Found paired fd-nd non ecc event")

            fd_nue_throws_passed[0] = 1

            #################################
            # Unpack info and store to output
            #################################
            print ("Saving...")

            # Save LArBath electron info
            electrondeps_trkID[:nElectronEdeps[0]] = np.array(electron_dep_trkID_list, dtype=np.int32)
            electrondeps_parentID[:nElectronEdeps[0]] = np.array(electron_dep_parentID_list, dtype=np.int32)
            electrondeps_pdg[:nElectronEdeps[0]] = np.array(electron_dep_pdg_list, dtype=np.int32)
            electrondeps_E_MeV[:nElectronEdeps[0]] = np.array(electron_edep_list, dtype=np.float32)
            electrondeps_start_t_us[:nElectronEdeps[0]] = np.array(electron_dep_starttime_list, dtype=np.float32)
            electrondeps_stop_t_us[:nElectronEdeps[0]] = np.array(electron_dep_stoptime_list, dtype=np.float32)
            larbath_electrondeps_start_x_cm[:nElectronEdeps[0]] = np.array(electron_dep_startpos_list[::3], dtype=np.float32) # every 3 element: x list
            larbath_electrondeps_start_y_cm[:nElectronEdeps[0]] = np.array(electron_dep_startpos_list[1::3], dtype=np.float32) # y list
            larbath_electrondeps_start_z_cm[:nElectronEdeps[0]] = np.array(electron_dep_startpos_list[2::3], dtype=np.float32) # z list
            larbath_electrondeps_stop_x_cm[:nElectronEdeps[0]] = np.array(electron_dep_stoppos_list[::3], dtype=np.float32)
            larbath_electrondeps_stop_y_cm[:nElectronEdeps[0]] = np.array(electron_dep_stoppos_list[1::3], dtype=np.float32)
            larbath_electrondeps_stop_z_cm[:nElectronEdeps[0]] = np.array(electron_dep_stoppos_list[2::3], dtype=np.float32)

            # Save selected FD nue event info
            fd_nuevtx_cm_pair_nd_nonecc[0] = fd_nuevtx_x_cm_pair_nd_nonecc[0]
            fd_nuevtx_cm_pair_nd_nonecc[1] = fd_nuevtx_y_cm_pair_nd_nonecc[0]
            fd_nuevtx_cm_pair_nd_nonecc[2] = fd_nuevtx_z_cm_pair_nd_nonecc[0]
            fd_nuedeps_start_x_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][0,:], dtype=np.float32)
            fd_nuedeps_start_y_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][1,:], dtype=np.float32)
            fd_nuedeps_start_z_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][2,:], dtype=np.float32)
            fd_nuedeps_stop_x_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][0,:], dtype=np.float32)
            fd_nuedeps_stop_y_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][1,:], dtype=np.float32)
            fd_nuedeps_stop_z_cm_pair_nd_nonecc[:nFDnueEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][2,:], dtype=np.float32)

            # Break the while loop, move on to next evt
            print ("Paired data saved, breaking fd throw loop")
            break
        else:
            print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc, "failed fd had veto!")

        # indentation is important!
        # if don't, put it in another random FD pos...until it passes FD veto
        tot_fd_throw_pair_nd_nonecc = tot_fd_throw_pair_nd_nonecc + 1

    if not fd_nue_throws_passed[0]:
        print("-- no selected fd throw found after max tries. Giving up!")
        nElectronEdeps[0] = 0
        nFDnueEdeps[0] = 0

    # event level
    myStage2Events.Fill()

myStage2Events.SetName("myStage2Events"); # update tree name
f_out.cd()
myStage2Events.Write()

print("\n")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
