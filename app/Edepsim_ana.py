#! /usr/bin/env python
"""
Runs over edepsim files and produce paired events in near and far detector.
Produce validation plots to see if the data is being manipulated correctly.
"""
import time
start_time = time.time()

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import os.path

import sys

import ROOT
from ROOT import TG4Event, TFile, TTree, TCanvas, TGraph, TH1D

from optparse import OptionParser
import xml.etree.ElementTree as ET

from array import array
from math import cos, sin
import random

with open("UserConfig.py") as infile:
    exec(infile.read())
import pyGeoEff

from ROOT import gROOT # for creating the output file

ROOT.gROOT.ProcessLine('#include<vector>') # for outputting data to ROOT file

#########
# INPUT
#########
a = str(sys.argv[1])
#print(a)
#print(sys.argv[1])

# this works interactively
#edep_file = TFile("edep.*.root")
edep_file = TFile(a, "OPEN")
if not edep_file:
  print ("Error: could not open file", a)
  sys.exit()
inputTree = edep_file.Get("EDepSimEvents")
event = TG4Event()
inputTree.SetBranchAddress("Event", event)
entries = inputTree.GetEntriesFast()
gDecayZ = TGraph(27, OffAxisPoints, meanPDPZ)
#print(entries)

#########
# OUTPUT
#########

# create directory for plots to be stored if it doesn't already exist
out_path = "/dune/app/users/weishi/testn2fd/DUNE_ND_GeoEff/app/output"
plotpath = out_path + "/plots"
rootpath = out_path + "/root_out" # for output ROOT files
if not os.path.exists( out_path):
    os.makedirs( out_path)
    print("out_path '" + out_path + "' did not exist. It has been created!")
if not os.path.exists( plotpath):
    os.makedirs( plotpath)
    print("plotpath '" + plotpath + "' did not exist. It has been created!")
if not os.path.exists( rootpath):
    os.makedirs( rootpath)
    print("rootpath '" + rootpath + "' did not exist. It has been created!")
print(" \n") # separate output

# Create the ROOT file that will hold the output of this script
# including but not limited to the output TH1D objects and raw data
f_out = TFile('{0}/n2fd_paired_out.root'.format(rootpath), 'RECREATE')
myEvents = TTree('myEvents', 'myEvents')
maxEdeps = 100000 # max number of edeps for initialize arrays

##################################
# LArBath evt info
##################################
nEdeps = array('i', [0])
myEvents.Branch('nEdeps', nEdeps, 'nEdeps/I')
deps_E_MeV = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_E_MeV', deps_E_MeV, 'deps_E_MeV[nEdeps]/F')
deps_start_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_start_t_us', deps_start_t_us, 'deps_start_t_us[nEdeps]/F')
deps_stop_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_stop_t_us', deps_stop_t_us, 'deps_stop_t_us[nEdeps]/F')
larbath_vtx_cm = array('f', 3*[0.0])
myEvents.Branch('larbath_vtx_cm', larbath_vtx_cm, 'larbath_vtx_cm[3]/F')
# edeps generated in LArBath (start point)
larbath_deps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_x_cm', larbath_deps_start_x_cm, 'larbath_deps_start_x_cm[nEdeps]/F') # larbath edeps x
larbath_deps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_y_cm', larbath_deps_start_y_cm, 'larbath_deps_start_y_cm[nEdeps]/F')
larbath_deps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_z_cm', larbath_deps_start_z_cm, 'larbath_deps_start_z_cm[nEdeps]/F')
# edeps generated in LArBath (stop point)
larbath_deps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_x_cm', larbath_deps_stop_x_cm, 'larbath_deps_stop_x_cm[nEdeps]/F')
larbath_deps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_y_cm', larbath_deps_stop_y_cm, 'larbath_deps_stop_y_cm[nEdeps]/F')
larbath_deps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_z_cm', larbath_deps_stop_z_cm, 'larbath_deps_stop_z_cm[nEdeps]/F')

##################################
# ND paired evt info
##################################
# Random thrown vertex in ND (paired evt)
throwVtx_nd_cm = array('f', 3*[0.0])
myEvents.Branch('throwVtx_nd_cm', throwVtx_nd_cm, 'throwVtx_nd_cm[3]/F')
# Edeps in ND random throw (start points)
ndthrow_deps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_x_cm', ndthrow_deps_start_x_cm, 'ndthrow_deps_start_x_cm[nEdeps]/F')
ndthrow_deps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_y_cm', ndthrow_deps_start_y_cm, 'ndthrow_deps_start_y_cm[nEdeps]/F')
ndthrow_deps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_z_cm', ndthrow_deps_start_z_cm, 'ndthrow_deps_start_z_cm[nEdeps]/F')
# Edeps in ND random throw (stop points)
ndthrow_deps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_x_cm', ndthrow_deps_stop_x_cm, 'ndthrow_deps_stop_x_cm[nEdeps]/F')
ndthrow_deps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_y_cm', ndthrow_deps_stop_y_cm, 'ndthrow_deps_stop_y_cm[nEdeps]/F')
ndthrow_deps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_z_cm', ndthrow_deps_stop_z_cm, 'ndthrow_deps_stop_z_cm[nEdeps]/F')

# Edeps in ND with earth curvature correction
throwVtx_nd_ecc_cm = array('f', 3*[0.0])
myEvents.Branch('throwVtx_nd_ecc_cm', throwVtx_nd_ecc_cm, 'throwVtx_nd_ecc_cm[3]/F')
ndthrow_ecc_deps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_x_cm', ndthrow_ecc_deps_start_x_cm, 'ndthrow_ecc_deps_start_x_cm[nEdeps]/F')
ndthrow_ecc_deps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_y_cm', ndthrow_ecc_deps_start_y_cm, 'ndthrow_ecc_deps_start_y_cm[nEdeps]/F')
ndthrow_ecc_deps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_z_cm', ndthrow_ecc_deps_start_z_cm, 'ndthrow_ecc_deps_start_z_cm[nEdeps]/F')
ndthrow_ecc_deps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_x_cm', ndthrow_ecc_deps_stop_x_cm, 'ndthrow_ecc_deps_stop_x_cm[nEdeps]/F')
ndthrow_ecc_deps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_y_cm', ndthrow_ecc_deps_stop_y_cm, 'ndthrow_ecc_deps_stop_y_cm[nEdeps]/F')
ndthrow_ecc_deps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_z_cm', ndthrow_ecc_deps_stop_z_cm, 'ndthrow_ecc_deps_stop_z_cm[nEdeps]/F')

##################################
# FD paired evt info
##################################
# Random thrown vertex in FD
throwVtx_fd_cm = array('f', 3*[0.0])
myEvents.Branch('throwVtx_fd_cm', throwVtx_fd_cm, 'throwVtx_fd_cm[3]/F')
# Edeps in FD random throw (start points)
fdthrow_deps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_x_cm', fdthrow_deps_start_x_cm, 'fdthrow_deps_start_x_cm[nEdeps]/F')
fdthrow_deps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_y_cm', fdthrow_deps_start_y_cm, 'fdthrow_deps_start_y_cm[nEdeps]/F')
fdthrow_deps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_z_cm', fdthrow_deps_start_z_cm, 'fdthrow_deps_start_z_cm[nEdeps]/F')
# Edeps in FD random throw (stop points)
fdthrow_deps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_x_cm', fdthrow_deps_stop_x_cm, 'fdthrow_deps_stop_x_cm[nEdeps]/F')
fdthrow_deps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_y_cm', fdthrow_deps_stop_y_cm, 'fdthrow_deps_stop_y_cm[nEdeps]/F')
fdthrow_deps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_z_cm', fdthrow_deps_stop_z_cm, 'fdthrow_deps_stop_z_cm[nEdeps]/F')

###########################
# Loop over edepsim events
##########################

for jentry in range(41):
    print("jentry = " + str(jentry))
    nb = inputTree.GetEntry(jentry) #number of bytes read
    #print("nb =" + str(nb))
    #print("event number: ", event.EventId)
    #print("number of primaries: ", event.Primaries.size())
    #print("number of trajectories: ", event.Trajectories.size())
    #print("number of segments: ", event.SegmentDetectors.size())
    all_dep_startpos_list = list()
    all_dep_stoppos_list = list()
    all_edep_list   = list()
    all_dep_starttime_list = list()
    all_dep_stoptime_list = list()
    had_dep_pos_list = list()
    had_edep_list   = list()

    # initialize LArBath vertex
    larbath_vtx_cm[0] = 0; larbath_vtx_cm[1] = 0; larbath_vtx_cm[2] = 0;
    nEdeps[0] = 0
    throwVtx_nd_cm[0] = 0; throwVtx_nd_cm[1] = 0; throwVtx_nd_cm[2] = 0;
    throwVtx_nd_ecc_cm[0] = 0; throwVtx_nd_ecc_cm[1] = 0; throwVtx_nd_ecc_cm[2] = 0;
    throwVtx_fd_cm[0] = 0; throwVtx_fd_cm[1] = 0; throwVtx_fd_cm[2] = 0;

    for primary in event.Primaries:
        #print("number of particles: ", primary.Particles.size())
        for ipart, particle in enumerate(primary.Particles):
            #print("ipart here = " +str(ipart))
            PDGCode = particle.GetPDGCode()
            #print("pdgcode: ", PDGCode)
            if abs(PDGCode) >= 11 and abs(PDGCode) <= 16:
                posx = primary.GetPosition().X() * edep2cm
                posy = primary.GetPosition().Y() * edep2cm
                posz = primary.GetPosition().Z() * edep2cm
                momx = particle.GetMomentum().X()
                momy = particle.GetMomentum().Y()
                momz = particle.GetMomentum().Z()
                trackIDforPrimary = particle.GetTrackId()
                """
                print("pos x: ", posx, " y: ", posy, " z: ", posz, " [cm]")
                print("mom x: ", momx, " y: ", momy, " z: ", momz, " [MeV]")
                print("trackIDforPrimary: ", trackIDforPrimary)
                """

    for containerName, hitSegments in event.SegmentDetectors:
        # iHit is the index, hitSEgment is the data stored at the index in the second item in event.SegementDectectors
        for iHit, hitSegment in enumerate(hitSegments):
            primaryID = hitSegment.GetPrimaryId()
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
            edep = hitSegment.GetEnergyDeposit()

            all_dep_startpos_list.append(edep_start_x)
            all_dep_startpos_list.append(edep_start_y)
            all_dep_startpos_list.append(edep_start_z)
            all_dep_starttime_list.append(edep_start_t)
            all_dep_stoppos_list.append(edep_stop_x)
            all_dep_stoppos_list.append(edep_stop_y)
            all_dep_stoppos_list.append(edep_stop_z)
            all_dep_stoptime_list.append(edep_stop_t)
            all_edep_list.append(edep)

            # Hadronic part of edepsim
            # Later need these hadronic edeps to evaluate hadronic veto
            if primaryID != trackIDforPrimary:
                had_dep_pos_list.append(edep_x)
                had_dep_pos_list.append(edep_y)
                had_dep_pos_list.append(edep_z)
                had_edep_list.append(edep)

    # for use in processing events before and after transformations
    nEdeps[0] = len(all_edep_list)

    ##################################################################
    # Initialize geometric efficiency module to manipulate energy deps
    # Only consider ND LAr on axis
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

    # Near detector active dimensions for hadronic veto
    geoEff.setActiveX(NDActiveVol_min[0], NDActiveVol_max[0])
    geoEff.setActiveY(NDActiveVol_min[1], NDActiveVol_max[1])
    geoEff.setActiveZ(NDActiveVol_min[2], NDActiveVol_max[2])

    # Far detector active dimensions for hadronic veto
    geoEff.setFDActiveX(FDActiveVol_min[0], FDActiveVol_max[0])
    geoEff.setFDActiveY(FDActiveVol_min[1], FDActiveVol_max[1])
    geoEff.setFDActiveZ(FDActiveVol_min[2], FDActiveVol_max[2])

    # Range for random translation throws in ND fiducial volume
    geoEff.setRangeX(ND_FV_min[0], ND_FV_max[0])
    geoEff.setRangeY(ND_FV_min[1], ND_FV_max[1])
    geoEff.setRangeZ(ND_FV_min[2], ND_FV_max[2])

    # Range for random translation throws in FD fiducial volume
    geoEff.setRangeXFD(FD_FV_min[0], FD_FV_max[0])
    geoEff.setRangeYFD(FD_FV_min[1], FD_FV_max[1])
    geoEff.setRangeZFD(FD_FV_min[2], FD_FV_max[2])

    # Set offset between ND MC coordinate system and volumes defined above.
    # Is this still needed for Alex's Genie Gen??? To be validated/discussed
    geoEff.setOffsetX(NDLAr_OnAxis_offset[0])
    geoEff.setOffsetY(NDLAr_OnAxis_offset[1])
    geoEff.setOffsetZ(NDLAr_OnAxis_offset[2])

    # Original Genie Gen event vtx, not guaranteed to be at (0,0,0)!
    # Later we need the vertex to calculate the rotation axis when event is randomly put at a new position in ND
    geoEff.setVertex(posx, posy, posz) #cm

    # Interpolate event neutrino production point (beam coordinate)
    # Input needs to be in unit of meters
    decayZbeamCoord = gDecayZ.Eval( posx/100. - NDLAr_OnAxis_offset[0]/100. - detRefBeamCoord[0] )

    # Calculate neutrino production point in detector coordinate
    decayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0]
    decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation)
    decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation)

    # Set production point in unit: cm
    # Later we need this decay position to calculate the rotation axis when event is randomly put at a new position in ND
    geoEff.setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.)

    # Randomly throw the Genie Gen event in ND
    tot_nd_throw = 0

    # If after max_nd_throws still don't pass at nd, stop and move to next event (otherwise too much computing resources)
    while tot_nd_throw < max_nd_throws:
        print ("-- tot nd throw:", tot_nd_throw)

        ####################################
        # Only do one throw in ND at a time
        ####################################
        geoEff.setNthrows(1)
        geoEff.throwTransforms() # this randomly generates new vtx position and a rotation angle w.r.t. the neutrino direction

        # Get the randomly generated vtx x, y, z, and the angle
        throwVtxX_nd  = geoEff.getCurrentThrowTranslationsX() # cm
        throwVtxY_nd  = geoEff.getCurrentThrowTranslationsY() # cm
        throwVtxZ_nd  = geoEff.getCurrentThrowTranslationsZ() # cm
        throwAngle = geoEff.getCurrentThrowRotations()

        # Require random thrown vtx pos outside in ND dead regions
        if ( IsInNDFV( throwVtxX_nd[0] - NDLAr_OnAxis_offset[0], throwVtxY_nd[0] - NDLAr_OnAxis_offset[1], throwVtxZ_nd[0] - NDLAr_OnAxis_offset[2]) ):
            print ("-- nd throw", tot_nd_throw, "is in nd FV")

            # Interpolate neutrino production point (beam coordinate) for the random throw, unit meter
            RandomthrowdecayZbeamCoord = gDecayZ.Eval( throwVtxX_nd[0]/100. - NDLAr_OnAxis_offset[0]/100. - detRefBeamCoord[0] )

            # Calculate neutrino production point in detector coordinate, unit meter
            RandomthrowdecayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0]
            RandomthrowdecayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( RandomthrowdecayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation)
            RandomthrowdecayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( RandomthrowdecayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation)

            # We have a new decay position because the above random throw changes x in ND
            # Need this to calculate the rotation axis for the random throw
            geoEff.setDecayPos4RandomThrowX(RandomthrowdecayXdetCoord*100., RandomthrowdecayYdetCoord*100., RandomthrowdecayZdetCoord*100.)

            # Set the hadronic part list of edeps to evaluate ND hadronic veto
            # NOTE: in the future may need to add lepton selection at ND as well (to be discussed)
            geoEff.setHitSegEdeps(had_edep_list)
            geoEff.setHitSegPoss(had_dep_pos_list)
            ndrandthrowresulthad = geoEff.getNDContainment4RandomThrowX() # returns a struct

            if (ndrandthrowresulthad.containresult[0][0][0] != 0):
                print ("-- nd throw", tot_nd_throw, "passed nd had veto")
                ###########################################################
                # Random throw passed ND hadronic veto !
                ###########################################################

                # Now change to the full list of edeps start points
                # the random thrown x/y/z/angle should remain the same because throw was done above already
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_dep_startpos_list)
                # And call the function again to get new transformed positions for all edeps
                ndrandthrowresultall_start = geoEff.getNDContainment4RandomThrowX()

                # Repeat for edepsim stop points !!!
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_dep_stoppos_list)
                ndrandthrowresultall_stop = geoEff.getNDContainment4RandomThrowX()

                ###########################################################
                # Now translate this event vertex to ND det coordinate (0,0,0) (and the edeps accordingly)
                # Do it for all edeps and also had part edeps (because later we need had part only to evaluate FD had veto)
                ###########################################################

                # First tell the module where is the random thrown vertex
                geoEff.setNDrandVertex(throwVtxX_nd[0], throwVtxY_nd[0], throwVtxZ_nd[0])
                print ("-- nd throw x: ", throwVtxX_nd[0], "y: ", throwVtxY_nd[0], ", z: ", throwVtxZ_nd[0])
                all_startposdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresultall_start.thrownEdepspos[0]) # returns Eigen::Matrix3Xf
                all_stopposdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresultall_stop.thrownEdepspos[0]) # repeat for stop points
                had_posdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresulthad.thrownEdepspos[0])

                ####################################################################################################################
                # Apply earth curvature correction to translate into FD coordinate system, vtx now at (0,0,0) in FD det coordinate,
                # this info will be used by both leg 1 and leg 2 transformations below
                ####################################################################################################################
                had_posdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(had_posdep_ndorig_matrix, beamLineRotation) # returns Eigen::Matrix3Xf
                all_startposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_startposdep_ndorig_matrix, beamLineRotation)
                all_stopposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_stopposdep_ndorig_matrix, beamLineRotation) # repeat for stop points

                ####################################################################################################################
                # Here is where things branch out, we have two kinds of pairs (two legs)
                #
                # leg 1 needed by Alex:
                #       this randomly rotated ECC event is now at FD origin,
                #       then randomly rotated it around the beam axis w.r.t. FD as if the event is at FD,
                #       we further translated it randomly in ND LAr and it needs to pass ND LAr hadronic veto.
                # leg 2 needed by Radi:
                #       same FD event obtained in leg 1 (after random rot), but ND event is without ECC (already obtained above)
                ####################################################################################################################

                ##################################################
                # leg 1: paired ND evt with ECC needed by Alex
                ##################################################
                # throw a random rotation angle around FD beam axis and a random point in ND LAr,
                # rotate the event and translate vtx from origin to that point and evaluate ND hadronic veto
                # if pass, continue to leg 1
                # if not, throw another one
                # if after max throws no one passed, continue to next nd throw, do not even bother the leg 1

                tot_nd_ecc_throw = 0

                while tot_nd_ecc_throw < max_nd_throws:
                    print ("---- tot nd ecc throw:", tot_nd_ecc_throw)
                    geoEff.setNthrowsNDECC(1)
                    geoEff.throwTransformsNDECC() # this randomly generates a rotation angle and a new vtx position in ND LAr

                    throwVtxX_nd_ecc = geoEff.getCurrentNDECCThrowTranslationsX()
                    throwVtxY_nd_ecc = geoEff.getCurrentNDECCThrowTranslationsY()
                    throwVtxZ_nd_ecc = geoEff.getCurrentNDECCThrowTranslationsZ()
                    throwAngle_nd_ecc = geoEff.getCurrentNDECCThrowRotations()

                    print ("---- nd ecc throw x: ", throwVtxX_nd_ecc[0], "y: ", throwVtxY_nd_ecc[0], ", z: ", throwVtxZ_nd_ecc[0], ", angle: ", throwAngle_nd_ecc[0])

                    geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list to evaluate hadronic veto
                    # Here no need to setHitSegPoss because it's feed in next to moveBack2ndVertex

                    # Rotate around fd beam axis and then translate the nd ecc event back to a randomly thrown position in ND LAr
                    ndeccrandthrowresulthad = geoEff.moveBack2ndVertex(had_posdep_fdorig_matrix, beamLineRotation)

                    if (ndeccrandthrowresulthad.containresult[0][0][0] != 0):
                        print ("---- nd ecc throw", tot_nd_ecc_throw, "passed nd had veto")

                        # Now change to the full list of edeps start points
                        # the random thrown x/y/z should remain the same because throw was done above already
                        geoEff.setHitSegEdeps(all_edep_list)
                        ndeccrandthrowresultall_start= geoEff.moveBack2ndVertex(all_startposdep_fdorig_matrix, beamLineRotation)

                        # Repeat for edepsim stop points !!!
                        geoEff.setHitSegEdeps(all_edep_list)
                        ndeccrandthrowresultall_stop = geoEff.moveBack2ndVertex(all_stopposdep_fdorig_matrix, beamLineRotation)

                        print ("Found paired ndecc-fd random thrown events")

                        print ("Breaking nd ecc throw loop")
                        break

                    else:
                        print ("-- nd ecc throw", tot_nd_ecc_throw, "failed fd had veto!")

                    # indentation is important!
                    tot_nd_ecc_throw = tot_nd_ecc_throw + 1

                if tot_nd_ecc_throw == max_nd_throws:
                    print ("Reached max nd ecc throw", max_nd_throws, ", continue to next nd throw")
                    tot_nd_throw = tot_nd_throw + 1
                    continue

                ######################################################
                # If it gets this far, the event already
                #   1) has a random rotation around FD beam axis,
                #   2) placed at a random position in ND LAr,
                #   3) passed ND LAr hadronic veto
                ######################################################

                ##################################################
                # leg 2: paired FD evt needed by both Radi and Alex
                #        here only involves a translation in fd
                ##################################################
                # Tell the module where the vertex is in FD
                # We pass the above random position in ND LAr, this doesn't matter as we are just doing translation
                geoEff.setVertexFD(throwVtxX_nd_ecc[0], throwVtxY_nd_ecc[0], throwVtxZ_nd_ecc[0])

                tot_fd_throw = 0

                while tot_fd_throw < max_fd_throws:
                    print ("---- tot fd throw:", tot_fd_throw)
                    ##########################################################################################
                    # Below do random throw (translate only) in FD similar to ND: only one throw in FD at a time
                    ##########################################################################################
                    geoEff.setNthrowsFD(1)
                    geoEff.throwTransformsFD() # this randomly generates new vtx position in FD FV

                    throwVtxX_fd = geoEff.getCurrentFDThrowTranslationsX()
                    throwVtxY_fd = geoEff.getCurrentFDThrowTranslationsY()
                    throwVtxZ_fd = geoEff.getCurrentFDThrowTranslationsZ()

                    # Check if it passes FD hadronic veto
                    geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list
                    fdrandthrowresulthad = geoEff.getFDContainment4RandomThrow(ndeccrandthrowresulthad.thrownEdepspos[0])

                    if (fdrandthrowresulthad.containresult[0][0][0] != 0):
                        print ("---- fd throw", tot_fd_throw, "passed fd had veto")
                        print ("---- fd throw x: ", throwVtxX_fd[0], "y: ", throwVtxY_fd[0], ", z: ", throwVtxZ_fd[0])
                        ###########################################################
                        # FD rand throw passes veto, write paired evt info
                        ###########################################################

                        # Now change to the full list of edeps
                        # the random thrown x/y/z should reamin the same because throw is done above already
                        geoEff.setHitSegEdeps(all_edep_list)
                        fdrandthrowresultall_start = geoEff.getFDContainment4RandomThrow(ndeccrandthrowresultall_start.thrownEdepspos[0])
                        # Repeat for edepsim stop points !!!
                        geoEff.setHitSegEdeps(all_edep_list)
                        fdrandthrowresultall_stop = geoEff.getFDContainment4RandomThrow(ndeccrandthrowresultall_stop.thrownEdepspos[0])

                        print ("Found paired nd-fd random thrown events")

                        #################################
                        # Unpack info and store to output
                        #################################
                        print ("Saving...")

                        # LArBath info
                        deps_E_MeV[:nEdeps[0]] = np.array(all_edep_list, dtype=np.float32)
                        deps_start_t_us[:nEdeps[0]] = np.array(all_dep_starttime_list, dtype=np.float32)
                        deps_stop_t_us[:nEdeps[0]] = np.array(all_dep_stoptime_list, dtype=np.float32)
                        larbath_vtx_cm[0] = posx
                        larbath_vtx_cm[1] = posy
                        larbath_vtx_cm[2] = posz
                        larbath_deps_start_x_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[::3], dtype=np.float32) # every 3 element: x list
                        larbath_deps_start_y_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[1::3], dtype=np.float32) # y list
                        larbath_deps_start_z_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[2::3], dtype=np.float32) # z list
                        larbath_deps_stop_x_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[::3], dtype=np.float32)
                        larbath_deps_stop_y_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[1::3], dtype=np.float32)
                        larbath_deps_stop_z_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[2::3], dtype=np.float32)

                        # Paired event in ND from random throw
                        throwVtx_nd_cm[0] = throwVtxX_nd[0]
                        throwVtx_nd_cm[1] = throwVtxY_nd[0]
                        throwVtx_nd_cm[2] = throwVtxZ_nd[0]
                        ndthrow_deps_start_x_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32)
                        ndthrow_deps_start_y_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_deps_start_z_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                        ndthrow_deps_stop_x_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                        ndthrow_deps_stop_y_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_deps_stop_z_cm[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                        # Paired event in ND from random throw with earth curvature corrected (as FD)
                        throwVtx_nd_ecc_cm[0] = throwVtxX_nd_ecc[0]
                        throwVtx_nd_ecc_cm[1] = throwVtxY_nd_ecc[0]
                        throwVtx_nd_ecc_cm[2] = throwVtxZ_nd_ecc[0]
                        ndthrow_ecc_deps_start_x_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32) # ecc: earth curvature corrected
                        ndthrow_ecc_deps_start_y_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_ecc_deps_start_z_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                        ndthrow_ecc_deps_stop_x_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                        ndthrow_ecc_deps_stop_y_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_ecc_deps_stop_z_cm[:nEdeps[0]] = np.array(ndeccrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                        # Paired event in FD from random throw
                        throwVtx_fd_cm[0] = throwVtxX_fd[0]
                        throwVtx_fd_cm[1] = throwVtxY_fd[0]
                        throwVtx_fd_cm[2] = throwVtxZ_fd[0]
                        fdthrow_deps_start_x_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32)
                        fdthrow_deps_start_y_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                        fdthrow_deps_start_z_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                        fdthrow_deps_stop_x_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                        fdthrow_deps_stop_y_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                        fdthrow_deps_stop_z_cm[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                        # event level
                        myEvents.Fill()

                        # Break the while loop, move on to next evt
                        print ("Paired data saved, breaking fd throw loop")
                        break
                    else:
                        print ("---- fd throw", tot_fd_throw, "failed fd had veto!")

                    # indentation is important!
                    # if don't, put it in another random FD pos...until it passes FD veto
                    tot_fd_throw = tot_fd_throw + 1

                # if reached max fd throw and still didn't pass FD veto, try next nd throw
                if tot_fd_throw == max_fd_throws:
                    print ("Reached max fd throw", max_fd_throws, ", continue to next nd throw")
                    tot_nd_throw = tot_nd_throw + 1
                    continue

                # if found paired fd evts, break nd throw loop
                print ("And breaking nd throw loop")
                break

            else:
                print ("-- nd throw", tot_nd_throw, "failed nd had veto!")

        else:
            print ("-- nd throw", tot_nd_throw, "outside nd FV!")

        # indentation is important!
        tot_nd_throw = tot_nd_throw + 1

f_out.cd()
myEvents.Write()

print("\n")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
