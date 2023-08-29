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
deps_E = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_E', deps_E, 'deps_E[nEdeps]/F')
larbath_vtx = array('f', 3*[0.0])
myEvents.Branch('larbath_vtx', larbath_vtx, 'larbath_vtx[3]/F')
# edeps generated in LArBath (start point)
larbath_deps_start_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_x', larbath_deps_start_x, 'larbath_deps_start_x[nEdeps]/F') # larbath edeps x
larbath_deps_start_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_y', larbath_deps_start_y, 'larbath_deps_start_y[nEdeps]/F')
larbath_deps_start_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_z', larbath_deps_start_z, 'larbath_deps_start_z[nEdeps]/F')
# edeps generated in LArBath (stop point)
larbath_deps_stop_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_x', larbath_deps_stop_x, 'larbath_deps_stop_x[nEdeps]/F')
larbath_deps_stop_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_y', larbath_deps_stop_y, 'larbath_deps_stop_y[nEdeps]/F')
larbath_deps_stop_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_z', larbath_deps_stop_z, 'larbath_deps_stop_z[nEdeps]/F')

##################################
# ND paired evt info
##################################
# Random thrown vertex in ND (paired evt)
throwVtx_nd = array('f', 3*[0.0])
myEvents.Branch('throwVtx_nd', throwVtx_nd, 'throwVtx_nd[3]/F')
# Edeps in ND random throw (start points)
ndthrow_deps_start_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_x', ndthrow_deps_start_x, 'ndthrow_deps_start_x[nEdeps]/F')
ndthrow_deps_start_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_y', ndthrow_deps_start_y, 'ndthrow_deps_start_y[nEdeps]/F')
ndthrow_deps_start_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_start_z', ndthrow_deps_start_z, 'ndthrow_deps_start_z[nEdeps]/F')
# Edeps in ND random throw (stop points)
ndthrow_deps_stop_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_x', ndthrow_deps_stop_x, 'ndthrow_deps_stop_x[nEdeps]/F')
ndthrow_deps_stop_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_y', ndthrow_deps_stop_y, 'ndthrow_deps_stop_y[nEdeps]/F')
ndthrow_deps_stop_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_deps_stop_z', ndthrow_deps_stop_z, 'ndthrow_deps_stop_z[nEdeps]/F')

# Edeps in ND with earth curvature correction
# this share the same vertex throwVtx_nd
ndthrow_ecc_deps_start_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_x', ndthrow_ecc_deps_start_x, 'ndthrow_ecc_deps_start_x[nEdeps]/F')
ndthrow_ecc_deps_start_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_y', ndthrow_ecc_deps_start_y, 'ndthrow_ecc_deps_start_y[nEdeps]/F')
ndthrow_ecc_deps_start_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_start_z', ndthrow_ecc_deps_start_z, 'ndthrow_ecc_deps_start_z[nEdeps]/F')
ndthrow_ecc_deps_stop_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_x', ndthrow_ecc_deps_stop_x, 'ndthrow_ecc_deps_stop_x[nEdeps]/F')
ndthrow_ecc_deps_stop_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_y', ndthrow_ecc_deps_stop_y, 'ndthrow_ecc_deps_stop_y[nEdeps]/F')
ndthrow_ecc_deps_stop_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('ndthrow_ecc_deps_stop_z', ndthrow_ecc_deps_stop_z, 'ndthrow_ecc_deps_stop_z[nEdeps]/F')

##################################
# FD paired evt info
##################################
# Random thrown vertex in FD
throwVtx_fd = array('f', 3*[0.0])
myEvents.Branch('throwVtx_fd', throwVtx_fd, 'throwVtx_fd[3]/F')
# Edeps in FD random throw (start points)
fdthrow_deps_start_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_x', fdthrow_deps_start_x, 'fdthrow_deps_start_x[nEdeps]/F')
fdthrow_deps_start_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_y', fdthrow_deps_start_y, 'fdthrow_deps_start_y[nEdeps]/F')
fdthrow_deps_start_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_start_z', fdthrow_deps_start_z, 'fdthrow_deps_start_z[nEdeps]/F')
# Edeps in FD random throw (stop points)
fdthrow_deps_stop_x = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_x', fdthrow_deps_stop_x, 'fdthrow_deps_stop_x[nEdeps]/F')
fdthrow_deps_stop_y = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_y', fdthrow_deps_stop_y, 'fdthrow_deps_stop_y[nEdeps]/F')
fdthrow_deps_stop_z = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fdthrow_deps_stop_z', fdthrow_deps_stop_z, 'fdthrow_deps_stop_z[nEdeps]/F')

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
    all_startposdep_list = list()
    all_stopposdep_list = list()
    all_edep_list   = list()
    had_posdep_list = list()
    had_edep_list   = list()

    # initialize LArBath vertex
    larbath_vtx[0] = 0; larbath_vtx[1] = 0; larbath_vtx[2] = 0;
    nEdeps[0] = 0
    throwVtx_nd[0] = 0; throwVtx_nd[1] = 0; throwVtx_nd[2] = 0;
    throwVtx_fd[0] = 0; throwVtx_fd[1] = 0; throwVtx_fd[2] = 0;

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
            edep_stop_x = hitSegment.GetStop().X() * edep2cm
            edep_stop_y = hitSegment.GetStop().Y() * edep2cm
            edep_stop_z = hitSegment.GetStop().Z() * edep2cm
            edep_x = (edep_start_x + edep_stop_x)/2 # use this for hadronic veto
            edep_y = (edep_start_y + edep_stop_y)/2
            edep_z = (edep_start_z + edep_stop_z)/2
            edep = hitSegment.GetEnergyDeposit()

            all_startposdep_list.append(edep_start_x)
            all_startposdep_list.append(edep_start_y)
            all_startposdep_list.append(edep_start_z)
            all_stopposdep_list.append(edep_stop_x)
            all_stopposdep_list.append(edep_stop_y)
            all_stopposdep_list.append(edep_stop_z)
            all_edep_list.append(edep)

            # Hadronic part of edepsim
            # Later need these hadronic edeps to evaluate hadronic veto
            if primaryID != trackIDforPrimary:
                had_posdep_list.append(edep_x)
                had_posdep_list.append(edep_y)
                had_posdep_list.append(edep_z)
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

    # If after max_throws still don't pass at nd, stop and move to next event (otherwise too much computing resources)
    while tot_nd_throw < max_throws:
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
            geoEff.setHitSegPoss(had_posdep_list)
            ndrandthrowresulthad = geoEff.getNDContainment4RandomThrowX() # returns a struct

            if (ndrandthrowresulthad.containresult[0][0][0] != 0):
                print ("-- nd throw", tot_nd_throw, "passed nd had veto")
                ###########################################################
                # Random throw passed ND hadronic veto !
                ###########################################################

                # Now change to the full list of edeps start points
                # the random thrown x/y/z/angle should remain the same because throw was done above already
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_startposdep_list)
                # And call the function again to get new transformed positions for all edeps
                ndrandthrowresultall_start = geoEff.getNDContainment4RandomThrowX()

                # Repeat for edepsim stop points !!!
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_stopposdep_list)
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

                # Apply earth curvature correction to translate into FD coordinate system, vtx now at (0,0,0) in FD det coordinate
                had_posdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(had_posdep_ndorig_matrix, beamLineRotation) # returns Eigen::Matrix3Xf
                all_startposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_startposdep_ndorig_matrix, beamLineRotation)
                all_stopposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_stopposdep_ndorig_matrix, beamLineRotation) # repeat for stop points

                # Translate the event vertex back to the randomly thrown ND position (this step doesn't interfere with next steps)
                # Do we care hadronic veto here???
                all_startposdep_nd_earthcurvature_matrix = geoEff.moveBack2ndVertex(all_startposdep_fdorig_matrix) # returns Eigen::Matrix3Xf
                all_stopposdep_nd_earthcurvature_matrix = geoEff.moveBack2ndVertex(all_stopposdep_fdorig_matrix) # returns Eigen::Matrix3Xf

                # Tell the module where the vertex is in FD
                geoEff.setVertexFD(0, 0, 0) # it's at FD origin because we moved it to ND origin and then just rotated at there

                tot_fd_throw = 0

                while tot_fd_throw < max_throws:
                    print ("---- tot fd throw:", tot_fd_throw)
                    # Below do random throw (translate only) in FD similar to ND:
                    geoEff.setNthrowsFD(1)
                    geoEff.throwTransformsFD() # this randomly generates new vtx position in FD FV

                    throwVtxX_fd = geoEff.getCurrentFDThrowTranslationsX()
                    throwVtxY_fd = geoEff.getCurrentFDThrowTranslationsY()
                    throwVtxZ_fd = geoEff.getCurrentFDThrowTranslationsZ()

                    # Check if it passes FD hadronic veto
                    geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list
                    fdrandthrowresulthad = geoEff.getFDContainment4RandomThrow(had_posdep_fdorig_matrix)

                    if (fdrandthrowresulthad.containresult[0][0][0] != 0):
                        print ("---- fd throw", tot_fd_throw, "passed fd had veto")
                        print ("---- fd throw x: ", throwVtxX_fd[0], "y: ", throwVtxY_fd[0], ", z: ", throwVtxZ_fd[0])
                        ###########################################################
                        # FD rand throw passes veto, write paired evt info
                        ###########################################################

                        # Now change to the full list of edeps
                        # the random thrown x/y/z should reamin the same because throw is done above already
                        geoEff.setHitSegEdeps(all_edep_list)
                        fdrandthrowresultall_start = geoEff.getFDContainment4RandomThrow(all_startposdep_fdorig_matrix)
                        # Repeat for edepsim stop points !!!
                        geoEff.setHitSegEdeps(all_edep_list)
                        fdrandthrowresultall_stop = geoEff.getFDContainment4RandomThrow(all_stopposdep_fdorig_matrix)

                        print ("Found paired nd-fd random thrown events, saving...")

                        """
                        nEdeps = len(all_edep_list)
                        print ("nEdeps:", nEdeps)

                        for iDep in range(len(all_edep_list)):
                            # Save edeps and their positions in nd and fd
                            print ("iDep #", iDep, ": E [MeV]:", all_edep_list[iDep], ", original genie gen pos[cm]:", all_startposdep_list[iDep*3:(iDep)*3+3], ", paired nd pos [cm]:", ndrandthrowresultall_start.thrownEdepspos[0][:, iDep], ", paired fd pos [cm]:", fdrandthrowresultall_start.thrownEdepspos[0][:, iDep])
                        """

                        #################################
                        # Unpack info and store to output
                        #################################

                        # LArBath info
                        larbath_vtx[0] = posx
                        larbath_vtx[1] = posy
                        larbath_vtx[2] = posz
                        deps_E[:nEdeps[0]] = np.array(all_edep_list, dtype=np.float32)
                        larbath_deps_start_x[:nEdeps[0]] = np.array(all_startposdep_list[::3], dtype=np.float32) # every 3 element: x list
                        larbath_deps_start_y[:nEdeps[0]] = np.array(all_startposdep_list[1::3], dtype=np.float32) # y list
                        larbath_deps_start_z[:nEdeps[0]] = np.array(all_startposdep_list[2::3], dtype=np.float32) # z list
                        larbath_deps_stop_x[:nEdeps[0]] = np.array(all_stopposdep_list[::3], dtype=np.float32)
                        larbath_deps_stop_y[:nEdeps[0]] = np.array(all_stopposdep_list[1::3], dtype=np.float32)
                        larbath_deps_stop_z[:nEdeps[0]] = np.array(all_stopposdep_list[2::3], dtype=np.float32)

                        # Paired event in ND from random throw
                        throwVtx_nd[0] = throwVtxX_nd[0]
                        throwVtx_nd[1] = throwVtxY_nd[0]
                        throwVtx_nd[2] = throwVtxZ_nd[0]
                        ndthrow_deps_start_x[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32)
                        ndthrow_deps_start_y[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_deps_start_z[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                        ndthrow_deps_stop_x[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                        ndthrow_deps_stop_y[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                        ndthrow_deps_stop_z[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                        # Paired event in ND from random throw with earth curvature corrected (as FD)
                        # vertx is the same as throwVtx_nd
                        ndthrow_ecc_deps_start_x[:nEdeps[0]] = np.array(all_startposdep_nd_earthcurvature_matrix[0,:], dtype=np.float32) # ecc: earth curvature corrected
                        ndthrow_ecc_deps_start_y[:nEdeps[0]] = np.array(all_startposdep_nd_earthcurvature_matrix[1,:], dtype=np.float32)
                        ndthrow_ecc_deps_start_z[:nEdeps[0]] = np.array(all_startposdep_nd_earthcurvature_matrix[2,:], dtype=np.float32)
                        ndthrow_ecc_deps_stop_x[:nEdeps[0]] = np.array(all_stopposdep_nd_earthcurvature_matrix[0,:], dtype=np.float32) # ecc: earth curvature corrected
                        ndthrow_ecc_deps_stop_y[:nEdeps[0]] = np.array(all_stopposdep_nd_earthcurvature_matrix[1,:], dtype=np.float32)
                        ndthrow_ecc_deps_stop_z[:nEdeps[0]] = np.array(all_stopposdep_nd_earthcurvature_matrix[2,:], dtype=np.float32)

                        # Paired event in FD from random throw
                        throwVtx_fd[0] = throwVtxX_fd[0]
                        throwVtx_fd[1] = throwVtxY_fd[0]
                        throwVtx_fd[2] = throwVtxZ_fd[0]
                        fdthrow_deps_start_x[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32)
                        fdthrow_deps_start_y[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                        fdthrow_deps_start_z[:nEdeps[0]] = np.array(fdrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                        fdthrow_deps_stop_x[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                        fdthrow_deps_stop_y[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                        fdthrow_deps_stop_z[:nEdeps[0]] = np.array(fdrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                        # event level
                        myEvents.Fill()


                        # had dep only
                        #for jDep in range(len(had_edep_list)):
                        #    print ("Had dep #", jDep, ": E [MeV]:", had_edep_list[jDep], ", paired nd pos [cm]:", ndrandthrowresulthad.thrownEdepspos[0][:, jDep], ", paired fd pos [cm]:", fdrandthrowresulthad.thrownEdepspos[0][:, jDep])

                        # these are the random throwed x y z in ND for all edeps
                        # for print out/debug purposes
                        #geoEff.getCurrentThrowDepsX(0) # we only have one throw (0th throw) ==> array x[ith edep]
                        #geoEff.getCurrentThrowDepsY(0)
                        #geoEff.getCurrentThrowDepsZ(0)

                        # Break the while loop, move on to next evt
                        print ("Paired data saved, breaking fd throw loop")
                        break
                    else:
                        print ("---- fd throw", tot_fd_throw, "failed fd had veto!")

                    # indentation is important!
                    # if don't, put it in another random FD pos...until it passes FD veto
                    tot_fd_throw = tot_fd_throw + 1

                print ("Breaking nd throw loop")
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
