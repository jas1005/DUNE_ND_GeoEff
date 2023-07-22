#! /usr/bin/env python
"""
Runs over ROOT files created by edepsim
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import os.path
import sys
import ROOT
from ROOT import TG4Event, TFile, TGraph
from optparse import OptionParser
import xml.etree.ElementTree as ET
from array import array
from math import cos, sin

with open("UserConfig.py") as infile:
    exec(infile.read())
import pyGeoEff

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
#print(entries)

trajectories_dtype = np.dtype([("pdgID", "i4"), ("trackID", "u4")]) #i4 - integer 4, u4 - unsigned integer(sort of abs value)

gDecayZ = ROOT.TGraph(27, OffAxisPoints, meanPDPZ)

# Loop over events
for jentry in range(entries):
    print("jentry = " +str(jentry))
    nb = inputTree.GetEntry(jentry) #number of bytes read
    #print("nb =" + str(nb))
    #print("event number: ", event.EventId)
    #print("number of primaries: ", event.Primaries.size())
    #print("number of trajectories: ", event.Trajectories.size())
    #print("number of segments: ", event.SegmentDetectors.size())
    all_edep_x_list=list()
    all_edep_y_list=list()
    all_edep_z_list=list()
    all_edep_list=list()
    all_posdep_list=list()
    had_edep_x_list=list()
    had_edep_y_list=list()
    had_edep_z_list=list()
    had_edep_list=list()
    had_posdep_list=list()

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
                #print("pos x: ", posx, " y: ", posy, " z: ", posz, " [cm]")
                #print("mom x: ", momx, " y: ", momy, " z: ", momz, " [MeV]")
                #print("trackIDforPrimary: ", trackIDforPrimary)

    trajectories = np.empty(len(event.Trajectories), dtype=trajectories_dtype)
    for iTraj, trajectory in enumerate(event.Trajectories):
        start_pt, end_pt = trajectory.Points[0], trajectory.Points[-1]
        trajectories[iTraj]["pdgID"] = trajectory.GetPDGCode()
        trajectories[iTraj]["trackID"] = trajectory.GetTrackId()

    for containerName, hitSegments in event.SegmentDetectors:
        # iHit is the index, hitSEgment is the data stored at the index in the second item in event.SegementDectectors
        for iHit, hitSegment in enumerate(hitSegments):
            primaryID = hitSegment.GetPrimaryId()
            # Energy deposit from primary particles
            all_edep_x = hitSegment.GetStop().X() * edep2cm
            all_edep_y = hitSegment.GetStop().Y() * edep2cm
            all_edep_z = hitSegment.GetStop().Z() * edep2cm
            # Later need these hadronic edeps to evaluate hadronic veto
            all_posdep_list.append(all_edep_x)
            all_posdep_list.append(all_edep_y)
            all_posdep_list.append(all_edep_z)
            all_edep = hitSegment.GetEnergyDeposit()
            all_edep_x_list.append(all_edep_x)
            all_edep_y_list.append(all_edep_y)
            all_edep_z_list.append(all_edep_z)
            all_edep_list.append(all_edep)
            # Hadronic part of edepsim
            # Later need these hadronic edeps to evaluate hadronic veto
            if primaryID != trackIDforPrimary:
                # Energy deposit from primary particles
                had_edep_x = hitSegment.GetStop().X() * edep2cm
                had_edep_y = hitSegment.GetStop().Y() * edep2cm
                had_edep_z = hitSegment.GetStop().Z() * edep2cm
                had_posdep_list.append(had_edep_x)
                had_posdep_list.append(had_edep_y)
                had_posdep_list.append(had_edep_z)
                had_edep = hitSegment.GetEnergyDeposit()
                had_edep_x_list.append(had_edep_x)
                had_edep_y_list.append(had_edep_y)
                had_edep_z_list.append(had_edep_z)
                had_edep_list.append(had_edep)

    #edep_array=np.array(had_edep_list)

    # Create plot of charge deposits per event
    plotpath = "plots"
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    plt.scatter(all_edep_x_list, all_edep_y_list, s=5, c=all_edep_list)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('[MeV]')
    plt.xlabel('x (cm)')
    plt.ylabel('y (cm)')
    plt.xlim([-dimension/2, dimension/2])
    plt.ylim([-dimension/2, dimension/2])
    plt.savefig("plots/evt_{0}_xy.png".format(jentry))
    plt.clf() # clear figure
    plt.close()
    plt.scatter(all_edep_z_list, all_edep_x_list, s=5, c=all_edep_list)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('[MeV]')
    plt.xlabel('z (cm)')
    plt.ylabel('x (cm)')
    plt.xlim([-dimension/2, dimension/2])
    plt.ylim([-dimension/2, dimension/2])
    plt.savefig("plots/evt_{0}_zx.png".format(jentry))
    plt.clf()
    plt.close()
    plt.scatter(all_edep_z_list, all_edep_y_list, s=5, c=all_edep_list)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('[MeV]')
    plt.xlabel('z (cm)')
    plt.ylabel('y (cm)')
    plt.xlim([-dimension/2, dimension/2])
    plt.ylim([-dimension/2, dimension/2])
    plt.savefig("plots/evt_{0}_zy.png".format(jentry))
    plt.clf()
    plt.close()

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    p = ax.scatter(all_edep_z_list, all_edep_x_list, all_edep_y_list, s=20, c=all_edep_list)
    cbar = fig.colorbar(p)
    cbar.set_label('[MeV]')
    ax.set_xlabel('z (cm)')
    ax.set_ylabel('x (cm)')
    ax.set_zlabel('y (cm)')
    ax.set_xlim([-dimension/2, dimension/2])
    ax.set_ylim([-dimension/2, dimension/2])
    ax.set_zlim([-dimension/2, dimension/2])
    fig.savefig("plots/evt_{0}_xyz.png".format(jentry))
    plt.clf()
    plt.close()

    # Initialize geometric efficiency module to manipulate energy deps
    # Only consider ND LAr on axis
    geoEff = pyGeoEff.geoEff(123, False)
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
    geoEff.setRangeXFD(FD_FV_min[0], FD_FV_min[0])
    geoEff.setRangeYFD(FD_FV_min[1], FD_FV_min[1])
    geoEff.setRangeZFD(FD_FV_min[2], FD_FV_min[2])
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
        # Only do one throw in ND
        geoEff.setNthrows(1)
        geoEff.throwTransforms() # this randomly generates new vtx position and a rotation angle w.r.t. the neutrino direction
        # Get the randomly generated vtx x, y, z, and the angle
        throwVtxX = geoEff.getCurrentThrowTranslationsX() #cm
        throwVtxY = geoEff.getCurrentThrowTranslationsY()
        throwVtxZ = geoEff.getCurrentThrowTranslationsZ()
        throwAngle = geoEff.getCurrentThrowRotations()

        # Require random thrown vtx pos outside in ND dead regions
        if ( IsInNDFV(throwVtxX[0] - NDLAr_OnAxis_offset[0], throwVtxY[0] - NDLAr_OnAxis_offset[1], throwVtxZ[0] - NDLAr_OnAxis_offset[2]) ):
            print ("-- nd throw", tot_nd_throw, "is in nd FV")
            # Interpolate neutrino production point (beam coordinate) for the random throw, unit meter
            RandomthrowdecayZbeamCoord = gDecayZ.Eval( throwVtxX[0]/100. - NDLAr_OnAxis_offset[0]/100. - detRefBeamCoord[0] )

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
                # Random throw passed ND hadronic veto

                # Now change to the full list of edeps (leptonic + hadronic)
                # the random thrown x/y/z/angle should reamin the same because throw is done above already
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_posdep_list)
                # And call the function again to get new transformed positions for all edeps
                ndrandthrowresultall = geoEff.getNDContainment4RandomThrowX()
                # Tell the module where is the random thrown vertex
                geoEff.setNDrandVertex(throwVtxX[0], throwVtxY[0], throwVtxZ[0])

                # Now translate this event vertex to ND det coordinate (0,0,0) (and the edeps accordingly)
                # Do it for all edeps and also had part edeps (because later we need had part only to evaluate FD had veto)
                all_posdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresultall.thrownEdepspos[0]) # returns Eigen::Matrix3Xf
                had_posdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresulthad.thrownEdepspos[0])

                # Apply earth curvature correction to translate into FD coordinate system, vtx now at (0,0,0) in FD det coordinate
                had_posdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(had_posdep_ndorig_matrix, beamLineRotation) #returns Eigen::Matrix3Xf
                all_posdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_posdep_ndorig_matrix, beamLineRotation)

                # Tell the module where the vertex is in FD
                geoEff.setVertexFD(0,0,0) # it's at 000 because we moved it to ND 000 and then just rotated at 000

                tot_fd_throw = 0

                while tot_fd_throw < max_throws:
                    print ("---- tot fd throw:", tot_fd_throw)
                    # Below do random throw (translate only) in FD similar to ND:
                    geoEff.setNthrowsFD(1)
                    geoEff.throwTransformsFD() # this randomly generates new vtx position in FD FV

                    # Check if it passes FD hadronic veto
                    geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list
                    fdrandthrowresulthad = geoEff.getFDContainment4RandomThrow(had_posdep_fdorig_matrix)

                    if (fdrandthrowresulthad.containresult[0][0][0] != 0):
                        print ("---- fd throw", tot_fd_throw, "passed fd had veto")
                        # FD rand throw passes veto, write paired evt info

                        # Now change to the full list of edeps (leptonic + hadronic)
                        # the random thrown x/y/z should reamin the same because throw is done above already
                        geoEff.setHitSegEdeps(all_edep_list)
                        fdrandthrowresultall = geoEff.getFDContainment4RandomThrow(all_posdep_fdorig_matrix)

                        print ("Found paired nd-fd random thrown events, saving...")

                        nEdeps = len(all_edep_list)
                        print ("nEdeps:", nEdeps)

                        for iDep in range(len(all_edep_list)):
                            # Save edeps and their positions in nd and fd
                            print ("iDep #", iDep, ": E [MeV]:", all_edep_list[iDep], ", original genie gen pos[cm]:", all_posdep_list[iDep*3:(iDep)*3+3], ", paired nd pos [cm]:", ndrandthrowresultall.thrownEdepspos[0][:, iDep], ", paired fd pos [cm]:", fdrandthrowresultall.thrownEdepspos[0][:, iDep])

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
