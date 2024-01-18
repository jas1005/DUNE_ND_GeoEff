#!/bin/bash
shopt -s expand_aliases
# need expand alias for those defined in edepsim below to work

################################################################################
# Script to run marley + edepsim on dunegpvm
################################################################################
# Options
GEOMETRY="LArBath.gdml"
EDEP_MAC="Gen_electron.mac"
################################################################################

echo "Run electron edepsim based on stage 1 lepton kinematics"

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

# Dump the current clean environment to to env.sh so it can be restored when needed
echo "Saving environment env.sh"
declare -px > env.sh

# Setting up for gensim stuff
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup geant4 v4_10_6_p01e -q e20:prof
setup edepsim v3_2_0 -q e20:prof

# this directory should exist already from stage 1
cp $EDEP_MAC output/nue_output
cp $GEOMETRY output/nue_output
cd output/nue_output
ls

echo 'Start read txt file:'
while read -r iEVENT LEP_E LEP_PX LEP_PY LEP_PZ; do
  # skip header at first line
  if [[ $iEVENT == "event_number" ]]; then
     continue
  fi
  echo "EVT #: ${iEVENT} E: $LEP_E PX: $LEP_PX PY: $LEP_PY PZ: $LEP_PZ"
  # Change lepton kinematics event by event
  cp $EDEP_MAC Gen_electron_${iEVENT}.mac
  sed -i "s/lep_e/$LEP_E/g" Gen_electron_${iEVENT}.mac
  sed -i "s/lep_px/$LEP_PX/g" Gen_electron_${iEVENT}.mac
  sed -i "s/lep_py/$LEP_PY/g" Gen_electron_${iEVENT}.mac
  sed -i "s/lep_pz/$LEP_PZ/g" Gen_electron_${iEVENT}.mac

  # For given energy electron, generate one event by default
  # Same energy electron may have variations, but hopefully these can be captured by large training data
  # Events don't pass nd throw/selection will have zero lep energy written in txt file, edepsim still runs on such events, but they shouldn't be used
  echo "Running: edep-sim -C -g $GEOMETRY -o edep_LArBath_electron_evtnb_${iEVENT}.root -u -e 1 Gen_electron_${iEVENT}.mac"
  edep-sim -C -g $GEOMETRY -o edep_LArBath_electron_evtnb_${iEVENT}.root -u -e 1 Gen_electron_${iEVENT}.mac

done < $1
