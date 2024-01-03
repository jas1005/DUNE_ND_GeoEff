# Set up

[First time only]
```
cd /dune/app/users/<your_username>
mkdir testn2fd
cd testn2fd

git clone --recurse-submodules -b N2FD https://github.com/weishi10141993/DUNE_ND_GeoEff.git      # Get geoEff library
# Note for git version (git --version) before 2.13, use: git clone --recursive -b N2FD https://github.com/weishi10141993/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
source setup.sh                                                                                  # Necessary setups for build
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j geoEff                                                                                   # Build geoEff
make -j pyGeoEff                                                                                 # Build pygeoEff
```

# Analyze edepsim data

```
source setup.sh

cd app
nohup python3 Edepsim_ana.py /dune/data/users/awilkins/extrapolation/edep.LArBath.NDGenieGen.root >& output.log &

# The first time you may need to install a few packages via pip install, depending on what it complains when you run, e.g.:
pip install --target=/dune/app/users/weishi/python3libs torch --upgrade
pip install --target=/dune/app/users/weishi/python3libs scipy --upgrade
```

## Numu to nue event translation: paired energy deposits

It reflects the [workflow](https://indico.fnal.gov/event/62304/contributions/280309/attachments/173208/234357/Numu2nue.pdf).

```
python3 Edepsim_ana_nueapp_stage1.py /dune/data/users/awilkins/extrapolation/edep.LArBath.NDGenieGen.root 
```

## Set up edep-sim

Once obtained a GDML file (below use ```LArBath.gdml``` as example), run edep-sim to simulate event energy deposits.

```
# Setup edep-sim
mkdir edepsim
cd edepsim

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup geant4 v4_10_6_p01e -q e20:prof
setup edepsim v3_2_0 -q e20:prof

# The following config file generates muons with energy of 10GeV
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/PRISM/N2FD/Gen_muon.mac

edep-sim \
    -C \
    -g LArBath.gdml \
    -o edep.LArBath.electron.root \
    -u \
    -e 1 \
    Gen_electron.mac

# Run options refer to https://github.com/ClarkMcGrew/edep-sim#running-the-detector-simulation
```
