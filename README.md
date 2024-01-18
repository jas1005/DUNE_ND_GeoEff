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

# ND numu to FD numu event pair

```
source setup.sh

cd app
python3 Edepsim_ana.py /dune/data/users/awilkins/extrapolation/edep.LArBath.NDGenieGen.root

# The first time you may need to install a few packages via pip install, depending on what it complains when you run, e.g.:
pip install --target=/dune/app/users/weishi/python3libs torch --upgrade
pip install --target=/dune/app/users/weishi/python3libs scipy --upgrade
```

# ND numu to FD nue event pair

It reflects the [workflow](https://indico.fnal.gov/event/62304/contributions/280309/attachments/173208/234357/Numu2nue.pdf).

```
source setup.sh

# Keeps non leptonic deposits at the earth curvature correction step and extracts lepton kinematics
python3 Edepsim_ana_nueapp_stage1.py /dune/data/users/awilkins/extrapolation/edep.LArBath.NDGenieGen.root

# Runs electron edepsim based on lepton kinematics from stage 1
# it assumes to look for the txt file in ./output/nue_output
./Electron_gen_edepsim.sh n2fd_nueapp_paired_stage1.txt

# Collate energy deposits of electron to those non-leptonic ones from stage 1
python3 Edepsim_ana_nueapp_stage2.py  n2fd_nueapp_paired_stage1.root
```
