# Instruction for running translation and rotations on FD n-tuples

[First time only]
```
cd /dune/app/users/<your_username>
mkdir NDEff (first time only)
cd NDEff
git clone --recurse-submodules -b FD_Wei https://github.com/weishi10141993/DUNE_ND_GeoEff.git      # Get geoEff library
cd DUNE_ND_GeoEff
source setup.sh                                                                                    # Necessary setups for build
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j geoEff                                                                                     # Build geoEff (can also use: make -j pyGeoEff)
```

To (re)compile and (re)run,
```
cd /dune/app/users/<your_username>/NDEff/DUNE_ND_GeoEff/
#
# In case you log out, need to source setup.sh to setup ROOT
#
source setup.sh                                                                           
cd app
make runGeoEffFDEvtSim                                                                             # Compile program
cd ../bin
./runGeoEffFDEvtSim                                                                                # Run program
```
this will produce a root file containing throws and the hadron throw result. 
