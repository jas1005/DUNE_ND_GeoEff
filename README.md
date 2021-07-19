# Instruction for running translation and rotations on FD n-tuples

```
cd /dune/app/users/<your_username>
mkdir NDEff (first time only)
cd NDEff
git clone --recurse-submodules -b FD_Wei https://github.com/weishi10141993/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
source setup.sh
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j pyGeoEff
```
