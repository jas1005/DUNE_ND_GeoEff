add_library(geoEff SHARED geoEff.cpp)

pybind11_add_module(pyGeoEff pyGeoEff.cpp)
target_link_libraries(pyGeoEff PRIVATE geoEff)
