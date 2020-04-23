#include_directories(${geoEff_SOURCE_DIR}/include)
add_executable(runGeoEff runGeoEff.cpp)
target_link_libraries(runGeoEff geoEff)
