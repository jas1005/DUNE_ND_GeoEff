find_package(ROOT)
find_package(Geant4)

if ( ${ROOT_FOUND} AND ${Geant4_FOUND} )
  # Use example EDepSim file to produce libraries
  execute_process(COMMAND root -b -q ${CMAKE_CURRENT_SOURCE_DIR}/EDepSimForMakeProject_20200424.root ${CMAKE_CURRENT_SOURCE_DIR}/makeEDepSimProject.C)
  file(GLOB EDEPSIMHEADERS ${geoEff_SOURCE_DIR}/EDepSimEvents/*.h)

  file(COPY ${geoEff_SOURCE_DIR}/EDepSimEvents/EDepSimEvents.so DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#  file(RENAME ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/EDepSimEvents.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libEDepSimEvents.so)
  add_library(EDepSimEvents SHARED IMPORTED)
  set_property(TARGET EDepSimEvents PROPERTY IMPORTED_LOCATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/EDepSimEvents.so)
  
  file(COPY ${EDEPSIMHEADERS} DESTINATION ${geoEff_SOURCE_DIR}/include)
  
  include_directories(${ROOT_INCLUDE_DIRS} ${Geant4_INCLUDE_DIRS})
  add_executable(runGeoEffEDepSimBeamDir runGeoEffEDepSimBeamDir.cpp)

  # Why is this necessary?
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  target_link_libraries(runGeoEffEDepSimBeamDir geoEff EDepSimEvents ${ROOT_LIBRARIES} ${Geant4_LIBRARIES})
else ()
  Message("Couldn't find either ROOT or Geant4. Skipping runGeoEffEDepSimBeamDir app. Library will still compile.")
endif()
