find_package(ROOT)

if ( ${ROOT_FOUND} )
  # Use example EDepSim file to produce libraries
  execute_process(COMMAND root -b -q ${CMAKE_CURRENT_SOURCE_DIR}/EDepSimForMakeProject_20200424.root ${CMAKE_CURRENT_SOURCE_DIR}/makeEDepSimProject.C)
  file(GLOB EDEPSIMHEADERS ${geoEff_SOURCE_DIR}/EDepSimEvents/*.h)

  file(COPY ${geoEff_SOURCE_DIR}/EDepSimEvents/EDepSimEvents.so DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  add_library(EDepSimEvents SHARED IMPORTED)
  set_property(TARGET EDepSimEvents PROPERTY IMPORTED_LOCATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/EDepSimEvents.so)
  
  file(COPY ${EDEPSIMHEADERS} DESTINATION ${geoEff_SOURCE_DIR}/include)
  
  include_directories(${ROOT_INCLUDE_DIRS})	
  add_executable(runGeoEffEDepSimBeamDir runGeoEffEDepSimBeamDir.cpp)

  target_link_libraries(runGeoEffEDepSimBeamDir geoEff EDepSimEvents ${ROOT_LIBRARIES})	
else ()
  Message("Couldn't find either ROOT. Skipping runGeoEffEDepSimBeamDir app. Library will still compile.")
endif()
