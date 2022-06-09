add_library(
  dc1aResponse 
  src/Aeff.cxx
  src/Edisp.cxx
  src/IrfLoader.cxx
  src/loadIrfs.cxx
  src/Psf.cxx
)

target_include_directories(
  dc1aResponse PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  dc1aResponse
  PUBLIC irfInterface
  PRIVATE st_facilities irfUtil astro tip CLHEP::CLHEP
)

###### Tests ######
add_executable(test_dc1aResponse src/test/test.cxx)
target_include_directories(
  test_dc1aResponse PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:>
)
target_link_libraries(
  test_dc1aResponse
  PRIVATE astro irfInterface dc1aResponse
)

###############################################################
# Installation
###############################################################
install(DIRECTORY dc1aResponse DESTINATION ${FERMI_INSTALL_INCLUDEDIR})

install(
  TARGETS dc1aResponse test_dc1aResponse
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
