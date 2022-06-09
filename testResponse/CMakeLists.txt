add_library(
  testResponse 
  src/Aeff.cxx
  src/Edisp.cxx
  src/IrfLoader.cxx
  src/loadIrfs.cxx
  src/Psf.cxx
)
add_library(Fermitools::testResponse ALIAS testResponse)

target_include_directories(
  testResponse PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  testResponse
  PRIVATE astro irfInterface irfUtil tip st_facilities CLHEP::CLHEP
)

###### Tests ######
add_executable(test_testResponse src/test/test.cxx)
target_link_libraries(
  test_testResponse
  PRIVATE testResponse irfInterface astro PkgConfig::cppunit
)

###############################################################
# Installation
###############################################################
install(DIRECTORY testResponse DESTINATION ${FERMI_INSTALL_INCLUDEDIR})

install(
  TARGETS testResponse test_testResponse
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
