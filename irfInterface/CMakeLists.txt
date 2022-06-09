add_library(
  irfInterface 
  src/IEdisp.cxx
  src/IEfficiencyFactor.cxx
  src/IPsf.cxx
  src/IrfRegistry.cxx
  src/IrfsFactory.cxx
)
add_library(Fermitools::irfInterface ALIAS irfInterface)

target_include_directories(
  irfInterface PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  irfInterface
  PUBLIC st_facilities
  PRIVATE astro tip CLHEP::CLHEP
)

###### Tests ######
add_executable(test_irfInterface src/test/Edisp.cxx src/test/Psf.cxx src/test/test.cxx)
target_link_libraries(test_irfInterface PRIVATE astro PkgConfig::cppunit irfInterface)

###############################################################
# Installation
###############################################################
install(DIRECTORY irfInterface DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_DATADIR})

install(
  TARGETS irfInterface test_irfInterface
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
