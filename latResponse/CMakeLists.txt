add_library(
  latResponse STATIC
  src/Aeff.cxx
  src/AeffEpochDep.cxx
  src/Bilinear.cxx
  src/Edisp.cxx
  src/Edisp2.cxx
  src/Edisp3.cxx
  src/EdispEpochDep.cxx
  src/EdispInterpolator.cxx
  src/EfficiencyFactor.cxx
  src/EfficiencyFactorEpochDep.cxx
  src/EpochDep.cxx
  src/FitsTable.cxx
  src/IrfLoader.cxx
  src/Irfs.cxx
  src/ParTables.cxx
  src/Psf.cxx
  src/Psf2.cxx
  src/Psf3.cxx
  src/PsfBase.cxx
  src/PsfEpochDep.cxx
  src/PsfIntegralCache.cxx
)
add_library(Fermitools::latResponse ALIAS latResponse)

target_include_directories(
  latResponse PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  latResponse PRIVATE
  astro
  CLHEP::CLHEP
  irfInterface
  irfUtil
  st_facilities
  st_stream
  tip
)

###### Tests ######
add_executable(test_latResponse src/test/test.cxx)
target_include_directories(
  test_latResponse PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:>
)
target_link_libraries(
  test_latResponse
  PRIVATE astro irfInterface CppUnit::CppUnit latResponse
)

###############################################################
# Installation
###############################################################
install(DIRECTORY latResponse DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_REFDATADIR}/latResponse)
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})

install(
  TARGETS latResponse test_latResponse
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
