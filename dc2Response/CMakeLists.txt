add_library(
  dc2Response SHARED
  src/Aeff.cxx
  src/DC2.cxx
  src/Edisp.cxx
  src/IrfLoader.cxx
  src/loadIrfs.cxx
  src/Psf.cxx
  src/PsfScaling.cxx
)

target_include_directories(
  dc2Response PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  dc2Response
  PRIVATE st_facilities
  irfInterface
  irfUtil
  astro
  tip
  CLHEP::CLHEP
)

###### Tests ######
add_executable(test_dc2Response src/test/test.cxx)
target_include_directories(
  test_dc2Response PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:>
)
target_link_libraries(
  test_dc2Response
  PRIVATE astro irfInterface dc2Response
)

###############################################################
# Installation
###############################################################
install(DIRECTORY dc2Response DESTINATION ${FERMI_INSTALL_INCLUDEDIR})

install(
  TARGETS dc2Response test_dc2Response
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
