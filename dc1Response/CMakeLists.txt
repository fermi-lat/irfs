add_library(
  dc1Response 
  src/AeffDC1.cxx
  src/DC1.cxx
  src/EdispDC1.cxx
  src/IrfLoader.cxx
  src/loadIrfs.cxx
  src/PsfDC1.cxx
)

target_include_directories(
  dc1Response PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  dc1Response
  PRIVATE st_facilities
  irfInterface
  irfUtil
  astro
  tip
  CLHEP::CLHEP
)

###### Tests ######
add_executable(test_dc1Response src/test/test.cxx)
target_include_directories(
  test_dc1Response PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/src
  $<INSTALL_INTERFACE:>
)
target_link_libraries(
  test_dc1Response
  PRIVATE astro irfInterface dc1Response
)

###############################################################
# Installation
###############################################################
install(DIRECTORY dc1Response DESTINATION ${FERMI_INSTALL_INCLUDEDIR})

install(
  TARGETS dc1Response test_dc1Response
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
