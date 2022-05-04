add_library(
  irfLoader SHARED
  src/Loader.cxx
  src/Registrar.cxx
)
add_library(Fermitools::irfLoader ALIAS irfLoader)

target_include_directories(
  irfLoader PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  irfLoader
  PUBLIC facilities st_app st_stream st_facilities irfInterface dc1aResponse
  dc1Response dc2Response latResponse testResponse tip
)

if(NOT APPLE)
  target_compile_definitions(irfLoader PUBLIC TRAP_FPE)
endif()

###### Executables ######
add_executable(gtirfs src/gtirfs/gtirfs.cxx)
target_link_libraries(gtirfs PRIVATE irfLoader)

###### Tests ######
add_executable(test_irfLoader src/test/test.cxx)
target_link_libraries(test_irfLoader PRIVATE PkgConfig::cppunit irfLoader)

###############################################################
# Installation
###############################################################
install(DIRECTORY irfLoader DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})

install(
  TARGETS irfLoader gtirfs test_irfLoader
  # EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
