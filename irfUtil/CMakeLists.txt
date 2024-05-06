add_library(
  irfUtil STATIC
  src/EventTypeMapper.cxx
  src/HdCaldb.cxx
  src/HDgtcalf.c
  src/HDgtcalf_internal.c
  src/HDgtcalf_standalone.c
  src/IrfHdus.cxx
  src/Util.cxx
)
add_library(Fermitools::irfUtil ALIAS irfUtil)

target_include_directories(
  irfUtil PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:>
)

target_link_libraries(
  irfUtil
  PRIVATE astro tip st_facilities st_facilities cfitsio::cfitsio
)

###### Tests ######
add_executable(test_irfUtil src/test/test.cxx)
add_executable(p8_test_irfUtil src/test/p8_test.cxx)
target_link_libraries(test_irfUtil PRIVATE irfUtil tip)
target_link_libraries(p8_test_irfUtil PRIVATE irfUtil tip)

###############################################################
# Installation
###############################################################
install(DIRECTORY irfUtil DESTINATION ${FERMI_INSTALL_INCLUDEDIR})

install(
  TARGETS irfUtil test_irfUtil p8_test_irfUtil
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
