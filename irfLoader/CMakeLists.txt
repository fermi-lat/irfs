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
  irfLoader PUBLIC
  facilities st_app st_stream st_facilities irfInterface latResponse testResponse tip
)

target_compile_definitions(irfLoader PRIVATE ST_USER)

###### Executables ######
add_executable(gtirfs src/gtirfs/gtirfs.cxx)
target_link_libraries(gtirfs PRIVATE irfLoader)

###### Tests ######
add_executable(test_irfLoader src/test/test.cxx)
target_link_libraries(test_irfLoader PRIVATE PkgConfig::cppunit irfLoader)
