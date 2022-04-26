####################################
# Create target and set properties #
####################################

include(UseSWIG)

set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
set(CMAKE_SWIG_FLAGS -w362 -w389 -w401 -w402 -w503 -w509)

### SWIG generates the python wrapper library
set_property(SOURCE src/pyIrfLoader.i PROPERTY CPLUSPLUS ON)
swig_add_library(pyIrfLoader LANGUAGE python SOURCES src/pyIrfLoader.i)
set_property(TARGET pyIrfLoader PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES TRUE)

target_include_directories(
  pyIrfLoader
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:>
)

target_link_libraries(pyIrfLoader PRIVATE irfLoader Python3::Python)

if(APPLE)
  target_link_options(pyIrfLoader PRIVATE -bundle -flat_namespace)
endif()

###############################################################
# Installation
###############################################################
install(DIRECTORY python/ DESTINATION ${FERMI_INSTALL_PYTHON})

install(
  TARGETS pyIrfLoader
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/pyIrfLoader.py
  DESTINATION ${FERMI_INSTALL_PYTHON}
)
