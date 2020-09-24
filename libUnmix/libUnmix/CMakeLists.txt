
# Find the QtWidgets library
find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

add_compile_definitions(QCUSTOMPLOT_USE_LIBRARY)
if (UNIX)
	set(LINK_LIB_PLOT debug "${CMAKE_CURRENT_SOURCE_DIR}/libqcustomplotd.so"
	optimized "${CMAKE_CURRENT_SOURCE_DIR}/libqcustomplot.so")
elseif (WIN32)
	set(LINK_LIB_PLOT debug "${CMAKE_CURRENT_SOURCE_DIR}/qcustomplotd2.lib"
	optimized "${CMAKE_CURRENT_SOURCE_DIR}/qcustomplot2.lib")
endif()

set (FIND_HDF_COMPONENTS C CXX HL HLCXX  static)
find_package (HDF5 NAMES "hdf5" COMPONENTS ${FIND_HDF_COMPONENTS})
if (HDF5_FOUND)
  if (HDF5_static_C_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})
  endif ()
  if (HDF5_static_CXX_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_CXX_STATIC_LIBRARY})
  endif ()
else ()
  message (FATAL_ERROR " HDF5 is Required")
endif ()

find_package(TIFF REQUIRED)

if (UNIX)
    find_package(Armadillo REQUIRED)
elseif (WIN32)
 # Please set path to armadillo include and link dirs (relative or absolute)
    set (ARMADILLO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include_arma")
    set (ARMADILLO_LIBRARIES "${PROJECT_SOURCE_DIR}/ext_lib/blas_win64_MT.lib" "${PROJECT_SOURCE_DIR}/ext_lib/lapack_win64_MT.lib" )
endif()
message(STATUS "${ARMADILLO_INCLUDE_DIRS}")
include_directories(${HDF5_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS})
add_library(libUnmix SHARED
  include/libUnmix_global.h
  src/libImage.cpp
  include/libImage.h
  src/libpixelsynth.cpp
  include/libpixelsynth.h
  src/libmcmc.cpp
  include/libmcmc.h
  src/libgraph.cpp
  include/libgraph.h
  include/qcustomplot.h
  include/unmixing.h
  include/readH5File.h
)

message(STATUS "QCustomPlot : ${LINK_LIB_PLOT}")
target_link_libraries (libUnmix PRIVATE Qt5::Widgets Qt5::PrintSupport ${LINK_LIBS} ${HDF5_CXX_STATIC_LIBRARY}
    ${TIFF_LIBRARIES} ${ARMADILLO_LIBRARIES} ${LINK_LIB_PLOT})
set_target_properties(libUnmix PROPERTIES PUBLIC_HEADER "include/qcustomplot.h;include/unmixing.h;include/libImage.h;include/libpixelsynth.h;include/libmcmc.h;include/libgraph.h;include/libUnmix_global.h;include/readH5File.h")

if (UNIX)
	install(TARGETS libUnmix
	  # IMPORTANT: Add the libUnmix library to the "export-set"
	  EXPORT libUnmixTargets
	  LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT shlib
	PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}
	COMPONENT dev)
elseif(WIN32)
	install(TARGETS libUnmix 
	EXPORT libUnmixTargets
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
	PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR})	
endif()


if (WIN32)
    add_definitions(-DLIBUNMIX_LIBRARY)
endif()

if (UNIX)
	install(FILES "libqcustomplotd.so" CONFIGURATIONS Debug DESTINATION ${INSTALL_LIB_DIR})
	install(FILES "libqcustomplot.so" CONFIGURATIONS Release DESTINATION ${INSTALL_LIB_DIR})
elseif(WIN32)
	install(FILES "qcustomplotd2.lib" CONFIGURATIONS Debug DESTINATION ${INSTALL_LIB_DIR})
	install(FILES "qcustomplot2.lib" CONFIGURATIONS Release DESTINATION ${INSTALL_LIB_DIR})
endif()


