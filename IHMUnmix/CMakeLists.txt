cmake_minimum_required(VERSION 3.5)

project(IHMUnmix LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 COMPONENTS Widgets PrintSupport REQUIRED)

find_package(libUnmix REQUIRED)

set(LIBUNMIX_DEBUG ${LIBUNMIX_LIBRARY})
string(REPLACE Debug Release LIBUNMIX_RELEASE ${LIBUNMIX_DEBUG})

# Use generator expressions, so variable can be used not only by target_link_libraries.
set(LIBUNMIX_LIBRARY
        "$<$<NOT:$<CONFIG:DEBUG>>:${LIBUNMIX_RELEASE}>"
        "$<$<CONFIG:DEBUG>:${LIBUNMIX_DEBUG}>"
)

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

if (UNIX)
	find_package(Armadillo REQUIRED)
elseif (WIN32)
    # Please set path to armadillo include and link dirs (relative or absolute)
	set (ARMADILLO_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/../libUnmix/include_arma")
    set (ARMADILLO_LIBRARIES "${PROJECT_SOURCE_DIR}/../libUnmix/ext_lib/blas_win64_MT.lib" "${PROJECT_SOURCE_DIR}/../libUnmix/ext_lib/lapack_win64_MT.lib" )
endif()
include_directories(${LIBUNMIX_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS})

add_executable(IHMUnmix
    main.cpp
    fenprimaire.h
    fenprimaire.cpp
    fensynthpixel.h
    fensynthpixel.cpp
    fenresult.h
    fenresult.cpp
    fenimagehyper.cpp
    fenimagehyper.h
    fensynthpixel.ui
    FenImageHyper.ui
  )

target_link_libraries(IHMUnmix PRIVATE Qt5::Widgets Qt5::PrintSupport ${LIBUNMIX_LIBRARY} ${LIB_QCUSTOMPLOT} ${LINK_LIBS}
    ${HDF5_CXX_STATIC_LIBRARY} ${ARMADILLO_LIBRARIES})

#qt5_create_translation(QM_FILES ${CMAKE_SOURCE_DIR} ${TS_FILES})
