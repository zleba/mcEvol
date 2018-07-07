cmake_minimum_required (VERSION 2.6)
#set (CMAKE_CXX_STANDARD 11)
project (mcEvolution C CXX Fortran)
#enable_language(Fortran)

cmake_policy(SET CMP0015 NEW)

execute_process(COMMAND root-config --incdir  OUTPUT_VARIABLE ROOTINC)
execute_process(COMMAND root-config --libdir  OUTPUT_VARIABLE ROOTLIB)

set(YODAdir YODA/install)

include_directories(inc ${YODAdir}/include/ ${ROOTINC})

link_directories(${YODAdir}/lib  ${ROOTLIB})

#find_library(MYYODA YODA)
find_library(MYYODA
    NAMES YODA
    HINTS "/afs/desy.de/user/z/zlebcr/h1/TMD/mcEvol/YODA/install/lib"
    )

if(MYYODA STREQUAL "MYYODA-NOTFOUND")
    message(FATAL_ERROR "Unable to find YODA libraries")
else()
    message(STATUS "YODA found ${MYYODA}")
endif()


#find_package(Boost REQUIRED iostreams)
#find_package(Root)


set(SOURCES src/alphaSpline.cpp  src/main.cpp    src/SplineGen.cpp   src/SplittingsLO.cpp   src/SplittingsNNLO.cpp  src/sudakovSpline.cpp        src/tmd.cpp src/integration.cpp  src/Spline.cpp  src/SplittingsNLO.cpp  src/Sudakov.cpp qcdnum/pij_nlo.f qcdnum/xpij2p.f qcdnum/xpns2p.f  qcdnum/ome.f qcdnum/wgplg.f )
add_executable(mcEvol ${SOURCES})
target_compile_features(mcEvol PRIVATE cxx_range_for cxx_override cxx_nullptr cxx_lambdas)

#find_library(yodaLIB libYODA.so PATHS ${YODAdir}/lib)

target_link_libraries(mcEvol  -lYODA)
