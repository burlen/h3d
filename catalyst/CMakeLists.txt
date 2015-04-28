cmake_minimum_required(VERSION 2.6)
project(UH3DAdaptor CXX Fortran)

find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)

INCLUDE(${PARAVIEW_USE_FILE})
if(NOT PARAVIEW_USE_MPI)
  message(SEND_ERROR "ParaView must be built with MPI enabled")
endif()
if(NOT MPI_Fortran_LIBRARIES)
  find_package(MPI)
endif()

configure_file(run.sh.in run.sh)

include_directories(${UH3DAdaptor_BINARY_DIR})
configure_file(insitu.py insitu.py COPYONLY)

set(UH3DAdaptor_SRCS
  UH3DAdaptor.cxx
)

add_library(UH3DAdaptor
  ${UH3DAdaptor_SRCS})
target_link_libraries(UH3DAdaptor
  vtkPVPythonCatalyst)
