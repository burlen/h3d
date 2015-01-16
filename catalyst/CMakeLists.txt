project(UH3DAdaptor)
cmake_minimum_required(VERSION 2.6)

include_directories(${FortranAdaptor_BINARY_DIR} ${FortranAdaptor_SOURCE_DIR}
  ${UH3DAdaptor_BINARY_DIR})

set(UH3DAdaptor_SRCS
  UH3DAdaptor.cxx
)
add_library(UH3DAdaptor ${UH3DAdaptor_SRCS})
target_link_libraries(
  UH3DAdaptor FortranAdaptor
)

