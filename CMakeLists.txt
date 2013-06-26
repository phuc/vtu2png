
cmake_minimum_required(VERSION 2.8)

PROJECT(vtu2png)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(vtu2png MACOSX_BUNDLE vtu2png)

if(VTK_LIBRARIES)
  target_link_libraries(vtu2png ${VTK_LIBRARIES})
else()
  target_link_libraries(vtu2png vtuHybrid)
endif()
