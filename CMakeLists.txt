cmake_minimum_required(VERSION 3.16)
project(bl4s_bismuth_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Geant4's "ui_all vis_all" components need Qt + OpenGL + X11, which conda-forge
# macOS builds do not provide. Build headless physics; opt in to the viewer.
option(BL4S_VIS "Build the interactive Qt/OpenGL viewer" OFF)

# CMake's FindEXPAT/FindZLIB prefer the macOS SDK stubs over the conda prefix,
# which makes Geant4's find_dependency(EXPAT 2.8.5) / (ZLIB 1.3.2) version
# checks fail. Point them at the active environment when one is present.
if(DEFINED ENV{CONDA_PREFIX})
  set(_conda "$ENV{CONDA_PREFIX}")
  list(PREPEND CMAKE_PREFIX_PATH "${_conda}")
  find_library(_expat NAMES expat PATHS "${_conda}/lib" NO_DEFAULT_PATH)
  find_path(_expat_inc NAMES expat.h PATHS "${_conda}/include" NO_DEFAULT_PATH)
  find_library(_zlib NAMES z PATHS "${_conda}/lib" NO_DEFAULT_PATH)
  find_path(_zlib_inc NAMES zlib.h PATHS "${_conda}/include" NO_DEFAULT_PATH)
  if(_expat AND _expat_inc)
    set(EXPAT_LIBRARY "${_expat}")
    set(EXPAT_INCLUDE_DIR "${_expat_inc}")
  endif()
  if(_zlib AND _zlib_inc)
    set(ZLIB_LIBRARY "${_zlib}")
    set(ZLIB_INCLUDE_DIR "${_zlib_inc}")
  endif()
endif()

if(BL4S_VIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})

add_executable(bl4s_sim src/sim.cc)
target_link_libraries(bl4s_sim ${Geant4_LIBRARIES})
if(BL4S_VIS)
  target_compile_definitions(bl4s_sim PRIVATE BL4S_VIS)
endif()
