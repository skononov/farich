#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(FARICH)

#----------------------------------------------------------------------------
# Where to find Geant4 and ROOT installations
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{G4INSTALL})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Redefine compiler flags as needed for ROOT
#
set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -Wno-non-virtual-dtor\
 -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual\ 
 -Wno-variadic-macros -Wshadow -pipe -DG4USE_STD14\
 -static-libstdc++ -pipe -m64 -fsigned-char -fPIC -pthread\ 
 -DR__HAVE_CONFIG -O3 -DNDEBUG -fno-trapping-math\
 -ftree-vectorize -fno-math-errno -std=c++14")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(aRichSim aRichSim.cc ${sources} ${headers})
target_link_libraries(aRichSim ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS aRichSim DESTINATION bin)
