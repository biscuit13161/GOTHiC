cmake_minimum_required (VERSION 3.4)
cmake_policy(SET CMP0074 NEW)

project (GOTHiC VERSION 0.1.3 HOMEPAGE_URL https://github.com/biscuit13161/GOTHiC LANGUAGES C CXX)

include(CTest) # note: this adds a BUILD_TESTING which defaults to ON

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp ")

set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

# to allow verbose output
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG_FLAG ")

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory (src) 
if(BUILD_TESTING)
  #message(STATUS "testing")
  add_subdirectory(tests)
endif()

## print list of cmake variables
get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()


