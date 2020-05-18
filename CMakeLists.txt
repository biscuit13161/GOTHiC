cmake_minimum_required (VERSION 3.4)
project (GOTHiC VERSION 0.0.2 HOMEPAGE_URL https://github.com/biscuit13161/GOTHiC LANGUAGES C CXX)

include(CTest) # note: this adds a BUILD_TESTING which defaults to ON

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory (src) 
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()


