cmake_minimum_required (VERSION 3.4)
project (GOTHiC CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory (src) 
#add_subdirectory (test)

#enable_testing ()
#add_test (NAME MyTest COMMAND test)
