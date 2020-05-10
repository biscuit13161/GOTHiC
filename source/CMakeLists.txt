cmake_minimum_required (VERSION 3.4)
project (GOTHiC CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
#-D_GLIBCXX_PARALLEL")
 
find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
#find_package(Boost 1.36.0 COMPONENTS iostreams)
#find_package (gzstream 1.5 REQUIRED "${CMAKE_SOURCE_DIR}/../gzstream")
 
include_directories(${Boost_INCLUDE_DIRS})

#add_executable (hello_world hello_world.cpp)

set(LIBS hicupData BinomData)

foreach(lib ${LIBS})
  add_library(${lib} ${lib}.cpp)
  target_link_libraries (${lib} ${Boost_LIBRARIES} seqan3::seqan3)
endforeach(lib)


add_executable (gothic gothic.cpp)
 
#target_link_libraries (hello_world seqan3::seqan3)
target_link_libraries (gothic seqan3::seqan3 ${LIBS} ${Boost_LIBRARIES})
#target_link_libraries (gothic libgzstream.a)
#arget_link_libraries (gothic -L${CMAKE_SOURCE_DIR}/../gzstream/)
