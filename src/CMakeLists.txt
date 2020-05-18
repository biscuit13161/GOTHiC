#find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(Boost 1.36.0 COMPONENTS iostreams)
#find_package (gzstream 1.5 REQUIRED "${CMAKE_SOURCE_DIR}/../gzstream")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(GOTH_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set(GOTH_MINOR_VERSION ${CMAKE_PROJECT_VERSION_MINOR})
set(GOTH_PATCH_VERSION ${CMAKE_PROJECT_VERSION_PATCH}) 
configure_file("version.h.in" "version.h" @ONLY)

set(CLIBS pbinom toms708)

foreach(lib ${CLIBS})
  add_library(${lib} ${lib}.c)
  target_link_libraries (${lib} ${Boost_LIBRARIES})
endforeach(lib)

set(LIBS Setup Utils hicupData BinomData binTest )

foreach(lib ${LIBS})
  add_library(${lib} ${lib}.cpp)
  target_link_libraries (${lib} ${CLIBS} ${Boost_LIBRARIES})
endforeach(lib)


add_executable (gothic gothic.cpp)
 
target_link_libraries (gothic ${CLIBS} ${LIBS} ${Boost_LIBRARIES})



