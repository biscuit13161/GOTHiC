#find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(Boost 1.36.0 COMPONENTS iostreams)
#find_package (gzstream 1.5 REQUIRED "${CMAKE_SOURCE_DIR}/../gzstream")

set(LIBS Setup Utils hicupData BinomData)

foreach(lib ${LIBS})
  add_library(${lib} ${lib}.cpp)
  target_link_libraries (${lib} ${Boost_LIBRARIES})# seqan3::seqan3)
endforeach(lib)


add_executable (gothic gothic.cpp)
 
target_link_libraries (gothic ${LIBS} ${Boost_LIBRARIES})


