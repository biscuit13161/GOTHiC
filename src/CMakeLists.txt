set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/modules")
#find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package(Boost 1.36.0 COMPONENTS iostreams)
#find_package (gzstream 1.5 REQUIRED "${CMAKE_SOURCE_DIR}/../gzstream")
find_package(TBB REQUIRED)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbbmalloc_proxy -ltbbmalloc")

message(STATUS "Source Path: ${CMAKE_SOURCE_DIR}")
message(STATUS "Module Path: ${CMAKE_MODULE_PATH}")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(GOTH_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set(GOTH_MINOR_VERSION ${CMAKE_PROJECT_VERSION_MINOR})
set(GOTH_PATCH_VERSION ${CMAKE_PROJECT_VERSION_PATCH}) 
configure_file("version.h.in" "version.h" @ONLY)

list (APPEND CLIBS dbinom pbinom toms708 bd0 lgamma stirlerr gamma)

foreach(lib ${CLIBS})
  list (SUBLIST CLIBS 0 -1 xlibs)
  list (REMOVE_ITEM xlibs ${lib}) 
  add_library(${lib} ${lib}.c)
  target_link_libraries (${lib} ${xlibs} ${Boost_LIBRARIES})
endforeach(lib)

list (APPEND LIBS Site padjust binomTest Utils UtilsComp)
list (APPEND LIBS Interactions SetupData SetupComp hicupData hicupDataComp)
list (APPEND LIBS BinomData BinomDataComp binTest Baits )

foreach(lib ${LIBS})
  ## copy list (zero-based, -1 refers to end) ##
  list (SUBLIST LIBS 0 -1 xlibs)
  ##  remove current library from list ##
  list (REMOVE_ITEM xlibs ${lib}) 
  #message(STATUS "xlibs = ${xlibs}")
  add_library(${lib} ${lib}.cpp)
  target_link_libraries (${lib} ${CLIBS} ${xlibs} ${Boost_LIBRARIES})
endforeach(lib)

list(APPEND EXE gothic gothicomp)

# Set LD_PRELOAD so that loader loads release version of proxy
set_property(GLOBAL PROPERTY ENVIRONMENT LD_PRELOAD=libtbbmalloc_proxy.so.2)

foreach(exe ${EXE})
  add_executable (${exe} ${exe}.cpp)
  target_link_libraries (${exe} ${CLIBS} ${LIBS} ${Boost_LIBRARIES} -flifetime-dse=1 -ltbb ${TBB_IMPORTED_TARGETS})
endforeach(exe)

add_custom_command(TARGET gothic POST_BUILD
                   COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/hicupToTable.sh ${CMAKE_CURRENT_BINARY_DIR}/. )
                   
install(TARGETS ${EXE}
    DESTINATION ${INSTALL_BIN_DIR}
    )
                   