find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)                      
                      
# Gather the header and source files
file(GLOB UNIT_TESTS_INC_FILES ${PROJECT_SOURCE_DIR}/test/*.hpp)
file(GLOB UNIT_TESTS_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)

# Include paths
set(UNIT_TESTS_INCLUDES ${PROJECT_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR})

# Assign the include directories
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${UNIT_TESTS_INCLUDES})

# Build unit tests
add_definitions(-DBOOST_TEST_DYN_LINK)

  message(STATUS "found boost")

set(LIBS sysTests hicupTest  overlapTest importHicupTest binInterTest calcFreqTest pBhAdjustTest)
set(LIBS ${LIBS} sortPositionsTest countDupsTest removeDupsTest findCisTest binomialTest binomTestTest)
set(LIBS ${LIBS} findTransTest fixChromoTest halfInteractionTest hicupCompInterTest)
set(LIBS ${LIBS} sumSquareTest IHWTest)

foreach(lib ${LIBS})
  add_executable(${lib} ${lib}.cpp)
  target_link_libraries (${lib} ${Boost_FILESYSTEM_LIBRARY}
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                       hicupData BinomData Utils Interactions)
  string(REPLACE "Test" "_Test" name ${lib})
  message(STATUS ${name})
  add_test(${name} ${lib})
endforeach(lib)
